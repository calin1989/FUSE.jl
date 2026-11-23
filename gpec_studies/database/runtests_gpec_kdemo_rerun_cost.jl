# Add ActorWholeFacility (engineering + costing) to all existing cases.
# Loads dd.equilibrium + dd.core_profiles from each dd.json, runs
# ActorWholeFacility with update_plasma=false, and writes dd.costing + dd.build
# back into the file.  Does NOT touch dd.mhd_linear (GPEC results are preserved).
#
# Run AFTER runtests_gpec_kdemo_cost_test.jl confirms the single-case test passes.
#
# Usage (on Feynman):
#   nohup julia --project=/mnt/homes_global/cap2297/FUSE.jl \
#     test/runtests_gpec_kdemo_rerun_cost.jl > ~/rerun_cost.log 2>&1 &

using Distributed
using FUSE
import IMAS
import JSON

const N_WORKERS   = 32
const SAVE_FOLDER = joinpath(homedir(), "gpec_kdemo_database")

FUSE.parallel_environment("localhost", N_WORKERS; memory_usage_fraction=0.5)

@everywhere begin
    using FUSE
    import IMAS
    import JSON

    function add_cost_one(folder::String)
        dd_file = joinpath(folder, "dd.json")
        !isfile(dd_file) && return (:no_file, NaN, "")

        d = open(JSON.parse, dd_file)

        # skip cases that never had a successful plasma solve
        !haskey(d, "equilibrium")   && return (:skip, NaN, "no equilibrium")
        !haskey(d, "core_profiles") && return (:skip, NaN, "no core_profiles")
        !haskey(d, "mhd_linear")    && return (:skip, NaN, "no mhd_linear")

        # skip if cost already computed
        if haskey(d, "costing")
            tmp = IMAS.dd()
            IMAS.dict2imas(d["costing"], tmp.costing)
            if !ismissing(tmp.costing.cost_direct_capital, :cost) &&
               tmp.costing.cost_direct_capital.cost > 0.0
                return (:already_done, tmp.costing.cost_direct_capital.cost, "")
            end
        end

        # load original ini/act from case folder, init fresh dd to get dd.build
        _, ini, act = FUSE.load(folder; load_dd=false)
        case_dd = IMAS.dd()
        FUSE.init(case_dd, ini, act)

        # overwrite plasma with converged solved values from dd.json
        IMAS.dict2imas(d["equilibrium"],   case_dd.equilibrium)
        IMAS.dict2imas(d["core_profiles"], case_dd.core_profiles)
        haskey(d, "wall") && IMAS.dict2imas(d["wall"], case_dd.wall)

        # run engineering + costing (no plasma re-solve)
        act.ActorWholeFacility.update_plasma = false

        try
            FUSE.ActorWholeFacility(case_dd, act)
            cost = case_dd.costing.cost_direct_capital.cost

            # write costing and build back into the file
            d["costing"] = IMAS.imas2dict(case_dd.costing)
            d["build"]   = IMAS.imas2dict(case_dd.build)
            open(dd_file, "w") do f
                JSON.print(f, d)
            end

            return (:ok, cost, "")
        catch e
            return (:crashed, NaN, sprint(showerror, e))
        end
    end
end

# ── collect all case folders ──────────────────────────────────────────────────
case_folders = sort(filter(
    f -> isdir(f) && occursin(r"^\d+__", basename(f)),
    readdir(SAVE_FOLDER; join=true)
))
println("Found $(length(case_folders)) case folders in $SAVE_FOLDER")
println("Adding cost data on $(N_WORKERS) workers...\n")

# ── parallel rerun ────────────────────────────────────────────────────────────
results = pmap(add_cost_one, case_folders)

# ── summary ───────────────────────────────────────────────────────────────────
n_ok           = count(r -> r[1] == :ok,           results)
n_already_done = count(r -> r[1] == :already_done, results)
n_skip         = count(r -> r[1] == :skip,         results)
n_crashed      = count(r -> r[1] == :crashed,      results)

println("\n=== Summary ===")
println("  Updated      : $n_ok")
println("  Already done : $n_already_done")
println("  Skipped      : $n_skip  (no plasma data)")
println("  Crashed      : $n_crashed")

costs_ok = [r[2] for r in results if r[1] ∈ (:ok, :already_done) && !isnan(r[2])]
if !isempty(costs_ok)
    println("\nCost statistics across successful cases:")
    println("  Min  : $(round(minimum(costs_ok); digits=1)) \$M")
    println("  Max  : $(round(maximum(costs_ok); digits=1)) \$M")
    println("  Mean : $(round(sum(costs_ok)/length(costs_ok); digits=1)) \$M")
end

if n_crashed > 0
    println("\nCrashed cases:")
    for (folder, r) in zip(case_folders, results)
        r[1] == :crashed && println("  $(basename(folder)): $(r[3])")
    end
end

println("\nDone. Run runtests_gpec_kdemo_extract.jl to rebuild CSV with cost column.")

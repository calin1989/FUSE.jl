# Re-run ActorGPEC with psihigh=0.99 on all existing successful cases.
# Loads only dd.equilibrium from each dd.json — no transport re-solve.
# Overwrites dd.mhd_linear in each dd.json with the corrected δW values.
#
# Usage (on Feynman):
#   nohup julia --project=/mnt/homes_global/cap2297/FUSE.jl \
#     test/runtests_gpec_kdemo_rerun_psihigh.jl > ~/rerun_psihigh.log 2>&1 &

using Distributed
using FUSE
import IMAS
import JSON

const PSIHIGH    = 0.99
const N_WORKERS  = 32
const SAVE_FOLDER = joinpath(homedir(), "gpec_kdemo_database")

FUSE.parallel_environment("localhost", N_WORKERS; memory_usage_fraction=0.5)

@everywhere begin
    using FUSE
    import IMAS
    import JSON

    function rerun_one(folder::String, psihigh::Float64)
        dd_file = joinpath(folder, "dd.json")
        !isfile(dd_file) && return (:no_file, NaN, NaN, "")

        d = open(JSON.parse, dd_file)

        # skip cases that never had GPEC output
        !haskey(d, "mhd_linear") && return (:skip, NaN, NaN, "")
        tmp = IMAS.dd()
        IMAS.dict2imas(d["mhd_linear"], tmp.mhd_linear)
        if isempty(tmp.mhd_linear.time_slice) ||
           isempty(tmp.mhd_linear.time_slice[1].toroidal_mode)
            return (:skip, NaN, NaN, "")
        end
        old_dW = real(tmp.mhd_linear.time_slice[1].toroidal_mode[1].energy_perturbed)

        # build dd with the fields GPEC needs
        case_dd = IMAS.dd()
        IMAS.dict2imas(d["equilibrium"], case_dd.equilibrium)
        haskey(d, "wall")          && IMAS.dict2imas(d["wall"],          case_dd.wall)
        haskey(d, "core_profiles") && IMAS.dict2imas(d["core_profiles"], case_dd.core_profiles)

        _, act = FUSE.case_parameters(:KDEMO)
        act.ActorGPEC.nn_low   = 1
        act.ActorGPEC.nn_high  = 1
        act.ActorGPEC.vac_flag = true
        act.ActorGPEC.verbose  = false
        act.ActorGPEC.psihigh  = psihigh

        try
            FUSE.ActorGPEC(case_dd, act)
            new_dW = real(case_dd.mhd_linear.time_slice[1].toroidal_mode[1].energy_perturbed)

            # overwrite mhd_linear in dd.json
            d["mhd_linear"] = IMAS.imas2dict(case_dd.mhd_linear)
            open(dd_file, "w") do f
                JSON.print(f, d)
            end

            return (:ok, old_dW, new_dW, "")
        catch e
            return (:crashed, old_dW, NaN, sprint(showerror, e))
        end
    end
end

# ── collect all case folders ──────────────────────────────────────────────────
case_folders = sort(filter(
    f -> isdir(f) && occursin(r"^\d+__", basename(f)),
    readdir(SAVE_FOLDER; join=true)
))
println("Found $(length(case_folders)) case folders in $SAVE_FOLDER")
println("Re-running ActorGPEC with psihigh=$PSIHIGH on $(N_WORKERS) workers...\n")

# ── parallel re-run ───────────────────────────────────────────────────────────
results = pmap(f -> rerun_one(f, PSIHIGH), case_folders)

# ── summarise ─────────────────────────────────────────────────────────────────
n_ok      = count(r -> r[1] == :ok,      results)
n_skip    = count(r -> r[1] == :skip,    results)
n_crashed = count(r -> r[1] == :crashed, results)

println("\n=== Summary ===")
println("  Updated : $n_ok")
println("  Skipped : $n_skip  (no prior GPEC output)")
println("  Crashed : $n_crashed")

if n_crashed > 0
    println("\nCrashed cases:")
    for (folder, r) in zip(case_folders, results)
        r[1] == :crashed && println("  $(basename(folder)): $(r[4])")
    end
end

println("\nDone. Run runtests_gpec_kdemo_extract.jl to rebuild the CSV and plots.")

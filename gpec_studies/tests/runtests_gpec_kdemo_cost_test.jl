# Test ActorWholeFacility (costing) on the first existing case that has mhd_linear.
# Loads dd.equilibrium + dd.core_profiles from dd.json, runs ActorWholeFacility
# with update_plasma=false (engineering only, no plasma re-solve), and prints the
# resulting capital cost.  Run this before the full rerun to confirm it works.
#
# Usage (on Feynman):
#   julia --project=/mnt/homes_global/cap2297/FUSE.jl \
#     test/runtests_gpec_kdemo_cost_test.jl

using FUSE
import IMAS
import JSON

const SAVE_FOLDER = joinpath(homedir(), "gpec_kdemo_database")

# ── find first case that has both equilibrium and mhd_linear ─────────────────
case_folders = sort(filter(
    f -> isdir(f) && occursin(r"^\d+__", basename(f)),
    readdir(SAVE_FOLDER; join=true)
))

test_folder = nothing
for folder in case_folders
    dd_file = joinpath(folder, "dd.json")
    !isfile(dd_file) && continue
    d = open(JSON.parse, dd_file)
    if haskey(d, "equilibrium") && haskey(d, "core_profiles") && haskey(d, "mhd_linear")
        global test_folder = folder
        break
    end
end

if test_folder === nothing
    error("No suitable case found in $SAVE_FOLDER")
end

println("Testing on: $(basename(test_folder))")

# ── load original ini/act from case folder + solved plasma from dd.json ───────
dd_file = joinpath(test_folder, "dd.json")
d = open(JSON.parse, dd_file)

# load the original ini and act that were used to generate this case
_, ini, act = FUSE.load(test_folder; load_dd=false)

# init a fresh dd from ini — this sets up dd.build with correct geometry
case_dd = IMAS.dd()
FUSE.init(case_dd, ini, act)

# overwrite plasma fields with the converged solved values from dd.json
IMAS.dict2imas(d["equilibrium"],   case_dd.equilibrium)
IMAS.dict2imas(d["core_profiles"], case_dd.core_profiles)
haskey(d, "wall") && IMAS.dict2imas(d["wall"], case_dd.wall)

println("Loaded ini from case folder, initialized build, restored solved plasma")
println("  R0  = $(round(case_dd.equilibrium.time_slice[1].global_quantities.magnetic_axis.r; digits=3)) m")
println("  B0  = $(round(case_dd.equilibrium.vacuum_toroidal_field.b0[1]; digits=3)) T")
println("  β_N = $(round(case_dd.equilibrium.time_slice[1].global_quantities.beta_normal; digits=3))")
println("  q95 = $(round(case_dd.equilibrium.time_slice[1].global_quantities.q_95; digits=3))")
println("  build layers: $(length(case_dd.build.layer))")

# ── run ActorWholeFacility (engineering + costing only) ───────────────────────
act.ActorWholeFacility.update_plasma = false  # plasma already solved — skip re-run

println("\nRunning ActorWholeFacility (engineering + costing)...")
try
    FUSE.ActorWholeFacility(case_dd, act)
    cost = case_dd.costing.cost_direct_capital.cost
    println("\n✓ ActorWholeFacility completed successfully")
    println("  Direct capital cost = $(round(cost; digits=2)) \$M")
    println("\nCost breakdown by system:")
    for sys in case_dd.costing.cost_direct_capital.system
        println("  $(sys.name): $(round(sys.cost; digits=2)) \$M")
    end
catch e
    println("\n✗ ActorWholeFacility crashed: $(sprint(showerror, e))")
    rethrow(e)
end

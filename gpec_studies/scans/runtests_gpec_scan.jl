# Parameter scan: δW(n=1) vs major radius R0
#
# Uses FUSE.StudyDatabaseGenerator to vary ini.equilibrium.R0 across
# a range of values, run ActorGPEC (n=1, no-wall) at each point, and
# plot δW vs R0.
#
# Expected result if IMAS path is working correctly:
#   - Smooth curve (R0 changes equilibrium geometry → changes stability)
#   - Noisy/jumping curve signals a bug in read_imas or the actor
#
# Usage (from FUSE.jl root):
#   julia --project=/Users/sophiapuscalau/Desktop/FUSE.jl test/runtests_gpec_scan.jl

using FUSE
import IMAS
import JSON
using Plots
using Printf
using Distributed

# ── base case ─────────────────────────────────────────────────────────────────
println("Loading D3D base case...")
ini, act = FUSE.case_parameters(:D3D, :default)

# Use scalar parameters (R0, a, κ, δ) to build the boundary instead of the ODS file,
# so that varying ini.equilibrium.R0 actually changes the equilibrium geometry
ini.equilibrium.boundary_from = :scalars

# GPEC: n=1 only, no-wall, quiet
act.ActorGPEC.nn_low      = 1
act.ActorGPEC.nn_high     = 1
act.ActorGPEC.vac_flag    = true
act.ActorGPEC.wall_shape  = "nowall"
act.ActorGPEC.verbose     = false
act.ActorGPEC.psihigh     = 0.95   # stay well within FUSE equilibrium domain
act.ActorGPEC.qlow        = 0.0    # no q threshold required

R0_base = ini.equilibrium.R0
println("Base major radius: R0 = $(round(R0_base; sigdigits=4)) m")

# ── parameter scan using ↔ notation ──────────────────────────────────────────
# Vary R0 ±25% around nominal in 11 uniform steps
R0_min = 0.88 * R0_base   # ±12%: stays within stable range for boundary_from=:scalars
R0_max = 1.12 * R0_base
ini.equilibrium.R0 = R0_base ↔ [R0_min, R0_max]

# ── study configuration ───────────────────────────────────────────────────────
sty = FUSE.study_parameters(:DatabaseGenerator)
sty.server                   = "localhost"
sty.n_workers                = 0          # 0 = no distributed, runs serially
sty.n_simulations            = 11
sty.save_folder              = let d = joinpath(homedir(), "gpec_R0_scan_results"); mkpath(d); d end
sty.file_save_mode           = :overwrite
sty.database_policy          = :separate_folders
sty.release_workers_after_run = true

println("Results will be saved to: $(sty.save_folder)")

# ── workflow: init + ActorGPEC ────────────────────────────────────────────────
@everywhere function workflow_DatabaseGenerator(dd, ini, act)
    FUSE.init(dd, ini, act)
    FUSE.ActorGPEC(dd, act)
    return nothing
end

# ── run ───────────────────────────────────────────────────────────────────────
println("\nRunning StudyDatabaseGenerator ($(sty.n_simulations) points)...")
study = FUSE.StudyDatabaseGenerator(sty, ini, act)
study.workflow = workflow_DatabaseGenerator
FUSE.run(study)

# ── extract results from dataframe ───────────────────────────────────────────
println("\nExtracting results...")

R0_values = Float64[]
dW_values = Float64[]

for (k, folder) in enumerate(sort(filter(isdir, readdir(sty.save_folder; join=true))))
    dd_file = joinpath(folder, "dd.json")
    !isfile(dd_file) && (println("  Case $k: dd.json not found in $folder"); continue)
    try
        # Parse JSON, then load only equilibrium + mhd_linear into fresh IMAS structures.
        # This uses the IMAS data model (validates write_imas paths) without triggering
        # the type-convert crash that json2imas hits on unrelated FUSE fields in the full dd.
        d = open(JSON.parse, dd_file)
        case_dd = IMAS.dd()
        IMAS.dict2imas(d["equilibrium"], case_dd.equilibrium)
        IMAS.dict2imas(d["mhd_linear"],  case_dd.mhd_linear)

        R0 = case_dd.equilibrium.time_slice[1].global_quantities.magnetic_axis.r
        if length(case_dd.mhd_linear.time_slice) > 0 &&
           length(case_dd.mhd_linear.time_slice[1].toroidal_mode) > 0
            δW = real(case_dd.mhd_linear.time_slice[1].toroidal_mode[1].energy_perturbed)
            push!(R0_values, Float64(R0))
            push!(dW_values, Float64(δW))
            @printf "  Case %2d: R0 = %.3f m  →  δW(n=1) = %+.5f\n" k R0 δW
        else
            println("  Case $k: no mhd_linear output (GPEC may have failed)")
        end
    catch e
        println("  Case $k: failed to load — $(sprint(showerror, e))")
    end
end

println("\n$(length(dW_values))/$(sty.n_simulations) points succeeded.")

if length(dW_values) < 2
    error("Too few successful points to plot — check ActorGPEC setup.")
end

# sort by R0 for clean plot
perm = sortperm(R0_values)
R0_sorted = R0_values[perm]
dW_sorted = dW_values[perm]

# ── plot ──────────────────────────────────────────────────────────────────────
p = plot(
    R0_sorted, dW_sorted;
    xlabel     = "Major radius  R0  [m]",
    ylabel     = "δW(n=1)  [normalized]",
    title      = "Ideal MHD stability vs major radius\nD3D :default  |  n=1  |  no-wall",
    marker     = :circle,
    markersize = 5,
    lw         = 2,
    label      = "δW(n=1)",
    legend     = :topright,
)

hline!(p, [0.0];
    linestyle = :dash,
    color     = :red,
    lw        = 1.5,
    label     = "marginal stability  (δW = 0)",
)

vline!(p, [R0_base];
    linestyle = :dot,
    color     = :gray,
    lw        = 1.5,
    label     = "nominal R0 = $(round(R0_base; sigdigits=4)) m",
)

display(p)

println("\nInterpret the plot:")
println("  Smooth curve  → IMAS path and ActorGPEC working correctly end-to-end")
println("  Noisy curve   → check read_imas spline construction")
println("  δW crosses 0  → marginal stability radius found")

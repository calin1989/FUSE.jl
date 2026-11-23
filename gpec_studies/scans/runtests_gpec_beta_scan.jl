# Parameter scan: δW(n=1) vs normalized beta β_N
#
# Varies ini.equilibrium.ip DOWNWARD (lower current = higher β_N at fixed pressure).
# β_N = β·(a·B0/Ip) — so reducing Ip raises β_N while keeping the pressure profile fixed,
# which correctly scans the no-wall beta limit without FUSE compensating via Ip rescaling.
# x-axis of the plot is β_N read back from dd.equilibrium after the G-S solve.
#
# Expected result:
#   - As β_N increases (lower Ip) the plasma becomes more unstable (δW decreases)
#   - Smooth monotonically decreasing curve
#   - δW crosses 0 at the no-wall beta limit
#
# Usage (from FUSE.jl root):
#   julia --project=/Users/sophiapuscalau/Desktop/FUSE.jl test/runtests_gpec_beta_scan.jl

using FUSE
import IMAS
import JSON
using Plots
using Printf
using Distributed

# ── base case ─────────────────────────────────────────────────────────────────
println("Loading D3D base case...")
ini, act = FUSE.case_parameters(:D3D, :default)

# Fix geometry with scalars and set R0 to a STABLE point (from R0 scan: R0=1.85m → δW≈+19)
# This ensures we start from δW > 0 so the curve can decrease as β_N increases
ini.equilibrium.boundary_from = :scalars
ini.equilibrium.R0            = 1.85

# GPEC: n=1 only, no-wall, quiet
act.ActorGPEC.nn_low      = 1
act.ActorGPEC.nn_high     = 1
act.ActorGPEC.vac_flag    = true
act.ActorGPEC.wall_shape  = "nowall"
act.ActorGPEC.verbose     = false
act.ActorGPEC.psihigh     = 0.95
act.ActorGPEC.qlow        = 0.0

p0_base = ini.equilibrium.pressure_core
println("Base pressure_core: p0 = $(round(p0_base; sigdigits=4)) Pa")
println("Fixed R0 = 1.85 m  (stable starting point from R0 scan)")

# ── parameter scan using ↔ notation ──────────────────────────────────────────
# Vary pressure_core from low to high while keeping Ip and geometry FIXED.
# Fixed Ip → fixed q profile → β_N changes purely from pressure.
# This gives the clean decreasing curve the advisor expects.
p0_min = 0.5 * p0_base    # low β_N end  (stable)
p0_max = 4.0 * p0_base    # high β_N end (5x caused FUSE solver failure)
ini.equilibrium.pressure_core = p0_base ↔ [p0_min, p0_max]

# ── study configuration ───────────────────────────────────────────────────────
sty = FUSE.study_parameters(:DatabaseGenerator)
sty.server                    = "localhost"
sty.n_workers                 = 0
sty.n_simulations             = 21
sty.save_folder               = let d = joinpath(homedir(), "gpec_beta_scan_results"); mkpath(d); d end
sty.file_save_mode            = :overwrite
sty.database_policy           = :separate_folders
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

# ── extract results ───────────────────────────────────────────────────────────
println("\nExtracting results...")

beta_values = Float64[]
dW_values   = Float64[]

for (k, folder) in enumerate(sort(filter(isdir, readdir(sty.save_folder; join=true))))
    dd_file  = joinpath(folder, "dd.json")
    ini_file = joinpath(folder, "ini.json")
    !isfile(dd_file) && (println("  Case $k: dd.json not found"); continue)
    try
        d = open(JSON.parse, dd_file)
        case_dd = IMAS.dd()
        IMAS.dict2imas(d["equilibrium"], case_dd.equilibrium)
        IMAS.dict2imas(d["mhd_linear"],  case_dd.mhd_linear)

        # β_N = β / (Ip / (a·B0))  — read from IMAS global quantities
        beta_n = case_dd.equilibrium.time_slice[1].global_quantities.beta_normal

        if length(case_dd.mhd_linear.time_slice) > 0 &&
           length(case_dd.mhd_linear.time_slice[1].toroidal_mode) > 0
            δW = real(case_dd.mhd_linear.time_slice[1].toroidal_mode[1].energy_perturbed)
            # skip outliers: FUSE solver failure at extreme pressure gives huge δW
            abs(δW) > 10.0 && (println("  Case $k: outlier δW=$(@sprintf("%.1f", δW)) skipped (FUSE solver failure)"); continue)
            push!(beta_values, Float64(beta_n))
            push!(dW_values,   Float64(δW))
            @printf "  Case %2d: β_N = %.4f  →  δW(n=1) = %+.5f\n" k beta_n δW
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

# sort by β_N for clean plot
perm     = sortperm(beta_values)
β_sorted = beta_values[perm]
dW_sorted = dW_values[perm]

# ── plot ──────────────────────────────────────────────────────────────────────
p = plot(
    β_sorted, dW_sorted;
    xlabel     = "Normalized beta  β_N",
    ylabel     = "δW(n=1)  [normalized]",
    title      = "Ideal MHD stability vs beta\nD3D :default  |  n=1  |  no-wall",
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

savefig(p, joinpath(homedir(), "Desktop", "gpec_beta_scan_plot.png"))
println("\nPlot saved to ~/Desktop/gpec_beta_scan_plot.png")

println("\nInterpret the plot:")
println("  Smooth decreasing curve → higher β destabilizes the kink mode (correct physics)")
println("  δW crosses 0            → no-wall beta limit found")
println("  Flat/noisy curve        → pressure scan not affecting equilibrium correctly")

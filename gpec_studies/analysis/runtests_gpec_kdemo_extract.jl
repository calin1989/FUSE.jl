# Extract results from test/runtests_gpec_kdemo_database.jl
#
# Iterates the per-case folders written by StudyDatabaseGenerator
# (database_policy = :separate_folders), pulls the *realized* operating-point
# quantities (β_N, q95, B0, R0, ne, Zeff — read back from the post-solve dd, not
# from ini) plus δW(n=1) from each dd.json, and writes them to a CSV. Also makes
# a quick δW(n=1) vs β_N scatter plot, colored by stable/unstable, following the
# same diagnostic used in runtests_gpec_beta_scan.jl.
#
# Usage (from FUSE.jl root):
#   julia --project=/Users/sophiapuscalau/Desktop/FUSE.jl test/runtests_gpec_kdemo_extract.jl

import JSON
using DataFrames
using CSV
using Plots
using Printf
using Statistics

save_folder = joinpath(homedir(), "gpec_kdemo_database")

# ── load each case's dd.json and pull out the quantities of interest ──────────
rows = NamedTuple[]

# StudyDatabaseGenerator case folders are named "<item>__<timestamp>__<pid>";
# skip its "checkpoints" bookkeeping folder, which has no dd.json.
case_folders = sort(filter(f -> isdir(f) && occursin(r"^\d+__", basename(f)), readdir(save_folder; join=true)))
for (k, folder) in enumerate(case_folders)
    dd_file = joinpath(folder, "dd.json")
    !isfile(dd_file) && continue
    try
        # read everything directly from JSON — no IMAS at all, avoids expression bugs
        d = open(JSON.parse, dd_file)
        !haskey(d, "equilibrium")   && continue
        !haskey(d, "core_profiles") && continue
        !haskey(d, "mhd_linear")    && continue

        mhd_ts = get(d["mhd_linear"], "time_slice", [])
        if isempty(mhd_ts) || isempty(get(mhd_ts[1], "toroidal_mode", []))
            println("  Case $k ($(basename(folder))): no mhd_linear output (GPEC may have failed)")
            continue
        end

        eq_gq  = d["equilibrium"]["time_slice"][1]["global_quantities"]
        eq_p1d = d["equilibrium"]["time_slice"][1]["profiles_1d"]
        cp1d_j = d["core_profiles"]["profiles_1d"][1]
        tm_j   = d["mhd_linear"]["time_slice"][1]["toroidal_mode"][1]

        beta_n    = Float64(eq_gq["beta_normal"])
        B0        = Float64(d["equilibrium"]["vacuum_toroidal_field"]["b0"][1])
        R0        = Float64(eq_gq["magnetic_axis"]["r"])
        ne_core   = Float64(cp1d_j["electrons"]["density_thermal"][1])
        zeff_core = Float64(cp1d_j["zeff"][1])
        n_tor     = Int(tm_j["n_tor"])
        ep_raw    = tm_j["energy_perturbed"]
        dW        = isa(ep_raw, AbstractVector) ? Float64(first(ep_raw)) : Float64(ep_raw)

        # q95: interpolate q profile at rho_tor_norm = 0.95
        q95 = NaN
        if haskey(eq_p1d, "rho_tor_norm") && haskey(eq_p1d, "q")
            rho_prof = Float64.(eq_p1d["rho_tor_norm"])
            q_prof   = Float64.(eq_p1d["q"])
            idx = searchsortedlast(rho_prof, 0.95)
            q95 = if idx <= 0
                Float64(q_prof[1])
            elseif idx >= length(rho_prof)
                Float64(q_prof[end])
            else
                q_prof[idx] + (q_prof[idx+1] - q_prof[idx]) *
                    (0.95 - rho_prof[idx]) / (rho_prof[idx+1] - rho_prof[idx])
            end
        end

        # delta: triangularity at LCFS — last value of profiles_1d
        delta = NaN
        if haskey(eq_p1d, "triangularity_upper") && haskey(eq_p1d, "triangularity_lower")
            delta = 0.5 * (Float64(last(eq_p1d["triangularity_upper"])) +
                           Float64(last(eq_p1d["triangularity_lower"])))
        end

        # capital cost — sum subsystem costs
        cost_M = NaN
        if haskey(d, "costing")
            try
                total = 0.0
                for sys in d["costing"]["cost_direct_capital"]["system"]
                    haskey(sys, "subsystem") || continue
                    for sub in sys["subsystem"]
                        haskey(sub, "cost") && (total += Float64(sub["cost"]))
                    end
                end
                cost_M = total > 0.0 ? total : NaN
            catch
                cost_M = NaN
            end
        end

        push!(rows, (
            case      = k,
            folder    = basename(folder),
            beta_n    = beta_n,
            q95       = q95,
            B0        = B0,
            R0        = R0,
            ne_core   = ne_core,
            zeff_core = zeff_core,
            delta     = delta,
            n_tor     = n_tor,
            dW        = dW,
            cost_M    = cost_M,
        ))

        @printf "  Case %2d (%s): β_N=%.3f q95=%.3f δ=%.3f B0=%.2f R0=%.3f → δW(n=%d)=%+.4f\n" k basename(folder) beta_n q95 delta B0 R0 n_tor dW
    catch e
        println("  Case $k ($(basename(folder))): failed to load — $(sprint(showerror, e))")
    end
end

df = DataFrame(rows)
println("\n$(nrow(df)) / $(length(case_folders)) cases succeeded.")

# ── write CSV ───────────────────────────────────────────────────────────────
csv_path = joinpath(save_folder, "gpec_kdemo_results.csv")
CSV.write(csv_path, df)
println("Wrote $csv_path")

if nrow(df) < 1
    error("No successful cases — nothing to plot.")
end

# ── plot: δW(n=1) vs β_N, colored by stable/unstable ──────────────────────────
stable = df.dW .>= 0

# clip y-axis to the 2nd–98th percentile so outliers don't squash the bulk
dw_lo = quantile(df.dW, 0.02)
dw_hi = quantile(df.dW, 0.98)

p = scatter(
    df.beta_n[stable], df.dW[stable];
    label      = "stable (δW ≥ 0)",
    color      = :blue,
    markersize = 4,
    markerstrokewidth = 0,
    xlabel     = "Normalized beta  β_N",
    ylabel     = "δW(n=1)  [normalized]",
    title      = "KDEMO operating-point scan — ideal MHD stability\n($(nrow(df)) cases)",
    legend     = :topright,
    ylims      = (dw_lo, dw_hi),
)
scatter!(
    p, df.beta_n[.!stable], df.dW[.!stable];
    label      = "unstable (δW < 0)",
    color      = :red,
    markersize = 4,
    markerstrokewidth = 0,
)
hline!(p, [0.0]; linestyle=:dash, color=:black, lw=1.5, label="marginal stability (δW = 0)")

plot_dir  = let d = joinpath(homedir(), "Desktop"); isdir(d) ? d : homedir() end
plot_path = joinpath(plot_dir, "gpec_kdemo_scan_plot.png")
savefig(p, plot_path)
println("Wrote $plot_path")

# ── plot 2: δW(n=1) vs β_N, colored by q95 ───────────────────────────────────
p2 = scatter(
    df.beta_n, df.dW;
    marker_z          = df.q95,
    color             = :viridis,
    colorbar_title    = "q₉₅",
    markersize        = 4,
    markerstrokewidth = 0,
    xlabel            = "Normalized beta  β_N",
    ylabel            = "δW(n=1)  [normalized]",
    title             = "KDEMO operating-point scan — colored by q₉₅\n($(nrow(df)) cases)",
    legend            = false,
    ylims             = (dw_lo, dw_hi),
)
hline!(p2, [0.0]; linestyle=:dash, color=:black, lw=1.5)

plot_path2 = joinpath(plot_dir, "gpec_kdemo_scan_plot_q95.png")
savefig(p2, plot_path2)
println("Wrote $plot_path2")

# ── plot 3: stability map — q95 vs β_N, colored by stable/unstable ────────────
p3 = scatter(
    df.beta_n[stable], df.q95[stable];
    label             = "stable (δW ≥ 0)",
    color             = :blue,
    markersize        = 3,
    markerstrokewidth = 0,
    markerstrokealpha = 0,
    alpha             = 0.5,
    xlabel            = "Normalized beta  β_N",
    ylabel            = "q₉₅",
    title             = "KDEMO stability map — n=1 ideal MHD\n($(nrow(df)) cases)",
    legend            = :topright,
)
scatter!(
    p3, df.beta_n[.!stable], df.q95[.!stable];
    label             = "unstable (δW < 0)",
    color             = :red,
    markersize        = 3,
    markerstrokewidth = 0,
    markerstrokealpha = 0,
    alpha             = 0.7,
)

plot_path3 = joinpath(plot_dir, "gpec_kdemo_stability_map.png")
savefig(p3, plot_path3)
println("Wrote $plot_path3")

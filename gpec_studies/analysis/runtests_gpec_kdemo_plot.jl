# Fast plotting from the pre-extracted CSV.
# Run this instead of runtests_gpec_kdemo_extract.jl whenever you want new plots —
# no dd.json parsing, starts in seconds.
#
# Usage:
#   julia --project=/mnt/homes_global/cap2297/FUSE.jl test/runtests_gpec_kdemo_plot.jl

using DataFrames
using CSV
using Plots
using Statistics

csv_path = joinpath(homedir(), "gpec_kdemo_database", "gpec_kdemo_results.csv")
df = CSV.read(csv_path, DataFrame)
println("Loaded $(nrow(df)) cases from $csv_path")

plot_dir = let d = joinpath(homedir(), "Desktop"); isdir(d) ? d : homedir() end
stable   = df.dW .>= 0
dw_lo    = quantile(df.dW, 0.02)
dw_hi    = quantile(df.dW, 0.98)

# ── plot 1: δW vs β_N, colored by stable/unstable ─────────────────────────────
p1 = scatter(
    df.beta_n[stable], df.dW[stable];
    label             = "stable (δW ≥ 0)",
    color             = :blue,
    markersize        = 4,
    markerstrokewidth = 0,
    alpha             = 0.5,
    xlabel            = "Normalized beta  β_N",
    ylabel            = "δW(n=1)  [normalized]",
    title             = "KDEMO operating-point scan — ideal MHD stability\n($(nrow(df)) cases)",
    legend            = :topright,
    ylims             = (dw_lo, dw_hi),
)
scatter!(p1, df.beta_n[.!stable], df.dW[.!stable];
    label="unstable (δW < 0)", color=:red, markersize=4, markerstrokewidth=0, alpha=0.7)
hline!(p1, [0.0]; linestyle=:dash, color=:black, lw=1.5, label="marginal stability (δW = 0)")
savefig(p1, joinpath(plot_dir, "gpec_kdemo_scan_plot.png"))
println("Wrote gpec_kdemo_scan_plot.png")

# ── plot 2: δW vs β_N, colored by q95 ────────────────────────────────────────
p2 = scatter(
    df.beta_n, df.dW;
    marker_z=df.q95, color=:viridis, colorbar_title="q₉₅",
    markersize=4, markerstrokewidth=0,
    xlabel="Normalized beta  β_N", ylabel="δW(n=1)  [normalized]",
    title="KDEMO operating-point scan — colored by q₉₅\n($(nrow(df)) cases)",
    legend=false, ylims=(dw_lo, dw_hi),
)
hline!(p2, [0.0]; linestyle=:dash, color=:black, lw=1.5)
savefig(p2, joinpath(plot_dir, "gpec_kdemo_scan_plot_q95.png"))
println("Wrote gpec_kdemo_scan_plot_q95.png")

# ── plot 3: stability map — q95 vs β_N, colored by stable/unstable ────────────
p3 = scatter(
    df.beta_n[stable], df.q95[stable];
    label="stable (δW ≥ 0)", color=:blue,
    markersize=3, markerstrokewidth=0, alpha=0.5,
    xlabel="Normalized beta  β_N", ylabel="q₉₅",
    title="KDEMO stability map — n=1 ideal MHD\n($(nrow(df)) cases)",
    legend=:topright,
)
scatter!(p3, df.beta_n[.!stable], df.q95[.!stable];
    label="unstable (δW < 0)", color=:red, markersize=3, markerstrokewidth=0, alpha=0.7)
savefig(p3, joinpath(plot_dir, "gpec_kdemo_stability_map.png"))
println("Wrote gpec_kdemo_stability_map.png")

# ── plot 4: cost vs q95, colored by stable/unstable ──────────────────────────
# only plot cases where cost was computed (not NaN)
has_cost = .!isnan.(df.cost_M)
if sum(has_cost) > 0
    df_cost = df[has_cost, :]
    stable_c = df_cost.dW .>= 0
    p4 = scatter(
        df_cost.q95[stable_c], df_cost.cost_M[stable_c];
        label="stable (δW ≥ 0)", color=:blue,
        markersize=4, markerstrokewidth=0, alpha=0.6,
        xlabel="q₉₅", ylabel="Direct capital cost  [\$M]",
        title="KDEMO stability vs cost — n=1 ideal MHD\n($(sum(has_cost)) cases with cost data)",
        legend=:topright,
    )
    scatter!(p4, df_cost.q95[.!stable_c], df_cost.cost_M[.!stable_c];
        label="unstable (δW < 0)", color=:red,
        markersize=4, markerstrokewidth=0, alpha=0.7)
    savefig(p4, joinpath(plot_dir, "gpec_kdemo_cost_stability.png"))
    println("Wrote gpec_kdemo_cost_stability.png")
else
    println("No cost data found — skipping cost plot (run runtests_gpec_kdemo_rerun_cost.jl first)")
end

# ── plot 5: stability map colored by δ (triangularity) ───────────────────────
if hasproperty(df, :delta)
    p5 = scatter(
        df.beta_n, df.q95;
        marker_z=df.delta, color=:plasma, colorbar_title="δ (triangularity)",
        markersize=3, markerstrokewidth=0, alpha=0.7,
        xlabel="Normalized beta  β_N", ylabel="q₉₅",
        title="KDEMO stability map — colored by triangularity δ\n($(nrow(df)) cases)",
        legend=false,
    )
    savefig(p5, joinpath(plot_dir, "gpec_kdemo_delta_map.png"))
    println("Wrote gpec_kdemo_delta_map.png")
else
    println("No delta column in CSV — re-run extract script first")
end

# ── plot 6: zoomed stability map with integer q95 lines ──────────────────────
# zoom into q95 = 1.5–6.5 where the resonant instability lives
zoom_mask = (df.q95 .>= 1.5) .& (df.q95 .<= 6.5)
df_zoom   = df[zoom_mask, :]
stable_z  = df_zoom.dW .>= 0

p6 = scatter(
    df_zoom.beta_n[stable_z], df_zoom.q95[stable_z];
    label="stable (δW ≥ 0)", color=:blue,
    markersize=3, markerstrokewidth=0, alpha=0.5,
    xlabel="Normalized beta  β_N", ylabel="q₉₅",
    title="KDEMO stability map — zoomed, integer q₉₅ resonances\n($(sum(zoom_mask)) cases with q₉₅ ∈ [1.5, 6.5])",
    legend=:topright,
    ylims=(1.5, 6.5),
)
scatter!(p6, df_zoom.beta_n[.!stable_z], df_zoom.q95[.!stable_z];
    label="unstable (δW < 0)", color=:red, markersize=3, markerstrokewidth=0, alpha=0.8)

# mark integer rational surfaces
for n in 2:6
    hline!(p6, [Float64(n)]; linestyle=:dash, color=:black, lw=1.2, label=(n==2 ? "integer q₉₅" : ""))
end

savefig(p6, joinpath(plot_dir, "gpec_kdemo_stability_map_zoomed.png"))
println("Wrote gpec_kdemo_stability_map_zoomed.png")

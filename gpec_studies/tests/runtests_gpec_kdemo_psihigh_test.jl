# Test psihigh=0.99 on the first 10 successful cases from the existing database.
# Loads only dd.equilibrium from each dd.json (no transport re-solve), re-runs
# ActorGPEC with psihigh=0.99, and compares new δW against the old psihigh=0.95 value.
#
# Usage:
#   julia --project=/mnt/homes_global/cap2297/FUSE.jl test/runtests_gpec_kdemo_psihigh_test.jl

using FUSE
import IMAS
import JSON
using Printf

save_folder = joinpath(homedir(), "gpec_kdemo_database")

# ── find first 10 successful cases ────────────────────────────────────────────
case_folders = sort(filter(f -> isdir(f) && occursin(r"^\d+__", basename(f)),
                           readdir(save_folder; join=true)))

test_cases = String[]
for folder in case_folders
    dd_file = joinpath(folder, "dd.json")
    !isfile(dd_file) && continue
    d = open(JSON.parse, dd_file)
    !haskey(d, "mhd_linear") && continue
    tmp = IMAS.dd()
    IMAS.dict2imas(d["mhd_linear"], tmp.mhd_linear)
    if !isempty(tmp.mhd_linear.time_slice) &&
       !isempty(tmp.mhd_linear.time_slice[1].toroidal_mode)
        push!(test_cases, folder)
    end
    length(test_cases) >= 10 && break
end

println("Found $(length(test_cases)) successful cases to test.\n")

# ── set up ActorGPEC parameters ───────────────────────────────────────────────
_, act = FUSE.case_parameters(:KDEMO)
act.ActorGPEC.nn_low   = 1
act.ActorGPEC.nn_high  = 1
act.ActorGPEC.vac_flag = true
act.ActorGPEC.verbose  = false
act.ActorGPEC.psihigh  = 0.99   # ← the value we are testing

println("Running ActorGPEC with psihigh=$(act.ActorGPEC.psihigh) ...\n")

# ── re-run GPEC on each case ──────────────────────────────────────────────────
n_ok = 0
n_fail = 0

for (k, folder) in enumerate(test_cases)
    dd_file = joinpath(folder, "dd.json")
    d = open(JSON.parse, dd_file)

    # read old δW (psihigh=0.95)
    tmp = IMAS.dd()
    IMAS.dict2imas(d["mhd_linear"], tmp.mhd_linear)
    old_dW = real(tmp.mhd_linear.time_slice[1].toroidal_mode[1].energy_perturbed)

    # build dd with just the fields GPEC needs
    case_dd = IMAS.dd()
    IMAS.dict2imas(d["equilibrium"], case_dd.equilibrium)
    haskey(d, "wall")          && IMAS.dict2imas(d["wall"],          case_dd.wall)
    haskey(d, "core_profiles") && IMAS.dict2imas(d["core_profiles"], case_dd.core_profiles)

    try
        FUSE.ActorGPEC(case_dd, act)
        new_dW = real(case_dd.mhd_linear.time_slice[1].toroidal_mode[1].energy_perturbed)
        @printf "  Case %2d: old δW(0.95)=%+8.4f  →  new δW(0.99)=%+8.4f  ✓\n" k old_dW new_dW
        global n_ok += 1
    catch e
        @printf "  Case %2d: CRASHED — %s\n" k sprint(showerror, e)
        global n_fail += 1
    end
end

println("\n$(n_ok)/$(length(test_cases)) passed, $(n_fail) crashed.")
if n_fail == 0
    println("psihigh=0.99 looks safe — ready to re-run the full database.")
else
    println("Crashes detected — check error messages above before proceeding.")
end

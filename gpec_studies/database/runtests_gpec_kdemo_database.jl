# ActorGPEC × KDEMO StudyDatabaseGenerator
#
# Couples ActorGPEC into a StudyDatabaseGenerator scan over a KDEMO operating-point
# space (B0, R0, Ip, shaping, density/Zeff, heating powers — see the ↔ ranges below).
#
# Phase A (n_simulations=3, n_workers=0, serial) confirmed the interface works:
# ActorGPEC runs cleanly right after FUSE.ActorStationaryPlasma on randomly-sampled
# KDEMO points, with no uncaught errors (after the read_imas fixes in JPEC).
#
# Phase C (this configuration): scale up to n_simulations=50, n_workers=4, still
# n=1 (external kink) only via act.ActorGPEC.nn_high=1, matching Phase A. Uses
# file_save_mode=:append so Phase A's 3 cases stay in save_folder (53 cases total
# once this finishes).
#
# Usage (from FUSE.jl root):
#   julia --project=/Users/sophiapuscalau/Desktop/FUSE.jl test/runtests_gpec_kdemo_database.jl

using FUSE
using Logging
using Distributed

FUSE.logging(Logging.Info; actors=Logging.Error)

ini, act = FUSE.case_parameters(:KDEMO)

### Ini settings — KDEMO operating-point scan space

resize!(ini.ec_launcher, 1)
ini.ec_launcher[1].power_launched = 1e7 ↔ [1e4, 5e7]
ini.ec_launcher[1].efficiency_conversion = 0.45
ini.ec_launcher[1].efficiency_transmission = 0.8

act.ActorSimpleEC.actuator[1].rho_0 = 0.2 ↔ [0.1, 0.9]

ini.ic_antenna[1].power_launched = 5.0e7 ↔ [1e4, 5e7]
ini.ic_antenna[1].efficiency_conversion = 0.6
ini.ic_antenna[1].efficiency_transmission = 0.8

act.ActorFluxMatcher.max_iterations = 500

ini.equilibrium.B0 = 10.0 ↔ [4.0, 15.0]
ini.equilibrium.ip = 11.0e6 ↔ [8.0e6, 22e6]
ini.equilibrium.R0 = 5.5 ↔ [4.0, 7.0]
ini.equilibrium.pressure_core = missing # don't change this

ini.equilibrium.δ = 0.4 ↔ [0.3, 0.7]
ini.equilibrium.κ = 1.8

ini.core_profiles.impurity = :Kr ↔ (:Kr, :Ar, :Ne, :Xe)
ini.core_profiles.ne_setting = :greenwald_fraction_ped
ini.core_profiles.ne_value = 0.5 ↔ [0.2, 1.0]
ini.core_profiles.zeff = 2.0 ↔ [1.5, 4.0]

### Act settings

act.ActorEquilibrium.model = :TEQUILA
act.ActorCoreTransport.model = :FluxMatcher
act.ActorFluxMatcher.algorithm = :simple

act.ActorPFdesign.model = :optimal
act.ActorFluxSwing.operate_oh_at_j_crit = true # this maximizes flattop inside the fluxswing actor

act.ActorWholeFacility.update_plasma = true

# This is handled by the constraint functions
act.ActorHFSsizing.error_on_performance = false
act.ActorHFSsizing.error_on_technology = false

act.ActorFluxMatcher.evolve_densities = :flux_match

act.ActorPedestal.model = :EPED
act.ActorFluxMatcher.rho_transport = 0.2:0.05:0.8
act.ActorTGLF.tglfnn_model = "sat0quench_em_d3d+mastu_azf+1"
act.ActorFluxMatcher.evolve_pedestal = true

### Study configuration — Phase C database run (parallel, 50 cases)

sty = FUSE.study_parameters(:DatabaseGenerator)
sty.server = "localhost" # change to a ClusterManager name to spread over nodes
sty.n_workers = 4 # parallel across 4 of this machine's 10 cores (Phase C)
sty.n_simulations = 50 # Phase C database (added to Phase A's 3 cases already in save_folder)
sty.save_folder = let d = joinpath(homedir(), "gpec_kdemo_database"); mkpath(d); d end
sty.file_save_mode = :append # :overwrite would rm -rf save_folder, deleting Phase A's 3 cases
sty.database_policy = :separate_folders # writes dd.json per case (needed to read back mhd_linear)

println("Results will be saved to: $(sty.save_folder)")

# The single-ini/act StudyDatabaseGenerator constructor (used below so each case gets
# its own rand(ini)/rand(act) draw from the ↔ ranges) does NOT call
# parallel_environment itself (only the Vector{ini}/Vector{act} constructor does) --
# provision the workers explicitly so Distributed.workers() has sty.n_workers entries
# before FUSE.run(study) asserts on it. memory_usage_fraction=0.2 keeps the per-worker
# --heap-size-hint (and thus 4 workers x hint) within this machine's 17GB RAM.
FUSE.parallel_environment(sty.server, sty.n_workers; memory_usage_fraction=0.2)

study = FUSE.StudyDatabaseGenerator(sty, ini, act)

@everywhere import FUSE
@everywhere ProgressMeter = FUSE.ProgressMeter

@everywhere function workflow_DatabaseGenerator(dd::FUSE.IMAS.dd, ini::FUSE.ParametersAllInits, act::FUSE.ParametersAllActors)
    FUSE.init(dd, ini, act)
    FUSE.ActorStationaryPlasma(dd, act) # consistent plasma profiles (transport, sources, equilibrium)

    # Couple ActorGPEC: evaluate ideal MHD stability at this operating point
    act.ActorGPEC.nn_low = 1
    act.ActorGPEC.nn_high = 1     # n=1 (external kink) only, consistent with Phase A
    act.ActorGPEC.vac_flag = true
    act.ActorGPEC.verbose = false
    act.ActorGPEC.psihigh = 0.95  # stay below GPEC's internal psiedge=0.99 (avoids findmax_dW_edge! crash)
    FUSE.ActorGPEC(dd, act)

    return nothing
end

study.workflow = workflow_DatabaseGenerator

println("\nRunning StudyDatabaseGenerator ($(sty.n_simulations) new points, $(sty.n_workers) workers)...")
FUSE.run(study)

println("\nDone. Inspect $(sty.save_folder) for per-case dd.json/log.txt/error.txt.")
println("Next: run test/runtests_gpec_kdemo_extract.jl to pull out δW(n=1) results (now ~53 cases).")

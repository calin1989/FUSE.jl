# Diagnostic: compare the FUSE-solved KDEMO equilibrium geometry (true magnetic
# axis location, R/Z grid extents) against what GPEC's read_imas sees, BEFORE
# calling the direct equilibrium solver. This is to understand why
# direct_position! ("Took too many iterations to get bz=0") fails for KDEMO.
#
# Usage (from FUSE.jl root):
#   julia --project=. test/runtests_gpec_kdemo_diag.jl

using FUSE
import GeneralizedPerturbedEquilibrium as GPEC
import Logging

FUSE.logging(Logging.Info; actors=Logging.Error)

ini, act = FUSE.case_parameters(:KDEMO)

# nominal (central) KDEMO operating point -- same as the email's scan, but
# without the ↔ ranges, so this single case is deterministic.
resize!(ini.ec_launcher, 1)
ini.ec_launcher[1].power_launched = 1e7
ini.ec_launcher[1].efficiency_conversion = 0.45
ini.ec_launcher[1].efficiency_transmission = 0.8
act.ActorSimpleEC.actuator[1].rho_0 = 0.2
ini.ic_antenna[1].power_launched = 5.0e7
ini.ic_antenna[1].efficiency_conversion = 0.6
ini.ic_antenna[1].efficiency_transmission = 0.8
act.ActorFluxMatcher.max_iterations = 500
ini.equilibrium.B0 = 10.0
ini.equilibrium.ip = 11.0e6
ini.equilibrium.R0 = 5.5
ini.equilibrium.δ = 0.4
ini.equilibrium.κ = 1.8
ini.core_profiles.impurity = :Kr
ini.core_profiles.ne_setting = :greenwald_fraction_ped
ini.core_profiles.ne_value = 0.5
ini.core_profiles.zeff = 2.0

act.ActorEquilibrium.model = :TEQUILA
act.ActorCoreTransport.model = :FluxMatcher
act.ActorFluxMatcher.algorithm = :simple
act.ActorPFdesign.model = :optimal
act.ActorFluxSwing.operate_oh_at_j_crit = true
act.ActorWholeFacility.update_plasma = true
act.ActorHFSsizing.error_on_performance = false
act.ActorHFSsizing.error_on_technology = false
act.ActorFluxMatcher.evolve_densities = :flux_match
act.ActorPedestal.model = :EPED
act.ActorFluxMatcher.rho_transport = 0.2:0.05:0.8
act.ActorTGLF.tglfnn_model = "sat0quench_em_d3d+mastu_azf+1"
act.ActorFluxMatcher.evolve_pedestal = true

println("boundary_from = ", ini.equilibrium.boundary_from)

dd = IMAS.dd()
FUSE.init(dd, ini, act)
FUSE.ActorStationaryPlasma(dd, act)

eqt = dd.equilibrium.time_slice[]
ro_true = eqt.global_quantities.magnetic_axis.r
zo_true = eqt.global_quantities.magnetic_axis.z
println("\n--- FUSE-solved equilibrium ---")
println("psi_axis     = ", eqt.global_quantities.psi_axis)
println("psi_boundary = ", eqt.global_quantities.psi_boundary)
println("magnetic_axis: R0 = $ro_true, Z0 = $zo_true")

println("\n--- GPEC read_imas (after grid_type.index==1 fix) ---")
eq_config = GPEC.Equilibrium.EquilibriumConfig(;
    eq_type="imas", eq_filename="from_fuse_dd", imas_cocos=11, psilow=0.01, psihigh=0.95)
run_input = GPEC.Equilibrium.read_imas(eq_config, dd)
println("DirectRunInput: rmin=$(run_input.rmin), rmax=$(run_input.rmax), zmin=$(run_input.zmin), zmax=$(run_input.zmax)")
println("psio = ", run_input.psio)
println("Ro_true within R grid? ", run_input.rmin <= ro_true <= run_input.rmax)
println("Zo_true within Z grid? ", run_input.zmin <= zo_true <= run_input.zmax)

# Walk the same +r search direct_position! uses, at z = (zmin+zmax)/2, and
# print bz along the way to see whether/where it crosses zero.
println("\n--- bz(r) scan at z = (zmin+zmax)/2 ---")
bfield = GPEC.Equilibrium.DirectBField()
sq_in_deriv = GPEC.Equilibrium.deriv1(run_input.sq_in)
z0 = (run_input.zmin + run_input.zmax) / 2
r0 = (run_input.rmin + run_input.rmax) / 2
dr = (run_input.rmax - run_input.rmin) / 20.0
for k in 0:10
    r = r0 + k * dr
    GPEC.Equilibrium.direct_get_bfield!(bfield, r, z0, run_input.psi_in, run_input.sq_in, sq_in_deriv, run_input.psio; derivs=1)
    println("  k=$k  r=$(round(r; sigdigits=6))  bz=$(bfield.bz)")
end

println("\n--- ActorGPEC ---")
act.ActorGPEC.nn_low = 1
act.ActorGPEC.nn_high = 1
act.ActorGPEC.vac_flag = true
act.ActorGPEC.verbose = false
act.ActorGPEC.psihigh = 0.95
FUSE.ActorGPEC(dd, act)
tm = dd.mhd_linear.time_slice[1].toroidal_mode[1]
println("n = $(tm.n), δW = $(real(tm.energy_perturbed))")

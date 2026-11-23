# GPEC studies with FUSE

Scripts used to build and analyse a database of n=1 ideal MHD stability for KDEMO
operating points, using `ActorGPEC` (FUSE PR #1101), plus the earlier validation and scan
scripts. They are shared for review and are not cleaned up for general use: paths point to
my home directories on my laptop and on Feynman, and the actor settings are the ones in use
over the summer (see "Settings used").

Run any script from the FUSE.jl root, e.g.

    julia --project=. gpec_studies/database/runtests_gpec_kdemo_database.jl

Developed with Julia 1.12.4. The summer runs used FUSE 1.1.0 (actor at commit `b2e28239a`);
this branch is on FUSE 1.2.1.

## KDEMO database workflow

| Step | Script | What it does | Output |
|---|---|---|---|
| 1 | `database/runtests_gpec_kdemo_database.jl` | `StudyDatabaseGenerator` scan of KDEMO (B0, Ip, R0, δ, density, Zeff, impurity, EC/IC power). Each case runs `ActorStationaryPlasma` (TEQUILA, FluxMatcher, EPED) and then `ActorGPEC` with `psihigh = 0.95`. | one folder per case in `~/gpec_kdemo_database/`, each with `dd.json` |
| 2 | `database/runtests_gpec_kdemo_rerun_psihigh.jl` | Re-runs `ActorGPEC` on the saved equilibria with `psihigh = 0.99` and overwrites `dd.mhd_linear`. No transport re-solve. Cases that crash at 0.99 keep their 0.95 result, so the database may mix the two. | updated `dd.json` |
| 3 | `database/runtests_gpec_kdemo_rerun_cost.jl` | Runs `ActorWholeFacility` with `update_plasma = false` on each case to add build and costing. Does not touch `dd.mhd_linear`. | `dd.costing`, `dd.build` added to `dd.json` |
| 4 | `analysis/runtests_gpec_kdemo_extract.jl` | Reads every `dd.json` and writes one row per case: β_N, q95, δ, B0, R0, core n_e, core Zeff, n, δW, capital cost. q95 is q interpolated at ρ_tor = 0.95 (not ψ_N). Also makes three quick plots. | `gpec_kdemo_results.csv` |
| 5 | `analysis/runtests_gpec_kdemo_plot.jl` | Plots from the CSV: δW vs β_N coloured by stability and by q95, the q95–β_N stability map, cost vs q95 coloured by stability, the δ-coloured map and the zoomed q95 map. | PNGs in `~/Desktop` (home directory if there is none) |

`tests/runtests_gpec_kdemo_psihigh_test.jl` and `tests/runtests_gpec_kdemo_cost_test.jl` are the
single-case checks run before steps 2 and 3. `tests/runtests_gpec_kdemo_diag.jl` compares the
FUSE (TEQUILA) equilibrium with what GPEC's `read_imas` sees; it was used to track down the
`direct_position!` failures on KDEMO.

The copy of step 1 here is the 50-case configuration. The full run (11,331 cases) used the copy
on Feynman with a larger `n_simulations`; I still need to sync that version.

## Settings used for the summer database

These are the `ActorGPEC` defaults of the time plus the overrides in the scripts. They are not
the current defaults in PR #1101.

| | Summer database | Current actor (PR #1101) |
|---|---|---|
| Modes | n = 1 | n = 1 |
| psihigh | 0.99 (0.95 for cases that crashed in step 2) | 0.994, edge scan from psiedge = 0.99 |
| psilow | 0.01 | 1e-4 |
| Wall | conformal, a = 0.2 minor radii (labelled metres in the old actor, but passed to GPEC as `a`) | none |
| Coordinates | Boozer | Hamada |
| GPEC | `calin1989/JPEC` branch `imas-v2`; the exact commit is in the Feynman Manifest | OpenFUSIONToolkit/GPEC (scripting API) |

## Data

On Feynman:

- `~/gpec_kdemo_database/`: 11,331 cases, each with `dd.json`
- `~/gpec_kdemo_database/gpec_kdemo_results.csv`: 10,549 cases. Step 4 skips cases with no
  `dd.json`, `equilibrium`, `core_profiles` or `mhd_linear`, and cases that fail to parse.

`summer_results/` has the plots from 22 July (δW vs β_N coloured by stability and by q95, the
stability map, cost vs q95). In that dataset 1,464 of the 10,549 cases have δW < 0. The 22 July
CSV predates the `delta` column, so the δ-coloured map is not among them.

## Before rerunning

- Results come from the summer actor settings above. The wall is stabilizing, so the unstable
  fraction is likely to go up with the current no-wall default.
- q95 spans 0.99 to 33 and δW goes down to about −3.4 × 10⁷. Both extremes are probably
  unconverged and need checking.
- The FUSE side of each case (equilibrium, profiles, cost) does not depend on the actor, so a
  rerun only needs to repeat the GPEC step on the saved `dd.json` files, as in step 2.

## Earlier scripts

| Script | Purpose |
|---|---|
| `scans/runtests_gpec_validation.jl` | Checks on the actor on D3D `:default`: it matches a direct `GPEC.main` run on the same `dd`, `write_imas` stores the eigenvalues correctly, and the n = 1 spectrum is stable. On the DIII-D-like g-file it gives et[1] ≈ 1.20 against the stored reference 1.707 (mpert = 34), so only the sign is checked there. |
| `scans/runtests_gpec_scan.jl` | δW(n = 1) vs R0 on D3D, 11 cases. |
| `scans/runtests_gpec_beta_scan.jl` | δW(n = 1) vs β_N on D3D by raising `pressure_core` from 0.5 to 4 times its base value at fixed Ip, 21 cases. Cases with \|δW\| > 10 are dropped as solver failures. |

Both scans use no wall and `psihigh = 0.95`, the workaround in place at the time. The header
comment of the β scan still describes an earlier version that lowered Ip.

The scan notebook is on the `gpec_stability_scan` branch of `calin1989/FuseExamples`.

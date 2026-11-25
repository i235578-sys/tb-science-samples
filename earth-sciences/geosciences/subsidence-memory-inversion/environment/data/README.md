# Input data and conventions

All task observations are synthetic, deterministic samples from a simplified aquifer-compaction and elastic-half-space model. They are not measurements of a real California basin and should not be interpreted as field estimates. The workflow is motivated by California DWR's public TRE ALTAMIRA InSAR subsidence monitoring and by published studies that combine InSAR deformation with groundwater-storage observations. The synthetic model makes the task self-contained and gives the verifier known, hidden observation labels.

## Forward model

For source block `j` and month `t`, `drawdown_m` is the water-level decline in metres relative to a reference head. Define the historical maximum drawdown

`peak_j(t) = max_{0 <= s <= t} drawdown_j(s)`.

The modeled compaction is

`C_j(t) = S_elastic_j * drawdown_j(t) + S_inelastic_j * peak_j(t)`.

Both coefficients are in metres of compaction per metre of drawdown (`m/m`) and must be in `[0, 0.004]`. The second term is a simple memory term: it remains after the head partly recovers. It is a task-specific approximation, not a complete constitutive law for an aquifer.

The surface vertical displacement at sensor `i` is

`u_i(t) = -1000 * sum_j G_ij * C_j(t) + b_i + a_t + noise_i(t)` millimetres,

where `b_i` and `a_t` are supplied site and epoch reference offsets. The Mogi point-volume-source gain is

`G_ij = (1 - nu) / pi * z_j * A_j / (r_ij^2 + z_j^2)^(3/2)`,

with Poisson ratio `nu = 0.25`, horizontal distance `r` and depth `z` in km, and block area `A` in km². This makes `G` dimensionless. Compression is negative vertical displacement.

## Files

- `source_blocks.csv`: source ID, centre coordinates, area, and depth.
- `sensor_locations.csv`: surface sensor ID and coordinates.
- `mogi_kernel.csv`: precomputed `G_ij` for every sensor/source pair.
- `head_history.csv`: monthly drawdown by source block.
- `epoch_reference_offsets.csv`: known epoch reference-frame offsets in mm.
- `site_reference_offsets.csv`: known sensor reference offsets in mm.
- `training_observations.csv`: observed vertical displacement in mm; use these rows for fitting.

Coordinates and distances use kilometres; time is an integer month index from 0 to 47. CSV files use UTF-8 and decimal values.

The underlying science is motivated by the [California DWR InSAR subsidence dataset](https://publish.data.ca.gov/dataset/tre-altamira-insar-subsidence-data) and the published [Sentinel-1/GRACE Tulare Basin study](https://doi.org/10.1038/s41598-022-07650-1). The DWR dataset-level metadata lists public use without restrictions. No DWR observations are copied into this synthetic task.

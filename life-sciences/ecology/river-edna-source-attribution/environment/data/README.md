# River eDNA task data

## Study provenance

These inputs derive from the Thur catchment fish eDNA and river-network analysis by Carraro, Blackman, and Altermatt (2023), *Modelling environmental DNA transport in rivers reveals highly resolved spatio-temporal biodiversity patterns*, *Scientific Reports* 13, 8854 ([doi:10.1038/s41598-023-35614-6](https://doi.org/10.1038/s41598-023-35614-6)). The authors' analysis repository is [lucarraro/eDITH_RiverDNA](https://github.com/lucarraro/eDITH_RiverDNA), pinned to commit `6ef61c76d69a739ffd508d66c89207ba8867ecb6`. The upstream license notice is included as `UPSTREAM_LICENSE.txt`. Inputs are local; no trial-time download is required.

This task estimates fish biodiversity from archived 12S marker read-abundance values. Values are continuous and intentionally retain their archived precision; they are not direct counts of fish. The three genera are `Barbus`, `Cottus`, and `Phoxinus`, each observed in spring, summer, and autumn. The original sample design has 73 stations overall and no summer sample at four stations. `training_observations.csv` contains visible measurements. A missing site-season row means no sample was collected; `read_abundance = 0` is an observed non-detection.

## Files

- `river_network.npz`: 1,908-reach network, standardized covariates, and seasonal hydraulic fields.
- `sample_sites.csv`: 73 eDNA stations and their mapped network node. Coordinates use the study's local projected coordinate system in metres.
- `training_observations.csv`: visible observations with columns `genus`, `season`, `site_id`, and `read_abundance`.
- `UPSTREAM_LICENSE.txt`: license notice distributed with the authors' analysis repository.

### `river_network.npz` arrays

| Array | Shape | Meaning |
|---|---:|---|
| `downstream_id` | `(1908,)` | Zero-based downstream reach ID; `-1` marks the outlet. |
| `upstream_to_downstream_order` | `(1907,)` | Each non-outlet reach once, in propagation order. |
| `outlet_id` | `(1,)` | Zero-based outlet reach ID. |
| `length_m` | `(1908,)` | Reach length in metres. |
| `width_m` | `(3, 1908)` | Reach width by season: spring, summer, autumn; metres. |
| `discharge_m3_s` | `(3, 1908)` | Discharge by season in cubic metres per second. |
| `velocity_m_s` | `(3, 1908)` | Water velocity by season in metres per second. |
| `covariates` | `(1908, 35)` | Standardized catchment/reach covariates. |
| `covariate_names` | `(35,)` | Names in the same order as `covariates` columns. |
| `reach_x_m`, `reach_y_m` | `(1908,)` each | Projected reach coordinates in metres. |
| `season_names` | `(3,)` | `spring`, `summer`, `autumn` in the same order as hydraulic arrays. |

Every reach array uses zero-based `node_id`; `sample_sites.csv.node_id` refers to that index. The network is directed and dendritic: every non-outlet reach has one downstream reach. Use `upstream_to_downstream_order` to propagate transport.

## Forward transport

For genus `g` and season `s`, let `p[r]` be local DNA production at reach `r`, `tau_hours` be the DNA decay time in hours, and `tau_seconds = 3600 * tau_hours`. Compute predicted concentration `C` as follows:

```text
C[r] = p[r] * width[r] * length[r] / discharge[r]
       * exp(-length[r] / (tau_seconds * velocity[r]))

for r in upstream_to_downstream_order:
    d = downstream_id[r]
    C[d] += C[r] * discharge[r] / discharge[d]
           * exp(-length[d] / (tau_seconds * velocity[d]))
```

The predicted mean read abundance at a station is `C[node_id]`. The reference solution parameterizes local production as `p[r] = p0 * exp(covariates[r] @ beta)` with shrinkage to limit overfitting. Other estimation methods are acceptable if their submitted maps predict held-out observations.

## Outcome score

The hidden outcome score uses natural-log `log1p` RMSE between predicted mean abundance and observed abundance, computed separately for each genus-season. The nine group scores are averaged equally, so a group with a larger raw abundance scale does not dominate. Both every group score and the macro-average must satisfy the hidden thresholds. The hidden sites and numeric thresholds are absent from the agent environment.

## Output schema

`/app/submission/production.csv` must have exactly these columns:

```text
genus,season,node_id,production_rate
```

It must contain one row for each of the 3 genera × 3 seasons × 1,908 reaches. `production_rate` is a nonnegative finite local source term in the scale implied by the read-abundance data and network inputs. `node_id` is the zero-based reach index.

`/app/submission/decay.csv` must have exactly these columns:

```text
genus,season,tau_hours
```

It must contain one row per genus-season. `tau_hours` is a positive finite decay time in hours.

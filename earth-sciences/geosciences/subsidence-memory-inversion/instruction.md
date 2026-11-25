# Infer spatial elastic and inelastic aquifer compaction

Estimate the two nonnegative compaction-coefficient maps that best explain the supplied synthetic InSAR-style deformation observations. Write `/app/submission/subsidence_parameters.csv` with one row for every source block and columns `source_id`, `elastic_m_per_m`, and `inelastic_m_per_m`. Coefficients must be finite and lie in `[0, 0.004]` metres of compaction per metre of drawdown.

The measurements combine reversible deformation with compaction that persists after water levels recover. The exact task model, units, input schemas, reference offsets, and Mogi half-space kernel are specified in `/app/data/README.md`. Fit the visible training observations and produce maps that predict deformation at held-out sensor blocks and later months. The verifier will recompute those predictions from your coefficient map. It scores held-out outcomes; it does not require the exact hidden coefficient values or a particular fitting method.

These are synthetic data generated from a simplified model inspired by satellite monitoring of groundwater subsidence. They are not observations or parameter estimates for a real basin.

The verifier computes RMSE in millimetres separately for held-out sensor locations and later months, then takes the equal-weight mean of those two RMSE values. Both split scores and their mean must meet the verifier's limits. You do not need to recover the generating coefficients exactly.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.

# Reconstruct fish eDNA sources across the Thur catchment

Estimate seasonal, reach-level local DNA production for *Barbus*, *Cottus*, and *Phoxinus* from the supplied 12S metabarcoding observations. A sample includes DNA produced at that reach and DNA transported from upstream. Account for the directed river network, seasonal discharge, width, reach length, and velocity when recovering each source map.

The authoritative input conventions and transport equation are in `/app/data/README.md`. Write `/app/submission/production.csv` with one finite nonnegative `production_rate` per `genus,season,node_id`, and `/app/submission/decay.csv` with one positive `tau_hours` per `genus,season`. Include all three genera, all three seasons, and all reaches.

The outcome is spatial prediction at held-out sampling sites. The verifier computes natural-log `log1p` RMSE between predicted and observed read abundance for each genus-season, then averages the nine scores equally. The task passes only if every genus-season score and the macro-average meet their hidden limits. Your maps should generalize to unsampled reaches; code style and the estimation algorithm are not graded.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.

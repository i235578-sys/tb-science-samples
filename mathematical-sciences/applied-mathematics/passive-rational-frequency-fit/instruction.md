# Fit compact passive rational models to frequency responses

Build a rational model for each supplied three-port frequency response and write `/app/submission/models.json`. Use the format and channel order in `/app/data/README.md`. Each case can have a different model order. Return finite coefficients and stable poles; your fitted responses must predict the hidden frequency grid, remain reciprocal, and satisfy the passivity check.

The input contains noisy samples over six frequency decades. For each case, the verifier computes complex-response RMSE and maximum absolute error on hidden frequencies, checks reciprocity by maximum matrix asymmetry, and checks passivity using the maximum singular value. It also reports the equal-weight mean of the three case RMSE values; every per-case limit and the mean limit must be met. The compact representation may contain at most 14 representative non-real poles (one from each conjugate pair) and at most two real poles. The verifier evaluates your model independently from its poles, residues, and direct terms. It does not require the exact generating poles or a particular fitting algorithm. All input samples are synthetic passive-system responses.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.

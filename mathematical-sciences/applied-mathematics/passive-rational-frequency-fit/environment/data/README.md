# Frequency-response inputs

The task contains three synthetic reciprocal three-port scattering responses. The data are generated from stable second-order all-pass factors mixed by a fixed real orthogonal matrix and scaled by gains below one. The construction produces a passive system with resonances over six frequency decades. Independent complex Gaussian noise with standard deviation `2e-5` is added to the agent-visible training samples.

This is a numerical rational-approximation and reduced-order modelling workflow, based on vector fitting for frequency-domain responses. The source method is Gustavsen and Semlyen, *IEEE Transactions on Power Delivery* 14(3), 1052–1061 (1999), DOI `10.1109/61.772353`.

## Input file

`training_response.csv` contains one row per frequency sample, with columns:

- `case_id`: `net1`, `net2`, or `net3`.
- `frequency_hz`: positive frequency in hertz.
- For each row-major channel `S11`, `S12`, `S13`, `S21`, `S22`, `S23`, `S31`, `S32`, and `S33`, columns `<channel>_re` and `<channel>_im` give the real and imaginary parts of the dimensionless scattering parameter.

All cases use a 50-ohm reference impedance and are reciprocal. The required model form is

`H_c(s) = d_c + e_c*s + sum_k r_ck/(s-p_k)`

with the conjugate term added for each non-real pole/residue pair so the overall model has real coefficients. For a real pole, include its term once. `s = 2*pi*i*f`. A valid stable model has `Re(p_k) < 0`; a passive scattering model has largest singular value of `H(i*2*pi*f)` no greater than one.

The three cases have different effective model orders. Do not assume a common order; select a compact model that predicts the independent verification frequencies while preserving stability, reciprocity, and passivity.

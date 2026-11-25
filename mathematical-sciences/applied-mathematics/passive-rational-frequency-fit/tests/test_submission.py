from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


HIDDEN = Path("/tests/hidden/heldout_response.csv")
THRESHOLDS = Path("/tests/hidden/score_thresholds.json")
SUBMISSION = Path("/app/submission/models.json")
CHANNELS = [f"S{i+1}{j+1}" for i in range(3) for j in range(3)]
EXPECTED_CASES = {"net1", "net2", "net3"}


def evaluate_model(model: dict[str, object], freq_hz: np.ndarray) -> np.ndarray:
    poles = np.asarray([complex(*p) for p in model["poles"]], dtype=complex)
    residues = np.asarray(
        [[complex(*c) for c in row] for row in model["residues"]], dtype=complex
    )
    direct = np.asarray(model["constant"], dtype=float)
    proportional = np.asarray(model["proportional"], dtype=float)
    assert residues.shape == (9, len(poles)), "residues must be 9 channels by number of poles"
    assert direct.shape == (9,) and proportional.shape == (9,)
    assert np.isfinite(poles.real).all() and np.isfinite(poles.imag).all()
    assert np.isfinite(residues.real).all() and np.isfinite(residues.imag).all()
    assert np.isfinite(direct).all() and np.isfinite(proportional).all()
    assert np.all(poles.real < -1e-8), "all poles must be strictly stable"
    assert np.all(poles.imag >= -1e-9), "use at most one representative from each conjugate pair"
    complex_poles = int(np.count_nonzero(poles.imag > 1e-8))
    real_poles = int(np.count_nonzero(np.abs(poles.imag) <= 1e-8))
    assert complex_poles <= 14 and real_poles <= 2
    assert np.all(np.abs(residues[:, np.abs(poles.imag) <= 1e-8].imag) <= 1e-7), (
        "residues for real poles must be real"
    )

    s = 2j * np.pi * freq_hz
    result = direct[None, :] + s[:, None] * proportional[None, :]
    for k, pole in enumerate(poles):
        residue = residues[:, k]
        result += residue[None, :] / (s[:, None] - pole)
        if pole.imag > 1e-8:
            result += np.conj(residue)[None, :] / (s[:, None] - np.conj(pole))
    return result.reshape((len(freq_hz), 3, 3))


def test_rational_models_predict_and_remain_passive():
    assert SUBMISSION.is_file(), "missing /app/submission/models.json"
    models = json.loads(SUBMISSION.read_text(encoding="utf-8"))
    assert set(models) == EXPECTED_CASES
    held = pd.read_csv(HIDDEN)
    limits = json.loads(THRESHOLDS.read_text(encoding="utf-8"))
    case_rmse = {}

    for case_id, frame in held.groupby("case_id", sort=True):
        assert case_id in models
        frame = frame.sort_values("frequency_hz")
        freq = frame["frequency_hz"].to_numpy(dtype=float)
        observed = np.empty((len(frame), 3, 3), dtype=complex)
        for c, channel in enumerate(CHANNELS):
            i, j = divmod(c, 3)
            observed[:, i, j] = frame[f"{channel}_re"].to_numpy() + 1j * frame[f"{channel}_im"].to_numpy()
        predicted = evaluate_model(models[case_id], freq)
        error = predicted - observed
        rmse = float(np.sqrt(np.mean(np.abs(error) ** 2)))
        max_abs = float(np.max(np.abs(error)))
        symmetry_error = float(np.max(np.abs(predicted - np.swapaxes(predicted, 1, 2))))
        sigma_max = np.linalg.svd(predicted, compute_uv=False)[:, 0]
        max_singular = float(np.max(sigma_max))
        case_rmse[case_id] = rmse
        assert rmse <= limits["case_rmse_max"], f"{case_id} response RMSE {rmse:.3e} exceeds limit"
        assert max_abs <= limits["max_abs_error"], f"{case_id} maximum response error exceeds limit"
        assert symmetry_error <= limits["symmetry_error_max"], f"{case_id} model is not reciprocal"
        assert max_singular <= limits["max_singular_value"], f"{case_id} model violates the passivity bound"

    macro = float(np.mean(list(case_rmse.values())))
    assert macro <= limits["macro_rmse_max"], f"macro response RMSE {macro:.3e} exceeds limit"

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import skrf as rf
from skrf.vectorFitting import VectorFitting


CHANNELS = [f"S{i+1}{j+1}" for i in range(3) for j in range(3)]


def as_network(frame: pd.DataFrame) -> rf.Network:
    frame = frame.sort_values("frequency_hz")
    freq = frame["frequency_hz"].to_numpy(dtype=float)
    s = np.empty((len(frame), 3, 3), dtype=complex)
    for c, name in enumerate(CHANNELS):
        i, j = divmod(c, 3)
        s[:, i, j] = frame[f"{name}_re"].to_numpy() + 1j * frame[f"{name}_im"].to_numpy()
    network = rf.Network(
        frequency=rf.Frequency.from_f(freq, unit="hz"),
        s=s,
        z0=50.0,
    )
    return network


def predict_fit(vf: VectorFitting, freq_hz: np.ndarray) -> np.ndarray:
    result = np.empty((len(freq_hz), 3, 3), dtype=complex)
    for i in range(3):
        for j in range(3):
            result[:, i, j] = vf.get_model_response(i, j, freq_hz)
    return result


def fit_network(network: rf.Network, n_poles_cmplx: int) -> VectorFitting:
    vf = VectorFitting(network)
    vf.vector_fit(
        n_poles_real=0,
        n_poles_cmplx=n_poles_cmplx,
        init_pole_spacing="log",
        parameter_type="s",
        fit_constant=True,
        fit_proportional=False,
    )
    return vf


def fit_case(frame: pd.DataFrame) -> dict[str, object]:
    frame = frame.sort_values("frequency_hz").reset_index(drop=True)
    validation_mask = np.arange(len(frame)) % 5 == 0
    train_frame = frame.loc[~validation_mask]
    validation_frame = frame.loc[validation_mask]
    validation_f = validation_frame["frequency_hz"].to_numpy(dtype=float)
    validation_s = as_network(validation_frame).s
    accepted: list[tuple[int, float]] = []
    for order in (4, 6, 8, 9, 10, 12, 14):
        candidate = fit_network(as_network(train_frame), order)
        prediction = predict_fit(candidate, validation_f)
        rmse = float(np.sqrt(np.mean(np.abs(prediction - validation_s) ** 2)))
        accepted.append((order, rmse))
        if rmse <= 3.5e-5:
            break
    selected_order = accepted[-1][0]
    selected_cv_rmse = accepted[-1][1]
    vf = fit_network(as_network(frame), selected_order)
    return {
        "cross_validation_rmse": selected_cv_rmse,
        "poles": [[float(p.real), float(p.imag)] for p in vf.poles],
        "residues": [
            [[float(r.real), float(r.imag)] for r in row]
            for row in vf.residues
        ],
        "constant": [float(x) for x in vf.constant_coeff],
        "proportional": [float(x) for x in vf.proportional_coeff],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=Path("/app/data/training_response.csv"))
    parser.add_argument("--output", type=Path, default=Path("/app/submission/models.json"))
    args = parser.parse_args()
    samples = pd.read_csv(args.data)
    result = {
        str(case_id): fit_case(group)
        for case_id, group in samples.groupby("case_id", sort=True)
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(",", ":")), encoding="utf-8")
    print(f"fitted {len(result)} networks")


if __name__ == "__main__":
    main()

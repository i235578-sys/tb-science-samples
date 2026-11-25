from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


DATA = Path("/app/data")
HIDDEN = Path("/tests/hidden")
SUBMISSION = Path("/app/submission/subsidence_parameters.csv")


def test_heldout_subsidence_predictions():
    assert SUBMISSION.is_file(), "missing subsidence_parameters.csv"
    parameters = pd.read_csv(SUBMISSION)
    assert list(parameters.columns) == ["source_id", "elastic_m_per_m", "inelastic_m_per_m"]
    assert parameters["source_id"].is_unique
    assert np.isfinite(parameters[["elastic_m_per_m", "inelastic_m_per_m"]].to_numpy()).all()
    assert (parameters[["elastic_m_per_m", "inelastic_m_per_m"]] >= 0.0).all().all()
    assert (parameters[["elastic_m_per_m", "inelastic_m_per_m"]] <= 0.004).all().all()

    sources = pd.read_csv(DATA / "source_blocks.csv")
    source_ids = sources["source_id"].tolist()
    assert set(parameters["source_id"]) == set(source_ids)
    parameters = parameters.set_index("source_id").loc[source_ids]
    elastic = parameters["elastic_m_per_m"].to_numpy()
    inelastic = parameters["inelastic_m_per_m"].to_numpy()

    kernel_long = pd.read_csv(DATA / "mogi_kernel.csv")
    kernel = kernel_long.pivot(index="sensor_id", columns="source_id", values="vertical_gain")[source_ids]
    sensor_ids = kernel.index.tolist()
    sensor_index = {name: i for i, name in enumerate(sensor_ids)}
    kernel_values = kernel.to_numpy()

    head = pd.read_csv(DATA / "head_history.csv")
    down = head.pivot(index="month", columns="source_id", values="drawdown_m")[source_ids].to_numpy()
    peak = np.maximum.accumulate(down, axis=0)

    site_ref = pd.read_csv(DATA / "site_reference_offsets.csv").set_index("sensor_id")["offset_mm"].to_dict()
    epoch_ref = pd.read_csv(DATA / "epoch_reference_offsets.csv").set_index("month")["offset_mm"].to_dict()
    heldout = pd.read_csv(HIDDEN / "heldout_observations.csv")
    predicted = np.empty(len(heldout), dtype=float)
    for k, row in enumerate(heldout.itertuples(index=False)):
        si = sensor_index[row.sensor_id]
        t = int(row.month)
        compaction = elastic * down[t] + inelastic * peak[t]
        predicted[k] = -1000.0 * float(kernel_values[si] @ compaction) + site_ref[row.sensor_id] + epoch_ref[t]
    heldout = heldout.copy()
    heldout["prediction"] = predicted
    heldout["squared_error"] = (heldout["prediction"] - heldout["displacement_mm"]) ** 2
    scores = heldout.groupby("split")["squared_error"].mean().pow(0.5).to_dict()

    limits = json.loads((HIDDEN / "score_thresholds.json").read_text(encoding="utf-8"))
    for split in ("spatial", "temporal"):
        assert split in scores
        assert scores[split] <= limits[f"{split}_rmse_max_mm"], (
            f"{split} held-out RMSE {scores[split]:.4f} mm exceeds its limit"
        )
    macro = float(np.mean([scores["spatial"], scores["temporal"]]))
    assert macro <= limits["macro_rmse_max_mm"], f"macro held-out RMSE {macro:.4f} mm exceeds its limit"

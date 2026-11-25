from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear


def build_arrays(data_dir: Path):
    sources = pd.read_csv(data_dir / "source_blocks.csv")
    sensors = pd.read_csv(data_dir / "sensor_locations.csv")
    source_ids = sources["source_id"].tolist()
    sensor_ids = sensors["sensor_id"].tolist()
    s_index = {key: i for i, key in enumerate(source_ids)}
    p_index = {key: i for i, key in enumerate(sensor_ids)}
    n_source = len(source_ids)

    kernel_long = pd.read_csv(data_dir / "mogi_kernel.csv")
    kernel = np.zeros((len(sensor_ids), n_source), dtype=float)
    for row in kernel_long.itertuples(index=False):
        kernel[p_index[row.sensor_id], s_index[row.source_id]] = row.vertical_gain

    head = pd.read_csv(data_dir / "head_history.csv")
    head_pivot = head.pivot(index="month", columns="source_id", values="drawdown_m")[source_ids]
    down = head_pivot.to_numpy()
    peak = np.maximum.accumulate(down, axis=0)
    return sources, source_ids, sensor_ids, p_index, kernel, down, peak


def smoothness_matrix(sources: pd.DataFrame) -> np.ndarray:
    xy = sources[["x_km", "y_km"]].to_numpy()
    edges: list[tuple[int, int]] = []
    for i in range(len(xy)):
        for j in range(i + 1, len(xy)):
            if abs(np.linalg.norm(xy[i] - xy[j]) - 6.0) < 0.15:
                edges.append((i, j))
    n = len(xy)
    result = np.zeros((2 * len(edges), 2 * n), dtype=float)
    for k, (i, j) in enumerate(edges):
        result[k, 2 * i] = 1.0
        result[k, 2 * j] = -1.0
        result[len(edges) + k, 2 * i + 1] = 1.0
        result[len(edges) + k, 2 * j + 1] = -1.0
    return result


def design_rows(obs: pd.DataFrame, p_index, kernel, down, peak, site_offset, epoch_offset):
    n_source = kernel.shape[1]
    a = np.empty((len(obs), 2 * n_source), dtype=float)
    b = np.empty(len(obs), dtype=float)
    sensor_indices = np.fromiter((p_index[k] for k in obs["sensor_id"]), dtype=int, count=len(obs))
    months = obs["month"].to_numpy(dtype=int)
    for r, (si, month) in enumerate(zip(sensor_indices, months)):
        k = kernel[si]
        a[r, 0::2] = -k * down[month]
        a[r, 1::2] = -k * peak[month]
    b[:] = obs["displacement_mm"].to_numpy(dtype=float)
    b -= np.array([site_offset[k] for k in obs["sensor_id"]])
    b -= epoch_offset[months]
    return a, b, sensor_indices


def fit_bounded(a, b, penalty, strength):
    augmented_a = np.vstack((a, np.sqrt(strength) * penalty))
    augmented_b = np.concatenate((b, np.zeros(penalty.shape[0])))
    return lsq_linear(
        augmented_a,
        augmented_b,
        bounds=(0.0, 4.0),
        method="trf",
        tol=1e-9,
        lsmr_tol="auto",
        max_iter=300,
    ).x


def solve(data_dir: Path, output_path: Path) -> None:
    sources, source_ids, sensor_ids, p_index, kernel, down, peak = build_arrays(data_dir)
    train = pd.read_csv(data_dir / "training_observations.csv")
    site_offset_df = pd.read_csv(data_dir / "site_reference_offsets.csv")
    site_offset = dict(zip(site_offset_df["sensor_id"], site_offset_df["offset_mm"]))
    epoch_offset_df = pd.read_csv(data_dir / "epoch_reference_offsets.csv")
    epoch_offset = np.zeros(int(epoch_offset_df["month"].max()) + 1)
    epoch_offset[epoch_offset_df["month"].to_numpy(dtype=int)] = epoch_offset_df["offset_mm"].to_numpy(dtype=float)
    a, b, sensor_indices = design_rows(train, p_index, kernel, down, peak, site_offset, epoch_offset)
    penalty = smoothness_matrix(sources)

    # Select spatial regularization with sensor-blocked cross-validation.
    strengths = [0.0, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]
    sensor_fold = {sensor_id: i % 4 for i, sensor_id in enumerate(sensor_ids)}
    folds = np.array([sensor_fold[k] for k in train["sensor_id"]], dtype=int)
    cv_rmse = []
    for strength in strengths:
        sum_sq = 0.0
        count = 0
        for fold in range(4):
            keep = folds != fold
            hold = ~keep
            q = fit_bounded(a[keep], b[keep], penalty, strength)
            residual = a[hold] @ q - b[hold]
            sum_sq += float(residual @ residual)
            count += int(hold.sum())
        cv_rmse.append(np.sqrt(sum_sq / count))
    best_strength = strengths[int(np.argmin(cv_rmse))]
    q = fit_bounded(a, b, penalty, best_strength)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "source_id": source_ids,
            "elastic_m_per_m": q[0::2] / 1000.0,
            "inelastic_m_per_m": q[1::2] / 1000.0,
        }
    ).to_csv(output_path, index=False)
    print(f"selected_spatial_penalty={best_strength:g}; blocked_cv_rmse_mm={min(cv_rmse):.4f}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=Path("/app/data"))
    parser.add_argument("--output", type=Path, default=Path("/app/submission/subsidence_parameters.csv"))
    args = parser.parse_args()
    solve(args.data, args.output)


if __name__ == "__main__":
    main()

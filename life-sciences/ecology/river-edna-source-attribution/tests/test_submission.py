"""Outcome checks for the submitted ecological production maps."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


APP = Path("/app")
TESTS = Path(__file__).resolve().parent
GENERA = ("Barbus", "Cottus", "Phoxinus")
SEASONS = ("spring", "summer", "autumn")
EXPECTED_PRODUCTION_COLUMNS = ["genus", "season", "node_id", "production_rate"]
EXPECTED_DECAY_COLUMNS = ["genus", "season", "tau_hours"]


def load_outputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    production_path = APP / "submission/production.csv"
    decay_path = APP / "submission/decay.csv"
    assert production_path.is_file(), "missing /app/submission/production.csv"
    assert decay_path.is_file(), "missing /app/submission/decay.csv"
    production = pd.read_csv(production_path)
    decay = pd.read_csv(decay_path)
    assert list(production.columns) == EXPECTED_PRODUCTION_COLUMNS
    assert list(decay.columns) == EXPECTED_DECAY_COLUMNS
    return production, decay


def validate_outputs(
    production: pd.DataFrame, decay: pd.DataFrame, n_nodes: int
) -> dict[tuple[str, str], tuple[np.ndarray, float]]:
    expected_groups = {(genus, season) for genus in GENERA for season in SEASONS}
    assert set(zip(decay["genus"], decay["season"])) == expected_groups
    assert len(decay) == len(expected_groups), "decay.csv has duplicate or extra rows"
    tau_values = decay["tau_hours"].to_numpy(dtype=np.float64)
    assert np.all(np.isfinite(tau_values)) and np.all(tau_values > 0.0)

    fits: dict[tuple[str, str], tuple[np.ndarray, float]] = {}
    expected_nodes = np.arange(n_nodes)
    for genus, season in expected_groups:
        subset = production[
            (production["genus"] == genus) & (production["season"] == season)
        ].sort_values("node_id")
        assert len(subset) == n_nodes, f"{genus}/{season} must include every reach"
        nodes = subset["node_id"].to_numpy()
        assert np.issubdtype(nodes.dtype, np.number)
        assert np.all(np.isfinite(nodes)) and np.array_equal(nodes, expected_nodes)
        values = subset["production_rate"].to_numpy(dtype=np.float64)
        assert np.all(np.isfinite(values)) and np.all(values >= 0.0)
        tau = float(
            decay.loc[
                (decay["genus"] == genus) & (decay["season"] == season), "tau_hours"
            ].iloc[0]
        )
        fits[(genus, season)] = (values, tau)
    assert len(production) == len(expected_groups) * n_nodes
    return fits


def forward_concentration(
    production: np.ndarray,
    tau_hours: float,
    network: dict[str, np.ndarray],
    season_index: int,
) -> np.ndarray:
    """Independent recurrence over reaches, matching the stated mass balance."""
    downstream = network["downstream_id"].astype(int)
    length = network["length_m"].astype(np.float64)
    width = network["width_m"][season_index].astype(np.float64)
    discharge = network["discharge_m3_s"][season_index].astype(np.float64)
    velocity = network["velocity_m_s"][season_index].astype(np.float64)
    tau_seconds = 3600.0 * tau_hours
    concentration = (
        production
        * width
        * length
        / discharge
        * np.exp(-length / (tau_seconds * velocity))
    )
    for reach in network["upstream_to_downstream_order"].astype(int):
        target = int(downstream[reach])
        concentration[target] += (
            concentration[reach]
            * discharge[reach]
            / discharge[target]
            * np.exp(-length[target] / (tau_seconds * velocity[target]))
        )
    return concentration


def test_artifact_schema_and_bounds() -> None:
    with np.load(APP / "data/river_network.npz") as data:
        n_nodes = len(data["downstream_id"])
    production, decay = load_outputs()
    validate_outputs(production, decay, n_nodes)


def test_heldout_spatial_prediction() -> None:
    with np.load(APP / "data/river_network.npz") as archive:
        network = {key: archive[key] for key in archive.files}
    production, decay = load_outputs()
    fits = validate_outputs(production, decay, len(network["downstream_id"]))
    heldout = pd.read_csv(TESTS / "hidden/heldout_observations.csv")
    thresholds = json.loads((TESTS / "hidden/score_thresholds.json").read_text())
    site_nodes = pd.read_csv(APP / "data/sample_sites.csv").set_index("site_id")["node_id"]

    expected_groups = {(genus, season) for genus in GENERA for season in SEASONS}
    assert set(zip(heldout["genus"], heldout["season"])) == expected_groups
    assert not heldout.duplicated(["genus", "season", "site_id"]).any()

    log_rmse: dict[str, float] = {}
    season_index = {season: i for i, season in enumerate(SEASONS)}
    for (genus, season), (source, tau) in fits.items():
        concentration = forward_concentration(
            source, tau, network, season_index[season]
        )
        assert np.all(np.isfinite(concentration)) and np.all(concentration >= 0.0)
        subset = heldout[
            (heldout["genus"] == genus) & (heldout["season"] == season)
        ]
        nodes = subset.site_id.map(site_nodes).to_numpy(dtype=int)
        observed = subset["read_abundance"].to_numpy(dtype=np.float64)
        predicted = concentration[nodes]
        key = f"{genus}/{season}"
        log_rmse[key] = float(
            np.sqrt(np.mean((np.log1p(predicted) - np.log1p(observed)) ** 2))
        )

    for group, score in log_rmse.items():
        limit = float(thresholds["log1p_rmse_per_group"][group])
        assert score <= limit, (
            f"held-out log1p RMSE for {group} is {score:.4f}; required <= {limit:.4f}"
        )
    macro_rmse = float(np.mean(list(log_rmse.values())))
    rmse_limit = float(thresholds["log1p_rmse_macro"])
    assert macro_rmse <= rmse_limit, (
        f"macro-averaged held-out log1p RMSE is {macro_rmse:.4f}; "
        f"required <= {rmse_limit:.4f}"
    )

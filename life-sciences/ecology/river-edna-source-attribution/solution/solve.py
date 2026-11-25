"""Reference MAP fit for the Thur fish-eDNA source reconstruction task."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize


DATA = Path("/app/data")
OUT = Path("/app/submission")
GENERA = ("Barbus", "Cottus", "Phoxinus")
SEASONS = ("spring", "summer", "autumn")
TAU_MU = np.log(5.0)
TAU_SIGMA = np.sqrt(np.log(5.0) - np.log(4.0))
BETA_PRIOR_SD = 0.25


def load_data() -> tuple[dict[str, np.ndarray], pd.DataFrame, pd.DataFrame]:
    with np.load(DATA / "river_network.npz") as source:
        network = {key: source[key] for key in source.files}
    obs = pd.read_csv(DATA / "training_observations.csv")
    sites = pd.read_csv(DATA / "sample_sites.csv")
    return network, obs, sites


def path_time_and_scale(
    network: dict[str, np.ndarray], nodes: np.ndarray, season_index: int
) -> tuple[np.ndarray, np.ndarray]:
    """Return travel time and source coefficient for each sample/source pair."""
    downstream = network["downstream_id"].astype(int)
    length = network["length_m"].astype(float)
    velocity = network["velocity_m_s"][season_index].astype(float)
    width = network["width_m"][season_index].astype(float)
    discharge = network["discharge_m3_s"][season_index].astype(float)
    n_nodes = len(downstream)

    upstream: list[list[int]] = [[] for _ in range(n_nodes)]
    for source, target in enumerate(downstream):
        if target >= 0:
            upstream[int(target)].append(source)

    travel_time = np.full((len(nodes), n_nodes), np.inf, dtype=np.float64)
    for row, target in enumerate(nodes.astype(int)):
        stack = [(int(target), 0.0)]
        while stack:
            node, time_after_node = stack.pop()
            total_time = time_after_node + length[node] / velocity[node]
            travel_time[row, node] = total_time
            stack.extend((source, total_time) for source in upstream[node])

    scale = np.broadcast_to(width * length, travel_time.shape) / discharge[nodes.astype(int), None]
    return travel_time, scale


def fit_one(
    x_all: np.ndarray,
    obs_nodes: np.ndarray,
    y: np.ndarray,
    network: dict[str, np.ndarray],
    season_index: int,
) -> tuple[np.ndarray, float, dict[str, float]]:
    """Fit tau, p0, and covariate effects against the visible observations."""
    x = x_all.astype(np.float64, copy=False)
    travel_time, source_scale = path_time_and_scale(network, obs_nodes, season_index)
    finite = np.isfinite(travel_time)
    travel_safe = np.where(finite, travel_time, 0.0)
    n_features = x.shape[1]

    def objective(theta: np.ndarray) -> tuple[float, np.ndarray]:
        log_tau = float(theta[0])
        tau_hours = float(np.exp(log_tau))
        p0 = float(theta[1])
        beta = theta[2:]

        linear = x @ beta
        linear_clipped = np.clip(linear, -35.0, 35.0)
        cov_mask = (linear > -35.0) & (linear < 35.0)
        production = p0 * np.exp(linear_clipped)

        exponent = np.exp(-travel_safe / (3600.0 * tau_hours))
        exponent[~finite] = 0.0
        kernel = source_scale * exponent
        mu = np.maximum(kernel @ production, 1e-200)
        log_den = np.log1p(mu)
        log_likelihood = np.sum(-log_den + y * (np.log(mu) - log_den))

        # Derivative of the geometric log likelihood with respect to its mean.
        d_ll_d_mu = (y - mu) / (mu * (1.0 + mu))
        d_mu_d_log_tau = (
            source_scale
            * exponent
            * (travel_safe / (3600.0 * tau_hours))
        ) @ production
        d_mu_d_log_tau[~np.isfinite(d_mu_d_log_tau)] = 0.0
        d_ll_d_log_tau = float(d_ll_d_mu @ d_mu_d_log_tau)

        d_mu_d_p0 = mu / p0
        d_ll_d_p0 = float(d_ll_d_mu @ d_mu_d_p0)
        d_mu_d_beta = (kernel * (production * cov_mask)[None, :]) @ x
        d_ll_d_beta = d_mu_d_beta.T @ d_ll_d_mu

        log_tau_prior = -np.log(tau_hours * TAU_SIGMA * np.sqrt(2.0 * np.pi))
        log_tau_prior -= (log_tau - TAU_MU) ** 2 / (2.0 * TAU_SIGMA**2)
        beta_prior = -0.5 * np.sum((beta / BETA_PRIOR_SD) ** 2)
        log_posterior = log_likelihood + log_tau_prior + beta_prior

        d_prior_d_log_tau = -(log_tau - TAU_MU) / (TAU_SIGMA**2) - 1.0
        gradient = np.concatenate(
            [
                np.asarray([d_ll_d_log_tau + d_prior_d_log_tau, d_ll_d_p0]),
                d_ll_d_beta - beta / (BETA_PRIOR_SD**2),
            ]
        )
        if not np.isfinite(log_posterior) or not np.all(np.isfinite(gradient)):
            return 1e300, np.zeros(2 + n_features, dtype=np.float64)
        return -float(log_posterior), -gradient

    bounds = [(np.log(0.05), np.log(24.0)), (1e-8, 1.0)] + [(-10.0, 10.0)] * n_features
    starts = []
    for tau0 in (2.0, 5.0, 12.0):
        decay = np.exp(-np.where(finite, travel_safe, 0.0) / (3600.0 * tau0))
        decay[~finite] = 0.0
        expected_transfer = np.sum(source_scale * decay, axis=1)
        scale_guess = float(np.mean(y) / max(float(np.mean(expected_transfer)), 1e-12))
        p0_guess = float(np.clip(scale_guess, 1e-5, 0.9))
        starts.append(np.concatenate([[np.log(tau0), p0_guess], np.zeros(n_features)]))

    candidates = [
        minimize(
            objective,
            start,
            method="L-BFGS-B",
            jac=True,
            bounds=bounds,
            options={"maxiter": 1200, "ftol": 2e-10, "gtol": 1e-6, "maxls": 40},
        )
        for start in starts
    ]
    best = min(candidates, key=lambda result: float(result.fun))
    theta = np.asarray(best.x, dtype=np.float64)
    tau_hours = float(np.exp(theta[0]))
    p0 = float(theta[1])
    beta = theta[2:]
    production = p0 * np.exp(np.clip(x_all @ beta, -35.0, 35.0))
    diagnostics = {
        "objective": float(best.fun),
        "iterations": float(getattr(best, "nit", -1)),
        "converged": float(bool(best.success)),
    }
    return production, tau_hours, diagnostics


def main() -> None:
    network, observations, sites = load_data()
    site_nodes = dict(zip(sites["site_id"], sites["node_id"].astype(int)))
    x_all = network["covariates"].astype(np.float64)
    n_nodes = len(network["downstream_id"])
    production_rows = []
    decay_rows = []

    for genus in GENERA:
        for season_index, season in enumerate(SEASONS):
            subset = observations[
                (observations["genus"] == genus) & (observations["season"] == season)
            ]
            nodes = subset["site_id"].map(site_nodes).to_numpy(dtype=int)
            y = subset["read_abundance"].to_numpy(dtype=np.float64)
            production, tau_hours, diagnostics = fit_one(
                x_all, nodes, y, network, season_index
            )
            print(genus, season, "tau_hours", round(tau_hours, 5), diagnostics, flush=True)
            decay_rows.append(
                {"genus": genus, "season": season, "tau_hours": tau_hours}
            )
            production_rows.extend(
                {
                    "genus": genus,
                    "season": season,
                    "node_id": node,
                    "production_rate": float(production[node]),
                }
                for node in range(n_nodes)
            )

    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(production_rows).to_csv(OUT / "production.csv", index=False, float_format="%.12g")
    pd.DataFrame(decay_rows).to_csv(OUT / "decay.csv", index=False, float_format="%.12g")


if __name__ == "__main__":
    main()

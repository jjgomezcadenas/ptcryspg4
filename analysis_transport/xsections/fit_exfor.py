"""Fit EXFOR excitation functions and generate correlated curve replicas.

The calculation uses published EXFOR central values without renormalizing
campaigns. Points lacking a positive quoted cross-section uncertainty remain
visible in plots but do not enter the weighted fit.
"""

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import BSpline

from .channels import CHANNELS
from .config import load


@dataclass
class PSpline:
    threshold: float
    energy_max: float
    threshold_power: float
    knots: np.ndarray
    degree: int = 3

    @classmethod
    def create(cls, threshold, energy_max, threshold_power, internal_knots):
        lower = threshold + 1.0e-3
        inside = np.asarray([value for value in internal_knots
                             if lower < value < energy_max], dtype=float)
        knots = np.concatenate((
            np.repeat(lower, 4), inside, np.repeat(energy_max, 4)))
        return cls(threshold, energy_max, threshold_power, knots)

    @property
    def coefficient_count(self):
        return len(self.knots) - self.degree - 1

    def design(self, energy):
        values = np.asarray(energy, dtype=float)
        clipped = np.clip(values, self.threshold + 1.0e-3, self.energy_max)
        return BSpline.design_matrix(
            clipped, self.knots, self.degree, extrapolate=True).toarray()

    def transform(self, energy, sigma):
        factor = np.maximum(np.asarray(energy) - self.threshold, 1.0e-12)
        return np.log(np.asarray(sigma)) - self.threshold_power * np.log(factor)

    def predict(self, coefficients, energy):
        values = np.asarray(energy, dtype=float)
        factor = np.maximum(values - self.threshold, 1.0e-12)
        log_sigma = (self.design(values) @ coefficients
                     + self.threshold_power * np.log(factor))
        sigma = np.exp(np.clip(log_sigma, -30.0, 20.0))
        return np.where(values > self.threshold, sigma, 0.0)


def _penalty(size):
    return np.diff(np.eye(size), n=2, axis=0)


def inverse_covariance(frame, log_uncertainty, campaign_log_spread,
                       campaign_spreads=None):
    """Point errors plus one fully correlated offset per campaign.

    A campaign with a documented normalization uncertainty (an entry in
    campaign_spreads) uses that value as its offset scale; the others use
    the fitted global spread. Anchored campaigns therefore pin the level of
    the curve to the precision their source published.
    """
    size = len(frame)
    inverse = np.zeros((size, size))
    campaigns = frame["campaign_id"].to_numpy()
    for campaign in np.unique(campaigns):
        spread = (campaign_spreads or {}).get(campaign, campaign_log_spread)
        indices = np.flatnonzero(campaigns == campaign)
        variances = np.maximum(log_uncertainty[indices], 1.0e-6) ** 2
        diagonal_inverse = np.diag(1.0 / variances)
        if spread > 0:
            column = (1.0 / variances)[:, None]
            denominator = (1.0 / spread ** 2
                           + np.sum(1.0 / variances))
            block = diagonal_inverse - column @ column.T / denominator
        else:
            block = diagonal_inverse
        inverse[np.ix_(indices, indices)] = block
    return inverse


def fit_coefficients(model, frame, log_sigma, log_uncertainty, smoothing,
                     campaign_log_spread, campaign_spreads=None):
    design = model.design(frame.energy_MeV)
    weight = inverse_covariance(frame, log_uncertainty, campaign_log_spread,
                                campaign_spreads)
    differences = _penalty(model.coefficient_count)
    system = (design.T @ weight @ design
              + smoothing * differences.T @ differences
              + 1.0e-10 * np.eye(model.coefficient_count))
    return np.linalg.solve(system, design.T @ weight @ log_sigma)


def choose_smoothing(model, frame, log_sigma, log_uncertainty, candidates,
                     campaign_log_spread, campaign_spreads=None):
    """Select smoothing by prediction of complete held-out campaigns."""
    campaigns = frame["campaign_id"].unique()
    scores = []
    for smoothing in candidates:
        campaign_scores = []
        for campaign in campaigns:
            test = frame["campaign_id"].to_numpy() == campaign
            train = ~test
            if train.sum() <= model.coefficient_count // 2:
                continue
            coefficients = fit_coefficients(
                model, frame.loc[train].reset_index(drop=True),
                log_sigma[train], log_uncertainty[train], smoothing,
                campaign_log_spread, campaign_spreads)
            prediction = model.predict(
                coefficients, frame.loc[test, "energy_MeV"].to_numpy())
            predicted_log = model.transform(
                frame.loc[test, "energy_MeV"].to_numpy(), prediction)
            residual = predicted_log - log_sigma[test]
            spread = (campaign_spreads or {}).get(campaign, campaign_log_spread)
            variance = log_uncertainty[test] ** 2 + spread ** 2
            campaign_scores.append(np.mean(
                residual ** 2 / variance + np.log(variance)))
        scores.append(np.mean(campaign_scores) if campaign_scores else np.inf)
    return float(candidates[int(np.argmin(scores))]), scores


def choose_campaign_spread(model, frame, log_sigma, log_uncertainty,
                           smoothing_candidates, spread_candidates,
                           campaign_spreads=None):
    results = []
    for spread in spread_candidates:
        smoothing, scores = choose_smoothing(
            model, frame, log_sigma, log_uncertainty,
            smoothing_candidates, float(spread), campaign_spreads)
        results.append((float(np.min(scores)), float(spread), smoothing, scores))
    return min(results, key=lambda result: result[0])


def _load_catalog_points(repo, selected, energy_min, energy_max):
    frames = []
    for metadata in selected.to_dict("records"):
        points = pd.read_csv(repo / metadata["point_file"])
        points["dataset_id"] = metadata["dataset_id"]
        points["campaign_id"] = metadata["dataset_id"].split("_")[1]
        points["label"] = metadata["label"]
        frames.append(points)
    if not frames:
        return pd.DataFrame()
    frame = pd.concat(frames, ignore_index=True)
    frame = frame[
        frame.energy_MeV.between(energy_min, energy_max)
        & (frame.sigma_mb > 0)
    ].copy()
    frame["sigma_unc_mb"] = frame[
        ["sigma_unc_minus_mb", "sigma_unc_plus_mb"]].mean(axis=1)
    return frame


def load_exfor_points(repo, channel, energy_min, energy_max,
                      extra_dataset_ids=()):
    catalog = pd.read_csv(repo / "data/xsections/normalized/datasets.csv")
    channel_catalog = catalog[
        (catalog.library == "EXFOR")
        & (catalog.target == channel.target)
        & (catalog.residual == channel.residual)
    ]
    curation_path = repo / "data/xsections/curation.csv"
    if not curation_path.exists():
        raise ValueError("Missing data/xsections/curation.csv; run curate_exfor first")
    curation = pd.read_csv(curation_path)
    channel_curation = curation[
        curation.dataset_id.isin(channel_catalog.dataset_id)]
    missing = set(channel_catalog.dataset_id) - set(channel_curation.dataset_id)
    if missing:
        raise ValueError(
            f"EXFOR datasets lack curation decisions: {sorted(missing)}")
    invalid_states = set(channel_curation.state) - {"accepted", "excluded", "pending"}
    if invalid_states:
        raise ValueError(f"Invalid EXFOR curation states: {sorted(invalid_states)}")
    accepted_ids = set(
        channel_curation.loc[channel_curation.state == "accepted", "dataset_id"])
    extra_ids = set(extra_dataset_ids)
    unknown = extra_ids - set(channel_catalog.dataset_id)
    if unknown:
        raise ValueError(
            f"Sensitivity datasets do not belong to {channel.channel_id}: "
            f"{sorted(unknown)}")
    extra_curation = channel_curation[
        channel_curation.dataset_id.isin(extra_ids)]
    invalid_extra = set(extra_curation.loc[
        extra_curation.state != "pending", "dataset_id"])
    if invalid_extra:
        raise ValueError(
            f"Sensitivity additions must be pending: {sorted(invalid_extra)}")
    selected = channel_catalog[
        channel_catalog.dataset_id.isin(accepted_ids | extra_ids)]
    return _load_catalog_points(repo, selected, energy_min, energy_max), selected


def load_pending_exfor_points(repo, channel, energy_min, energy_max):
    """Load pending series for display without admitting them to the fit."""
    catalog = pd.read_csv(repo / "data/xsections/normalized/datasets.csv")
    curation = pd.read_csv(repo / "data/xsections/curation.csv")
    pending_ids = set(curation.loc[
        (curation.channel_id == channel.channel_id)
        & (curation.state == "pending"), "dataset_id"])
    selected = catalog[
        (catalog.library == "EXFOR")
        & catalog.dataset_id.isin(pending_ids)]
    points = _load_catalog_points(repo, selected, energy_min, energy_max)
    if points.empty:
        return points
    return points.merge(
        curation[["dataset_id", "reason_code"]], on="dataset_id", how="left")


def _apply_normalizations(log_values, frame, config, rng, skip=()):
    """Correlated normalization draws for shared groups; datasets whose
    campaign offset already uses the documented value are skipped."""
    campaign_uncertainties = config.get("campaign_normalization_fraction", {})
    shared_groups = config.get("shared_normalization_group", {})
    shifts = {}
    for dataset_id, fraction in campaign_uncertainties.items():
        if dataset_id in skip:
            continue
        shifts[dataset_id] = rng.normal(0.0, float(fraction))
    grouped_draws = {}
    for dataset_id, group in shared_groups.items():
        if dataset_id not in campaign_uncertainties:
            continue
        if group not in grouped_draws:
            grouped_draws[group] = rng.normal(
                0.0, float(campaign_uncertainties[dataset_id]))
        shifts[dataset_id] = grouped_draws[group]
    adjusted = log_values.copy()
    for dataset_id, shift in shifts.items():
        mask = frame.dataset_id.to_numpy() == dataset_id
        adjusted[mask] += math.log(max(1.0 + shift, 1.0e-6))
    return adjusted


def _write_csv(path, frame):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, lineterminator="\n")


def fit_channel(repo, channel, config, rng):
    threshold = float(config["threshold_MeV"][channel.channel_id])
    all_points, campaigns = load_exfor_points(
        repo, channel, config["energy_min_MeV"], config["energy_max_MeV"])
    usable = all_points[
        (all_points.energy_MeV > threshold)
        & all_points.sigma_unc_mb.notna()
        & (all_points.sigma_unc_mb > 0)
    ].copy().reset_index(drop=True)
    if len(usable) < 10:
        raise ValueError(f"Too few usable EXFOR points for {channel.channel_id}")

    model = PSpline.create(
        threshold, config["energy_max_MeV"], config["threshold_power"],
        config["internal_knots_MeV"])
    log_sigma = model.transform(usable.energy_MeV, usable.sigma_mb)
    relative_uncertainty = usable.sigma_unc_mb.to_numpy() / usable.sigma_mb.to_numpy()
    log_uncertainty = np.sqrt(np.log1p(relative_uncertainty ** 2))

    # Campaigns whose source documents a normalization uncertainty anchor the
    # curve level at that documented precision (config table, with citations).
    documented = config.get("campaign_normalization_fraction", {})
    anchors = {}
    anchored_datasets = set()
    for dataset_id, fraction in documented.items():
        rows = usable[usable.dataset_id == dataset_id]
        if rows.empty:
            continue
        campaign = rows.campaign_id.iloc[0]
        spread = float(np.log1p(float(fraction)))
        anchors[campaign] = min(anchors.get(campaign, spread), spread)
        anchored_datasets.add(dataset_id)
    if anchors:
        print(f"{channel.channel_id}: anchored campaigns "
              f"{sorted(anchors)} at documented normalization")

    # Optional two-segment fit: the campaign-offset model applied separately
    # below and above a configured boundary, blended over 2 MeV. Prevents a
    # resonance region with large inter-campaign spread from setting the
    # plateau level (docs/xsection_fit.tex).
    boundary = config.get("segment_boundary_MeV", {}).get(channel.channel_id)
    if boundary is not None:
        boundary = float(boundary)
        segment_masks = [
            usable.energy_MeV.to_numpy() <= boundary + 2.0,
            usable.energy_MeV.to_numpy() > boundary,
        ]
    else:
        segment_masks = [np.ones(len(usable), dtype=bool)]

    segment_selection = []
    for mask in segment_masks:
        segment = usable.loc[mask].reset_index(drop=True)
        _, seg_spread, seg_smoothing, cv_scores = choose_campaign_spread(
            model, segment, log_sigma[mask], log_uncertainty[mask],
            config["lambda_candidates"],
            config["campaign_log_spread_candidates"], anchors)
        segment_selection.append((mask, seg_spread, seg_smoothing))
    # The reported spread/smoothing are those of the (last) plateau segment.
    campaign_log_spread = segment_selection[-1][1]
    smoothing = segment_selection[-1][2]

    def blended_predict(coefficient_sets, energies):
        values = np.asarray(energies, dtype=float)
        if boundary is None:
            return model.predict(coefficient_sets[0], values)
        low = model.predict(coefficient_sets[0], values)
        high = model.predict(coefficient_sets[1], values)
        weight = np.clip((values - boundary) / 2.0, 0.0, 1.0)
        with np.errstate(divide="ignore"):
            log_low = np.where(low > 0, np.log(np.maximum(low, 1e-300)), -690.0)
            log_high = np.where(high > 0, np.log(np.maximum(high, 1e-300)), -690.0)
        blended = np.exp((1.0 - weight) * log_low + weight * log_high)
        return np.where((low <= 0) & (weight < 1.0), high * weight,
                        np.where((high <= 0) & (weight > 0.0),
                                 low * (1.0 - weight), blended))

    # Per-campaign draw scale: documented anchor if present, else the fitted
    # spread of the segment holding the campaign's median energy.
    draw_spread = {}
    for campaign in usable.campaign_id.unique():
        cmask = usable.campaign_id.to_numpy() == campaign
        median_energy = float(np.median(usable.energy_MeV.to_numpy()[cmask]))
        if boundary is not None and median_energy <= boundary:
            draw_spread[campaign] = segment_selection[0][1]
        else:
            draw_spread[campaign] = campaign_log_spread
    campaign_scales = {**draw_spread, **anchors}

    dense_energy = np.linspace(
        config["energy_min_MeV"], config["energy_max_MeV"],
        config["fit_grid_points"])
    table_energy = np.asarray(config["table_energies_MeV"], dtype=float)
    replica_count = int(config["n_replicas"])
    dense_replicas = np.empty((replica_count, len(dense_energy)))
    table_replicas = np.empty((replica_count, len(table_energy)))
    smoothing_values = np.empty(replica_count)

    # Energy errors enter only after a dataset has been explicitly classified
    # as reporting an energy uncertainty rather than an energy-bin width.
    energy_uncertainty_datasets = set(
        config.get("energy_uncertainty_dataset_ids", []))
    energy_uncertainty = usable[
        ["energy_unc_minus_MeV", "energy_unc_plus_MeV"]].mean(axis=1).fillna(0.0)
    energy_uncertainty = energy_uncertainty.where(
        usable.dataset_id.isin(energy_uncertainty_datasets), 0.0).to_numpy()

    # Reselect smoothing for each replica from the same fixed candidate set.
    for replica_id in range(replica_count):
        varied = usable.copy()
        varied_energy = (usable.energy_MeV.to_numpy()
                         + rng.normal(size=len(usable)) * energy_uncertainty)
        varied["energy_MeV"] = np.clip(
            varied_energy, threshold + 1.0e-3, config["energy_max_MeV"])
        varied_log_sigma = log_sigma + rng.normal(size=len(usable)) * log_uncertainty
        for campaign in usable.campaign_id.unique():
            mask = usable.campaign_id.to_numpy() == campaign
            varied_log_sigma[mask] += rng.normal(
                0.0, campaign_scales[campaign])
        varied_log_sigma = _apply_normalizations(
            varied_log_sigma, usable, config, rng, skip=anchored_datasets)
        # Keep the sampled physical cross section fixed when the sampled
        # energy changes, then re-express it in the threshold transform.
        varied_sigma = np.exp(
            varied_log_sigma + config["threshold_power"]
            * np.log(usable.energy_MeV.to_numpy() - threshold))
        varied_log_sigma = model.transform(varied.energy_MeV, varied_sigma)
        varied["sigma_mb"] = varied_sigma
        coefficient_sets = []
        replica_smoothing = None
        for mask, seg_spread, _ in segment_selection:
            segment = varied.loc[mask].reset_index(drop=True)
            seg_smoothing, _ = choose_smoothing(
                model, segment, varied_log_sigma[mask], log_uncertainty[mask],
                config["lambda_candidates"], seg_spread, anchors)
            coefficient_sets.append(fit_coefficients(
                model, segment, varied_log_sigma[mask], log_uncertainty[mask],
                seg_smoothing, seg_spread, anchors))
            replica_smoothing = seg_smoothing
        dense_replicas[replica_id] = blended_predict(coefficient_sets, dense_energy)
        table_replicas[replica_id] = blended_predict(coefficient_sets, table_energy)
        smoothing_values[replica_id] = replica_smoothing
        if (replica_id + 1) % 100 == 0:
            print(f"{channel.channel_id}: fitted {replica_id + 1}/{replica_count} replicas")

    dense_nominal = np.median(dense_replicas, axis=0)
    dense_lower = np.quantile(dense_replicas, 0.16, axis=0)
    dense_upper = np.quantile(dense_replicas, 0.84, axis=0)
    table_nominal = np.median(table_replicas, axis=0)
    table_lower = np.quantile(table_replicas, 0.16, axis=0)
    table_upper = np.quantile(table_replicas, 0.84, axis=0)

    nominal_at_data = np.interp(usable.energy_MeV, dense_energy, dense_nominal)
    replica_at_data = np.asarray([
        np.interp(usable.energy_MeV, dense_energy, curve)
        for curve in dense_replicas
    ])
    nominal_log_at_data = model.transform(
        usable.energy_MeV.to_numpy(), nominal_at_data)
    replica_log_at_data = np.asarray([
        model.transform(usable.energy_MeV.to_numpy(), values)
        for values in replica_at_data
    ])
    covariance_inverse = inverse_covariance(
        usable, log_uncertainty, campaign_log_spread, campaign_scales)
    replica_residual = replica_log_at_data - nominal_log_at_data[None, :]
    distances = np.einsum(
        "ri,ij,rj->r", replica_residual, covariance_inverse,
        replica_residual) / len(usable)
    quantiles = config["representative_distance_quantiles"]
    representative_indices = []
    for quantile in quantiles:
        target = np.quantile(distances, quantile)
        ordered = np.argsort(np.abs(distances - target))
        representative_indices.append(
            next(int(index) for index in ordered if int(index) not in representative_indices))

    output = repo / "data/xsections/fits"
    curve = pd.DataFrame({
        "energy_MeV": dense_energy,
        "sigma_nominal_mb": dense_nominal,
        "sigma_lower_16_mb": dense_lower,
        "sigma_upper_84_mb": dense_upper,
        "relative_half_width": np.divide(
            0.5 * (dense_upper - dense_lower), dense_nominal,
            out=np.zeros_like(dense_nominal), where=dense_nominal > 0),
    })
    _write_csv(output / f"{channel.channel_id}_curve.csv", curve)
    table = pd.DataFrame({
        "energy_MeV": table_energy,
        "sigma_nominal_mb": table_nominal,
        "sigma_lower_16_mb": table_lower,
        "sigma_upper_84_mb": table_upper,
    })
    _write_csv(output / f"{channel.channel_id}_table.csv", table)

    replica_columns = {"replica_id": np.arange(replica_count),
                       "distance_D": distances,
                       "smoothing_lambda": smoothing_values}
    for index, energy in enumerate(table_energy):
        replica_columns[f"sigma_{energy:g}_MeV_mb"] = table_replicas[:, index]
    _write_csv(output / f"{channel.channel_id}_replicas.csv",
               pd.DataFrame(replica_columns))

    representative_rows = []
    for rank, replica_index in enumerate(representative_indices):
        for energy, sigma in zip(dense_energy, dense_replicas[replica_index]):
            representative_rows.append({
                "representative_rank": rank,
                "distance_quantile": quantiles[rank],
                "replica_id": replica_index,
                "distance_D": distances[replica_index],
                "energy_MeV": energy,
                "sigma_mb": sigma,
            })
    _write_csv(output / f"{channel.channel_id}_representatives.csv",
               pd.DataFrame(representative_rows))

    histogram_count, histogram_edges = np.histogram(
        distances, bins=int(config["histogram_bins"]))
    histogram = pd.DataFrame({
        "bin_left": histogram_edges[:-1],
        "bin_right": histogram_edges[1:],
        "count": histogram_count,
    })
    _write_csv(output / f"{channel.channel_id}_distance_histogram.csv", histogram)

    residual_log = log_sigma - nominal_log_at_data
    chi2 = float(residual_log @ covariance_inverse @ residual_log)
    design = model.design(usable.energy_MeV)
    weight = covariance_inverse
    difference = _penalty(model.coefficient_count)
    inverse = np.linalg.inv(
        design.T @ weight @ design
        + smoothing * difference.T @ difference
        + 1.0e-10 * np.eye(model.coefficient_count))
    effective_parameters = np.trace(
        design @ inverse @ design.T @ weight)
    dof = max(len(usable) - effective_parameters, 1.0)
    peak_index = int(np.argmax(dense_nominal))
    production_region = dense_nominal >= 0.05 * dense_nominal[peak_index]
    relative_half_width = np.divide(
        0.5 * (dense_upper - dense_lower), dense_nominal,
        out=np.zeros_like(dense_nominal), where=dense_nominal > 0)
    return {
        "channel_id": channel.channel_id,
        "title": channel.title,
        "threshold_MeV": threshold,
        "campaigns_available": len(campaigns),
        "independent_campaigns_used": usable.campaign_id.nunique(),
        "points_in_range": len(all_points),
        "points_used": len(usable),
        "points_with_energy_uncertainty_varied": int(
            np.count_nonzero(energy_uncertainty)),
        "points_without_positive_uncertainty": int(
            ((all_points.sigma_unc_mb.isna()) | (all_points.sigma_unc_mb <= 0)).sum()),
        "subthreshold_points_excluded": int((all_points.energy_MeV <= threshold).sum()),
        "nominal_smoothing_lambda": smoothing,
        "campaign_log_spread": campaign_log_spread,
        "campaign_fractional_spread": math.exp(segment_selection[0][1]) - 1.0,
        "chi2": chi2,
        "effective_dof": dof,
        "chi2_per_dof": chi2 / dof,
        "peak_energy_MeV": dense_energy[peak_index],
        "peak_sigma_mb": dense_nominal[peak_index],
        "peak_lower_16_mb": dense_lower[peak_index],
        "peak_upper_84_mb": dense_upper[peak_index],
        "median_relative_half_width_production": float(
            np.median(relative_half_width[production_region])),
        "p90_relative_half_width_production": float(
            np.quantile(relative_half_width[production_region], 0.90)),
        "median_replica_distance_D": float(np.median(distances)),
        "cv_scores": json.dumps(dict(zip(config["lambda_candidates"], cv_scores))),
        "parent_dataset_ids": ";".join(campaigns.dataset_id),
    }


def fit_sensitivity_curves(repo, config):
    """Fit named pending datasets together with the nominal accepted set."""
    requested = set(config.get("sensitivity_include_dataset_ids", []))
    if not requested:
        return pd.DataFrame()
    catalog = pd.read_csv(repo / "data/xsections/normalized/datasets.csv")
    unknown = requested - set(catalog.dataset_id)
    if unknown:
        raise ValueError(f"Unknown sensitivity datasets: {sorted(unknown)}")
    rows = []
    dense_energy = np.linspace(
        config["energy_min_MeV"], config["energy_max_MeV"],
        config["fit_grid_points"])
    for channel in CHANNELS:
        channel_ids = requested & set(catalog.loc[
            (catalog.target == channel.target)
            & (catalog.residual == channel.residual), "dataset_id"])
        if not channel_ids:
            continue
        threshold = float(config["threshold_MeV"][channel.channel_id])
        points, selected = load_exfor_points(
            repo, channel, config["energy_min_MeV"], config["energy_max_MeV"],
            extra_dataset_ids=channel_ids)
        usable = points[
            (points.energy_MeV > threshold)
            & points.sigma_unc_mb.notna()
            & (points.sigma_unc_mb > 0)
        ].copy().reset_index(drop=True)
        model = PSpline.create(
            threshold, config["energy_max_MeV"], config["threshold_power"],
            config["internal_knots_MeV"])
        log_sigma = model.transform(usable.energy_MeV, usable.sigma_mb)
        relative_uncertainty = (
            usable.sigma_unc_mb.to_numpy() / usable.sigma_mb.to_numpy())
        log_uncertainty = np.sqrt(np.log1p(relative_uncertainty ** 2))
        _, spread, smoothing, _ = choose_campaign_spread(
            model, usable, log_sigma, log_uncertainty,
            config["lambda_candidates"],
            config["campaign_log_spread_candidates"])
        coefficients = fit_coefficients(
            model, usable, log_sigma, log_uncertainty, smoothing, spread)
        sensitivity = model.predict(coefficients, dense_energy)
        curve_path = (repo / "data/xsections/fits"
                      / f"{channel.channel_id}_sensitivity.csv")
        _write_csv(curve_path, pd.DataFrame({
            "energy_MeV": dense_energy,
            "sigma_sensitivity_mb": sensitivity,
        }))
        nominal = pd.read_csv(
            repo / "data/xsections/fits" / f"{channel.channel_id}_curve.csv")
        production = nominal.sigma_nominal_mb >= 0.05 * nominal.sigma_nominal_mb.max()
        fractional = np.abs(
            sensitivity[production] / nominal.loc[production, "sigma_nominal_mb"] - 1.0)
        peak = int(np.argmax(sensitivity))
        rows.append({
            "channel_id": channel.channel_id,
            "included_pending_dataset_ids": ";".join(sorted(channel_ids)),
            "points_used": len(usable),
            "campaign_fractional_spread": math.exp(spread) - 1.0,
            "smoothing_lambda": smoothing,
            "peak_energy_MeV": dense_energy[peak],
            "peak_sigma_mb": sensitivity[peak],
            "median_absolute_fractional_difference": float(np.median(fractional)),
            "maximum_absolute_fractional_difference": float(np.max(fractional)),
            "parent_dataset_ids": ";".join(selected.dataset_id),
        })
    result = pd.DataFrame(rows)
    _write_csv(repo / "data/xsections/fits/sensitivity_summary.csv", result)
    return result


def update_model_catalog(repo, summaries):
    path = repo / "data/xsections/models/models.csv"
    fields = ("model_id", "channel_id", "curve_file", "parent_dataset_ids",
              "construction", "energy_min_MeV", "energy_max_MeV",
              "interpolation", "extrapolation", "version")
    existing = []
    if path.exists():
        with path.open(newline="", encoding="utf-8") as stream:
            existing = [row for row in csv.DictReader(stream)
                        if not row["model_id"].startswith("exfor_fit_")]
    for summary in summaries:
        existing.append({
            "model_id": f"exfor_fit_{summary['channel_id']}",
            "channel_id": summary["channel_id"],
            "curve_file": f"data/xsections/fits/{summary['channel_id']}_curve.csv",
            "parent_dataset_ids": summary["parent_dataset_ids"],
            "construction": "threshold-aware penalized cubic spline; EXFOR replicas",
            "energy_min_MeV": summary.get("energy_min_MeV", 5.0),
            "energy_max_MeV": summary.get("energy_max_MeV", 150.0),
            "interpolation": "linear on stored dense grid",
            "extrapolation": "none",
            "version": "1",
        })
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(existing)


def run(repo, config_path):
    config = load(config_path)
    output = repo / "data/xsections/fits"
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(int(config["random_seed"]))
    summaries = [fit_channel(repo, channel, config, rng) for channel in CHANNELS]
    summary_frame = pd.DataFrame(summaries)
    summary_frame["energy_min_MeV"] = config["energy_min_MeV"]
    summary_frame["energy_max_MeV"] = config["energy_max_MeV"]
    _write_csv(output / "fit_summary.csv", summary_frame)
    fit_sensitivity_curves(repo, config)
    metadata = {key: value for key, value in config.items()
                if key not in ("config_path",)}
    metadata["config_sha256"] = hashlib.sha256(config_path.read_bytes()).hexdigest()
    curation_path = repo / "data/xsections/curation.csv"
    metadata["curation_sha256"] = hashlib.sha256(
        curation_path.read_bytes()).hexdigest()
    (output / "fit_meta.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    update_model_catalog(repo, summaries)
    return summary_frame


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--config", type=Path, default=Path("config/xsection_fit.toml"))
    args = parser.parse_args()
    summary = run(args.repo.resolve(), args.config.resolve())
    print(f"Wrote EXFOR fits and replicas for {len(summary)} channels")


if __name__ == "__main__":
    main()

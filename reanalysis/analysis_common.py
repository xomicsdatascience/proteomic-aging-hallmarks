"""Shared, side-effect-free operations for the frozen HepG2 reanalysis."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, polygamma

ROOT = Path(__file__).resolve().parents[1]
PRIMARY_RHO = 0.01
MIN_TREATED_OBS = 3
MIN_CONTROL_OBS = 3
PROGRAM_NAMES = {
    "PROG_DNA_REPAIR": "DNA repair",
    "PROG_TELOMERE_MAINT": "Telomere maintenance",
    "PROG_CHROMATIN_ORG": "Chromatin organization",
    "PROG_PROTEOSTASIS_UPR": "Unfolded protein response",
    "PROG_PROTEASOME": "Ubiquitination / proteasome",
    "PROG_MACROAUTOPHAGY": "Macroautophagy machinery",
    "PROG_MTORC1_SIGNALING": "mTORC1 response",
    "PROG_OXPHOS": "Oxidative phosphorylation",
    "PROG_SENESCENCE_CORE": "Senescence-associated machinery",
    "PROG_CELL_CELL_COMM": "Cell-cell communication",
    "PROG_INFLAMMATORY": "Inflammatory response",
    "PROG_INNATE_TLR": "Host innate TLR sensing",
    "PRIOR_DEMAGALHAES_UP": "Aging-up RNA prior",
    "PRIOR_DEMAGALHAES_DN": "Aging-down RNA prior",
    "NUIS_E2F_PROLIF": "E2F proliferation response",
    "NUIS_G2M_PROLIF": "G2M checkpoint response",
    "NUIS_APOPTOSIS": "Apoptosis-associated response",
}
PRIMARY_IDS = [k for k in PROGRAM_NAMES if not k.startswith("NUIS_")]
SENTINELS = {
    "Mitochondrial proteins": ["NDUFS1", "IDH2", "ACO2", "SDHA", "COX4I1", "ATP5F1A"],
    "Inflammatory annotation and signaling proteins": ["ICAM1", "LDLR", "NFKB1", "RELA", "STAT3", "ATP2B1"],
}
KINASE_PANEL = ["Binimetinib", "Zanubrutinib", "Almonertinib (hydrochloride)",
                "Erlotinib (Hydrochloride)", "Mobocertinib (succinate)", "Lapatinib (ditosylate)"]
HIGHLIGHTS = ["Binimetinib", "Baicalein", "Cyclothiazide", "Zanubrutinib",
              "Orelabrutinib", "Almonertinib (hydrochloride)", "Osimertinib",
              "Erlotinib (Hydrochloride)", "Mobocertinib (succinate)",
              "Lapatinib (ditosylate)", "Icotinib", "Ferulic acid", "alpha-Lipoic Acid",
              "Lipoic acid", "Nordihydroguaiaretic acid", "Lonafarnib",
              "Atorvastatin (hemicalcium salt)", "Simvastatin acid (ammonium)",
              "Desloratadine", "Amiodarone (hydrochloride)", "DFO"]


def canonical_compound(value):
    return {"\u00ce\u00b1-Lipoic Acid": "alpha-Lipoic Acid", "\u03b1-Lipoic Acid": "alpha-Lipoic Acid"}.get(str(value), str(value))


def short_drug(value):
    return (str(value).split(" (")[0].replace("alpha-Lipoic Acid", "alpha-Lipoic acid")
            .replace("\u00ce\u00b3-Oryzanol", "gamma-Oryzanol"))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_table(df, path):
    """Stable text and gzip outputs, including a fixed gzip timestamp."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        import io
        with path.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            with io.TextIOWrapper(gz, encoding="utf-8", newline="") as f:
                df.to_csv(f, sep="\t", index=False, float_format="%.9g", na_rep="NA", lineterminator="\n")
    else:
        df.to_csv(path, sep="\t", index=False, float_format="%.9g", na_rep="NA", lineterminator="\n")


def read_table(path):
    target = ROOT / path
    if not target.exists() and str(path) == 'reanalysis/A02/protein_effects.tsv.gz':
        directory = ROOT / 'reanalysis/A02/protein_effects'
        manifest = pd.read_csv(directory / 'manifest.tsv', sep='\t')
        return pd.concat([pd.read_csv(directory / name, sep='\t') for name in manifest.file_name.drop_duplicates()], ignore_index=True)
    return pd.read_csv(target, sep="\t", keep_default_na=True)


def write_protein_partitions(frame):
    """Store complete results in bounded partitions, with a per-contrast index."""
    directory = ROOT / 'reanalysis/A02/protein_effects'
    directory.mkdir(exist_ok=True)
    for pattern in ['contrast_*.tsv.gz', 'part_*.tsv.gz']:
        for old in directory.glob(pattern): old.unlink()
    records = []
    groups = list(frame.groupby('compound_id', sort=True))
    for start in range(0, len(groups), 4):
        batch = groups[start:start+4]
        name = f'part_{start // 4 + 1:03d}.tsv.gz'
        write_table(pd.concat([group for _, group in batch], ignore_index=True), directory / name)
        digest = sha256(directory / name)
        for compound, group in batch:
            records.append({'compound_id': compound, 'file_name': name, 'n_rows': len(group), 'sha256': digest})
    write_table(pd.DataFrame(records), directory / 'manifest.tsv')


def load_raw():
    with h5py.File(ROOT / "data/raw (1).h5ad", "r") as f:
        X = f["X"][:]
        features = f["var/protein"].asstr()[:]
        obs = pd.DataFrame({"sample_id": f["obs/_index"].asstr()[:]})
        for key in ["batch", "drug", "sample_type"]:
            node = f["obs"][key]
            obs[key] = node["categories"].asstr()[:][node["codes"][:]]
    obs["original_compound_id"] = obs["drug"]
    obs["drug"] = obs["drug"].map(canonical_compound)
    if len(set(features)) != len(features):
        raise ValueError("Protein feature identifiers must be unique")
    return X, obs, features


def normalized_data():
    X, obs, all_features = load_raw()
    mapping = read_table("reanalysis/A01/protein_mapping.tsv")
    retained = mapping["retained_50pct_filter"].to_numpy(bool)
    if not np.array_equal(mapping.feature_id.to_numpy(), all_features):
        raise ValueError("Raw feature order does not match the frozen mapping")
    X = X[:, retained]
    Y = np.full(X.shape, np.nan)
    np.log2(X, out=Y, where=X > 0)
    medians = np.nanmedian(Y, axis=1, keepdims=True)
    valid_samples = obs.drug.ne("Drug").to_numpy()
    Y = Y - medians + np.nanmedian(medians[valid_samples])
    return Y, obs, all_features[retained]


def load_gmt():
    return {p[0]: set(p[2:]) for p in (line.rstrip().split("\t") for line in (ROOT / "reanalysis/A03/programs.gmt").read_text().splitlines())}


def bh(p):
    """BH across the declared family, conservatively padding unavailable tests."""
    p = np.asarray(p, dtype=float)
    valid = np.isfinite(p)
    out = stats.false_discovery_control(np.where(valid, p, 1.0), method="bh")
    out[~valid] = np.nan
    return out


def make_design(samples, compounds):
    plates = sorted(samples.plate.unique())
    names = ["Intercept"] + list(compounds) + ["plate:" + p for p in plates[1:]]
    D = np.column_stack([np.ones(len(samples))] +
                        [(samples.compound_id == c).to_numpy(float) for c in compounds] +
                        [(samples.plate == p).to_numpy(float) for p in plates[1:]])
    if np.linalg.matrix_rank(D) != D.shape[1]:
        raise ValueError("The primary design is not full rank")
    return D, names


def fit_observed_model(Y, D, contrasts):
    """Fit each observed subset; do not impute or treat absent coefficients as zero."""
    n_genes, n_contrasts = Y.shape[1], contrasts.shape[1]
    effects = np.full((n_genes, n_contrasts), np.nan)
    unscaled = effects.copy()
    variances = np.full(n_genes, np.nan)
    dfs = np.zeros(n_genes)
    residuals = np.full(Y.shape, np.nan)
    complete = np.isfinite(Y).all(axis=0)
    inverse = np.linalg.inv(D.T @ D)
    coef = inverse @ D.T @ Y[:, complete]
    res = Y[:, complete] - D @ coef
    effects[complete] = (contrasts.T @ coef).T
    unscaled[complete] = np.sum(contrasts * (inverse @ contrasts), axis=0)
    residuals[:, complete] = res
    dfs[complete] = len(D) - D.shape[1]
    variances[complete] = np.sum(res**2, axis=0) / (len(D) - D.shape[1])
    for j in np.where(~complete)[0]:
        observed = np.isfinite(Y[:, j])
        Dj = D[observed]
        columns = np.any(Dj != 0, axis=0)
        Dj = Dj[:, columns]
        if len(Dj) <= Dj.shape[1] + 5 or np.linalg.matrix_rank(Dj) != Dj.shape[1]:
            continue
        inv = np.linalg.inv(Dj.T @ Dj)
        bj = inv @ (Dj.T @ Y[observed, j])
        rj = Y[observed, j] - Dj @ bj
        residuals[observed, j] = rj
        dfs[j] = len(Dj) - Dj.shape[1]
        variances[j] = rj @ rj / dfs[j]
        estimable = np.all((contrasts == 0) | columns[:, None], axis=0)
        Cj = contrasts[columns]
        effects[j, estimable] = (Cj.T @ bj)[estimable]
        unscaled[j, estimable] = np.sum(Cj * (inv @ Cj), axis=0)[estimable]
    return effects, unscaled, variances, dfs, residuals


def moderate_variances(variances, dfs, fit_mask=None):
    valid = np.isfinite(variances) & (variances > 0) & (dfs > 5)
    if fit_mask is not None:
        valid &= np.asarray(fit_mask, bool)
    e = np.log(variances[valid]) - digamma(dfs[valid] / 2) + np.log(dfs[valid] / 2)
    target = max(np.var(e, ddof=1) - np.mean(polygamma(1, dfs[valid] / 2)), 1e-4)
    lo, hi = 0.1, 200.0
    for _ in range(50):
        mid = (lo + hi) / 2
        if polygamma(1, mid / 2) > target:
            lo = mid
        else:
            hi = mid
    d0 = (lo + hi) / 2
    s0 = np.exp(e.mean() + digamma(d0 / 2) - np.log(d0 / 2))
    post = (dfs * variances + d0 * s0) / (dfs + d0)
    return post, float(d0), float(s0)


def competitive_test(statistic, membership, rho=PRIMARY_RHO):
    """Non-rank CAMERA-PR two-sample statistic for the actual finite universe."""
    values = np.asarray(statistic, float)
    inside = np.asarray(membership, bool) & np.isfinite(values)
    outside = (~np.asarray(membership, bool)) & np.isfinite(values)
    m, mo = int(inside.sum()), int(outside.sum())
    if m < 5 or mo < 2 or not np.isfinite(rho):
        return np.nan, np.nan, m, np.nan
    vif = max(1.0, 1.0 + (m - 1) * float(rho))
    pool = ((m - 1) * values[inside].var(ddof=1) + (mo - 1) * values[outside].var(ddof=1)) / (m + mo - 2)
    se = np.sqrt(pool * (vif / m + 1 / mo))
    if not np.isfinite(se) or se <= 0:
        return np.nan, np.nan, m, vif
    t = (values[inside].mean() - values[outside].mean()) / se
    return float(t), float(2 * stats.t.sf(abs(t), df=m + mo - 2)), m, vif


def residual_correlations(residuals, features, gene_sets):
    """Pairwise complete observations; failed estimates remain unavailable."""
    output = {}
    for pid, genes in gene_sets.items():
        idx = np.where(np.isin(features, list(genes)))[0]
        corr = pd.DataFrame(residuals[:, idx]).corr(min_periods=20).to_numpy()
        output[pid] = (idx, corr)
    return output


def empirical_rho(correlation, selected):
    sub = correlation[np.ix_(selected, selected)]
    pairs = sub[np.triu_indices(len(sub), 1)]
    valid = np.isfinite(pairs)
    fraction = float(valid.mean()) if len(pairs) else 0.0
    if len(pairs) == 0 or fraction < 0.8:
        return np.nan, fraction
    return max(0.0, float(pairs[valid].mean())), fraction


def primary_mask(samples):
    return (~samples.excluded) & samples.plate.ne("RedMix") & samples.condition_role.ne("analytical_qc")


def json_write(value, path):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def plate_program_comparisons(Y, samples, features, tested_features, comparisons, gene_sets):
    """Competitive mean log2FC per independent plate, followed by a one-sample t test.

    Proteins within a plate are not treated as independent biological replicates.
    Membership and observation thresholds are fixed from the full contrast. Within
    each plate, treated and reference means use the same available feature IDs.
    """
    rows, unit_rows = [], []
    for label, treated, reference in comparisons:
        allowed = np.asarray(tested_features[label], bool)
        common = sorted(set(samples.loc[samples.compound_id.eq(treated), 'plate']) &
                        set(samples.loc[samples.compound_id.eq(reference), 'plate']))
        fc_by_plate = []
        for plate in common:
            a = samples.compound_id.eq(treated) & samples.plate.eq(plate)
            b = samples.compound_id.eq(reference) & samples.plate.eq(plate)
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                fc = np.nanmean(Y[a], axis=0) - np.nanmean(Y[b], axis=0)
            fc_by_plate.append((plate, fc))
        for pid, members in gene_sets.items():
            membership = np.isin(features, list(members))
            scores, abundance = [], []
            for plate, fc in fc_by_plate:
                measured = allowed & np.isfinite(fc)
                inside, outside = measured & membership, measured & ~membership
                if inside.sum() < 5 or outside.sum() < 100:
                    continue
                raw_effect = float(fc[inside].mean())
                score = raw_effect - float(fc[outside].mean())
                scores.append(score); abundance.append(raw_effect)
                unit_rows.append({'compound_id': label, 'program_id': pid, 'plate': plate,
                                  'treated': treated, 'reference': reference,
                                  'n_program_features': int(inside.sum()),
                                  'n_background_features': int(outside.sum()),
                                  'program_mean_log2fc': raw_effect,
                                  'competitive_mean_log2fc': score})
            n = len(scores)
            effect = float(np.mean(scores)) if n else np.nan
            se = float(np.std(scores, ddof=1) / np.sqrt(n)) if n >= 3 else np.nan
            statistic = effect / se if np.isfinite(se) and se > 0 else np.nan
            pvalue = float(2 * stats.t.sf(abs(statistic), n - 1)) if np.isfinite(statistic) else np.nan
            margin = stats.t.ppf(.975, n - 1) * se if n >= 3 else np.nan
            rows.append({'compound_id': label, 'program_id': pid, 'n_effective_plates': n,
                         'plate_t': statistic, 'p_plate': pvalue,
                         'plate_competitive_log2fc': effect,
                         'plate_ci95_low': effect - margin, 'plate_ci95_high': effect + margin,
                         'plate_program_mean_log2fc': float(np.mean(abundance)) if n else np.nan,
                         'plate_test_status': 'estimable' if np.isfinite(pvalue) else 'insufficient_matched_plates'})
    result = pd.DataFrame(rows)
    for primary in [True, False]:
        use = result.program_id.isin(PRIMARY_IDS) == primary
        if use.any(): result.loc[use, 'q_BH_global_plate'] = bh(result.loc[use, 'p_plate'])
    return result, pd.DataFrame(unit_rows)

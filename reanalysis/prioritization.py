"""Exploratory response indices, with no score-level null or significance test."""
import numpy as np
import pandas as pd

# Original 14-program target hypotheses; these do not establish benefit.
COMBINED_DIRECTIONS = {
    'PROG_OXPHOS': 1, 'PRIOR_DEMAGALHAES_DN': 1,
    'PROG_MACROAUTOPHAGY': 1, 'PROG_PROTEOSTASIS_UPR': 1,
    'PROG_PROTEASOME': 1, 'PROG_DNA_REPAIR': 1,
    'PROG_TELOMERE_MAINT': 1, 'PROG_CHROMATIN_ORG': 1,
    'PROG_INFLAMMATORY': -1, 'PROG_INNATE_TLR': -1,
    'PROG_SENESCENCE_CORE': -1, 'PROG_MTORC1_SIGNALING': -1,
    'PRIOR_DEMAGALHAES_UP': -1, 'PROG_CELL_CELL_COMM': -1,
}
ASSAY_DIRECTIONS = {
    'PROG_OXPHOS': 1, 'PROG_MACROAUTOPHAGY': 1,
    'PROG_PROTEASOME': 1, 'PROG_SENESCENCE_CORE': -1,
}
GROWTH_PROXIES = ['NUIS_E2F_PROLIF', 'NUIS_G2M_PROLIF']
EFFECT = 'plate_competitive_log2fc'


def library_rank(values, eligible=None):
    """Descending library ranks; DFO never enters the ranking."""
    keep = values.index != 'DFO'
    if eligible is not None:
        keep &= eligible
    return values.where(keep).rank(ascending=False, method='min')


def indices(programs):
    """Average signed primary effects divided by library-only program SD.

    No mean is subtracted: zero remains no competitive abundance response.
    SDs equalize observed program scales, not sampling errors or null variances.
    Missing components propagate; the axis count is never silently reduced.
    """
    effect = programs.pivot(index='compound_id', columns='program_id', values=EFFECT)
    required = list(COMBINED_DIRECTIONS) + GROWTH_PROXIES
    effect = effect[required].replace([np.inf, -np.inf], np.nan)
    library = effect.drop(index='DFO', errors='ignore')
    scales = library.std(ddof=1).where(lambda x: x > 0)
    scaled = effect / scales
    result = pd.DataFrame(index=effect.index)
    for name, directions in [('combined_response_index', COMBINED_DIRECTIONS),
                             ('four_axis_index', ASSAY_DIRECTIONS)]:
        signed = scaled[list(directions)] * pd.Series(directions)
        result[name] = signed.mean(axis=1, skipna=False)
        result[name + '_rank'] = library_rank(result[name])
        result[name + '_n_requested_directions'] = (signed > 0).sum(axis=1)
        result[name + '_complete'] = signed.notna().all(axis=1)
        if name == 'four_axis_index':
            all_four = (signed > 0).all(axis=1) & signed.notna().all(axis=1)
            result['all_four_requested_directions'] = all_four
            result['four_axis_eligible_rank'] = library_rank(result[name], all_four)
            result['four_axis_weakest_signed_component'] = signed.min(axis=1, skipna=False)

    # Post-treatment E2F/G2M sensitivity, fitted on library entries only.
    X = np.column_stack([np.ones(len(scaled)), scaled[GROWTH_PROXIES].to_numpy()])
    for name in ['combined_response_index', 'four_axis_index']:
        y = result[name].to_numpy()
        fit = np.isfinite(y) & np.isfinite(X).all(axis=1) & (effect.index != 'DFO')
        residual = np.full(len(y), np.nan)
        if fit.sum() > X.shape[1] and np.linalg.matrix_rank(X[fit]) == X.shape[1]:
            beta = np.linalg.lstsq(X[fit], y[fit], rcond=None)[0]
            predict = np.isfinite(y) & np.isfinite(X).all(axis=1)
            residual[predict] = y[predict] - X[predict] @ beta
        result[name + '_growth_proxy_residual'] = residual
        result[name + '_growth_proxy_residual_rank'] = library_rank(pd.Series(residual, index=effect.index))

    definitions = pd.DataFrame({
        'target_direction': pd.Series(COMBINED_DIRECTIONS),
        'four_axis_direction': pd.Series(ASSAY_DIRECTIONS),
    }).rename_axis('program_id').reset_index()
    definitions['combined_weight'] = 1 / len(COMBINED_DIRECTIONS)
    definitions['four_axis_weight'] = definitions.four_axis_direction.notna().astype(float) / len(ASSAY_DIRECTIONS)
    scale_table = pd.DataFrame({'library_effect_sd': scales, 'n_library_effects': library.count()})
    definitions = definitions.merge(scale_table, on='program_id', validate='one_to_one')
    definitions['effect_field'] = EFFECT
    definitions['centering'] = 'none; zero is no competitive response'
    components = programs[programs.program_id.isin(COMBINED_DIRECTIONS)].merge(
        definitions, on='program_id', validate='many_to_one')
    components['signed_scaled_effect'] = components[EFFECT] / components.library_effect_sd * components.target_direction
    components['four_axis_component'] = components.program_id.isin(ASSAY_DIRECTIONS)
    columns = ['compound_id', 'program_id', EFFECT, 'plate_ci95_low', 'plate_ci95_high',
               'plate_t', 'q_BH_global_plate', 'n_effective_plates', 'plate_test_status',
               'nes_gsea', 'target_direction', 'signed_scaled_effect', 'four_axis_component']
    return (result.reset_index(), definitions,
            components[columns].sort_values(['compound_id', 'program_id']).reset_index(drop=True))

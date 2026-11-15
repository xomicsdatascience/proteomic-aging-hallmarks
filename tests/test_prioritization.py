"""Regression checks for exploratory rankings and their scientific boundaries."""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reanalysis'))
from prioritization import indices, COMBINED_DIRECTIONS, ASSAY_DIRECTIONS, GROWTH_PROXIES


class PrioritizationTests(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(42)
        names=[f'drug{i}' for i in range(12)]+['DFO']
        rows=[]
        for drug in names:
            for pid in list(COMBINED_DIRECTIONS)+GROWTH_PROXIES:
                effect=float(rng.normal())
                rows.append({'compound_id':drug,'program_id':pid,
                             'plate_competitive_log2fc':effect,'plate_ci95_low':effect-.2,
                             'plate_ci95_high':effect+.2,'plate_t':effect*2,
                             'q_BH_global_plate':.5,'n_effective_plates':6,
                             'plate_test_status':'estimable','nes_gsea':effect})
        self.programs=pd.DataFrame(rows)

    def test_assay_control_cannot_change_library_scales_scores_or_residuals(self):
        a,da,_=indices(self.programs)
        changed=self.programs.copy()
        changed.loc[changed.compound_id.eq('DFO'),'plate_competitive_log2fc']=1e9
        b,db,_=indices(changed)
        pd.testing.assert_frame_equal(da,db)
        pd.testing.assert_frame_equal(a[a.compound_id.ne('DFO')],b[b.compound_id.ne('DFO')])
        self.assertTrue(a.loc[a.compound_id.eq('DFO')].filter(regex='_rank$').isna().all(axis=None))

    def test_zero_effect_keeps_zero_anchor(self):
        self.programs.loc[self.programs.compound_id.eq('drug0'),'plate_competitive_log2fc']=0
        result,_,_=indices(self.programs)
        row=result.set_index('compound_id').loc['drug0']
        self.assertEqual(row.combined_response_index,0)
        self.assertEqual(row.four_axis_index,0)
        self.assertFalse(row.all_four_requested_directions)

    def test_missing_axis_does_not_silently_reduce_denominator(self):
        m=self.programs.compound_id.eq('drug0')&self.programs.program_id.eq('PROG_OXPHOS')
        self.programs.loc[m,'plate_competitive_log2fc']=np.nan
        result,_,_=indices(self.programs)
        row=result.set_index('compound_id').loc['drug0']
        self.assertTrue(np.isnan(row.combined_response_index))
        self.assertTrue(np.isnan(row.four_axis_index))
        self.assertTrue(np.isnan(row.four_axis_eligible_rank))
        self.assertFalse(row.all_four_requested_directions)

    def test_high_compensatory_mean_cannot_enter_four_direction_shortlist(self):
        for pid,direction in ASSAY_DIRECTIONS.items():
            m=self.programs.compound_id.eq('drug0')&self.programs.program_id.eq(pid)
            self.programs.loc[m,'plate_competitive_log2fc']=100*direction if pid!='PROG_PROTEASOME' else -.01
        result,_,_=indices(self.programs)
        row=result.set_index('compound_id').loc['drug0']
        self.assertGreater(row.four_axis_index,0)
        self.assertEqual(row.four_axis_index_n_requested_directions,3)
        self.assertFalse(row.all_four_requested_directions)
        self.assertTrue(np.isnan(row.four_axis_eligible_rank))

    def test_sign_eligibility_does_not_require_significance(self):
        for pid,direction in ASSAY_DIRECTIONS.items():
            m=self.programs.compound_id.eq('drug0')&self.programs.program_id.eq(pid)
            self.programs.loc[m,'plate_competitive_log2fc']=.01*direction
        result,_,_=indices(self.programs)
        row=result.set_index('compound_id').loc['drug0']
        self.assertTrue(row.all_four_requested_directions)
        self.assertTrue(np.isfinite(row.four_axis_eligible_rank))
        self.assertFalse(any(c.startswith(('p_','q_')) for c in result.columns))

    def test_program_unit_rescaling_does_not_change_indices(self):
        a,_,_=indices(self.programs)
        changed=self.programs.copy()
        changed.loc[changed.program_id.eq('PROG_OXPHOS'),'plate_competitive_log2fc']*=100
        b,_,_=indices(changed)
        np.testing.assert_allclose(a.combined_response_index,b.combined_response_index,atol=1e-12)
        np.testing.assert_allclose(a.four_axis_index,b.four_axis_index,atol=1e-12)


if __name__=='__main__':unittest.main()

import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reanalysis'))
from analysis_common import bh, competitive_test, canonical_compound, fit_observed_model, plate_program_comparisons, write_table


class AnalysisRegressionTests(unittest.TestCase):
    def test_missing_pvalues_do_not_erase_real_discoveries(self):
        adjusted=bh([1e-7,.02,np.nan,.8])
        self.assertLess(adjusted[0],.001)
        self.assertTrue(np.isnan(adjusted[2]))
        self.assertAlmostEqual(adjusted[1],.04)

    def test_competitive_test_uses_actual_available_background(self):
        values=np.array([2.,3.,4.,5.,6., -2.,-1.,0.,1.,2.,np.nan])
        members=np.array([True]*5+[False]*6)
        t,p,n,_=competitive_test(values,members,rho=0)
        expected=stats.ttest_ind(values[:5],values[5:10],equal_var=True)
        self.assertAlmostEqual(t,expected.statistic,places=12)
        self.assertAlmostEqual(p,expected.pvalue,places=12)
        self.assertEqual(n,5)

    def test_correlated_proteins_do_not_inflate_biological_replicates(self):
        changes=np.array([.1,.3,.2,.4,.5,.25])
        sample=pd.DataFrame({'compound_id':['drug']*6+['DMSO']*6,'plate':[f'P{i}' for i in range(6)]*2})
        results=[]
        for n_members in [5,50]:
            features=np.array([f'g{i}' for i in range(n_members+100)])
            Y=np.zeros((12,len(features)));Y[:6,:n_members]=changes[:,None]
            table,_=plate_program_comparisons(Y,sample,features,{'drug':np.ones(len(features),bool)},[('drug','drug','DMSO')],{'module':set(features[:n_members])})
            results.append(table.iloc[0])
        self.assertEqual(results[0].n_effective_plates,6)
        self.assertEqual(results[1].n_effective_plates,6)
        self.assertAlmostEqual(results[0].plate_t,results[1].plate_t,places=12)
        self.assertAlmostEqual(results[0].p_plate,stats.ttest_1samp(changes,0).pvalue,places=12)

    def test_unmatched_plates_are_not_treated_as_pairs(self):
        features=np.array([f'g{i}' for i in range(105)])
        sample=pd.DataFrame({'compound_id':['drug']*3+['ref']*3,'plate':['A','B','C','D','E','F']})
        table,_=plate_program_comparisons(np.zeros((6,105)),sample,features,{'comparison':np.ones(105,bool)},[('comparison','drug','ref')],{'module':set(features[:5])})
        self.assertEqual(table.iloc[0].n_effective_plates,0)
        self.assertTrue(np.isnan(table.iloc[0].p_plate))

    def test_missing_treatment_is_not_a_zero_effect(self):
        D=np.column_stack([np.ones(12),[0]*6+[1]*6])
        Y=np.column_stack([np.arange(12,dtype=float),np.r_[np.arange(6,dtype=float),[np.nan]*6]])
        B,_,_,_,_=fit_observed_model(Y,D,np.array([[0.],[1.]]))
        self.assertTrue(np.isfinite(B[0,0]))
        self.assertTrue(np.isnan(B[1,0]))

    def test_alias_resolution_is_explicit(self):
        self.assertEqual(canonical_compound('\u00ce\u00b1-Lipoic Acid'),'alpha-Lipoic Acid')
        self.assertEqual(canonical_compound('\u03b1-Lipoic Acid'),'alpha-Lipoic Acid')
        self.assertEqual(canonical_compound('Lipoic acid'),'Lipoic acid')

    def test_gzip_tables_have_deterministic_bytes(self):
        frame=pd.DataFrame({'compound':['drug'],'effect':[.12345]})
        with tempfile.TemporaryDirectory() as directory:
            a=Path(directory)/'a.tsv.gz';b=Path(directory)/'b.tsv.gz'
            write_table(frame,a);write_table(frame,b)
            self.assertEqual(a.read_bytes(),b.read_bytes())

if __name__=='__main__':unittest.main()

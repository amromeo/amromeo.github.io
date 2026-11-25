from pathlib import Path
import numpy as np,pandas as pd
from inputs import folders
scores=pd.read_csv('dataset-scores.csv.gz').set_index(['study','j','method'])
for folder in folders:
 r=Path(folder)
 fit=pd.read_csv(r/'evidence/fits.csv.gz',low_memory=False)
 check=pd.read_csv(r/'evidence/checks.csv.gz').query("panel=='guidance' and seed_kind=='primary'")
 for method,arm in [('B','t0.25'),('C','baseline')]:
  q=fit[fit.arm==arm].merge(check[check.arm==arm],on=['j','arm'],validate='one_to_one').set_index('j')
  one=(q.Analyte=='CRP') if 'Analyte' in q else pd.Series(False,index=q.index)
  keep=q.status.eq('keep') & np.isfinite(q.released_native_upper) & (one | np.isfinite(q.released_native_lower))
  z=q.native_z_upper.copy()
  if 'Analyte' in q:
   zero=q.released_native_upper.where(keep).eq(0).groupby(q.Analyte).transform('any')
   z=q.native_z_upper_before_group.where(~zero,q.native_z_upper_if_group_zero)
  error=z.abs().clip(upper=5).where(one,(q.native_z_lower.abs().clip(upper=5)+z.abs().clip(upper=5))/2)
  study='RIbench' if 'Analyte' in q else 'confirmation'
  d=scores.xs((study,method),level=('study','method')).loc[q.index]
  assert np.array_equal(d.returned,keep)
  np.testing.assert_allclose(d.error,error.where(keep),atol=1e-9,rtol=1e-10,equal_nan=True)
print('Standalone B and C reproduce the frozen native endpoint scores on every released dataset.')

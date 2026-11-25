from pathlib import Path
import numpy as np,pandas as pd
def read(pat):return pd.concat([pd.read_csv(p).assign(job=p.parent.name) for p in Path('art').glob(pat)]).set_index('j').sort_index()
a=read('*/audit-*.csv');p=read('*/primary-*.csv').loc[a.index]
assert len(p)==749
for c in ['routine_hash','anchor_hash','seed','h']:assert a[c].equals(p[c]),c
diff={}
for c in ['family','components','status','reason','lrt','lower_flag','optimizer']:
 diff[c]=~(a[c].eq(p[c])|(a[c].isna()&p[c].isna()))
for c in ['lower','upper','raw_lower','raw_upper']:
 diff[c]=~np.isclose(a[c],p[c],atol=1e-10,rtol=0,equal_nan=True)
bad=pd.DataFrame(diff,index=a.index)
details=p.add_prefix('primary_').join(a.add_prefix('audit_')).join(bad.add_suffix('_diff'))
details[bad.any(axis=1)].to_csv('first-audit-discrepancies.csv')
print('Original audit comparisons',len(a),'discrepant cases',bad.any(axis=1).sum());print(bad.sum())

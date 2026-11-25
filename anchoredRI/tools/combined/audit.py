from pathlib import Path
import pandas as pd,numpy as np
from inputs import folders
def read(pat):
 return pd.concat([pd.read_csv(p,keep_default_na=False) for p in sorted(Path('art').glob(pat))],ignore_index=True)
expected=pd.read_csv('audit.csv').j
for p in Path('art').glob('*/audit.csv'):
 assert pd.read_csv(p).j.equals(expected)
p=read('*/primary-*.csv').set_index('j').sort_index()
a=read('*/audit-*.csv').set_index('j').sort_index()
assert len(p)==7488 and p.index.is_unique and set(a.index)==set(expected) and len(a)==749
for c in p:
 x=p.loc[a.index,c];y=a[c]
 if c in ['lower','upper','raw_lower','raw_upper']:
  np.testing.assert_allclose(pd.to_numeric(x.replace('NA',np.nan)),pd.to_numeric(y.replace('NA',np.nan)),atol=1e-10,rtol=0,equal_nan=True)
 else:assert x.astype(str).eq(y.astype(str)).all(),c
keys=['j','arm']
chk=read('*/checks-*.csv').set_index(keys).sort_index()
official=pd.read_csv(Path(folders[0],'evidence/checks.csv.gz'),keep_default_na=False)
official=official[official.j.isin(expected)&official.arm.isin(['baseline','w1','t0.25'])&(official.panel=='guidance')&(official.seed_kind=='primary')].set_index(keys).sort_index()
assert len(chk)==2247 and chk.index.equals(official.index)
for c in ['status','reason','seed','released_native_upper','released_native_lower','statistic','p']:
 x=chk[c].replace({'NA':'','NaN':''});y=official[c].replace({'NA':'','NaN':''})
 if c in ['status','reason']:assert x.equals(y),c
 else:np.testing.assert_allclose(pd.to_numeric(x),pd.to_numeric(y),rtol=1e-12,atol=1e-10,equal_nan=True)
print('All 749 independent repeats and 2,247 archived decisions matched.')

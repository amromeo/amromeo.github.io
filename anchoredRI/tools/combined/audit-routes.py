import pandas as pd,numpy as np
from pathlib import Path
scores=pd.read_csv('confirmation-scored.csv',low_memory=False).set_index(['j','method'])
def read(p):return pd.concat([pd.read_csv(f) for f in Path('art').glob(p)])
a=read('*/audit-*.csv').set_index('j')
c=read('*/checks-*.csv').set_index(['j','arm'])
for j,v in a.iterrows():
 for method,order in [('A',['anchored','t0.25']),('R',['t0.25','anchored']),('weighted',['anchored','w1']),('limited',['anchored','t0.25']),('B',['t0.25']),('C',['baseline'])]:
  r=scores.loc[j,method];route='none';limits=[np.nan,np.nan]
  for k in order:
   if k=='anchored':
    if v.status!='reported':continue
    lim=[v.lower,v.upper]
   else:
    if method=='limited' and r.N<5000:continue
    s=c.loc[j,k]
    if s.status!='keep':continue
    lim=[s.released_native_lower,s.released_native_upper]
   if np.isfinite(lim).all(): route=k;limits=lim;break
  assert route==r.route,(j,method)
  np.testing.assert_allclose(limits,[r.lower,r.upper],atol=1e-10,rtol=0,equal_nan=True)
print('All six routes and released limits reproduced on all 749 audit cases.')

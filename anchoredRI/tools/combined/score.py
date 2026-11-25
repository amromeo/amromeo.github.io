from pathlib import Path
import subprocess,sys
import numpy as np,pandas as pd
from inputs import folders
methods=['A','B','C','limited','weighted','R']
def run(folder):
 r=Path(folder)
 sim='confirmation' in folder
 idx=pd.read_csv(r/'task-g-index.csv').set_index('j').sort_index()
 fits=pd.read_csv(r/'evidence/fits.csv.gz',low_memory=False)
 checks=pd.read_csv(r/'evidence/checks.csv.gz')
 checks=checks.query("panel=='guidance' and seed_kind=='primary'")
 f=fits.query("arm=='baseline'").set_index('j').sort_index()
 arms={}
 for key in ['baseline','t0.25','w1']:
  d=checks.query('arm==@key').set_index('j').sort_index()
  assert d.index.equals(idx.index)
  arms[key]=d[['released_native_lower','released_native_upper']].set_axis(['lower','upper'],axis=1).where(d.status=='keep')
 one=pd.Series(False,index=f.index) if sim else f.Analyte=='CRP'
 if sim:
  b=pd.concat([pd.read_csv(p) for p in Path('art').glob('*/primary-*.csv')]).set_index('j').sort_index()
  assert b.index.equals(idx.index)
  arm=b[['lower','upper']].where(b.status=='reported')
  n=idx.N
 else:
  b=pd.read_csv(r/'task-g-arm-b.csv.gz').set_index('Index').sort_index()
  original=pd.read_csv(r/'task-g-arm-b-original.csv.gz')
  assert set(original.Index)==set(idx.index) and len(original)==len(idx)
  arm=b[['reported_lower','reported_upper']].set_axis(['lower','upper'],axis=1).where(b.reported_status=='returned')
  n=idx.routine_n
 for name in ['routine_hash','anchor_hash']:
  assert np.array_equal(b[name],idx[name]) and np.array_equal(f[name],idx[name])
 arms['anchored']=arm
 frames=[]
 for m,order in zip(methods,[['anchored','t0.25'],['t0.25'],['baseline'],['anchored','t0.25'],['anchored','w1'],['t0.25','anchored']]):
  d=f.copy()
  d['route']='none'
  d['method']=m
  d['lower']=d['upper']=np.nan
  for i,k in enumerate(order):
   q=arms[k]
   keep=(d.route=='none') & np.isfinite(q.upper) & (one|np.isfinite(q.lower))
   if m=='limited' and i:keep &= n>=5000
   d.loc[keep,['lower','upper']]=q.loc[keep]
   d.loc[keep,'route']=k
  if not sim:
   for e in ['lower','upper']:
    s=d['nonp_shift'];lam=d.nonp_lambda
    x=d[e].mask(d[e]<=s,s+1e-20)
    if e=='upper':
     zero=(d.upper==0).groupby(d.Analyte).transform('any')
     x=x.mask(zero,d.upper.mask(d.upper==0,1e-20))
    def bc(x):
     return np.where(abs(lam)<1e-20,np.log(x), (x**lam-1)/lam)
    d['z_'+e]=(bc(d['truth_'+e]-s)-bc(x-s))/d.nonp_sigma
   d['fraction']=d.fractionPathol;d['disease']=d.Scenario;d['cell_id']=d.Analyte.astype(str)+'/'+d.Scenario.astype(str)
  frames.append(d.reset_index())
 d=pd.concat(frames,ignore_index=True)
 if sim:
  d.to_csv('confirmation-limits.csv',index=False)
  subprocess.run(['R','--vanilla','--slave','-f','confirmation-score.R'],check=True)
  d=pd.read_csv('confirmation-scored.csv',low_memory=False)
 kept=d.route!='none'
 upper=d.z_upper.abs().clip(upper=5);lower=d.z_lower.abs().clip(upper=5)
 d['error']=(upper+lower)/2
 if not sim:d.loc[d.Analyte=='CRP','error']=upper[d.Analyte=='CRP']
 assert d.error[kept].notna().all()
 d['loss']=d.error.where(kept,5);d['returned']=kept
 d['study']='confirmation' if sim else 'RIbench'
 return d[['study','j','method','route','N','fraction','disease','cell_id','error','loss','returned']]
if __name__=='__main__':
 d=pd.concat([run(f) for f in folders if 'confirmation' not in f or '--bench' not in sys.argv])
 d.to_csv('dataset-scores.csv.gz',index=False)

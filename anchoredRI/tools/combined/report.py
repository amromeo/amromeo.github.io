import pandas as pd,numpy as np
from pathlib import Path
d=pd.read_csv('dataset-scores.csv.gz')
out=Path('results');out.mkdir(exist_ok=True)
d['pathology']=np.where(d.fraction<=.2,'<=20%','>20%')
summ=[]
for level in [[],['pathology'],['fraction'],['disease'],['fraction','disease']]:
 s=d.groupby(['study','N','method']+level).agg(assigned=('j','size'),returned=('returned','sum'),share=('returned','mean'),returned_error=('error','mean'),loss=('loss','mean')).reset_index()
 s['stratum']='+'.join(level) or 'overall';summ.append(s)
pd.concat(summ).to_csv(out/'accuracy.csv',index=False)
d.groupby(['study','N','fraction','disease','method','route']).size().rename('count').to_csv(out/'routes.csv')
rows=[]
def ci(g,diff,sim):
 rng=np.random.default_rng(37000000)
 num=np.zeros(2000);den=num.copy()
 clusters=list(g.groupby('cell_id').indices.values())
 if sim:
  for ids in clusters:
   x=diff[ids];v=rng.choice(x,(2000,len(x)));num+=np.nansum(v,axis=1);den+=np.isfinite(v).sum(axis=1)
 else:
  a=np.array([(np.nansum(diff[i]),np.isfinite(diff[i]).sum()) for i in clusters])
  draws=rng.integers(len(a),size=(2000,len(a)));num,den=a[draws].sum(axis=1).T
 return np.nanpercentile(np.divide(num,den,out=np.full(2000,np.nan),where=den>0),[2.5,97.5])
for (study,n),data in d.groupby(['study','N']):
 for subset in ['overall','<=20%']:
  g=data if subset=='overall' else data[data.fraction<=.2]
  w=g.pivot(index='j',columns='method',values=['loss','error'])
  meta=g[g.method=='A'].set_index('j').loc[w.index]
  for other in ['B','R']:
   for metric in ['loss','error']:
    x=w[metric]['A'];y=w[metric][other];diff=(x-y).to_numpy();common=np.isfinite(diff)
    delta=np.mean(diff[common]);threshold=0 if metric=='loss' and subset=='overall' else .02
    l,u=ci(meta,diff,study=='confirmation') if other=='R' else (np.nan,np.nan)
    test=study=='confirmation' and other=='B'
    rows.append(dict(study=study,N=n,subset=subset,other=other,metric=metric,count=common.sum(),A=x[common].mean(),comparator=y[common].mean(),difference=delta,lower95=l,upper95=u,passes=(delta<0 if threshold==0 else delta<=threshold) if test else None,order_flag=(delta>0) if other=='R' and metric=='loss' else None))
pd.DataFrame(rows).sort_values(['N','other','subset','metric','study']).to_csv(out/'headline.csv',index=False)

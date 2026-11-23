from pathlib import Path
import pandas as pd,numpy as np,hashlib,zipfile,json
p=Path('assets/anchors-only');a=pd.read_csv(p/'estimates.csv');b=pd.read_csv(p/'anchoredRI-matched.csv')
assert len(a)==1800 and len(b)==900 and not a.duplicated(['case','replicate','method']).any()
for m in ['anchors_auto','anchors_oracle']:
 for _,r in a[a.method==m].groupby(['family','replicate']):
  assert len(r)==15 and r.lower.nunique(dropna=False)==1 and r.upper.nunique(dropna=False)==1
x=a[a.method=='anchors_auto'].merge(b,on=['case','replicate'],suffixes=('_auto','_anchored'),validate='one_to_one')
assert len(x)==900
for field in ['seed','true_lower','true_upper','actual_affected']:
 assert np.allclose(x[field+'_auto'],x[field+'_anchored'],rtol=0,atol=1e-12)
x['difference']=x.width_error_anchored-x.width_error_auto
x.to_csv(p/'paired-datasets.csv',index=False)
rows=[]
for case,r in x.groupby('case'):
 valid=(r.status_auto=='ok')&(r.status_anchored=='ok')&np.isfinite(r.difference)
 d=r.loc[valid,'difference'];n=len(d);mu=d.mean();se=d.std(ddof=1)/np.sqrt(n) if n>1 else np.nan
 overlap=r.overlap_auto.iloc[0]
 decision='descriptive only' if overlap=='heavy' else 'indeterminate' if n!=20 else 'adds value' if mu+2*se<0 else 'not demonstrated'
 oracle=a[(a.case==case)&(a.method=='anchors_oracle')]
 rows.append(dict(case=case,family=r.family_auto.iloc[0],overlap=overlap,affected=r.affected_auto.iloc[0],attempted_pairs=20,valid_pairs=n,
  mean_anchored=r.width_error_anchored.mean(),mean_auto=r.width_error_auto.mean(),mean_oracle=oracle.width_error.mean(),
  mean_difference=mu,median_difference=d.median(),mcse=se,two_mcse=2*se,wins_anchored=int((d<0).sum()),ties=int((d==0).sum()),wins_auto=int((d>0).sum()),
  failed_auto=int((r.status_auto!='ok').sum()),failed_anchored=int((r.status_anchored!='ok').sum()),failed_oracle=int((oracle.status!='ok').sum()),
  auto_at_bound=int((r.at_bound_auto==True).sum()),decision=decision))
s=pd.DataFrame(rows);s.to_csv(p/'paired-scenarios.csv',index=False)
primary=s[s.overlap!='heavy'];lines=['# Anchors-only comparator results','',f"Completed 45 scenarios × 20 paired datasets. Routine data met the specified rule in {sum(primary.decision=='adds value')} of 30 separated/moderate cells.",'',
'Negative paired differences favor anchoredRI. The rule is mean difference + 2 × MCSE < 0, applied separately to each separated or moderate cell. Heavy-overlap cells are descriptive only.','',
'| Family | Separated: cells meeting rule | Moderate: cells meeting rule |','|---|---:|---:|']
for fam in ['normal','lognormal','shifted']:
 rr=primary[primary.family==fam];lines.append(f"| {fam} | {sum((rr.overlap=='separated')&(rr.decision=='adds value'))}/5 | {sum((rr.overlap=='moderate')&(rr.decision=='adds value'))}/5 |")
lines+=['','Each table cell reports **paired mean difference ± 2 MCSE; anchoredRI wins / 20; decision**. Differences are percentage points of error normalized by the true reference interval width. The ±2 MCSE display describes simulation precision; it is not a clinical confidence interval.','']
for fam in ['normal','lognormal','shifted']:
 lines += ['## '+fam,'','| Overlap | 10% affected | 30% affected | 50% affected | 70% affected | 90% affected |','|---|---|---|---|---|---|']
 for ov in ['separated','moderate','heavy']:
  vals=[]
  for frac in [.1,.3,.5,.7,.9]:
   r=s[(s.family==fam)&(s.overlap==ov)&np.isclose(s.affected,frac)].iloc[0]
   vals.append(f"{r.mean_difference:+.2f} ± {r.two_mcse:.2f}; {r.wins_anchored}/20; {r.decision}")
  lines.append('| '+ov+' | '+' | '.join(vals)+' |')
 lines+=['']
lines+=['## Anchors-only context','','These means use the 20 distinct anchor samples for each family, without counting their reuse across scenarios as additional independent samples.','',
'| Family | Auto mean error | Oracle mean error | Auto scale choices: identity / log / shifted-log | Auto offset bound hits |','|---|---:|---:|---|---:|']
for fam in ['normal','lognormal','shifted']:
 u=a[a.family==fam].drop_duplicates(['replicate','method']);au=u[u.method=='anchors_auto'];o=u[u.method=='anchors_oracle'];counts=[int((au.model==f).sum()) for f in ['identity','log','shifted_log']]
 lines.append(f"| {fam} | {au.width_error.mean():.2f}% | {o.width_error.mean():.2f}% | {' / '.join(map(str,counts))} | {int((au.at_bound==True).sum())}/20 |")
lines += ['',f"Failed anchors-only estimates: {sum(a.status!='ok')} of 1,800 evaluated rows. Failed archived anchoredRI estimates: {sum(b.status!='ok')} of 900.",'',
'The oracle knows the family and true shifted-log offset, but estimates its mean and sample SD from the same 40 anchors. Auto uses maximum-likelihood SD for BIC and sample SD for limits, exactly as specified in the protocol. Both use 1.96.','',
'No changes were made to R/ or tools/stress-design.R. AnchoredRI estimates were reused from the frozen completed comparison. This is a follow-up on the existing grid, not independent validation. The 30 primary decisions are not a multiplicity-adjusted global superiority test.','']
(p/'REPORT.md').write_text('\n'.join(lines))
provenance={'input_archive_sha256':hashlib.sha256(Path('assets/comparison-study/summary/comparison-results.zip').read_bytes()).hexdigest(),'estimator_rows':1800,'paired_datasets':900,'scenarios':45,'unique_anchor_samples':60,'source_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in Path('tools/anchors-only').glob('*') if f.is_file()}}
(p/'provenance.json').write_text(json.dumps(provenance,indent=2))
with zipfile.ZipFile(p/'anchors-only-results.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in p.iterdir():
  if f.suffix in ['.csv','.md','.txt','.json']:z.write(f,f.name)
 for f in Path('tools/anchors-only').glob('*'):
  if f.is_file():z.write(f,'source/'+f.name)
print('\n'.join(lines))

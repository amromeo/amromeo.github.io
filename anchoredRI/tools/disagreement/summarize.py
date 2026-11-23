from pathlib import Path
import pandas as pd, numpy as np, hashlib, json, zipfile
out=Path('assets/disagreement')
d=pd.read_csv(out/'datasets.csv')
assert len(d)==900 and not d[['case','replicate']].duplicated().any()
assert d.loc[d.primary_declined,['reported_lower','reported_upper','reported_error']].isna().all().all()
assert np.allclose(d.loc[~d.primary_declined,'reported_lower'],d.loc[~d.primary_declined,'diagnostic_anchored_lower'])
assert np.allclose(d.loc[~d.primary_declined,'reported_upper'],d.loc[~d.primary_declined,'diagnostic_anchored_upper'])
def stats(g,check):
 flag=g[check+'_flag'].fillna(False).astype(bool);decl=g[check+'_declined'];e=g.original_error;keep=e[~decl];bad=e>10;good=e<5
 z=dict(check=check,n=len(g),flagged=int(flag.sum()),unavailable=int(g[check+'_flag'].isna().sum()),declined=int(decl.sum()),reported=len(keep),error_gt10_n=int(bad.sum()),error_gt10_flagged=int((bad&flag).sum()),error_lt5_n=int(good.sum()),error_lt5_flagged=int((good&flag).sum()))
 z.update(flag_pct=100*flag.mean(),decline_pct=100*decl.mean(),error_gt10_flag_pct=100*flag[bad].mean(),error_lt5_flag_pct=100*flag[good].mean(),reported_mean_error=keep.mean(),reported_median_error=keep.median(),reported_p90_error=keep.quantile(.9),reported_max_error=keep.max(),reported_gt10_pct=100*(keep>10).mean())
 return z
rows=[]
for check in ['primary','secondary']:
 rows.append(dict(family='all',overlap='all',affected='all',**stats(d,check)))
 for (f,o,a),g in d.groupby(['family','overlap','affected']):rows.append(dict(family=f,overlap=o,affected=a,**stats(g,check)))
 for (o,a),g in d.groupby(['overlap','affected']):rows.append(dict(family='all',overlap=o,affected=a,**stats(g,check)))
s=pd.DataFrame(rows);s.to_csv(out/'summary.csv',index=False)
lines=['# Step B: disagreement check','', 'All 900 existing datasets were evaluated. Each of 60 distinct anchor samples was bootstrapped 1,000 times, with family and offset selection repeated. Shared bootstrap calculations were checked against identical anchors in all 15 corresponding scenarios. No thresholds were changed; replicates 21–40 were therefore not used.','', 'Primary rule: report the original anchoredRI interval only if both limits fall within their respective marginal 95% bootstrap intervals. Otherwise decline and recommend at least 120 screened patients or verification of a suitable published interval with 20. Anchors-only estimates are diagnostic, never a fallback. This is a study evaluation; production package functions have not been changed.','']
for check in ['primary','secondary']:
 z=stats(d,check)
 lines += [f'## {check.capitalize()} check','',f"Flagged {z['flagged']}/900 ({z['flag_pct']:.1f}%); unavailable {z['unavailable']}; declined {z['declined']}/900 ({z['decline_pct']:.1f}%).",'',f"Flagged {z['error_gt10_flagged']}/{z['error_gt10_n']} datasets with original error >10% ({z['error_gt10_flag_pct']:.1f}%), and {z['error_lt5_flagged']}/{z['error_lt5_n']} with error <5% ({z['error_lt5_flag_pct']:.1f}%).",'',f"Retained {z['reported']} intervals: mean error {z['reported_mean_error']:.2f}%, median {z['reported_median_error']:.2f}%, 90th percentile {z['reported_p90_error']:.2f}%, maximum {z['reported_max_error']:.2f}%; {z['reported_gt10_pct']:.1f}% still had error >10%.",'']
 if check=='secondary':lines+=['Secondary declines and retained errors are hypothetical applications of the LRT alone. They do not change the primary output. The 2-df chi-squared reference is exploratory because the transformation was selected and anchors contributed to the original fit.','']
 for family in ['all','normal','lognormal','shifted']:
  lines += [f'### {family}: flag counts and rates','', '| Overlap | 10% affected | 30% | 50% | 70% | 90% |','|---|---:|---:|---:|---:|---:|']
  for o in ['separated','moderate','heavy']:
   cells=[]
   for a in [.1,.3,.5,.7,.9]:
    z=s[(s.check==check)&(s.family==family)&(s.overlap==o)&(s.affected==a)].iloc[0];cells.append(f"{z.flagged}/{z.n} ({z.flag_pct:.0f}%)")
   lines.append('| '+o+' | '+' | '.join(cells)+' |')
  lines+=['']
lines+=['## Interpretation limits','','These percentages describe this simulation, with shared anchors across scenarios. The two marginal CIs do not imply a 5% overall flag rate. Performance among reported intervals is conditional on passing; it must be read together with the decline rate and the errors that remain undetected. No thresholds were tuned after inspection. The original anchoredRI error is measured against simulated truth; that truth is never used by the disagreement rule.','','Full per-scenario counts, denominators and retained-error summaries are in summary.csv. datasets.csv separates diagnostic estimates from reported limits and includes recommendations. bootstrap-estimates.csv contains all 60,000 bootstrap fits.']
(out/'REPORT.md').write_text('\n'.join(lines)+'\n')
paths=sorted(Path('R').glob('*'))+[Path('tools/stress-design.R'),Path('tools/anchors-only/PROTOCOL.md'),Path('tools/anchors-only/estimators.R'),Path('tools/disagreement/run.R'),Path('tools/disagreement/summarize.py'),Path('assets/comparison-study/summary/comparison-results.zip')]
(out/'provenance.json').write_text(json.dumps({'protocol_commit':'dea1da81da6344e57b5c82361e682b75e5329df4','threshold_changes':False,'bootstrap_replicates':1000,'datasets':900,'unique_anchor_samples':60,'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()}},indent=2)+'\n')
with zipfile.ZipFile(out/'disagreement-results.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in sorted(out.glob('*')):
  if p.suffix in ['.csv','.json','.md','.txt']:z.write(p,p.name)
 for p in paths:
  if str(p).startswith('tools/'):z.write(p,'source/'+str(p))
print(s[(s.family=='all')&(s.overlap=='all')].to_string(index=False))

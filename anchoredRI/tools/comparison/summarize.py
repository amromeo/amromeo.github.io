"""Validate the complete paired grid and summarize without dropping failures."""
from pathlib import Path
import csv, math, statistics as st, zipfile, json
from collections import defaultdict
ROOT=Path('assets/comparison-study'); OUT=ROOT/'summary';OUT.mkdir(parents=True,exist_ok=True)
METHODS=['refineR','gmm_largest','gmm_lowest','pooled40','anchored40','direct40','direct120','direct400']
def read(p):
 with p.open() as f:return list(csv.DictReader(f))
def write(name,rows):
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def num(r,k):
 try:return float(r[k])
 except (ValueError,TypeError):return float('nan')
def mean(x):return st.mean(x) if x else float('nan')
def q(x,p):
 if not x:return float('nan')
 a=sorted(x);i=(len(a)-1)*p;lo=int(i);hi=math.ceil(i);return a[lo]+(a[hi]-a[lo])*(i-lo)
files=sorted(ROOT.glob('results-*.csv'))
rows=[r for p in files for r in read(p)]
expected={(c,r,m) for c in range(1,46) for r in range(1,21) for m in METHODS}
keys=[(int(r['case']),int(r['replicate']),r['method']) for r in rows]
missing=expected-set(keys)
(OUT/'completeness.json').write_text(json.dumps({'case_files':len(files),'rows':len(rows),'expected_rows':7200,'missing_rows':len(missing),'duplicates':len(keys)-len(set(keys))},indent=2))
assert not missing and len(keys)==len(set(keys)) and set(keys)==expected, 'Incomplete or duplicated grid; see completeness.json'
hashes=[read(p) for p in sorted(ROOT.glob('hashes-*.csv'))]
assert len(hashes)==45 and all(h==hashes[0] for h in hashes), 'Mixed estimator or runner versions'
index=dict(zip(keys,rows));group=defaultdict(list)
for r in rows:
 c=int(r['case']);rep=int(r['replicate'])
 for field in ['seed','actual_affected','actual_a','actual_b','true_lower','true_upper']:
  assert r[field]==index[c,rep,'anchored40'][field], (c,rep,field)
 assert int(r['n_screened'])==({'pooled40':40,'anchored40':40,'direct40':40,'direct120':120,'direct400':400}.get(r['method'],0))
 group[c,r['method']].append(r)
summary=[]
for (c,m),rr in sorted(group.items()):
 good=[r for r in rr if r['status']=='ok' and math.isfinite(num(r,'width_error'))];e=[num(r,'width_error') for r in good];base=rr[0]
 summary.append(dict(case=c,family=base['family'],overlap=base['overlap'],affected=base['affected'],method=m,
  attempted=len(rr),successful=len(good),failed=len(rr)-len(good),mean_error=mean(e),median_error=q(e,.5),p90_error=q(e,.9),
  mcse_mean=st.stdev(e)/math.sqrt(len(e)) if len(e)>1 else float('nan'),
  large_error_count=sum(v>10 for v in e),large_or_failed_fraction=(len(rr)-len(good)+sum(v>10 for v in e))/len(rr),
  boundary_count=sum(r['at_bound']=='TRUE' for r in good),mean_lower_error=mean([num(r,'lower_error') for r in good]),
  mean_upper_error=mean([num(r,'upper_error') for r in good]),
  mean_below_reference=mean([num(r,'below_reference') for r in good]),mean_above_reference=mean([num(r,'above_reference') for r in good])))
contrasts=[('anchored40','refineR'),('pooled40','gmm_largest'),('anchored40','pooled40'),('pooled40','gmm_lowest'),('anchored40','direct40'),('anchored40','direct120'),('anchored40','direct400')]
paired=[]
for c in range(1,46):
 for a,b in contrasts:
  diff=[];fail_a=fail_b=0
  for rep in range(1,21):
   ra,rb=index[c,rep,a],index[c,rep,b];ea,eb=num(ra,'width_error'),num(rb,'width_error')
   oka=ra['status']=='ok' and math.isfinite(ea);okb=rb['status']=='ok' and math.isfinite(eb)
   fail_a+=not oka;fail_b+=not okb
   if oka and okb:diff.append(ea-eb)
  base=index[c,1,a]
  paired.append(dict(case=c,family=base['family'],overlap=base['overlap'],affected=base['affected'],method_a=a,method_b=b,
   attempted_pairs=20,finite_pairs=len(diff),failed_a=fail_a,failed_b=fail_b,
   mean_difference=mean(diff),mcse_difference=st.stdev(diff)/math.sqrt(len(diff)) if len(diff)>1 else float('nan'),
   median_difference=q(diff,.5),wins_a=sum(d<0 for d in diff),ties=sum(d==0 for d in diff),wins_b=sum(d>0 for d in diff)))
write('all-results.csv',sorted(rows,key=lambda r:(int(r['case']),int(r['replicate']),r['method'])))
write('scenario-summary.csv',summary);write('paired-comparisons.csv',paired)
lines=['# Repeated comparison results','', 'All 45 scenarios and 20 replicates completed: 900 matched datasets. Errors are percentages of true healthy interval width. refineR uses fullDataEst, without bootstrap aggregation.','', '| Method | Mean error | Failed / 900 | Error >10% or failed / 900 |','|---|---:|---:|---:|']
for m in METHODS:
 ss=[r for r in summary if r['method']==m]
 lines.append(f"| {m} | {mean([r['mean_error'] for r in ss]):.2f} | {sum(r['failed'] for r in ss)} | {sum(r['failed']+r['large_error_count'] for r in ss)} |")
lines+=['','Mean errors average scenario means; failure counts are reported separately. The 10% threshold is exploratory, not a clinical acceptance threshold.','', '| Paired comparison (A minus B) | Mean difference | Scenarios favoring A / 45 | Finite paired wins for A |','|---|---:|---:|---:|']
for a,b in contrasts:
 pp=[r for r in paired if r['method_a']==a and r['method_b']==b]
 lines.append(f"| {a} minus {b} | {mean([r['mean_difference'] for r in pp]):.2f} | {sum(r['mean_difference']<0 for r in pp)} | {sum(r['wins_a'] for r in pp)} / {sum(r['finite_pairs'] for r in pp)} |")
lines+=['','Negative differences favor A. Shared anchors across scenarios prevent treating all 900 paired differences as independent for a pooled confidence interval. Scenario-specific Monte Carlo SEs are in paired-comparisons.csv.','', 'These are follow-up simulations within the model family. They do not establish clinical validity, novelty, or performance under misspecification. The lowest-component rule exploits the direction of disease shifts in this design.']
(OUT/'REPORT.md').write_text('\n'.join(lines)+'\n')
with zipfile.ZipFile(OUT/'comparison-results.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in OUT.iterdir():
  if p.suffix in ['.csv','.md','.json']:z.write(p,p.name)
 for p in ROOT.glob('hashes-*.csv'):z.write(p,'provenance/'+p.name)
 for p in ROOT.glob('versions-*.txt'):z.write(p,'provenance/'+p.name)
print('\n'.join(lines))

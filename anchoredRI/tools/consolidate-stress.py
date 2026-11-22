"""Validate downloaded study artifacts and build published data. Run from anchoredRI/."""
import csv, json, math, statistics as st, pathlib, sys, zipfile, hashlib
src=pathlib.Path(sys.argv[1]); dest=pathlib.Path('assets/stress-grid');dest.mkdir(exist_ok=True,parents=True)
examples=[]; repeats=[]
for case in range(1,46):
 root=src/f'stress-case-{case}'
 hashes=list(csv.DictReader((root/'source-hashes.csv').open()))
 for h in hashes: assert hashlib.md5(pathlib.Path(h['path']).read_bytes()).hexdigest()==h['md5'],h['path']
 for prefix,target,n in [('example',examples,5),('repeat',repeats,80)]:
  rr=list(csv.DictReader((root/'cases'/f'{prefix}-{case:02}.csv').open()));assert len(rr)==n
  assert all(int(r['case'])==case for r in rr)
  for r in rr:
   if r['status']=='ok':
    e=100*max(abs(float(r['lower'])-float(r['true_lower'])),abs(float(r['upper'])-float(r['true_upper'])))/(float(r['true_upper'])-float(r['true_lower']))
    assert math.isclose(e,float(r['width_error']),rel_tol=1e-9,abs_tol=1e-9)
  target.extend(rr)
assert len({(r['case'],r['replicate'],r['method']) for r in repeats})==3600
for c in range(1,46):
 rr=[r for r in examples if int(r['case'])==c]
 assert len({(r['seed'],r['true_lower'],r['true_upper'],r['actual_affected']) for r in rr})==1
 assert all(r['bootstrap_requested']=='200' for r in rr)
 for m in ['anchored40','direct40','direct120','direct400']:
  a=next(r for r in rr if r['method']==m);b=next(r for r in repeats if int(r['case'])==c and r['method']==m and r['replicate']=='1')
  assert math.isclose(float(a['lower']),float(b['lower']),rel_tol=1e-7)
  assert math.isclose(float(a['upper']),float(b['upper']),rel_tol=1e-7)
summary=[]
for c in range(1,46):
 for m in ['anchored40','direct40','direct120','direct400']:
  rr=[r for r in repeats if int(r['case'])==c and r['method']==m];good=[r for r in rr if r['status']=='ok'];v=[float(r['width_error']) for r in good]
  summary.append(dict(case=c,family=rr[0]['family'],overlap=rr[0]['overlap'],affected=float(rr[0]['affected']),method=m,attempted=len(rr),successful=len(good),failed=len(rr)-len(good),mean_error=st.mean(v),mcse=st.stdev(v)/math.sqrt(len(v)),large=sum(x>10 for x in v),mean_mpe=st.mean(float(r['mpe']) for r in good)))
def csvtext(rows):
 import io
 s=io.StringIO();w=csv.DictWriter(s,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows);return s.getvalue()
(dest/'examples.csv').write_text(csvtext(examples));(dest/'summary.csv').write_text(csvtext(summary))
(dest/'results.json').write_text(json.dumps(dict(examples=examples,summary=summary),separators=(',',':')))
with zipfile.ZipFile(dest/'all-results.zip','w',zipfile.ZIP_DEFLATED) as z:
 z.writestr('examples.csv',csvtext(examples));z.writestr('repeated-fits.csv',csvtext(repeats));z.writestr('summary.csv',csvtext(summary))
 z.write('STRESS_PROTOCOL.md','STRESS_PROTOCOL.md')
 for name in ['design.csv','source-hashes.csv']:z.write(src/'stress-case-1'/name,name)
 for c in range(1,46):
  for p in (src/f'stress-case-{c}'/'cases').glob('*'):
   if not p.name.startswith(('example-','repeat-')):z.write(p,'diagnostics/'+p.name)
print('Verified 45 matched comparisons, 900 repeated datasets, source hashes and paired point estimates.')
print('Archive bytes:',(dest/'all-results.zip').stat().st_size)

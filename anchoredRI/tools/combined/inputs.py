from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from urllib.request import urlopen
import hashlib,json,time
root='https://raw.githubusercontent.com/amromeo/anchorRefine/638e2d3500faa01d4e2816895a95d7a6685567a0/'
def fetch(folder,p):
 path=Path(folder,p);path.parent.mkdir(parents=True,exist_ok=True)
 for t in range(4):
  try:
   b=urlopen(root+'results/'+folder+'/'+p,timeout=100).read();path.write_bytes(b);return
  except Exception:
   if t==3: raise
folders=['stage3-confirmation-45463873162','stage3-ribench-small-45530368829','stage3-ribench-large-45541935930']
if __name__=='__main__':
 with ThreadPoolExecutor(12) as pool:
  jobs=[pool.submit(fetch,f,p) for f in folders for p in ['task-g-index.csv','task-g-manifest.json','evidence/fits.csv.gz','evidence/checks.csv.gz']+(['task-g-arm-b.csv.gz','task-g-arm-b-original.csv.gz'] if 'ribench' in f else [])]
  for j in jobs:j.result()
 for f in folders:
  md=json.loads(Path(f,'task-g-manifest.json').read_text())
  for p,h in md['sha256'].items():
   path=Path(f,p)
   if path.exists():assert hashlib.sha256(path.read_bytes()).hexdigest()==h,path

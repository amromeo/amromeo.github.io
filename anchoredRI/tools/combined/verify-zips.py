from concurrent.futures import ThreadPoolExecutor
import hashlib
from urllib.request import urlopen
from inputs import root,folders
def check(line):
 h,name=line.split(maxsplit=1); name=name.strip().lstrip('*')
 for attempt in range(3):
  try:
   dig=hashlib.sha256()
   with urlopen(root+'results/'+folders[0]+'/'+name,timeout=90) as r:
    while b:=r.read(1<<20):dig.update(b)
   assert dig.hexdigest()==h,name
   return name
  except Exception:
   if attempt==2:raise
lines=[x for x in open('source-sums') if 'artifacts/' in x and x.strip().endswith('.zip')]
assert len(lines)>=96
with ThreadPoolExecutor(12) as p: names=list(p.map(check,lines))
print(len(names),'Stage 3 confirmation ZIPs verified')

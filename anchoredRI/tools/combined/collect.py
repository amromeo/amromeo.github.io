from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen
from io import BytesIO
import json,zipfile,hashlib,time
run='45614043039'
api='https://api.github.com/repos/amromeo/amromeo.github.io/actions/runs/'+run
def get(a):
 name=a['name'];p=Path('art',name)
 if p.exists():return
 url='https://nightly.link/amromeo/amromeo.github.io/actions/runs/'+run+'/'+name+'.zip'
 try: data=urlopen(url,timeout=45).read()
 except Exception:
  data=urlopen('https://nightly.link/amromeo/amromeo.github.io/actions/artifacts/'+str(a['id'])+'.zip',timeout=45).read()
 if a.get('digest'):assert 'sha256:'+hashlib.sha256(data).hexdigest()==a['digest']
 with zipfile.ZipFile(BytesIO(data)) as z:
  assert all('/' not in n and '..' not in n for n in z.namelist())
  z.extractall(p)
for i in range(240):
 try:
  state=json.load(urlopen(api,timeout=30))
  a=json.load(urlopen(api+'/artifacts?per_page=100',timeout=30))['artifacts']
  a=[x for x in a if x['name'].startswith(('taskg-audit-','taskg-primary-'))]
  Path('artifacts.json').write_text(json.dumps(a))
  print(state['status'],state['conclusion'],len(a),flush=True)
  if state['conclusion'] not in (None,'success'):raise SystemExit(2)
  with ThreadPoolExecutor(8) as pool:list(pool.map(get,a))
  if len(a)==64 and state['conclusion']=='success':break
 except Exception as e:print(str(e),flush=True)
 time.sleep(30)

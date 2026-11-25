from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
import hashlib,json
p=Path('results');target=p/'task-g-evidence.zip'
with ZipFile(target,'w',ZIP_DEFLATED,compresslevel=6) as z:
 for f in sorted(Path('art').glob('*/*')):
  if f.is_file():z.write(f)
 for f in ['dataset-scores.csv.gz','audit.csv','source-sums','zip.log','artifacts.json']:
  z.write(f)
 for f in sorted(Path('.').glob('stage3-*/task-g-manifest.json')):z.write(f)
 for f in ['headline.csv','accuracy.csv','routes.csv']:z.write(p/f,'tables/'+f)
with ZipFile(target) as z:assert z.testzip() is None
print(hashlib.sha256(target.read_bytes()).hexdigest(),target)

from pathlib import Path
import hashlib
files=list(Path('art').glob('*/runtime.sha256'))
assert len(files)==64 and len({p.read_text().strip() for p in files})==1
print('All 64 fit jobs used the same compiled R library:',files[0].read_text().split()[0])

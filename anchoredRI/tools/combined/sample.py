import hashlib
import pandas as pd
import sys
d = pd.read_csv(sys.argv[1])
assert len(d)==7488 and d.j.is_unique
keys=['N','fraction','disease']
groups=d.groupby(keys)
count=groups.size().rename('n').reset_index()
count['quota']=count.n*749//7488
count['rem']=count.n*749%7488
count=count.sort_values(['rem']+keys,ascending=[False,True,True,True])
count.loc[count.index[:749-count.quota.sum()],'quota']+=1
d['rank']=d.j.map(lambda j:hashlib.sha256(f'TaskG-audit-v1:{j}'.encode()).hexdigest())
d=d.merge(count,on=keys).sort_values('rank')
audit=d[d.groupby(keys).cumcount()<d.quota]
assert len(audit)==749
audit[['j']].sort_values('j').to_csv('audit.csv',index=False)

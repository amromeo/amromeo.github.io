import csv,json,math,statistics,sys
from collections import Counter
from pathlib import Path
def read(p):
    with p.open() as f:return list(csv.DictReader(f))
def num(s):return float(s) if s not in (None, '', 'NA') else float('nan')
def group(rows,*fields):
    return {tuple(r[k] for k in fields):r for r in rows}
def q(v,p):
    if not v:return None
    v=sorted(v)
    a=(len(v)-1)*p
    i=int(a)
    return v[i]+(a-i)*(v[min(i+1,len(v)-1)]-v[i])
def pct(v):
    return {k:q(v,p) for k,p in [('median',.5),('q90',.9),('q95',.95),('max',1)]}
base,score_path,out=map(Path,sys.argv[1:])
roster=read(Path('anchoredRI/tools/subsample/pilot-roster.csv'))
dirs={p.parent.name:p.parent for p in base.rglob('COMPLETE')}
assert set(dirs)=={r['Index'] for r in roster}, (len(dirs), list(dirs))
score=group(read(score_path),'Index','method')
assert len(score)==240, len(score)
dataset={}
for row in roster:
    idx=row['Index'];d=dirs[idx]
    med=group(read(d/'medians.csv'),'cap','seed_set','definition')
    fits=read(d/'fits.csv');fs=next(f for f in fits if f['fit_type']=='full' and f['cap']=='5')
    dec=group(read(d/'decisions.csv'),'cap','seed_set','definition','estimator','rule')
    dataset[idx]=dict(m=med,f=fs,d=dec,analyte=row['Analyte'])
summary={}
movement=[]
for cap in ('5','3'):
    for definition in ('all','passing'):
        cell={}
        summary[f'{cap}_{definition}']=cell
        changed=[]; gained=[];lost=[];both=[]
        extremes=[];dm=[];zs=[];contrib=[];zero=[];rel_extreme=[];raw_extreme_method=Counter()
        for idx,obj in dataset.items():
            rz={}
            ds=obj['d'];m=obj['m'];analyte=obj['analyte']
            for bank in ('A','B'):
                item=m[(cap,bank,definition)]
                count=int(item['contributors']);contrib.append(count)
                if count==0:zero.append((idx,bank))
                rz[bank]=ds[(cap,bank,definition,'median','majority')]['status']
                scr=score[(idx,f'median_{cap}_{bank}_{definition}')]
                if num(scr['d'])>5:dm.append((idx,bank,round(num(scr['d']),3),rz[bank]))
                for side in ('lower','upper'):
                    z=num(scr['z_'+side])
                    if math.isfinite(z) and abs(z)>5:
                        extremes.append((idx,bank,side,round(z,3),rz[bank]))
                        if rz[bank]=='reported':rel_extreme.append((idx,bank,side))
            if rz['A']!=rz['B']:changed.append((idx,obj['analyte'],rz,int(m[(cap,'A',definition)]['votes_decline']),int(m[(cap,'B',definition)]['votes_decline'])))
            if rz['A']=='reported' and rz['B']=='reported':both.append(idx)
            if rz['A']!='reported' and rz['B']=='reported':gained.append(idx)
            if rz['A']=='reported' and rz['B']!='reported':lost.append(idx)
            for side in ('lower','upper'):
                a=num(m[(cap,'A',definition)][side]);b=num(m[(cap,'B',definition)][side])
                za=num(score[(idx,f'median_{cap}_A_{definition}')]['z_'+side]);zb=num(score[(idx,f'median_{cap}_B_{definition}')]['z_'+side])
                if math.isfinite(a+b):
                    entry=dict(cap=cap,definition=definition,Index=idx,analyte=analyte,side=side,
                               delta_units=b-a,delta_z=zb-za,both=rz['A']==rz['B']=='reported')
                    movement.append(entry)
        cell.update(changed=changed, gained=gained,lost=lost,both_release=len(both),available=sum(x>0 for x in contrib),min_contributors=min(contrib),max_contributors=max(contrib),zero=zero,extremes=extremes,d_gt5=dm,released_extremes=rel_extreme)
        for pop in ('raw','both'):
            for side in ('lower','upper'):
                pairs=[e for e in movement if e['cap']==cap and e['definition']==definition and (pop=='raw' or e['both']) and e['side']==side and math.isfinite(e['delta_z'])]
                vals=[abs(e['delta_z']) for e in pairs]
                cell[pop+'_'+side] = dict(n=len(vals), **pct(vals), largest=sorted(pairs,key=lambda e:-abs(e['delta_z']))[:2])
fallback=[(idx, d['analyte'], d['f']['first_status'],d['f']['status'],d['f']['fallback_status'],d['f']['method']) for idx,d in dataset.items() if d['f']['fallback_status']]
def score_ext(method):
    table=[]
    for idx in dataset:
        s=score[(idx,method)]
        for side in ('lower','upper'):
            z=num(s['z_'+side])
            if abs(z)>5: table.append((idx,side,round(z,3)))
    return table
full=dict(results=Counter((o['f']['status'] for o in dataset.values())), fallback=fallback,
          extremes={method:score_ext(method) for method in ('standard_raw','package','W5','W3')})
out.mkdir(exist_ok=True)
with (out/'readback.json').open('w') as f:json.dump(dict(full=full,summary=summary),f,indent=2)
with (out/'movement.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(movement[0]));writer.writeheader();writer.writerows(movement)

def save_rows(name, rows):
    with (out/name).open('w') as outf:
        wr=csv.DictWriter(outf, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)
grouped=[];stratified=[];gates=[]
for cap in ('5','3'):
    for definition in ('all','passing'):
        for pop in ('raw','joint_majority_releases'):
            pairs=[e for e in movement if e['cap']==cap and e['definition']==definition and (pop=='raw' or e['both'])]
            for side in ('lower','upper'):
                sel=[e for e in pairs if e['side']==side]
                # Measurement units are pooled only within a single analyte.
                for analyte in sorted(set(r['Analyte'] for r in roster)):
                    if analyte=='CRP' and side=='lower':continue
                    analyte_sel=[e for e in sel if e['analyte']==analyte]
                    for mode in ('signed','absolute'):
                        v=[e['delta_units'] if mode=='signed' else abs(e['delta_units']) for e in analyte_sel]
                        grouped.append(dict(cap=cap, definition=definition, population=pop, side=side, analyte=analyte, movement=mode, n=len(v), **pct(v), largest_abs_Index=max(analyte_sel,key=lambda e:abs(e['delta_units']))['Index'] if v else None))
                sel=[e for e in sel if math.isfinite(e['delta_z'])]
                vals=[abs(e['delta_z']) for e in sel]
                stratified.append(dict(cap=cap,definition=definition,population=pop,side=side,n=len(vals),**pct(vals),largest_abs_Index=max(sel,key=lambda e:abs(e['delta_z']))['Index'] if vals else None))
        for idx,o in dataset.items():
            counts={}
            for bank in ('A','B'):
                m=o['m'][cap,bank,definition]
                d=o['d'][cap,bank,definition,'median','majority']
                for field in ('contributors','votes_decline','votes_technical','votes_pass'):counts[bank+'_'+field]=m[field]
                counts[bank+'_decision']=d['status']
            gates.append(dict(Index=idx,analyte=o['analyte'],cap=cap,definition=definition,**counts))
save_rows('seed-units.csv',grouped)
save_rows('seed-z.csv',stratified)
save_rows('seed-decisions.csv',gates)
print(json.dumps(full))
for key,s in summary.items():
 print(key,{k:v for k,v in s.items() if k not in ('raw_lower','raw_upper','both_lower','both_upper')})
 for side in ['raw_lower','raw_upper','both_lower','both_upper']:
    v=s[side]; print(side, {k:round(val,4) if isinstance(val,float) else val for k,val in v.items() if k!='largest'},v['largest'][:1])

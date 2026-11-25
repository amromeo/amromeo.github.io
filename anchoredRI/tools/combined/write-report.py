from pathlib import Path
import pandas as pd
p=Path('results')
h=pd.read_csv(p/'headline.csv')
g=h[h.passes.notna()]
assert len(g)==16
ng=int(g.passes.sum())
parts=["# Task G confirmation results",f"**{ng} of 16 conditions passed.** "+("The combined method met the prespecified success rule." if ng==16 else "The combined method missed the prespecified success rule.")]
flags=h[(h.study=='confirmation')&(h.order_flag==True)]
if len(flags):parts.append("**Order review required.** The blend-first order had lower mean loss in "+", ".join(f"N = {r.N:,} ({r.subset})" for r in flags.itertuples())+". The study stops here for review.")
else:parts.append("The blend-first order did not have lower confirmation loss in any prespecified stratum. The study stops here for review.")
rbflags=h[(h.study=='RIbench')&(h.order_flag==True)]
if len(rbflags):parts.append('Order flag on descriptive RIbench: R had lower loss at '+', '.join(f'N = {r.N:,} ({r.subset})' for r in rbflags.itertuples())+'.')
parts.append("A tries safeguarded anchoredRI first. B is the checked Stage 3 post-fit blend alone. R tries B first and uses anchoredRI as fallback. A and R can return on the same datasets. Their order matters when both component methods return different intervals. C is checked refineR. All losses count a nonreturn as 5. Accuracy uses only paired common returns.")
for comp in ['B','R']:
 parts.append('## '+('A versus B: success conditions' if comp=='B' else 'A versus R: co-primary order comparison'))
 parts.append('Differences are A minus the comparator. Negative values favor A. Parentheses give the paired denominator.'+(' Brackets give paired 95% bootstrap intervals.' if comp=='R' else ' Pass/fail uses unrounded confirmation values.'))
 tab=['| N | Pathology | Metric | Confirmation | RIbench (descriptive) |','|---:|---|---|---|---|']
 for (n,sub,metric),grp in h[h.other==comp].groupby(['N','subset','metric'],sort=True):
  cells=[]
  for study in ['confirmation','RIbench']:
   t=grp[grp.study==study]
   if t.empty:cells.append('Unavailable');continue
   r=t.iloc[0]
   v=f'{r.difference:+.4f} ({int(r["count"]):,})'
   if comp=='R':v+=f' [{r.lower95:+.4f}, {r.upper95:+.4f}]'
   if comp=='B' and study=='confirmation':v+=' pass' if r.passes else ' FAIL'
   cells.append(v)
  tab.append(f'| {n:,} | {sub} | {metric} | '+' | '.join(cells)+' |')
 parts.append('\n'.join(tab))
parts.append("## Study population and reproducibility")
parts.append("Confirmation comprised 7,488 Stage 3 synthetic datasets. The healthy families follow anchorRefine’s Box-Cox design, including lambda = 0.5, which anchoredRI’s transforms do not represent exactly. Extra replication of high-pathology shoulder and overlapping patterns contributes to the overall means. RIbench comprised 4,320 datasets and informed the original hypothesis. It has no N = 1,500 stratum. Its results are descriptive.")
parts.append("The first run used different Ubuntu image revisions across some paired jobs and failed its audit. All primary and audit jobs were repeated using the same immutable R 4.5.0 container and the same compiled library. No fits were selected between runs. All 749 prespecified audit repeats in the pinned run matched the primary anchoredRI fits and decisions at the protocol tolerance. Independent reconstruction checked all six routes against the shard releases. Standalone B and C were checked against the exported native scores. The paired intervals resample confirmation datasets within the original design cells and RIbench parent scenarios as clusters. [Actions run](https://github.com/amromeo/amromeo.github.io/actions/runs/45614043039) · [First-run discrepancies](first-audit-discrepancies.csv) · [Audit IDs](../audit.csv) · [Protocol](../PROTOCOL.md)")
parts.append("The tables give each method’s loss, returned-only error, return count and share by size, pathological fraction and disease pattern. Route counts identify which method supplied each interval. [All six methods by stratum](accuracy.csv) · [Route counts](routes.csv) · [Paired headline comparisons and denominators](headline.csv)")
parts.append("[Download the per-dataset scores, fit records, audit repeats and source manifests](task-g-evidence.zip).")
(p/'REPORT.md').write_text('\n\n'.join(parts)+'\n')

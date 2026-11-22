'use strict';
const labels={anchored40:'anchoredRI + 40 anchors + services',pooled40:'Pooled mixture + 40 anchors',refineR:'refineR',gmm_largest:'Mixture without anchors (largest component)',gmm_lowest:'Mixture without anchors (lowest-component sensitivity)',direct40:'40 screened alone',direct120:'120 screened alone',direct400:'400 screened alone'};
const methods=['anchored40','pooled40','refineR','gmm_largest','gmm_lowest','direct40','direct120','direct400'];
const overlaps=['separated','moderate','heavy'],names=['Separated','Moderate','Heavy'],fractions=[.1,.3,.5,.7,.9];
const el=id=>document.getElementById(id),fmt=x=>x===''||x==null||!Number.isFinite(Number(x))?'Unavailable':Number(x).toFixed(2);
const get=url=>fetch(url).then(r=>{if(!r.ok)throw Error('Could not load results');return r.json();});
Promise.all([get('assets/stress-grid/results.json'),get('assets/comparison-study/summary/results.json')]).then(([old,data])=>{
 if(old.examples.length!==225||data.summary.length!==360||data.paired.length!==315||data.overall.length!==8)throw Error('Incomplete result files');
 el('overall-rows').innerHTML=data.overall.map(r=>'<tr><th scope="row">'+labels[r.method]+'</th><td>'+fmt(r.mean_error)+'%</td><td>'+r.large+'/'+r.attempted+'</td><td>'+r.failed+'/'+r.attempted+'</td></tr>').join('');
 function render(){
 const family=el('family').value,repeat=el('evidence').value==='repeat';
 el('explanation').textContent=repeat?'Mean error across 20 fresh datasets per scenario, now including refineR and both pooled mixture comparisons. The small number counts estimates with error >10%. All methods use the same datasets; refineR uses fullDataEst without bootstrapping.':'Earlier study: one matched dataset per scenario with 200 bootstrap resamples per method. refineR uses medianBS. These values are individual examples, not the repeated comparison above; pooled mixtures were not bootstrapped.';
 const shown=repeat?methods:['anchored40','refineR','direct40','direct120','direct400'];
 el('maps').innerHTML=shown.map(m=>'<article><h3>'+labels[m]+'</h3><table><thead><tr><th>Overlap</th>'+fractions.map(p=>'<th scope="col">'+p*100+'%</th>').join('')+'</tr></thead><tbody>'+overlaps.map((o,i)=>'<tr><th scope="row">'+names[i]+'</th>'+fractions.map(p=>{
 const r=(repeat?data.summary:old.examples).find(r=>r.family===family&&r.overlap===o&&+r.affected===p&&r.method===m);
 const e=Number(repeat?r.mean_error:r.width_error),color=e<=5?'#d5eadc':e<=10?'#f1ebba':e<=25?'#f4cfb3':'#efb4b0';
 return '<td style="background:'+color+'">'+e.toFixed(1)+'%'+(repeat?'<small>'+r.large_error_count+'/'+r.attempted+' large</small>':'')+'</td>';
 }).join('')+'</tr>').join('')+'</tbody></table></article>').join('');
 const selected=r=>r.family===family&&r.overlap===el('overlap').value&&+r.affected===+el('affected').value;
 const rr=old.examples.filter(selected),first=rr[0];
 el('truth').textContent='True interval: '+fmt(first.true_lower)+'–'+fmt(first.true_upper)+'. Overall affected fraction: '+(+first.affected*100)+'%; simulated services: '+Math.round((+first.affected-.05)*100)+'% / '+Math.round((+first.affected+.05)*100)+'%.';
 el('scenario-rows').innerHTML=methods.map(m=>data.summary.find(r=>selected(r)&&r.method===m)).map(r=>'<tr><th scope="row">'+labels[r.method]+'</th><td>'+fmt(r.mean_error)+'%</td><td>'+fmt(r.median_error)+'%</td><td>'+fmt(r.p90_error)+'%</td><td>'+fmt(r.mean_lower_error)+' / '+fmt(r.mean_upper_error)+'</td><td>'+r.large_error_count+'/'+r.attempted+'</td><td>'+r.failed+'</td></tr>').join('');
 el('paired-rows').innerHTML=data.paired.filter(selected).map(r=>'<tr><th scope="row">'+labels[r.method_a]+'<br>vs '+labels[r.method_b]+'</th><td>'+fmt(r.mean_difference)+'</td><td>'+fmt(r.mcse_difference)+'</td><td>'+r.wins_a+'/'+r.finite_pairs+'</td></tr>').join('');
 el('rows').innerHTML=rr.map(r=>'<tr><th scope="row">'+labels[r.method]+'</th><td>'+fmt(r.lower)+'<small>'+fmt(r.lower_ci_low)+'–'+fmt(r.lower_ci_high)+'</small></td><td>'+fmt(r.upper)+'<small>'+fmt(r.upper_ci_low)+'–'+fmt(r.upper_ci_high)+'</small></td><td>'+fmt(r.width_error)+'%</td><td>'+r.bootstrap_successful+'/'+r.bootstrap_requested+'</td></tr>').join('');
 const s=data.summary.find(r=>selected(r)&&r.method==='anchored40');el('repeat-detail').textContent='For this scenario, the new anchoredRI 20-dataset mean error is '+fmt(s.mean_error)+'% (Monte Carlo SE '+fmt(s.mcse_mean)+' percentage points); '+s.large_error_count+'/20 exceeded 10%.';
 el('status').textContent='45 scenarios complete. All four approaches evaluated on 900 matched datasets. Earlier bootstrap examples shown separately.';
 }
 ['family','evidence','overlap','affected'].forEach(id=>el(id).addEventListener('change',render));render();
}).catch(e=>{el('maps').textContent=e.message;el('status').textContent='Results could not be displayed. Download the CSV files below.';});

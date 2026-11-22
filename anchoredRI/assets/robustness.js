'use strict';
const labels={anchored40:'anchoredRI + 40',refineR:'refineR',direct40:'40 screened alone',direct120:'120 screened alone',direct400:'400 screened alone'};
const overlaps=['separated','moderate','heavy'];const names=['Separated','Moderate','Heavy'];const fractions=[.1,.3,.5,.7,.9];
const el=id=>document.getElementById(id);const fmt=x=>x===''||x==null?'Unavailable':Number(x).toFixed(2);
fetch('assets/stress-grid/results.json').then(r=>{if(!r.ok)throw Error('Could not load results');return r.json()}).then(data=>{
 if(data.examples.length!==225||data.summary.length!==180)throw Error('Incomplete result files');
 function render(){const family=el('family').value,repeat=el('evidence').value==='repeat';
 el('explanation').textContent=repeat?'Mean error across 20 fresh datasets per scenario. The small number is the count with error >10%. refineR was not rerun on these 20 replicates; use the matched example view for that comparison.':'One matched dataset per scenario, with 200 bootstrap resamples per method. These are individual examples, not 20-dataset averages. anchoredRI also receives healthy anchors and service labels.';
 const methods=repeat?['anchored40','direct40','direct120','direct400']:['anchored40','refineR','direct40','direct120','direct400'];
 el('maps').innerHTML=methods.map(m=>'<article><h3>'+labels[m]+'</h3><table><thead><tr><th>Overlap</th>'+fractions.map(p=>'<th scope="col">'+p*100+'%</th>').join('')+'</tr></thead><tbody>'+overlaps.map((o,i)=>'<tr><th scope="row">'+names[i]+'</th>'+fractions.map(p=>{const r=(repeat?data.summary:data.examples).find(r=>r.family===family&&r.overlap===o&&+r.affected===p&&r.method===m);const e=+(repeat?r.mean_error:r.width_error);const color=e<=5?'#d5eadc':e<=10?'#f1ebba':e<=25?'#f4cfb3':'#efb4b0';return '<td style="background:'+color+'">'+e.toFixed(1)+'%'+(repeat?'<small>'+r.large+'/'+r.attempted+' large</small>':'')+'</td>';}).join('')+'</tr>').join('')+'</tbody></table></article>').join('');
 const rr=data.examples.filter(r=>r.family===family&&r.overlap===el('overlap').value&&+r.affected===+el('affected').value);const first=rr[0];
 el('truth').textContent='True interval: '+fmt(first.true_lower)+'–'+fmt(first.true_upper)+'. Overall affected fraction: '+(+first.affected*100)+'%; simulated services: '+Math.round((+first.affected-.05)*100)+'% / '+Math.round((+first.affected+.05)*100)+'%.';
 el('rows').innerHTML=rr.map(r=>'<tr><th scope="row">'+labels[r.method]+'</th><td>'+fmt(r.lower)+'<small>'+fmt(r.lower_ci_low)+'–'+fmt(r.lower_ci_high)+'</small></td><td>'+fmt(r.upper)+'<small>'+fmt(r.upper_ci_low)+'–'+fmt(r.upper_ci_high)+'</small></td><td>'+fmt(r.width_error)+'%</td><td>'+r.bootstrap_successful+'/'+r.bootstrap_requested+'</td></tr>').join('');
 const s=data.summary.find(r=>r.case===+first.case&&r.method==='anchored40');el('repeat-detail').textContent='For this same scenario, anchoredRI’s 20-dataset mean error was '+fmt(s.mean_error)+'% (Monte Carlo standard error '+fmt(s.mcse)+' percentage points); '+s.large+'/20 exceeded 10%, and '+s.failed+' fits failed. The standard error describes simulation precision, not uncertainty in a patient’s result.';
 el('status').textContent='45 scenarios complete. 900 repeated anchoredRI estimates; 45 matched refineR examples.';
 }
 ['family','evidence','overlap','affected'].forEach(id=>el(id).addEventListener('change',render));render();
}).catch(e=>{el('maps').textContent=e.message;el('status').textContent='Results could not be displayed. Download the CSV files below.';});

const el = id => document.getElementById(id);
const defs = {
  normal: {note:'A normal reference population with disease shifted upward. Expected affected fraction: 10% overall (5% and 15% by service).', model:'BoxCox'},
  overlap: {note:'Disease overlaps the reference distribution. Expected affected fraction: 25% overall (20% and 30% by service). The clean-region assumption is more challenging here.', model:'BoxCox'},
  lognormal: {note:'A lognormal reference population with an upper disease component. Expected affected fraction: 10%. The updated anchoredRI chooses its transformation from the data.', model:'BoxCox'},
  shifted: {note:'A lognormal reference population shifted 10 units from zero. Expected affected fraction: 10%. The updated anchoredRI estimates the offset; it is not supplied the true shift.', model:'modBoxCox'}
};
const labels = {refineR:'refineR', anchored40:'anchoredRI 0.1', anchored_auto40:'anchoredRI 0.2', direct120:'Screened alone · 120', direct400:'Screened alone · 400'};
const colors = {refineR:'#426da3', anchored40:'#7c8982', anchored_auto40:'#006b60', direct120:'#9a692c', direct400:'#81558e'};
const order = Object.keys(labels);
let rows = [], models = [];

function csvParse(text) {
  const records = []; let record = [], field = '', quoted = false;
  for (let i=0; i<text.length; i++) {
    const c=text[i];
    if (c==='"') { if (quoted && text[i+1]==='"') {field+='"';i++;} else quoted=!quoted; }
    else if (c===',' && !quoted) {record.push(field);field='';}
    else if ((c==='\n'||c==='\r')&&!quoted) {
      if (c==='\r' && text[i+1]==='\n') i++;
      record.push(field); if(record.some(v=>v!=='')) records.push(record);record=[];field='';
    } else field+=c;
  }
  if(field!==''||record.length){record.push(field);records.push(record);}
  if(quoted) throw new Error('Incomplete CSV field.');
  const columns=records.shift();
  return records.map(r=>Object.fromEntries(columns.map((c,i)=>[c,r[i]??''])));
}
function fmt(v,d=2){return Number.isFinite(Number(v)) && v!=='' ? Number(v).toFixed(d) : 'Unavailable';}
function add(parent,tag,text,cls){const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(cls)n.className=cls;parent.append(n);return n;}
function svgNode(tag, attrs={}, text){const n=document.createElementNS('http://www.w3.org/2000/svg',tag);for(const [k,v] of Object.entries(attrs))n.setAttribute(k,v);if(text!==undefined)n.textContent=text;return n;}
function drawChart(selected) {
  const valid=selected.filter(r=>r.status==='ok');
  if(!valid.length){el('guided-chart').textContent='No successful fits for this dataset.';return;}
  const t=[Number(valid[0].true_lower),Number(valid[0].true_upper)];
  const values=[...t,...valid.flatMap(r=>['lower','upper','lower_ci_low','lower_ci_high','upper_ci_low','upper_ci_high'].map(k=>r[k]===''?NaN:Number(r[k])))].filter(Number.isFinite);
  const span=Math.max(...values)-Math.min(...values), step=span>80?20:span>35?10:span>15?5:2;
  const lo=Math.floor((Math.min(...values)-step*.4)/step)*step,hi=Math.ceil((Math.max(...values)+step*.4)/step)*step;
  const X=v=>190+(v-lo)/(hi-lo)*560;
  const bottom=selected.length*62+8;
  const svg=svgNode('svg',{viewBox:`0 0 790 ${bottom+49}`,role:'img','aria-label':'Reference intervals and 95 percent confidence intervals'});
  svg.append(svgNode('rect',{x:X(t[0]),y:20,width:X(t[1])-X(t[0]),height:bottom-20,fill:'#e5eee8'}));
  for(let v=lo;v<=hi+1e-8;v+=step){svg.append(svgNode('line',{x1:X(v),x2:X(v),y1:20,y2:bottom,stroke:'#dce3de','stroke-width':1}));svg.append(svgNode('text',{x:X(v),y:bottom+23,'text-anchor':'middle','font-size':14,fill:'#4c5e55'},String(v)));}
  selected.forEach((r,i)=>{
    const y=45+i*62,color=colors[r.method];
    svg.append(svgNode('text',{x:8,y:y+5,'font-size':15,fill:'#233c31'},labels[r.method]));
    if(r.status!=='ok'){svg.append(svgNode('text',{x:200,y:y+5,fill:'#a32626','font-size':15},'Fit failed'));return;}
    svg.append(svgNode('line',{x1:X(Number(r.lower)),x2:X(Number(r.upper)),y1:y,y2:y,stroke:color,'stroke-width':5}));
    for(const side of ['lower','upper']){
      const a=Number(r[side+'_ci_low']),b=Number(r[side+'_ci_high']);
      if(r[side+'_ci_low']!==''&&r[side+'_ci_high']!==''){
        svg.append(svgNode('line',{x1:X(a),x2:X(b),y1:y,y2:y,stroke:color,'stroke-width':1.5}));
        for(const v of [a,b])svg.append(svgNode('line',{x1:X(v),x2:X(v),y1:y-8,y2:y+8,stroke:color,'stroke-width':1.5}));
      }
      svg.append(svgNode('circle',{cx:X(Number(r[side])),cy:y,r:4,fill:color}));
    }
  });
  el('guided-chart').replaceChildren(svg);
}
function render(){
  const profile=el('profile').value,rep=Number(el('dataset').value);
  el('profile-note').textContent=defs[profile].note;
  const selected=order.map(m=>rows.find(r=>r.profile===profile&&Number(r.replicate)===rep&&r.method===m)).filter(Boolean);
  if(selected.length!==order.length){el('guided-chart').textContent='Results for this selection are unavailable.';return;}
  el('truth-note').textContent='Known reference interval: '+fmt(selected[0].true_lower)+' to '+fmt(selected[0].true_upper)+' units';
  drawChart(selected);el('result-rows').replaceChildren();
  for(const r of selected){
    const tr=add(el('result-rows'),'tr');const name=add(tr,'td',labels[r.method]);
    add(name,'small',r.method==='refineR'?'5,000 routine':r.method.startsWith('anchored')?'5,000 routine + 40 screened':r.n_screened+' screened');
    if(r.status!=='ok'){const td=add(tr,'td','Fit failed: '+r.error,'failed');td.colSpan=4;continue;}
    for(const side of ['lower','upper']){const td=add(tr,'td',fmt(r[side]));add(td,'small',fmt(r[side+'_ci_low'])+' to '+fmt(r[side+'_ci_high']));}
    const errorCell=add(tr,'td',fmt(r.worst_error));
    add(errorCell,'small','MPE '+fmt(r.mpe)+'%');
    add(tr,'td',fmt(Number(r.reference_coverage)*100,1)+'%');
  }
  const successes=selected.map(r=>labels[r.method]+': '+r.bootstrap_successful+'/'+r.bootstrap_requested+' bootstrap fits');
  const bounds=selected.filter(r=>r.at_bound==='TRUE').map(r=>labels[r.method]);
  el('fit-status').textContent=successes.join('; ')+'.'+(bounds.length?' Numerical boundary reached: '+bounds.join(', ')+'.':'');
  const model=models.find(m=>m.profile===profile);
  const familyNames={identity:'normal',log:'lognormal',shifted_log:'shifted lognormal'};
  el('anchor-note').textContent=model?'anchoredRI 0.2 selected '+familyNames[model.selected_transform]+' with '+model.components+' components.'+(model.selected_transform==='shifted_log'?' Estimated offset: '+fmt(model.estimated_shift)+' units.':'')+' Bootstrap choices: '+model.bootstrap_identity+' normal, '+model.bootstrap_log+' lognormal, '+model.bootstrap_shifted_log+' shifted lognormal. '+model.bootstrap_boundary+' selected bootstrap fits reached a numerical bound.':'Model details unavailable.';
  el('anchor-selection-link').href='assets/offset/'+profile+'-selection.csv';
  el('anchor-bootstrap-link').href='assets/offset/'+profile+'-bootstrap.csv';
  const ref=selected.find(r=>r.method==='refineR');
  el('fit-note').textContent='Released refineR '+defs[profile].model+' model, medianBS point estimate, 200 bootstrap resamples. Estimated non-pathological fraction: '+fmt(Number(ref.np_fraction)*100,1)+'%. This estimate is not proof that the input meets the method’s assumptions.';
  const path='assets/guided/'+profile+'-'+String(rep).padStart(2,'0')+'.svg';
  el('fit-image').src=path;el('fit-link').href=path;el('default-fit').hidden=profile!=='shifted';
  el('default-fit-link').href='assets/guided/shifted-default-'+String(rep).padStart(2,'0')+'.svg';
}
async function load(){
 try{
   const responses=await Promise.all(['fits.csv','models.csv'].map(p=>fetch('assets/offset/'+p+'?v=offset2')));
   if(responses.some(r=>!r.ok))throw new Error('Unavailable results');
   [rows,models]=await Promise.all(responses.map(async r=>csvParse(await r.text())));
   if(!rows.length||rows.some(r=>!labels[r.method]||!defs[r.profile]||!['ok','failed'].includes(r.status)))throw new Error('Invalid results');
   const reps=[...new Set(rows.map(r=>Number(r.replicate)))].sort((a,b)=>a-b);
   el('dataset').replaceChildren();for(const rep of reps){const o=add(el('dataset'),'option','Dataset '+rep);o.value=String(rep);}
   render();
 }catch(e){el('guided-chart').textContent='The fitted results could not be loaded. Reload the page or download the CSV.';}
}
async function loadValidation(){
 try{
  const response=await fetch('assets/offset/validation-summary.csv?v=offset2');
  if(!response.ok)throw new Error('Unavailable validation');
  const summary=csvParse(await response.text());
  for(const r of summary){
    const tr=add(el('validation-rows'),'tr');add(tr,'td',r.title);
    add(tr,'td',fmt(r.original_mpe)+'%');add(tr,'td',fmt(r.auto_mpe)+'%');
    add(tr,'td',r.auto_better+' / '+r.ties+' / '+r.auto_worse);
    add(tr,'td',r.paired_successful+'/'+r.attempted);
  }
 }catch(e){el('validation-status').textContent='The validation results could not be loaded. Reload or download the CSV.';}
}
const requestedProfile=new URLSearchParams(window.location.search).get('profile');
if(Object.prototype.hasOwnProperty.call(defs,requestedProfile))el('profile').value=requestedProfile;
el('profile').addEventListener('change',render);el('dataset').addEventListener('change',render);load();loadValidation();

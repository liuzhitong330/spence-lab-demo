(function () {
  'use strict';
  const DATA = window.SPENCE_DATA;
  if (!DATA || !DATA.traits || !DATA.traits.length) {
    document.getElementById('readout').textContent = 'The research data did not load. Please reload this page; no synthetic values are substituted.';
    return;
  }
  const $ = id => document.getElementById(id);
  const PAGE = 6, Z = 1.959963984540054;
  const state = { trait: DATA.meta.defaultTrait || DATA.traits[0].id, locus: DATA.meta.defaultLocus, gene: DATA.meta.defaultGene, page:0, delta:.1, mode:'pointwise', count:20 };
  const esc = v => String(v ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const num = v => Number.isFinite(v) ? v.toFixed(3) : 'unavailable';
  const trait = () => DATA.traits.find(t => t.id === state.trait) || DATA.traits[0];
  const locus = () => trait().loci.find(l => l.id === state.locus) || trait().loci[0];
  const zValue = () => state.mode === 'family' ? trait().zFamily : Z;
  const region = l => `chr${l.chr}: ${(l.start/1e6).toFixed(2)}–${(l.end/1e6).toFixed(2)} Mb`;
  function classify(g, delta=state.delta, z=zValue()) {
    if (!Number.isFinite(g.beta) || !Number.isFinite(g.se) || g.se <= 0) return {kind:'unavailable',lo:null,hi:null};
    const lo=g.beta-z*g.se, hi=g.beta+z*g.se;
    return {lo,hi,kind:lo>=-delta && hi<=delta ? 'small' : lo>delta || hi<-delta ? 'large' : 'unresolved'};
  }
  const labels = {small:'within the selected range',large:'entirely beyond the selected range',unresolved:'unresolved at this threshold',unavailable:'no usable estimate'};
  const colors = {small:'#1f7a8c',large:'#a56e31',unresolved:'#666',unavailable:'#999'};
  function counts(l) { const c={small:0,large:0,unresolved:0,unavailable:0}; l.genes.forEach(g=>c[classify(g).kind]++); return c; }
  function option(value,text) { const el=document.createElement('option'); el.value=value; el.textContent=text; return el; }
  function choices(reset=false) {
    const t=trait();
    $('locus-select').replaceChildren(...t.loci.map(l=>option(l.id,`#${l.rank} · ${region(l)}`)));
    if(reset || !t.loci.some(l=>l.id===state.locus)) { state.locus=t.loci[0].id; state.gene=null; state.page=0; }
    const l=locus(); state.locus=l.id;
    $('locus-select').value=l.id;
    $('gene-select').replaceChildren(...l.genes.map(g=>option(g.id,`${g.symbol || g.id}${classify(g).kind==='unavailable'?' · no estimate':''}`)));
    if (!l.genes.some(g=>g.id===state.gene)) state.gene=l.genes[0]?.id;
    $('gene-select').value=state.gene || '';
    $('gene-select').disabled=!l.genes.length;
  }
  function draw() {
    const l=locus(), z=zValue(), rows=l.genes.slice(state.page*PAGE,(state.page+1)*PAGE);
    const intervals=rows.map(g=>classify(g));
    const extent=Math.max(state.delta*1.35,...intervals.flatMap(b=>b.lo===null?[]:[Math.abs(b.lo)*1.1,Math.abs(b.hi)*1.1]));
    const x=v=>108+(v+extent)/(extent*2)*472;
    let svg=`<title>Loss-of-function burden effect intervals at ${esc(region(l))}</title><desc>Points are gene effect estimates. Horizontal lines are uncertainty intervals; the shaded band marks plus or minus ${state.delta.toFixed(2)}. Use the gene selector for an accessible textual result.</desc><rect x="${x(-state.delta)}" y="16" width="${x(state.delta)-x(-state.delta)}" height="234" fill="#e5eeee"/><line x1="${x(0)}" y1="16" x2="${x(0)}" y2="250" stroke="#aab6b6" stroke-dasharray="3 3"/>`;
    rows.forEach((g,i)=>{
      const b=intervals[i], y=38+i*34, selected=g.id===state.gene;
      svg+=`<g class="gene-row" role="button" tabindex="0" data-gene="${esc(g.id)}" aria-label="Inspect ${esc(g.symbol || g.id)}"><rect x="0" y="${y-16}" width="600" height="30" fill="${selected?'#eaf0f0':'transparent'}"/><text x="99" y="${y+4}" text-anchor="end" font-family="system-ui,sans-serif" font-size="12" fill="#333">${esc(g.symbol || g.id)}</text>`;
      if(b.lo!==null) svg+=`<line x1="${x(b.lo)}" x2="${x(b.hi)}" y1="${y}" y2="${y}" stroke="${colors[b.kind]}" stroke-width="2"/><line x1="${x(b.lo)}" x2="${x(b.lo)}" y1="${y-4}" y2="${y+4}" stroke="${colors[b.kind]}"/><line x1="${x(b.hi)}" x2="${x(b.hi)}" y1="${y-4}" y2="${y+4}" stroke="${colors[b.kind]}"/><circle cx="${x(g.beta)}" cy="${y}" r="4" fill="${colors[b.kind]}"/>`;
      else svg+=`<text x="116" y="${y+4}" font-family="system-ui,sans-serif" font-size="12" fill="#888">No usable estimate</text>`;
      svg+='</g>';
    });
    if(!rows.length) svg+='<text x="300" y="120" text-anchor="middle" font-size="15">No annotated genes in this region</text>';
    svg+='<line x1="108" x2="580" y1="265" y2="265" stroke="#aaa"/>';
    [-extent,-extent/2,0,extent/2,extent].forEach(v=>{svg+=`<line x1="${x(v)}" x2="${x(v)}" y1="265" y2="270" stroke="#aaa"/><text x="${x(v)}" y="285" text-anchor="middle" font-family="system-ui,sans-serif" font-size="11" fill="#666">${v.toFixed(2)}</text>`;});
    svg+='<text x="344" y="307" text-anchor="middle" font-family="system-ui,sans-serif" font-size="11" fill="#666">Burden effect, inverse-rank-normalized trait units</text>';
    $('hero-viz').innerHTML=svg;
    $('hero-viz').querySelectorAll('[data-gene]').forEach(el=>{
      const select=()=>{state.gene=el.dataset.gene; $('gene-select').value=state.gene; render();};
      el.addEventListener('click',select); el.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();select();}});
    });
    $('previous-page').disabled=state.page===0;
    $('next-page').disabled=(state.page+1)*PAGE>=l.genes.length;
    $('page-label').textContent=l.genes.length?`${state.page*PAGE+1}–${Math.min((state.page+1)*PAGE,l.genes.length)} of ${l.genes.length}`:'0 genes';
    const g=l.genes.find(g=>g.id===state.gene);
    if(g){
      const b=classify(g), detail=b.kind==='unavailable' ? 'This row is retained for coverage transparency, not treated as a zero effect.' : `β ${num(g.beta)}; interval [${num(b.lo)}, ${num(b.hi)}]. ${b.kind==='small'?'This bounds this burden effect within the chosen range; it does not establish biological unimportance.':b.kind==='large'?'The interval lies on one side of the chosen range. This does not identify the gene as the causal mediator of the GWAS signal.':'The interval crosses a range boundary. These data do not settle whether the effect is inside or outside the chosen range.'}`;
      const pNote=Number.isFinite(g.burdenP) ? ` Conventional burden p ${g.burdenP===0?'below floating-point precision':g.burdenP.toExponential(2)}${g.burdenP < .05/trait().nAssayed?' (below':' (above'} the within-trait Bonferroni threshold). Significance and the chosen effect range answer different questions.` : '';
      $('readout').innerHTML=`<strong>${esc(g.symbol || g.id)}: ${labels[b.kind]}.</strong> ${esc(detail)}<span class="p-note">${esc(pNote)}</span>`;
    } else $('readout').textContent='No gene-level estimates can be displayed in this interval.';
    $('figure-caption').textContent=`All ${l.genes.length} annotated genes are retained; six are shown at a time. Teal: within range; ochre: beyond; grey: unresolved. The horizontal scale follows the visible intervals, so compare axis labels between pages. ${state.mode==='family'?'Bonferroni-adjusted':'Pointwise 95%'} intervals; individual-level observations are not available here.`;
  }
  function renderTable() {
    const ls=trait().loci.slice(0,state.count), totals={small:0,large:0,unresolved:0,unavailable:0};
    $('region-table').innerHTML=ls.map(l=>{const c=counts(l); Object.keys(totals).forEach(k=>totals[k]+=c[k]); return `<tr class="${l.id===state.locus?'current':''}"><td><button type="button" data-locus="${esc(l.id)}">#${l.rank} · ${esc(region(l))}</button></td><td>${c.small}</td><td>${c.large}</td><td>${c.unresolved}</td><td>${c.unavailable}</td></tr>`;}).join('');
    $('region-table').querySelectorAll('[data-locus]').forEach(b=>b.addEventListener('click',()=>{state.locus=b.dataset.locus;state.gene=null;state.page=0;choices();render();$('a1').scrollIntoView();}));
    $('cohort-note').textContent=`${ls.length} regions; ${totals.small+totals.large+totals.unresolved} gene–region entries with usable estimates and ${totals.unavailable} without. These entries are not independent participants or independent tests. Threshold ±${state.delta.toFixed(2)}; ${state.mode==='family'?'Bonferroni-adjusted':'pointwise'} intervals.`;
  }
  function render() {
    const t=trait(),l=locus(),c=counts(l);
    $('threshold-label').textContent=`±${state.delta.toFixed(2)}`;
    $('interval-note').textContent=`z = ${zValue().toFixed(4)}. Family correction uses all ${t.nAssayed.toLocaleString()} valid assayed genes for this trait, not only the displayed regions. Pointwise intervals do not control simultaneous coverage.`;
    $('scope-note').textContent=`${region(l)} · lead variant ${l.snp} · ${l.genes.length} annotated genes. Region boundaries come from the authors’ deposited GWAS output.`;
    $('metric-one').textContent=c.small; $('metric-two').textContent=c.large; $('metric-three').textContent=c.unresolved;
    $('missing-note').textContent=`Counts cover the entire selected region. ${c.unavailable} additional annotated ${c.unavailable===1?'gene has':'genes have'} no usable estimate and ${c.unavailable===1?'is':'are'} not counted as a small effect.`;
    draw(); renderTable();
  }
  function download(name,headers,rows) {
    const field=v=>'"'+String(v??'').replace(/"/g,'""')+'"';
    const csv=[headers,...rows].map(r=>r.map(field).join(',')).join('\r\n')+'\r\n';
    const url=URL.createObjectURL(new Blob([csv],{type:'text/csv;charset=utf-8'}));
    const a=document.createElement('a');a.href=url;a.download=name;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);
  }
  $('trait-select').replaceChildren(...DATA.traits.map(t=>option(t.id,t.name))); $('trait-select').value=trait().id;
  document.querySelectorAll('[data-example]').forEach(button=>button.addEventListener('click',()=>{
    const symbol=button.dataset.example, t=DATA.traits.find(t=>t.id===(symbol==='HHIP'?'50_irnt':'30780_irnt'));
    const l=t.loci.find(l=>l.genes.some(g=>g.symbol===symbol)), g=l.genes.find(g=>g.symbol===symbol);
    state.trait=t.id;state.locus=l.id;state.gene=g.id;state.page=Math.floor(l.genes.indexOf(g)/PAGE);
    $('trait-select').value=t.id;choices();render();
  }));
  $('trait-select').addEventListener('change',e=>{state.trait=e.target.value;choices(true);render();});
  $('locus-select').addEventListener('change',e=>{state.locus=e.target.value;state.gene=null;state.page=0;choices();render();});
  $('gene-select').addEventListener('change',e=>{state.gene=e.target.value;state.page=Math.floor(locus().genes.findIndex(g=>g.id===state.gene)/PAGE);render();});
  $('effect-threshold').addEventListener('input',e=>{state.delta=Number(e.target.value);render();});
  $('interval-mode').addEventListener('change',e=>{state.mode=e.target.value;render();});
  $('previous-page').addEventListener('click',()=>{state.page--;state.gene=locus().genes[state.page*PAGE].id;$('gene-select').value=state.gene;render();});
  $('next-page').addEventListener('click',()=>{state.page++;state.gene=locus().genes[state.page*PAGE].id;$('gene-select').value=state.gene;render();});
  $('region-count').addEventListener('change',e=>{state.count=Number(e.target.value);renderTable();});
  $('download-locus').addEventListener('click',()=>download(`spence-${trait().id}-region-${locus().rank}.csv`,['trait_id','locus_id','gene_id','symbol','beta_irnt','se_irnt','interval_low','interval_high','interval_mode','z','delta_irnt','classification','gene_start','gene_end','burden_p_normal_approx','availability'],locus().genes.map(g=>{const b=classify(g);return[trait().id,locus().id,g.id,g.symbol,g.beta,g.se,b.lo,b.hi,state.mode,zValue(),state.delta,b.kind,g.start,g.end,g.burdenP,g.availability];})));
  $('download-summary').addEventListener('click',()=>download(`spence-${trait().id}-regions.csv`,['trait_id','locus_id','rank','chr','start','end','within','beyond','unresolved','no_estimate','delta_irnt','interval_mode','z'],trait().loci.slice(0,state.count).map(l=>{const c=counts(l);return[trait().id,l.id,l.rank,l.chr,l.start,l.end,c.small,c.large,c.unresolved,c.unavailable,state.delta,state.mode,zValue()];})));
  $('provenance-note').textContent=DATA.meta.coverageNote || 'The data download and project README document source hashes, coordinate mapping, trait coverage and exclusions.';
  choices();
  const initialIndex=locus().genes.findIndex(g=>g.id===state.gene);if(initialIndex>=0)state.page=Math.floor(initialIndex/PAGE);
  render();
})();

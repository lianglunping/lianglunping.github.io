'use strict';
const $=id=>document.getElementById(id),esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const numKeys=Object.keys(RPS.DEFAULT).filter(k=>typeof RPS.DEFAULT[k]==='number'&&$(k));
let mode='sanger',inputMode='single',refFile=null,entries=[],refSha=null,refRevision=0,demoInfo=null,results=[],runManifest=null,busy=false,activeWorker=null,activeReject=null,cancelRequested=false;
let gffFile=null,gffIndex=null,gffSha=null,gffAudit=null,gffSelectedGene=null;
let indexAudit=null;let workerURL=null,nativeInfo=null,nativeAbort=null;let bgRows=[],maskText='',repeatRows=[],repeatFileMeta=null;
function workerCode(){return '('+RPSEngineFactory.toString()+')(self);\n('+RPSIOFactory.toString()+')(self);'+`\nself.onmessage=async ev=>{const {task,payload}=ev.data;try{let result;const progress=p=>self.postMessage({type:'progress',value:p});if(task==='design')result=RPS.design(payload.context,payload.mode,payload.params);else if(task==='index')result=await RPSIO.buildIndex(payload.file,progress);else if(task==='hash')result=await RPSIO.hashFile(payload.file,progress);else if(task==='approx')result=await RPSIO.scanApprox(payload.file,payload.queries,payload.maxMismatches,progress);else if(task==='scan')result=await RPSIO.scanExact(payload.file,payload.queries,progress);else throw Error('Unknown task');self.postMessage({type:'result',result});}catch(e){self.postMessage({type:'error',error:e.message});}};`;}
function runWorker(task,payload,progress=()=>{}){if(!workerURL)workerURL=URL.createObjectURL(new Blob([workerCode()],{type:'text/javascript'}));return new Promise((resolve,reject)=>{const w=new Worker(workerURL);activeWorker=w;activeReject=reject;const stop=()=>{w.terminate();if(activeWorker===w){activeWorker=null;activeReject=null;}};w.onmessage=ev=>{const m=ev.data;if(m.type==='progress')progress(m.value);if(m.type==='result'){stop();resolve(m.result);}if(m.type==='error'){stop();reject(Error(m.error));}};w.onerror=ev=>{stop();reject(Error(ev.message||'Worker 执行失败。请使用支持本地文件与 Web Worker 的桌面浏览器。'));};w.postMessage({task,payload});});}
function toast(s){$('toast').textContent=s;$('toast').style.display='block';clearTimeout(toast.timer);toast.timer=setTimeout(()=>$('toast').style.display='none',4200);}
function error(s){$('errorBox').textContent=s;$('errorBox').classList.remove('hidden');}
function clearError(){$('errorBox').classList.add('hidden');}
function showProgress(label,p=0){$('progressText').textContent=label;$('progressBar').style.width=Math.max(0,Math.min(100,p*100))+'%';}
function setBusy(value){busy=value;document.querySelectorAll('input,select,textarea,.needs-idle').forEach(x=>x.disabled=value);document.querySelectorAll('[data-scan]').forEach(x=>x.disabled=value||!refFile||results[Number(x.dataset.row)]?.refRevision!==refRevision||results[Number(x.dataset.row)]?.context.source!=='LOCAL_REFERENCE');$('progress').classList.toggle('visible',value);if(!value)activeWorker=null;}
function dirty(){if(results.length)$('staleNote').classList.remove('hidden');}
function save(name,text,type='text/plain;charset=utf-8'){const url=URL.createObjectURL(new Blob([text],{type})),a=document.createElement('a');a.href=url;a.download=name;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),3000);}
function shaText(t){return new RPSIO.SHA256().update(new TextEncoder().encode(t)).digest();}
function refMeta(){return {assemblyClaim:'Os-Nipponbare-Reference-IRGSP-1.0',annotationClaim:'MSU/RGAP7',sourceURL:'https://rice.uga.edu/download_osa1r7.shtml',identityStatus:'USER_SUPPLIED_NOT_VERIFIED_AGAINST_PUBLISHED_DIGEST',filename:refFile?.name??null,size:refFile?.size??null,lastModified:refFile?.lastModified??null,sha256:refSha,indexAudit,index:entries,indexSha256:entries.length?shaText(faiText()):null,gff:gffFile?{filename:gffFile.name,size:gffFile.size,lastModified:gffFile.lastModified,sha256:gffSha,indexedGenes:gffIndex?.stats?.genes??null,audit:gffAudit}:null};}
function faiText(){return entries.map(x=>[x.name,x.length,x.offset,x.lineBases,x.lineWidth].join('\t')).join('\n')+'\n';}
function renderRef(){$('refStatus').textContent=refFile?`${refFile.name} · ${entries.length?entries.length+' 个 contig，索引已载入':'等待 FAI 或构建索引'} · ${refSha?'已计算文件 SHA-256':'尚未计算文件 SHA-256'}`:'尚未选择本地 FASTA · 也可先运行合成序列示例';$('refMeta').textContent=`FASTA: ${refFile?.name||'未选择'}${refFile?' ('+refFile.size.toLocaleString()+' bytes)':''}\nFAI: ${entries.length?entries.length+' contigs':'未载入'}\nSHA-256: ${refSha||'未计算'}\n身份状态: 用户声明；未对照官方校验值鉴定`; $('saveFai').disabled=!entries.length;}
function setInput(v){inputMode=v;document.querySelectorAll('[data-input]').forEach(b=>{const yes=b.dataset.input===v;b.classList.toggle('active',yes);b.setAttribute('aria-selected',yes);});for(const m of ['single','batch','sequence'])$('input-'+m).classList.toggle('hidden',m!==v);dirty();}
function resetParams(){for(const k of numKeys)$(k).value=RPS.DEFAULT[k];if(mode==='kasp'){$('minProduct').value=70;$('maxProduct').value=180;$('rescueMaxProduct').value=400;$('rescueStep').value=50;}else{$('rescueMaxProduct').value=1600;$('rescueStep').value=200;}$('autoRescue').checked=true;$('rescueTrigger').value='repeat';$('direction').value='both';$('famAllele').value='REF';dirty();}
function setMode(v){mode=v;document.querySelectorAll('[data-mode]').forEach(x=>x.classList.toggle('active',x.dataset.mode===v));$('bufferField').classList.toggle('hidden',v==='kasp');$('kaspOptions').classList.toggle('hidden',v!=='kasp');$('modeNote').textContent=v==='sanger'?'设计位于变异两侧的 PCR 引物。目标与引物 3′ 端之间预留读序距离；产物范围同时检查 REF 与 ALT。':'分别构建 REF / ALT 并比较两个方向。SNP 与短纯 InDel 输出三引物候选；InDel 始终要求人工复核。尾序不代表寡核苷酸需直接加荧光修饰。';resetParams();}
function params(){const p={...RPS.DEFAULT};for(const k of numKeys)p[k]=Number($(k).value);p.direction=$('direction').value;p.famAllele=$('famAllele').value;p.regionGcHigh=Number($('gcRegionHigh').value);p.regionGcLow=Number($('gcRegionLow').value);p.sangerMaxReadSpan=Number($('sangerMaxReadSpan').value);if(!(p.regionGcLow>=0&&p.regionGcHigh<=100&&p.regionGcLow<p.regionGcHigh))throw Error('区域 GC 提示线必须满足 0 ≤ low < high ≤ 100。');if(!Number.isSafeInteger(p.sangerMaxReadSpan)||p.sangerMaxReadSpan<100||p.sangerMaxReadSpan>2000)throw Error('读长规划值须为 100–2000 bp 整数；不是实测读长。');p.thermoTemp=Number($('thermoTemp').value);if(!Number.isFinite(p.thermoTemp)||p.thermoTemp<0||p.thermoTemp>100)throw Error('热力学温度须在 0–100°C。');return RPS.validateParams(p,mode);}
function loadDemo(kind){demoInfo=RPS.demo(kind);setInput('sequence');const d=demoInfo;$('sequenceText').value=d.seq.slice(0,d.target0)+`[${d.ref}/${d.alt}]`+d.seq.slice(d.target0+d.ref.length);$('demoNotice').classList.remove('hidden');dirty();toast('已载入合成序列；点击“生成候选”进行真实计算。');}
function resolveNameOnly(name){return RPS.resolveContig(entries,name).name;}
function gffMappingAudit(){if(!gffIndex||!entries.length)return {status:'NOT_AUDITED',mapped:0,unmapped:[],ambiguous:[],seqids:gffIndex?[...gffIndex.seqids].length:0};const unmapped=[],mapped=[];for(const seqid of gffIndex.seqids){try{const e=RPS.resolveContig(entries,seqid);if((gffIndex.maxEndBySeq.get(seqid)||0)>e.length)throw Error('GFF3 coordinate out of range');mapped.push([seqid,e.name]);}catch(e){unmapped.push(seqid);}}return {status:unmapped.length?'PARTIAL':'MAPPED',mapped:mapped.length,unmapped:unmapped.slice(0,50),unmappedCount:unmapped.length,seqids:gffIndex.seqids.size,examples:mapped.slice(0,20)};}
function updateGffMeta(){const m=$('gffMeta');if(!m)return;m.textContent=`GFF3: ${gffFile?gffFile.name+' ('+gffFile.size.toLocaleString()+' bytes)':'未选择'}\n基因索引: ${gffIndex?gffIndex.stats.genes+' genes / '+gffIndex.stats.features+' features':'未构建'}\nSHA-256: ${gffSha||'未计算'}\n坐标映射: ${gffAudit?gffAudit.status+' · mapped '+gffAudit.mapped+'/'+gffAudit.seqids+(gffAudit.unmappedCount?' · unmapped '+gffAudit.unmappedCount:''):'未审计'}`;}
async function indexGffFile(){if(!gffFile)return;clearError();cancelRequested=false;setBusy(true);try{const idx=RPSGFF.newIndex(),reader=gffFile.stream().getReader(),dec=new TextDecoder('utf-8',{fatal:true}),hasher=new RPSIO.SHA256(),total=gffFile.size;let read=0,buf='';while(true){const {done,value}=await reader.read();if(done)break;if(cancelRequested)throw Error('CANCELLED');read+=value.byteLength;hasher.update(value);buf+=dec.decode(value,{stream:true});if(buf.length>4*1024*1024)throw Error('GFF3 line exceeds 4 MB');let cut;while((cut=buf.indexOf('\n'))>=0){const line=buf.slice(0,cut).replace(/\r$/,'');buf=buf.slice(cut+1);idx.stats.lines++;const f=RPSGFF.parseLine(line);if(f&&['gene','mrna','transcript','exon','cds','five_prime_UTR','three_prime_UTR','utr'].includes(f.type))RPSGFF.add(idx,f);else if(line&&line[0]!=='#'&&!f)idx.stats.skipped++;}showProgress(`索引 GFF3 · ${(100*read/Math.max(1,total)).toFixed(0)}%`,read/Math.max(1,total));}buf+=dec.decode();if(buf.trim()){idx.stats.lines++;const f=RPSGFF.parseLine(buf);if(f)RPSGFF.add(idx,f);}gffIndex=RPSGFF.finalize(idx);gffSha=hasher.digest();gffAudit=gffMappingAudit();updateGffMeta();$('geneQuery').disabled=false;$('geneSearch').disabled=false;$('clearGff').disabled=false;toast(`GFF3 索引完成：${gffIndex.stats.genes} genes。`);}catch(e){if(e.message!=='CANCELLED')error(e.message);}finally{setBusy(false);}}
function escAttr(s){return esc(s==null?'':String(s));}
function geneStructureSvg(g){const fs=RPSGFF.descendants(gffIndex,g.id),tx=fs.filter(x=>x.type==='mrna'||x.type==='transcript'),pick=(tx[0]||g);let feats=pick===g?fs:RPSGFF.descendants(gffIndex,pick.id);feats=feats.filter(x=>['exon','cds','five_prime_UTR','three_prime_UTR','utr'].includes(x.type));const lo=g.start,hi=g.end,w=760,sc=x=>25+(x-lo)/Math.max(1,hi-lo)*(w-50);let z=`<line x1="25" y1="50" x2="${w-25}" y2="50" stroke="#aebcaf" stroke-width="2"/>`;for(const f of feats){const x=sc(f.start),x2=sc(f.end),ww=Math.max(2,x2-x);const y=f.type==='cds'?39:43,h=f.type==='cds'?22:14;z+=`<rect x="${x}" y="${y}" width="${ww}" height="${h}" rx="2" fill="${f.type==='cds'?'#367b59':'#a9c6b0'}"/>`;}z+=`<text x="25" y="88" fill="#6d7d76" font-size="10">${escAttr(g.seqid)}:${g.start.toLocaleString()}–${g.end.toLocaleString()} · ${escAttr(g.strand)} · ${escAttr(pick.id||'gene')}</text>`;return `<div class="gene-graphic"><svg viewBox="0 0 ${w} 105" role="img" aria-label="gene structure">${z}</svg></div>`;}
function renderGeneDetail(g){gffSelectedGene=g;let mapped=null,mapMsg='';try{const e=RPS.resolveContig(entries,g.seqid);mapped=e.name;if(g.end>e.length)mapMsg='基因坐标超过当前参考 contig 长度；禁止按此注释定位。';else mapMsg=`映射到当前 FASTA: ${e.name} · 坐标在 contig 长度内`; }catch(e){mapMsg='无法将该 GFF3 seqid 唯一映射到当前 FASTA；仅显示注释，不建议用于坐标设计。';}
const fs=RPSGFF.descendants(gffIndex,g.id),counts={};for(const f of fs)counts[f.type]=(counts[f.type]||0)+1;const chips=Object.entries(counts).map(([k,v])=>`<span class="feature-chip">${escAttr(k)} ${v}</span>`).join('');$('geneDetail').innerHTML=`<div class="gene-card"><h3>${escAttr(g.id||g.name)}</h3><div class="hint">${escAttr(g.seqid)}:${g.start.toLocaleString()}–${g.end.toLocaleString()} · strand ${escAttr(g.strand)} · ${chips}</div>${geneStructureSvg(g)}<div class="gene-audit">${escAttr(mapMsg)}</div><div class="button-row" style="margin-top:10px"><button class="btn small" id="useGeneStart" ${mapped&&!mapMsg.startsWith('基因坐标超过')?'':'disabled'}>定位到 gene start</button><button class="btn small" id="useGeneCenter" ${mapped&&!mapMsg.startsWith('基因坐标超过')?'':'disabled'}>定位到 gene center</button><button class="btn small" id="useGeneEnd" ${mapped&&!mapMsg.startsWith('基因坐标超过')?'':'disabled'}>定位到 gene end</button></div><div class="hint">定位按钮只填写 Chr/POS，并从 FASTA 读取 1 bp REF；ALT 必须由用户根据真实变异填写。</div></div>`;if(mapped&&!mapMsg.startsWith('基因坐标超过')){const use=async pos=>{setInput('single');$('chrom').value=mapped;$('pos').value=pos;$('refFetchLen').value=1;if(!(await fetchReferenceAtTarget()))return;$('alt').value='';$('targetId').value=(g.id||'gene')+'_'+pos;dirty();$('workbench').scrollIntoView({behavior:'smooth'});};$('useGeneStart').onclick=()=>use(g.start);$('useGeneCenter').onclick=()=>use(Math.floor((g.start+g.end)/2));$('useGeneEnd').onclick=()=>use(g.end);}}
function searchGenes(){if(!gffIndex)return;const q=$('geneQuery').value.trim(),hits=RPSGFF.search(gffIndex,q,30),box=$('geneResults');box.classList.remove('hidden');box.innerHTML=hits.length?hits.map((g,i)=>`<div class="gene-row"><strong>${escAttr(g.id||g.name)}</strong><small>${escAttr(g.seqid)}:${g.start.toLocaleString()}–${g.end.toLocaleString()} · ${escAttr(g.strand)}</small><button class="btn small" data-gene-hit="${i}">查看</button></div>`).join(''):`<div class="gene-empty">未找到匹配的 gene。建议输入完整 LOC_Os ID；搜索只针对已加载 GFF3 的 gene 索引。</div>`;box.querySelectorAll('[data-gene-hit]').forEach(b=>b.onclick=()=>renderGeneDetail(hits[Number(b.dataset.geneHit)]));if(hits.length===1)renderGeneDetail(hits[0]);}
// A file read may complete after the user changes the target. Commit only the
// newest response against exactly the file, index and form state it requested.
let referenceLookupRequest = 0;
async function fetchReferenceAtTarget() {
 const request = ++referenceLookupRequest;
 if (busy) return false;
 $('refFetchStatus').textContent = '';
 if (!refFile || !entries.length) {
  error('请先选择 FASTA 并载入 FAI。');
  return false;
 }
 const snapshot = {
  file: refFile, entries, revision: refRevision, inputMode, mode,
  chrom: $('chrom').value.trim(), pos: Number($('pos').value),
  len: Number($('refFetchLen').value), refValue: $('ref').value
 };
 const {chrom,pos,len} = snapshot;
 if (!chrom || !Number.isSafeInteger(pos) || pos < 1 ||
     !Number.isSafeInteger(len) || len < 1 || len > 200) {
  error('请先填写有效 Chr/POS；REF 读取长度须为 1–200 bp。');
  return false;
 }
 const stillCurrent = () => request === referenceLookupRequest && !busy &&
  refFile === snapshot.file && entries === snapshot.entries &&
  refRevision === snapshot.revision && inputMode === snapshot.inputMode &&
  mode === snapshot.mode && $('chrom').value.trim() === chrom &&
  Number($('pos').value) === pos && Number($('refFetchLen').value) === len &&
  $('ref').value === snapshot.refValue;
 const discard = () => {
  if (request === referenceLookupRequest) {
   $('refFetchStatus').textContent = '输入或参考已变化，未写入旧 REF；请重新读取。';
  }
  return false;
 };
 try {
  const entry = RPS.resolveContig(snapshot.entries, chrom);
  if (pos + len - 1 > entry.length) throw Error('读取区间超过 contig 长度。');
  $('refFetchStatus').textContent = `正在读取 ${entry.name}:${pos}–${pos+len-1}…`;
  const seq = await RPS.fetchSeq(snapshot.file,entry,pos-1,pos-1+len);
  if (!stillCurrent()) return discard();
  $('ref').value = seq;
  $('chrom').value = entry.name;
  $('refFetchStatus').textContent = `${entry.name}:${pos}–${pos+len-1} · ${seq}`;
  dirty();
  toast('已从当前 FASTA 读取 REF；请根据真实变异填写 ALT。');
  return true;
 } catch (e) {
  if (!stillCurrent()) return discard();
  $('refFetchStatus').textContent = '读取失败，REF 未更新。';
  error(e.message);
  return false;
 }
}
function featureOverlap(a,b,f){return f.start<=b&&f.end>=a;}
function txFeatures(g){const all=RPSGFF.descendants(gffIndex,g.id),txs=all.filter(x=>x.type==='mrna'||x.type==='transcript');return (txs.length?txs:[g]).map(tx=>{const fs=tx===g?all:RPSGFF.descendants(gffIndex,tx.id);return {tx,features:fs};});}
function boundaryDistance(start,end,coord){return coord<start?start-coord:coord>end?coord-end:0;}
function featureContextForGene(g,start,end){const txOut=[];for(const {tx,features} of txFeatures(g)){const exons=features.filter(f=>f.type==='exon').sort((a,b)=>a.start-b.start);const cds=features.filter(f=>f.type==='cds');const utr=features.filter(f=>['five_prime_UTR','three_prime_UTR','utr'].includes(f.type));const hitEx=exons.filter(f=>featureOverlap(start,end,f)),hitCds=cds.filter(f=>featureOverlap(start,end,f)),hitUtr=utr.filter(f=>featureOverlap(start,end,f));let region=hitCds.length?'CDS':hitUtr.length?'UTR':hitEx.length?'EXON_NONCDS':(start>=g.start&&end<=g.end?'INTRON_OR_OTHER_GENE_REGION':'GENE_OVERLAP');const txStrand=['+','-'].includes(tx.strand)?tx.strand:g.strand;const order=txStrand==='-'?[...exons].sort((a,b)=>b.start-a.start):exons;const donors=[],acceptors=[];for(let i=0;i<order.length-1;i++){const donor=txStrand==='-'?order[i].start:order[i].end;const acceptor=txStrand==='-'?order[i+1].end:order[i+1].start;if(Math.abs(acceptor-donor)>2&&['+','-'].includes(txStrand)){donors.push(donor);acceptors.push(acceptor);}}const nearest=(arr,type)=>arr.length?arr.map(coord=>({type,coord,distance_bp:boundaryDistance(start,end,coord)})).sort((a,b)=>a.distance_bp-b.distance_bp||a.coord-b.coord)[0]:null;if(!exons.length&&!hitCds.length&&!hitUtr.length)region='EXON_STRUCTURE_UNAVAILABLE';txOut.push({transcript_id:tx.id||tx.name||g.id,strand:txStrand,region,exon_hits:hitEx.map(f=>({id:f.id,start:f.start,end:f.end})),cds_hits:hitCds.map(f=>({id:f.id,start:f.start,end:f.end,phase:f.phase})),utr_hits:hitUtr.map(f=>({id:f.id,type:f.type,start:f.start,end:f.end})),donor_boundaries:donors,acceptor_boundaries:acceptors,nearest_donor_boundary:nearest(donors,'DONOR_BOUNDARY'),nearest_acceptor_boundary:nearest(acceptors,'ACCEPTOR_BOUNDARY')});}return {gene_id:g.id||g.name,seqid:g.seqid,start:g.start,end:g.end,strand:g.strand,transcripts:txOut};}
function targetFeatureContext(context){if(!gffIndex||!context?.chrom||!Number.isInteger(context.pos))return [];const st=context.pos,en=st+Math.max(1,(context.ref||'').length)-1;let genes=[];try{genes=RPSGFF.overlapGenes(gffIndex,context.chrom,st,en,n=>RPS.resolveContig(entries,n).name);}catch{return [];}return genes.slice(0,20).map(g=>featureContextForGene(g,st,en));}
function compactFeatureLabel(fc){if(!fc?.length)return '无 gene overlap';const out=[];for(const g of fc){const regs=[...new Set(g.transcripts.map(t=>t.region))];out.push(`${g.gene_id}:${regs.join('|')}`);}return out.join('; ');}
function nearestBoundarySummary(fc,key){let best=null;for(const g of fc||[])for(const t of g.transcripts||[]){const x=t[key];if(x&&(!best||x.distance_bp<best.distance_bp))best={...x,gene:g.gene_id,tx:t.transcript_id};}return best;}
function variantStructureSvg(context,fc){if(!fc?.length)return '';const g=fc[0],gene=gffIndex.byId.get(g.gene_id);if(!gene)return '';const all=RPSGFF.descendants(gffIndex,gene.id),tx=all.find(x=>x.type==='mrna'||x.type==='transcript')||gene,feats=(tx===gene?all:RPSGFF.descendants(gffIndex,tx.id)).filter(x=>['exon','cds','five_prime_UTR','three_prime_UTR','utr'].includes(x.type));const lo=gene.start,hi=gene.end,w=820,sc=x=>30+(x-lo)/Math.max(1,hi-lo)*(w-60),v1=sc(Math.max(lo,context.pos)),v2=sc(Math.min(hi,context.pos+Math.max(1,context.ref.length)-1));let z=`<line x1="30" y1="56" x2="${w-30}" y2="56" stroke="#aebcaf" stroke-width="2"/>`;for(const f of feats){const x=sc(f.start),x2=sc(f.end),ww=Math.max(2,x2-x),isC=f.type==='cds';z+=`<rect x="${x}" y="${isC?43:48}" width="${ww}" height="${isC?26:16}" rx="2" fill="${isC?'#367b59':'#a9c6b0'}"/>`;}z+=`<line x1="${v1}" y1="25" x2="${v1}" y2="84" stroke="#ac4540" stroke-width="2"/><line x1="${Math.max(v1+1,v2)}" y1="25" x2="${Math.max(v1+1,v2)}" y2="84" stroke="#ac4540" stroke-width="2"/><text x="${Math.min(w-130,v1+4)}" y="20" fill="#ac4540" font-size="10">variant ${context.pos}</text><text x="30" y="105" fill="#6d7d76" font-size="10">${escAttr(g.gene_id)} · ${escAttr(tx.id||'gene')} · ${escAttr(g.strand)}</text>`;return `<div class="variant-track"><svg viewBox="0 0 ${w} 120" role="img" aria-label="variant on gene structure">${z}</svg></div>`;}
function rc2(s){return RPS.rc(String(s||''));}
function intervalOverlaps1(a,b,c,d){return a<=d&&b>=c;}
async function spliceMotifAudit(context,fc,maxDistance=20){
 if(!refFile||!entries.length||!fc?.length||context.source!=='LOCAL_REFERENCE')return {status:'NOT_RUN',reason:'LOCAL_FASTA_AND_GFF_REQUIRED',records:[]};
 const e=RPS.resolveContig(entries,context.chrom);
 const observed=await RPS.fetchSeq(refFile,e,context.pos-1,context.pos-1+context.ref.length);
 if(observed!==context.ref)throw Error(`结构检查 REF 不匹配：输入 ${context.ref}，当前 FASTA ${observed}。`);
 const edit=RPSAudit.normalizeEdit(context),vStart=context.pos,vEnd=context.pos+context.ref.length-1,records=[];
 for(const g of fc)for(const tx of g.transcripts||[]){
  if(!['+','-'].includes(tx.strand||g.strand))continue;
  const strand=tx.strand||g.strand;
  for(const kind of ['DONOR','ACCEPTOR'])for(const boundary of (kind==='DONOR'?tx.donor_boundaries:tx.acceptor_boundaries)||[]){
   const dist=boundaryDistance(vStart,vEnd,boundary);if(dist>maxDistance)continue;
   const mStart=strand==='-'?(kind==='DONOR'?boundary-2:boundary+1):(kind==='DONOR'?boundary+1:boundary-2),mEnd=mStart+1;
   if(mStart<1||mEnd>e.length)continue;
   const genomic=await RPS.fetchSeq(refFile,e,mStart-1,mEnd),oriented=strand==='-'?RPS.rc(genomic):genomic;
   const relation=RPSAudit.motifRelationship(context,mStart,mEnd);
   records.push({gene:g.gene_id,transcript:tx.transcript_id,strand,kind,boundary_coord:boundary,distance_to_boundary_bp:dist,
    motif_genomic_start:mStart,motif_genomic_end:mEnd,motif_genomic:genomic,motif_transcript_5to3:oriented,
    canonical:oriented===(kind==='DONOR'?'GT':'AG'),motifContainsAmbiguity:!/^[ACGT]{2}$/.test(oriented),
    record_overlaps_motif:relation.recordOverlaps,changed_reference_overlaps_motif:relation.changedReferenceOverlaps,
    insertion_inside_motif:relation.insertionInside,insertion_at_motif_edge:relation.insertionAtEdge,
    variant_overlaps_motif:relation.editDirectlyIntersectsMotif});
  }
 }
 return {status:'PERFORMED_LOCAL_FASTA',refCheck:'MATCH',definitionVersion:'EDIT_FOOTPRINT_v1',editFootprint:edit,max_distance_bp:maxDistance,records,
  variant_overlaps_any_motif:records.some(x=>x.variant_overlaps_motif),noncanonical_nearby:records.some(x=>!x.canonical&&!x.motifContainsAmbiguity),
  interpretation:'REFERENCE_MOTIF_GEOMETRY_ONLY_NOT_ALTERNATIVE_SPLICING_OR_FUNCTION'};
}
function spliceAuditHTML(a){if(!a||a.status!=='PERFORMED_LOCAL_FASTA')return '<div class="hint">剪接二核苷酸尚未从本地 FASTA 检查。</div>';if(!a.records.length)return '<div class="hint">当前变异 ±20 bp 内未发现已索引 transcript 的 exon junction boundary。</div>';return `<div class="table-wrap"><table><thead><tr><th>Transcript</th><th>类型</th><th>边界距离</th><th>转录方向 motif</th><th>canonical</th><th>改变区间 / 插入断点直接涉及 motif</th></tr></thead><tbody>${a.records.map(x=>`<tr><td class="mono">${esc(x.transcript)}</td><td>${x.kind}</td><td>${x.distance_to_boundary_bp} bp</td><td class="mono">${esc(x.motif_transcript_5to3)}</td><td>${x.canonical?'YES':'NO'}</td><td>${x.variant_overlaps_motif?'YES':x.insertion_at_motif_edge?'EDGE_ONLY':'NO'}</td></tr>`).join('')}</tbody></table></div><div class="hint">先核对 REF，再去掉 REF/ALT 公共前后缀；未改变的 anchor 不算 motif 改变。插入用断点，EDGE_ONLY 仅指相邻边缘。canonical 仅比较常见 GT/AG，不否定 GC/AG、AT/AC 等其他类型；不预测 ALT 剪接后果。</div>`;}
function featureDiscordant(fc){for(const g of fc||[]){const x=[...new Set((g.transcripts||[]).map(t=>t.region))];if(x.length>1)return true;}return false;}

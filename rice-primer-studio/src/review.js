/* Auditable input masks, amplicon extraction and review-only JSON import.
 * Variant POS: 1-based. Mask coordinates: 0-based, half-open.
 * No genotype/sample inference; every supplied background row is treated as an exclusion.
 */
(function(root){'use strict';
const R=root.RPS||(typeof require!=='undefined'?require('./engine.js'):null);
const IO=root.RPSIO||(typeof require!=='undefined'?require('./io.js'):null);
function assert(x,msg){if(!x)throw Error(msg);}
function hash(s){return new IO.SHA256().update(new TextEncoder().encode(s)).digest();}
function alias(s){return String(s).replace(/^chr/i,'').replace(/^0+(?=\d)/,'');}
function parseBackground(text){
 text=String(text).replace(/^\uFEFF/,'');if(!text.trim())return [];
 assert(text.length<=20*1024*1024,'背景差异文本超过 20 MB；请提供经过筛选的相关清单。');
 const lines=text.split(/\r?\n/).filter(x=>x.trim()&&!x.startsWith('##'));
 const h=lines[0].replace(/^#/,'').trim().split(/\t| +/).map(x=>x.toLowerCase());
 assert(h.includes('chrom')&&h.includes('pos')&&h.includes('ref')&&h.includes('alt'),'背景差异需带 chrom、pos、ref、alt 表头的 TSV 或标准 VCF。');
 const at=Object.fromEntries(h.map((v,i)=>[v,i])),out=[];
 for(let i=1;i<lines.length;i++){
  if(lines[i].startsWith('#'))continue;
  const f=lines[i].trim().split(/\t| +/),chrom=f[at.chrom],pos=Number(f[at.pos]);
  const ref=(f[at.ref]||'').toUpperCase(),alts=(f[at.alt]||'').toUpperCase().split(',');
  assert(chrom&&Number.isSafeInteger(pos)&&pos>0,`背景第 ${i+1} 行坐标无效。`);
  assert(/^[ACGT]+$/.test(ref)&&ref.length<=10000&&alts.length<=20&&alts.every(a=>/^[ACGT]+$/.test(a)&&a!==ref&&a.length<=10000),`背景第 ${i+1} 行不是支持的序列等位；符号等位须先另行处理。`);
  out.push({id:f[at.id]||`bg_${i}`,chrom,pos,ref,alts});
  assert(out.length<=50000,'背景清单最多 50,000 行；请先按目标区域筛选。');
 }
 return out;
}
function parseLocalRegions(text,n){
 if(!String(text).trim())return [];
 return String(text).split(/[;,\n，；]+/).filter(s=>s.trim()).map(s=>{
  const m=s.trim().match(/^(\d+)\s*[-:]\s*(\d+)$/);assert(m,'手工排除区间请写作 100-150，每行一段；坐标相对于本次提取的局部模板。');
  const st=Number(m[1])-1,en=Number(m[2]);assert(Number.isSafeInteger(st)&&Number.isSafeInteger(en)&&st>=0&&en>st&&en<=n,'手工排除区间越界；输入为局部模板的 1-based 闭区间。');
  return {start0:st,end0:en,source:'MANUAL_LOCAL_INTERVAL'};
 });
}
function applyExclusions(context,background=[],localText=''){
 const c=R.validateContext({...context}),regions=[...(c.excludedRegions||[]),...parseLocalRegions(localText,c.seq.length)];
 const audit={status:background.length?'PROVIDED_LIST_APPLIED_NOT_GENOME_WIDE_CERTIFICATION':'NOT_PROVIDED',providedRecords:background.length,nearbyRecords:0,excludedTargetRecords:0,refChecked:0,boundaryPartialChecks:0,unknownGenotypes:true};
 if(background.length){assert(c.chrom&&Number.isInteger(c.start1),'背景差异使用基因组坐标，不能套用于未定位的侧翼序列模式。');
  for(const v of background){
   if(v.resolvedChrom?v.resolvedChrom!==c.chrom:alias(v.chrom)!==alias(c.chrom))continue;
   const st=v.pos-c.start1,en=st+v.ref.length;
   if(en<=0||st>=c.seq.length)continue;
   audit.nearbyRecords++;
   // Check the REF overlap even when the full background record crosses the window.
   const a=Math.max(0,st),b=Math.min(c.seq.length,en);
   assert(c.seq.slice(a,b)===v.ref.slice(a-st,b-st),`背景 REF 不匹配：${v.chrom}:${v.pos}；请核对参考版本、文件与坐标。`);
   if(st>=0&&en<=c.seq.length)audit.refChecked++;else audit.boundaryPartialChecks++;
   // Only the exact same biallelic target is exempt. Additional alternate alleles remain masked.
   if(v.pos===c.pos&&v.ref===c.ref&&v.alts.length===1&&v.alts[0]===c.alt){audit.excludedTargetRecords++;continue;}
   const insertion=v.alts.some(a=>a.length>v.ref.length);
   regions.push({start0:a,end0:Math.min(c.seq.length,b+(insertion?1:0)),source:'BACKGROUND_VARIANT',id:v.id,pos1:v.pos,ref:v.ref,alts:v.alts});
  }
 }
 regions.sort((a,b)=>a.start0-b.start0||a.end0-b.end0);
 return {...c,excludedRegions:regions,backgroundAudit:{...audit,maskedIntervals:regions.length,manualInput:localText}};
}
function amplicons(c,candidate){
 c=R.validateContext(c);const delta=c.alt.length-c.ref.length,altSeq=c.seq.slice(0,c.target0)+c.alt+c.seq.slice(c.target0+c.ref.length);
 if(candidate.mode==='sanger'){
  const f=candidate.primers.find(x=>x.role==='F'),r=candidate.primers.find(x=>x.role==='R');
  assert(f&&r,'缺少 F/R。');return {ref:c.seq.slice(f.start0,r.end0),alt:altSeq.slice(f.start0,r.end0+delta),refStart0:f.start0,altStart0:f.start0,orientation:'+',includesTail:false,targetRef0:c.target0-f.start0,targetAlt0:c.target0-f.start0};
 }
 const x=candidate.primers.find(x=>x.role==='AS_REF'),y=candidate.primers.find(x=>x.role==='AS_ALT'),z=candidate.primers.find(x=>x.role==='COMMON');
 assert(x&&y&&z,'缺少 KASP 三引物。');const dir=candidate.direction,ref=dir==='+'?c.seq:R.rc(c.seq),alt=dir==='+'?altSeq:R.rc(altSeq);
 const xs=dir==='+'?x.start0:c.seq.length-x.end0,ze=dir==='+'?z.end0:c.seq.length-z.start0;
 const ys=y.haplotypeStart0;
 assert(Number.isInteger(ys),'缺少 ALT 单倍型坐标。');
 return {ref:ref.slice(xs,ze),alt:alt.slice(ys,ze+delta),refStart0:xs,altStart0:ys,orientation:dir,includesTail:false};
}
function validateImport(raw){
 assert(raw&&typeof raw==='object'&&raw.manifest&&Array.isArray(raw.results),'不是支持的运行 JSON。');
 assert(['rps.run/2','rps.run/3','rps.run/4','rps.run/5'].includes(raw.schemaVersion),'不支持的运行快照 schema；不会猜测其结构。');
 assert(raw.results.length<=100,'运行快照位点数超过 100。');
 const data=JSON.parse(JSON.stringify(raw));
 for(const r of data.results){
  assert(r&&r.context&&typeof r.context.id==='string'&&r.context.id.length<=500,'快照缺少位点标识。');
  assert(['sanger','kasp'].includes(r.mode)&&Array.isArray(r.candidates)&&r.candidates.length<=10,'快照候选格式无效。');
  if(r.error){assert(typeof r.error==='string'&&r.candidates.length===0,'失败位点格式异常。');continue;}
  assert(r.inputContext&&r.parameters,'快照缺少输入模板或参数。');const c=R.validateContext(r.inputContext);R.validateParams(r.parameters,r.mode);
  for(const key of ['ref','alt','target0','start1','chrom','pos','source'])assert((r.context[key]??null)===(c[key]??null),'快照输入与上下文不一致：'+key);
  assert(!r.context.templateSha256||hash(c.seq)===r.context.templateSha256,'模板 SHA-256 不匹配，不能载入。');
  for(const k of r.candidates){
   assert(['+','-'].includes(k.direction),'候选方向无效。');
   assert(k.mode===r.mode&&Array.isArray(k.primers)&&k.primers.length===(r.mode==='sanger'?2:3),'候选模式或引物数量错误。');
   assert(Number.isFinite(k.score)&&Number.isFinite(k.crossEnd)&&Number.isInteger(k.productRef)&&Number.isInteger(k.productAlt),'候选统计字段无效。');
   const roles=k.primers.map(o=>o.role).sort().join(',');assert(roles===(r.mode==='sanger'?'F,R':'AS_ALT,AS_REF,COMMON'),'引物角色错误或重复。');
   for(const checkName of ['exactCheck','approxCheck']){const v=k[checkName];if(v){assert(v&&Array.isArray(v.results)&&v.results.length===k.primers.length,'扫描快照结构无效。');for(const x of v.results){assert(typeof x.id==='string'&&Number.isSafeInteger(x.totalHits)&&x.totalHits>=0&&Array.isArray(x.hits)&&x.hits.length<=10000,'扫描命中字段无效。');if(checkName==='approxCheck')assert(Array.isArray(x.byMismatches)&&Number.isInteger(v.contigsScanned)&&Number.isInteger(v.basesScanned)&&Number.isInteger(v.ambiguousBases),'近似扫描统计无效。');}}}
   for(const o of k.primers){
    assert(typeof o.core==='string'&&/^[ACGT]{16,32}$/.test(o.core)&&typeof o.tail==='string'&&/^[ACGT]*$/.test(o.tail)&&o.tail.length<=30&&o.sequence===o.tail+o.core,'候选序列或尾序不一致。');
    assert([o.tm,o.gc,o.selfEnd,o.hairpin].every(Number.isFinite)&&o.length===o.core.length&&o.fullLength===o.sequence.length,'引物统计字段无效。');
    assert(['+','-'].includes(o.strand),'引物方向无效。');
    if(o.role==='AS_ALT'){const as=c.seq.slice(0,c.target0)+c.alt+c.seq.slice(c.target0+c.ref.length),oriented=k.direction==='+'?as:R.rc(as);assert(Number.isInteger(o.haplotypeStart0)&&Number.isInteger(o.haplotypeEnd0)&&o.haplotypeStart0>=0&&o.haplotypeEnd0-o.haplotypeStart0===o.core.length&&oriented.slice(o.haplotypeStart0,o.haplotypeEnd0)===o.core,'ALT 引物与构建单倍型不一致。');}
    if(o.role!=='AS_ALT'){
     assert(Number.isInteger(o.start0)&&Number.isInteger(o.end0)&&o.start0>=0&&o.end0<=c.seq.length&&o.end0-o.start0===o.core.length,'参考结合区坐标无效。');
     const s=c.seq.slice(o.start0,o.end0);assert((o.strand==='+'?s:R.rc(s))===o.core,'引物结合区与快照模板不一致。');
    }
   }
   const a=amplicons(c,k);assert(a.ref.length===k.productRef&&a.alt.length===k.productAlt,'扩增子长度与快照不一致。');
   k.originalImportedStatus=k.status;k.status=k.variantClass==='INDEL'?'INDEL_REQUIRES_EXPERT_REVIEW':'CANDIDATE_UNVALIDATED';
   k.selected=false;
  }
  r.refRevision=-1;r.importedReviewOnly=true;
  if(r.context.spliceMotifAudit){r.context.spliceMotifAudit.originalStatus=r.context.spliceMotifAudit.status;r.context.spliceMotifAudit.status='HISTORICAL_NOT_RECOMPUTED';}
  r.importAudit={status:'ARCHIVED_REVIEW_ONLY',checksNotAuthenticated:true};
 }
 data.manifest.importStatus='REVIEW_ONLY_NOT_RECOMPUTED';return data;
}
const api={hash,parseBackground,parseLocalRegions,applyExclusions,amplicons,validateImport};
if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.RPSReview=api;
})(typeof globalThis!=='undefined'?globalThis:this);

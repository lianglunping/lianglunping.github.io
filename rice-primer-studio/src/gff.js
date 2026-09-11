/* GFF3 structural navigation. Strict coordinates; multipart feature IDs retained.
 * This parser does not identify an assembly, infer RNA products or score consequences.
 */
(function(root){'use strict';
function decoded(s){try{return decodeURIComponent(s);}catch{throw Error('Malformed GFF3 percent encoding');}}
function attrs(text){const out=Object.create(null);for(const p of String(text||'').split(';')){const i=p.indexOf('=');if(i<1)continue;const k=p.slice(0,i).trim();out[k]=decoded(p.slice(i+1).trim());}return out;}
function normType(t){t=String(t||'').toLowerCase();return t==='five_prime_utr'||t==='5utr'?'five_prime_UTR':t==='three_prime_utr'||t==='3utr'?'three_prime_UTR':t;}
function parseLine(line){
 if(!line||line.startsWith('#'))return null;
 const f=line.split('\t');if(f.length!==9)throw Error('GFF3 must contain 9 tab-separated columns');
 const start=Number(f[3]),end=Number(f[4]);
 if(!f[0]||!Number.isSafeInteger(start)||!Number.isSafeInteger(end)||start<1||end<start||!['+','-','.','?'].includes(f[6]))throw Error('Invalid GFF3 coordinate or strand');
 const a=attrs(f[8]),rawParent=f[8].split(';').find(x=>x.startsWith('Parent='));
 // Split literal list commas before decoding. An encoded comma belongs to an ID.
 const parent=rawParent?rawParent.slice(7).split(',').filter(Boolean).map(decoded):[];
 const type=normType(f[2]);
 if(['gene','mrna','transcript'].includes(type)&&!a.ID)throw Error('gene/transcript requires a GFF3 ID');
 return {seqid:f[0],source:f[1],type,start,end,strand:f[6],phase:f[7],attrs:a,id:a.ID||null,parent,name:a.Name||a.ID||null};
}
function newIndex(){return {genes:[],byId:new Map(),segments:new Map(),children:new Map(),seqids:new Set(),bySeq:new Map(),maxEndBySeq:new Map(),
 stats:{lines:0,features:0,genes:0,transcripts:0,exons:0,cds:0,skipped:0,multipartSegments:0,orphanParents:0},warnings:[]};}
function add(idx,f){
 if(!f)return;
 if(f.id&&idx.byId.has(f.id)){
  const old=idx.byId.get(f.id);
  if(['gene','mrna','transcript'].includes(f.type)||old.type!==f.type||old.seqid!==f.seqid||old.strand!==f.strand||old.parent.join(',')!==f.parent.join(','))throw Error('Conflicting or unsupported repeated structural ID: '+f.id);
  idx.stats.multipartSegments++;
 } else if(f.id)idx.byId.set(f.id,f);
 if(f.id){if(!idx.segments.has(f.id))idx.segments.set(f.id,[]);idx.segments.get(f.id).push(f);}
 idx.stats.features++;idx.seqids.add(f.seqid);idx.maxEndBySeq.set(f.seqid,Math.max(idx.maxEndBySeq.get(f.seqid)||0,f.end));
 for(const p of f.parent){if(p===f.id)throw Error('Self-parent GFF3 feature: '+p);if(!idx.children.has(p))idx.children.set(p,[]);idx.children.get(p).push(f);}
 if(f.type==='gene'){idx.genes.push(f);idx.stats.genes++;}
 else if(f.type==='mrna'||f.type==='transcript')idx.stats.transcripts++;
 else if(f.type==='exon')idx.stats.exons++;
 else if(f.type==='cds')idx.stats.cds++;
}
function finalize(idx){
 idx.genes.sort((a,b)=>a.seqid.localeCompare(b.seqid,undefined,{numeric:true})||a.start-b.start||a.end-b.end);
 for(const p of idx.children.keys())if(!idx.byId.has(p)){idx.stats.orphanParents++;if(idx.warnings.length<100)idx.warnings.push('Missing parent: '+p);}
 for(const f of idx.byId.values())if(f.type==='mrna'||f.type==='transcript'){
  for(const p of f.parent){const parent=idx.byId.get(p);if(parent&&(parent.seqid!==f.seqid||parent.start>f.start||parent.end<f.end))throw Error('Transcript outside parent bounds: '+f.id);}
 }
 idx.bySeq=new Map();
 for(const g of idx.genes){if(!idx.bySeq.has(g.seqid))idx.bySeq.set(g.seqid,{genes:[],prefixMaxEnd:[]});const b=idx.bySeq.get(g.seqid);b.genes.push(g);b.prefixMaxEnd.push(Math.max(b.prefixMaxEnd.at(-1)||0,g.end));}
 return idx;
}
function geneOf(idx,f){if(!f)return null;if(f.type==='gene')return f;const seen=new Set(),queue=[...(f.parent||[])];while(queue.length){const id=queue.shift();if(seen.has(id))continue;seen.add(id);const p=idx.byId.get(id);if(!p)continue;if(p.type==='gene')return p;queue.push(...p.parent);}return null;}
function search(idx,query,limit=30){const raw=String(query||'').trim(),q=raw.toLowerCase();if(!q)return [];const exact=idx.byId.get(raw),out=[];if(exact){const g=geneOf(idx,exact);if(g)out.push(g);}for(const g of idx.genes){const hay=[g.id,g.name,g.attrs.Alias,g.attrs.Note,g.attrs.description].filter(Boolean).join(' ').toLowerCase();if(hay.includes(q)&&!out.includes(g)){out.push(g);if(out.length>=limit)break;}}return out;}
function directChildren(idx,id){return (idx.children.get(id)||[]).slice().sort((a,b)=>a.start-b.start||a.end-b.end);}
function descendants(idx,id,max=10000){const out=[],queue=[id],seenId=new Set([id]),seenObjects=new Set();while(queue.length){for(const c of idx.children.get(queue.shift())||[]){if(!seenObjects.has(c)){if(out.length>=max)throw Error('Feature traversal limit exceeded; refusing silent truncation');seenObjects.add(c);out.push(c);}if(c.id&&!seenId.has(c.id)){seenId.add(c.id);queue.push(c.id);}}}return out;}
function overlapGenes(idx,chrom,start,end,resolve){
 let target;try{target=resolve?resolve(chrom):chrom;}catch{return [];}
 const out=[];
 for(const [seqid,b] of idx.bySeq){let resolved;try{resolved=resolve?resolve(seqid):seqid;}catch{continue;}if(resolved!==target)continue;
  let lo=0,hi=b.genes.length;
  while(lo<hi){const m=(lo+hi)>>1;if(b.genes[m].start<=end)lo=m+1;else hi=m;}
  for(let i=lo-1;i>=0&&b.prefixMaxEnd[i]>=start;i--){const g=b.genes[i];if(g.end>=start)out.push(g);}
 }
 return out.sort((a,b)=>a.start-b.start||a.end-b.end);
}
const api={attrs,parseLine,newIndex,add,finalize,search,geneOf,directChildren,descendants,overlapGenes};
if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.RPSGFF=api;
})(typeof globalThis!=='undefined'?globalThis:this);

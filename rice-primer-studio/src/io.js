/* Streaming reference operations. No network requests. */
function RPSIOFactory(root){'use strict';
const K=new Uint32Array([0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]);
const rotr=(x,n)=>(x>>>n)|(x<<(32-n));
class SHA256 {
 constructor(){this.h=new Uint32Array([0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]);this.buffer=new Uint8Array(64);this.used=0;this.length=0;this.w=new Uint32Array(64);this.done=false;}
 block(b,off=0){const w=this.w;for(let i=0;i<16;i++){const k=off+4*i;w[i]=((b[k]<<24)|(b[k+1]<<16)|(b[k+2]<<8)|b[k+3])>>>0;}for(let i=16;i<64;i++){const x=w[i-15],y=w[i-2];w[i]=(w[i-16]+(rotr(x,7)^rotr(x,18)^(x>>>3))+w[i-7]+(rotr(y,17)^rotr(y,19)^(y>>>10)))>>>0;}
 let [a,b0,c,d,e,f,g,h]=this.h;
 for(let i=0;i<64;i++){const t1=(h+(rotr(e,6)^rotr(e,11)^rotr(e,25))+((e&f)^(~e&g))+K[i]+w[i])>>>0;const t2=((rotr(a,2)^rotr(a,13)^rotr(a,22))+((a&b0)^(a&c)^(b0&c)))>>>0;h=g;g=f;f=e;e=(d+t1)>>>0;d=c;c=b0;b0=a;a=(t1+t2)>>>0;}
 [a,b0,c,d,e,f,g,h].forEach((v,i)=>this.h[i]=(this.h[i]+v)>>>0);}
 update(bytes){if(this.done)throw Error('SHA256 already finalized');this.length+=bytes.length;let i=0;if(this.used){const n=Math.min(64-this.used,bytes.length);this.buffer.set(bytes.subarray(0,n),this.used);this.used+=n;i=n;if(this.used===64){this.block(this.buffer);this.used=0;}}
 for(;i+64<=bytes.length;i+=64)this.block(bytes,i);if(i<bytes.length){this.buffer.set(bytes.subarray(i),0);this.used=bytes.length-i;}return this;}
 digest(){if(this.done)return this.result;const bits=BigInt(this.length)*8n;this.buffer[this.used++]=0x80;if(this.used>56){this.buffer.fill(0,this.used);this.block(this.buffer);this.used=0;}this.buffer.fill(0,this.used,56);for(let i=0;i<8;i++)this.buffer[63-i]=Number((bits>>BigInt(i*8))&255n);this.block(this.buffer);this.done=true;this.result=Array.from(this.h,x=>x.toString(16).padStart(8,'0')).join('');return this.result;}
}
function assert(ok,msg){if(!ok)throw Error(msg);}
async function eachLine(file,callback,progress=()=>{},hash=null){let carry='',offset=0;const decoder=new TextDecoder('utf-8',{fatal:true,ignoreBOM:true});const chunkSize=2**20;for(let start=0;start<file.size;start+=chunkSize){const bytes=new Uint8Array(await file.slice(start,start+chunkSize).arrayBuffer());if(hash)hash.update(bytes);const s=carry+decoder.decode(bytes,{stream:true});let from=0;while(true){const n=s.indexOf('\n',from);if(n<0)break;const raw=s.slice(from,n+1);callback(raw.replace(/\r?\n$/,''),offset,raw.length);offset+=raw.length;from=n+1;}carry=s.slice(from);assert(carry.length<32*1024*1024,'单行 FASTA 超过 32 MB；请先用 seqkit seq -w 60 等工具规则换行。');progress(Math.min(1,(start+bytes.length)/file.size));}carry+=decoder.decode();if(carry)callback(carry,offset,carry.length);}
async function buildIndex(file,progress=()=>{}){const entries=[],sha=new SHA256();let cur=null,prev=null;
 function finish(){if(cur){assert(cur.length>0,'FASTA 存在空 contig。');entries.push(cur);}cur=null;prev=null;}
 await eachLine(file,(line,offset,width)=>{assert(/^[\x00-\x7F]*$/.test(line),'FASTA 需使用 ASCII 标题及序列，不能包含 BOM 或非 ASCII 字符。');if(line.startsWith('>')){finish();const name=line.slice(1).trim().split(/\s+/)[0];assert(name,'FASTA contig 名为空。');cur={name,length:0,offset:null,lineBases:0,lineWidth:0};}
 else{assert(cur,'FASTA 第一行须为 >contig。');assert(/^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$/.test(line),'FASTA 序列含空行、空格或非法碱基。');if(prev)assert(prev.bases===cur.lineBases&&prev.width===cur.lineWidth,'FASTA 换行不规则；短行只能位于 contig 最后一行。');if(cur.offset===null){cur.offset=offset;cur.lineBases=line.length;cur.lineWidth=width;}cur.length+=line.length;prev={bases:line.length,width};}
 },progress,sha);finish();assert(entries.length>0,'未找到 FASTA 序列。');assert(new Set(entries.map(x=>x.name)).size===entries.length,'重复 contig 名。');return {entries,sha256:sha.digest()};}
async function hashFile(file,progress=()=>{}){const s=new SHA256(),n=2**20;for(let i=0;i<file.size;i+=n){const b=new Uint8Array(await file.slice(i,i+n).arrayBuffer());s.update(b);progress(Math.min(1,(i+b.length)/file.size));}return s.digest();}
async function scanExact(file,queries,progress=()=>{}){const RPS=root.RPS||(typeof require!=='undefined'?require('./engine.js'):null);assert(Array.isArray(queries)&&queries.length>0,'扫描引物列表为空。');const normalized=queries.map(q=>({...q,sequence:q.sequence.toUpperCase()}));assert(new Set(normalized.map(q=>q.id)).size===normalized.length,'扫描引物 ID 不得重复。');normalized.forEach(q=>assert(/^[ACGT]+$/.test(q.sequence),'比对引物只能包含 A/C/G/T。'));
 const results=normalized.map(q=>({id:q.id,sequence:q.sequence,totalHits:0,hits:[],positionsTruncated:false}));const maxLen=Math.max(...normalized.map(q=>q.sequence.length));const patterns=normalized.flatMap((q,i)=>[{s:q.sequence,strand:'+',i},{s:RPS.rc(q.sequence),strand:'-',i}]);
 let chrom=null,buffer='',total=0,keep='',basesScanned=0,ambiguousBases=0;const seen=new Set();
 function flush(){if(!buffer)return;const text=(keep+buffer).toUpperCase(),base=total-keep.length,old=keep.length;for(const q of patterns){let at=text.indexOf(q.s);while(at>=0){if(at+q.s.length>old){const r=results[q.i];r.totalHits++;if(r.hits.length<1000)r.hits.push({chrom,start0:base+at,end0:base+at+q.s.length,strand:q.strand});else r.positionsTruncated=true;}at=text.indexOf(q.s,at+1);}}total+=buffer.length;keep=text.slice(-(maxLen-1));buffer='';}
 await eachLine(file,(line)=>{assert(/^[\x00-\x7F]*$/.test(line),'FASTA 不能包含 BOM 或非 ASCII 字符。');if(line.startsWith('>')){assert(chrom===null||total+buffer.length>0,'FASTA 存在空 contig。');flush();chrom=line.slice(1).trim().split(/\s+/)[0];assert(chrom&&!seen.has(chrom),'FASTA contig 名为空或重复。');seen.add(chrom);total=0;keep='';}else{assert(chrom&&/^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$/.test(line),'FASTA 格式错误。');basesScanned+=line.length;ambiguousBases+=(line.match(/[^ACGTacgt]/g)||[]).length;buffer+=line;if(buffer.length>=262144)flush();}},progress);assert(seen.size>0&&total+buffer.length>0,'FASTA 为空或最后一个 contig 无序列。');flush();return {method:'FULL_LENGTH_EXACT_MATCH_BOTH_STRANDS',referenceScope:'UPLOADED_FASTA_ONLY',mismatchesAllowed:0,gapsAllowed:false,tailIncluded:false,contigsScanned:seen.size,basesScanned,ambiguousBases,completed:true,results,warning:'只统计完整结合区的零错配命中；未检查近似匹配、错配延伸或真实 PCR。ALT 引物在 MSU 上无命中可能是预期现象。'};}
function exactProducts(a,b,min=40,max=3000){const found=[],seen=new Set();let total=0;for(const x of a.hits)for(const y of b.hits){if(x.chrom!==y.chrom||x.strand===y.strand)continue;const plus=x.strand==='+'?x:y,minus=x.strand==='-'?x:y;if(plus.end0>minus.start0)continue;const n=minus.end0-plus.start0;if(n<min||n>max)continue;const key=`${x.chrom}:${plus.start0}:${minus.end0}`;if(seen.has(key))continue;seen.add(key);total++;if(found.length<200)found.push({chrom:x.chrom,start1:plus.start0+1,end1:minus.end0,length:n});}
 return {count:total,products:found,truncated:total>200||a.positionsTruncated||b.positionsTruncated,minProduct:min,maxProduct:max,scope:'UPLOADED_REFERENCE_EXACT_ONLY'};}

// Exhaustive ungapped Hamming search by k+1 disjoint exact seeds.
// For <=k substitutions, at least one of the k+1 seeds must match exactly.
// Ambiguous reference windows are left unresolved, not treated as verified negatives.
async function scanApprox(file,queries,maxMismatches=2,progress=()=>{},options={}){
 const RPS=root.RPS||(typeof require!=='undefined'?require('./engine.js'):null);
 assert(Number.isInteger(maxMismatches)&&maxMismatches>=0&&maxMismatches<=2,'允许错配数须为 0、1 或 2。');
 assert(Array.isArray(queries)&&queries.length>=1&&queries.length<=30,'一次扫描须为 1–30 条引物。');
 const normalized=queries.map(q=>({...q,sequence:q.sequence.toUpperCase()}));
 normalized.forEach(q=>assert(typeof q.id==='string'&&/^[ACGT]{16,32}$/.test(q.sequence),'比对结合区须为 16–32 nt 的 A/C/G/T。'));
 assert(new Set(normalized.map(q=>q.id)).size===normalized.length,'扫描引物 ID 不得重复。');
 const limit=options.positionLimit??1000,blockSize=options.blockSize??262144;
 assert(Number.isInteger(limit)&&limit>=1&&limit<=10000&&Number.isInteger(blockSize)&&blockSize>=64,'扫描资源参数无效。');
 const maxLen=Math.max(...normalized.map(q=>q.sequence.length)),patterns=[];
 const results=normalized.map(q=>({id:q.id,sequence:q.sequence,totalHits:0,hits:[],positionsTruncated:false,byMismatches:Array(maxMismatches+1).fill(0)}));
 for(let i=0;i<normalized.length;i++)for(const strand of ['+','-']){
  const s=strand==='+'?normalized[i].sequence:RPS.rc(normalized[i].sequence),seeds=[];
  for(let k=0;k<=maxMismatches;k++){const start=Math.floor(k*s.length/(maxMismatches+1)),end=Math.floor((k+1)*s.length/(maxMismatches+1));seeds.push({start,s:s.slice(start,end)});}
  patterns.push({s,strand,i,seeds});
 }
 let chrom=null,buffer='',total=0,keep='',ambiguousBases=0,basesScanned=0;const seen=new Set();
 function flush(final=false){
  if(!buffer&&!keep)return;
  const text=(keep+buffer).toUpperCase(),base=total-keep.length;
  // Defer starts whose full longest query would cross a non-final block boundary.
  const safeStart=final?text.length:text.length-maxLen+1;
  for(const q of patterns){const visited=new Set();
   for(const seed of q.seeds){let at=text.indexOf(seed.s);
    while(at>=0){const st=at-seed.start;at=text.indexOf(seed.s,at+1);
     if(st<0||st+q.s.length>text.length||st>=safeStart||visited.has(st))continue;
     visited.add(st);let differences=[],valid=true;
     for(let j=0;j<q.s.length;j++){const b=text[st+j];if(!'ACGT'.includes(b)){valid=false;break;}if(b!==q.s[j]){differences.push(j);if(differences.length>maxMismatches){valid=false;break;}}}
     if(!valid)continue;
     const r=results[q.i];r.totalHits++;r.byMismatches[differences.length]++;
     if(r.hits.length<limit){const from3= differences.map(j=>q.strand==='+'?q.s.length-j:j+1).sort((a,b)=>a-b);r.hits.push({chrom,start0:base+st,end0:base+st+q.s.length,strand:q.strand,mismatches:differences.length,mismatchFrom3prime1:from3,last5Mismatches:from3.filter(x=>x<=5).length});}
     else r.positionsTruncated=true;
    }
   }
  }
  total+=buffer.length;keep=final?'':text.slice(-(maxLen-1));buffer='';
 }
 await eachLine(file,line=>{
  assert(/^[\x00-\x7F]*$/.test(line),'FASTA 不能包含 BOM 或非 ASCII 字符。');
  if(line.startsWith('>')){assert(chrom===null||total+buffer.length>0,'FASTA 存在空 contig。');flush(true);chrom=line.slice(1).trim().split(/\s+/)[0];assert(chrom&&!seen.has(chrom),'FASTA contig 名为空或重复。');seen.add(chrom);total=0;keep='';}
  else{assert(chrom&&/^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$/.test(line),'FASTA 格式错误。');ambiguousBases+=(line.match(/[^ACGTacgt]/g)||[]).length;basesScanned+=line.length;buffer+=line;if(buffer.length>=blockSize)flush(false);}
 },progress);
 assert(seen.size>0&&total+buffer.length>0,'FASTA 为空或最后一个 contig 无序列。');
 flush(true);
 results.forEach(r=>r.hits.sort((a,b)=>a.chrom.localeCompare(b.chrom)||a.start0-b.start0||a.strand.localeCompare(b.strand)));
 return {method:'DISJOINT_SEED_EXHAUSTIVE_UNGAPPED_HAMMING_v1',referenceScope:'UPLOADED_FASTA_ONLY',mismatchesAllowed:maxMismatches,gapsAllowed:false,tailIncluded:false,threePrimePerfectMatchRequired:false,contigsScanned:seen.size,basesScanned,ambiguousBases,ambiguousReferencePolicy:'WINDOWS_WITH_NON_ACGT_UNRESOLVED',results,completed:true,warning:'仅检查上传参考内、完整结合区至多 '+maxMismatches+' 处替换，不允许插缺；含模糊碱基的窗口未判定。3′ 错配逐项记录，不按错配自动判定能否扩增。KASP 的 REF/ALT 在目标位点可能出现预期交叉命中；不是实验放行。'};
}



async function verifyIndexBoundaries(file,entries){
 assert(entries.length>0,'FAI 为空。');const sorted=[...entries].sort((a,b)=>a.offset-b.offset);
 const first=(await file.slice(0,sorted[0].offset).text());
 assert(/^>[^\r\n]+\r?\n$/.test(first),'FAI 首条序列不是 FASTA 首个 contig，或标题超过支持范围。');
 for(let i=0;i<sorted.length;i++){
  const e=sorted[i];assert(e.offset>0&&e.lineBases<=e.length,'FAI 字段与序列长度不一致。');
  const before=await file.slice(Math.max(0,e.offset-8192),e.offset).text();
  const m=before.match(/(?:^|\n)>([^\s\r\n]+)[^\r\n]*\r?\n$/);
  assert(m&&m[1]===e.name,'FAI contig 名或首碱基偏移与 FASTA 标题不一致：'+e.name);
  const head=await file.slice(e.offset,e.offset+e.lineWidth).text();
  assert(head.replace(/\r?\n$/,'').length===e.lineBases&&/^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+(?:\r?\n)?$/.test(head),'FAI 首行长度或换行宽度错误：'+e.name);
  const last=e.offset+Math.floor((e.length-1)/e.lineBases)*e.lineWidth+(e.length-1)%e.lineBases+1;
  assert(last<=file.size&&(!sorted[i+1]||last<sorted[i+1].offset),'FAI 末端越过文件或下一 contig。');
  const suffix=await file.slice(last,Math.min(file.size,last+8192)).text();
  if(i===sorted.length-1)assert(suffix===''||suffix==='\n'||suffix==='\r\n','FAI 最后 contig 长度不符，或漏列参考序列。');
  else{const next=suffix.match(/^\r?\n>([^\s\r\n]+)/);assert(next&&next[1]===sorted[i+1].name,'FAI 长度、顺序或下一 contig 与 FASTA 不一致：'+e.name);}
 }
 return {method:'HEADER_FIRST_LINE_AND_END_BOUNDARY_CHECK',contigsChecked:sorted.length,fullFileHashed:false,identityVerified:false};
}

const api={verifyIndexBoundaries,scanApprox,SHA256,eachLine,buildIndex,hashFile,scanExact,exactProducts};if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.RPSIO=api;
}
RPSIOFactory(typeof globalThis!=='undefined'?globalThis:this);

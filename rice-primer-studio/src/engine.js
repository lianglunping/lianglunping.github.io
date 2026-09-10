/* Rice Primer Studio 1.1.0-rc3 — deterministic candidate screening, NOT Primer3.
 * Internal sequence coordinates: 0-based, half-open. Oligos: 5' -> 3'.
 * Tm: SantaLucia 1998 DNA_NN3, entropy salt correction, sodium-equivalent Mg.
 */
function RPSEngineFactory(root){
'use strict';
const VERSION='1.1.0-rc3';
const TAILS={FAM:'GAAGGTGACCAAGTTCATGCT',HEX:'GAAGGTCGGAGTCAACGGATT'};
const NN={AA:[-7.9,-22.2],TT:[-7.9,-22.2],AT:[-7.2,-20.4],TA:[-7.2,-21.3],CA:[-8.5,-22.7],TG:[-8.5,-22.7],GT:[-8.4,-22.4],AC:[-8.4,-22.4],CT:[-7.8,-21.0],AG:[-7.8,-21.0],GA:[-8.2,-22.2],TC:[-8.2,-22.2],CG:[-10.6,-27.2],GC:[-9.8,-24.4],GG:[-8.0,-19.9],CC:[-8.0,-19.9]};
const DEFAULT={minLen:18,maxLen:28,tmMin:57,tmOpt:60,tmMax:64,gcMin:30,gcMax:70,maxPoly:4,maxTmDiff:3,maxSelfEnd:4,maxHairpin:5,maxPairEnd:4,minProduct:400,maxProduct:800,readBuffer:80,numReturn:3,na:50,mg:1.5,dntp:0.2,dna:250,indelSuffix:12,direction:'both',famAllele:'REF'};
function assert(ok,msg){if(!ok)throw Error(msg);}
function rc(s){return s.toUpperCase().split('').reverse().map(x=>({A:'T',C:'G',G:'C',T:'A',N:'N',R:'Y',Y:'R',S:'S',W:'W',K:'M',M:'K',B:'V',V:'B',D:'H',H:'D'}[x]||'N')).join('');}
function tm(s,p={}){p={...DEFAULT,...p};s=s.toUpperCase();assert(/^[ACGT]{2,60}$/.test(s),'Tm 仅接受 2–60 nt 的 A/C/G/T 序列。');
 let h=0,e=0;for(const b of [s[0],s.at(-1)]){if(b==='A'||b==='T'){h+=2.3;e+=4.1;}else{h+=0.1;e-=2.8;}}
 for(let i=0;i<s.length-1;i++){const n=NN[s.slice(i,i+2)];h+=n[0];e+=n[1];}
 const self=s===rc(s);if(self)e-=1.4;
 const mon=(p.na+(p.mg>p.dntp?120*Math.sqrt(p.mg-p.dntp):0))*1e-3;
 assert(mon>0&&p.dna>0,'离子浓度与寡核苷酸浓度必须为正。');
 e+=0.368*(s.length-1)*Math.log(mon);
 const k=p.dna*1e-9/(self?1:4);
 return 1000*h/(e+1.987*Math.log(k))-273.15;
}
function gc(s){return 100*(s.match(/[GC]/g)||[]).length/s.length;}
function poly(s){let best=0,run=0,last='';for(const b of s){run=b===last?run+1:1;best=Math.max(best,run);last=b;}return best;}
// Ungapped contiguous complements only. These are NOT thermodynamic dimer energies.
function complement(a,b){const c=rc(b);let any=0,end=0;for(let shift=1-c.length;shift<a.length;shift++){let run=0;for(let i=Math.max(0,shift);i<Math.min(a.length,shift+c.length);i++){const j=i-shift;if(a[i]===c[j]){run++;any=Math.max(any,run);if(i===a.length-1||j-run+1===0)end=Math.max(end,run);}else run=0;}}return {any,end};}
function hairpin(s){let best=0;const comp={A:'T',T:'A',G:'C',C:'G'};for(let i=0;i<s.length;i++)for(let j=i+4;j<s.length;j++){let k=0;while(j-k-(i+k)>3&&s[i+k]===comp[s[j-k]])k++;best=Math.max(best,k);}return best;}
function suffixMatch(s,other){for(let n=Math.min(20,s.length);n>0;n--)if(other.includes(s.slice(-n)))return n;return 0;}
function oligo(core,role,start,end,strand,p,tail='',allele=null){const full=tail+core;return {role,allele,core,tail,sequence:full,length:core.length,fullLength:full.length,tm:tm(core,p),gc:gc(core),poly:poly(core),selfEnd:complement(full,full).end,hairpin:hairpin(full),start0:start,end0:end,strand};}
function basic(s,p){if(!/^[ACGT]+$/.test(s)||poly(s)>p.maxPoly)return false;const g=gc(s);if(g<p.gcMin||g>p.gcMax)return false;const t=tm(s,p);return t>=p.tmMin&&t<=p.tmMax;}
function structureOK(o,p){return o.selfEnd<=p.maxSelfEnd&&o.hairpin<=p.maxHairpin;}
function merit(o,p){return Math.abs(o.tm-p.tmOpt)+0.10*Math.abs(o.length-22)+0.012*Math.abs(o.gc-50)+0.10*o.hairpin+0.10*o.selfEnd;}
function validateParams(p,mode){p={...DEFAULT,...p};for(const key of Object.keys(DEFAULT).filter(k=>typeof DEFAULT[k]==='number'))assert(Number.isFinite(p[key]),`参数 ${key} 必须为数值。`);
 assert(Number.isInteger(p.minLen)&&Number.isInteger(p.maxLen)&&p.minLen>=16&&p.maxLen<=32&&p.minLen<=p.maxLen,'引物长度须为 16–32 nt，且最小值不大于最大值。');
 assert(p.tmMin<=p.tmOpt&&p.tmOpt<=p.tmMax&&p.tmMin>=30&&p.tmMax<=85,'Tm 范围或目标值无效。');
 assert(p.gcMin>=0&&p.gcMax<=100&&p.gcMin<=p.gcMax,'GC 范围无效。');
 assert(p.minProduct>=40&&p.maxProduct<=3000&&p.minProduct<=p.maxProduct&&Number.isInteger(p.minProduct)&&Number.isInteger(p.maxProduct),'产物长度须为 40–3000 bp 的整数。');
 assert(p.readBuffer>=0&&p.readBuffer<=1000&&Number.isInteger(p.readBuffer),'读序间隔须为 0–1000 的整数。');
 assert(p.numReturn>=1&&p.numReturn<=10&&Number.isInteger(p.numReturn),'候选数须为 1–10 的整数。');
 for(const k of ['maxPoly','maxSelfEnd','maxHairpin','maxPairEnd'])assert(Number.isInteger(p[k])&&p[k]>=0&&p[k]<=32,`${k} 须为 0–32 的整数。`);
 assert(p.maxTmDiff>=0&&p.maxTmDiff<=20,'最大 Tm 差值须为 0–20。');
 assert(p.na>0&&p.na<=1000&&p.mg>=0&&p.mg<=100&&p.dntp>=0&&p.dntp<=100&&p.dna>0&&p.dna<=10000,'溶液计算参数超出允许范围。');
 assert(['both','+','-'].includes(p.direction),'方向参数无效。');assert(['REF','ALT'].includes(p.famAllele),'FAM 对应等位无效。');
 assert(Number.isInteger(p.indelSuffix)&&p.indelSuffix>=8&&p.indelSuffix<=20,'KASP 另一等位 3′ 后缀阈值须为 8–20 nt。');
 if(mode==='kasp')assert(p.maxProduct<=500,'KASP 候选模式将产物上限限制为 500 bp。');return p;
}
function validateContext(c){assert(c&&typeof c.seq==='string'&&typeof c.ref==='string'&&typeof c.alt==='string','缺少参考序列或 REF/ALT。');c={...c,seq:c.seq.toUpperCase(),ref:c.ref.toUpperCase(),alt:c.alt.toUpperCase()};
 assert(/^[ACGTRYSWKMBDHVN]+$/.test(c.seq),'模板含有不支持的字符。');assert(/^[ACGT]+$/.test(c.ref)&&/^[ACGT]+$/.test(c.alt),'REF/ALT 必须是 A/C/G/T；InDel 保留 VCF anchor，不用 - 或 <DEL>。');
 assert(c.ref!==c.alt,'REF 与 ALT 相同。');assert(c.ref.length<=201&&c.alt.length<=201,'此版本仅处理长度不超过 201 nt 的序列等位。');
 assert(Number.isInteger(c.target0)&&c.target0>=0&&c.target0+c.ref.length<=c.seq.length,'目标坐标超出序列。');
 assert(c.seq.length<=20000,'局部模板不得超过 20,000 nt。');
 if(c.start1!=null){assert(Number.isSafeInteger(c.start1)&&c.start1>=1,'局部窗口起点无效。');if(c.pos!=null)assert(Number.isSafeInteger(c.pos)&&c.pos===c.start1+c.target0,'POS / start1 / target0 坐标不一致。');}
 c.excludedRegions=c.excludedRegions||[];assert(Array.isArray(c.excludedRegions)&&c.excludedRegions.length<=20000,'排除区间格式或数量无效。');
 for(const r of c.excludedRegions)assert(r&&Number.isInteger(r.start0)&&Number.isInteger(r.end0)&&r.start0>=0&&r.end0>r.start0&&r.end0<=c.seq.length,'排除区间越界；内部区间须为 0-based、左闭右开。');
 assert(c.seq.slice(c.target0,c.target0+c.ref.length)===c.ref,`REF 不匹配：输入 ${c.ref}，模板为 ${c.seq.slice(c.target0,c.target0+c.ref.length)}。`);
 return c;
}
function trimChange(c){let r=c.ref,a=c.alt,left=c.seq.slice(0,c.target0),right=c.seq.slice(c.target0+r.length);while(r&&a&&r[0]===a[0]){left+=r[0];r=r.slice(1);a=a.slice(1);}while(r&&a&&r.at(-1)===a.at(-1)){right=r.at(-1)+right;r=r.slice(0,-1);a=a.slice(0,-1);}return {left,r,a,right};}
function mapRef(o,c,dir,n){if(o.start0===null)return o;const s=dir==='+'?o.start0:n-o.end0,e=dir==='+'?o.end0:n-o.start0;return {...o,start0:s,end0:e,strand:o.strand==='+'?dir:(dir==='+'?'-':'+'),genomicStart1:c.start1!=null?c.start1+s:null,genomicEnd1:c.start1!=null?c.start1+e-1:null};}
function excluded(c,start,end){return (c.excludedRegions||[]).some(r=>r.start0<end&&r.end0>start);}
function alleleMasks(c){const ref=new Uint8Array(c.seq.length);for(const r of c.excludedRegions||[])ref.fill(1,r.start0,r.end0);const alt=new Uint8Array(c.seq.length+c.alt.length-c.ref.length),t=c.target0;alt.set(ref.slice(0,t));if(ref.slice(t,t+c.ref.length).some(Boolean))alt.fill(1,t,t+c.alt.length);alt.set(ref.slice(t+c.ref.length),t+c.alt.length);return {ref,alt};}
function makeSanger(c,p){const t=c.target0,z=t+c.ref.length,delta=c.alt.length-c.ref.length,seq=c.seq,F=[],R=[];
 for(let end=Math.max(p.minLen,t-p.maxProduct);end<=t-p.readBuffer;end++)for(let l=p.minLen;l<=p.maxLen;l++){const start=end-l;if(start<0||excluded(c,start,end))continue;const s=seq.slice(start,end);if(!basic(s,p))continue;const o=oligo(s,'F',start,end,'+',p);if(structureOK(o,p))F.push(o);}
 for(let start=z+p.readBuffer;start<=Math.min(seq.length-p.minLen,z+p.maxProduct);start++)for(let l=p.minLen;l<=p.maxLen;l++){if(start+l>seq.length||excluded(c,start,start+l))continue;const s=rc(seq.slice(start,start+l));if(!basic(s,p))continue;const o=oligo(s,'R',start,start+l,'-',p);if(structureOK(o,p))R.push(o);}
 // Deterministic spatial buckets retain diversity without unbounded pair enumeration.
 function keep(xs){const bins=new Map();for(const o of xs){const key=Math.floor(o.start0/20);if(!bins.has(key))bins.set(key,[]);bins.get(key).push(o);}return [...bins.values()].flatMap(x=>x.sort((a,b)=>merit(a,p)-merit(b,p)||a.start0-b.start0||a.length-b.length).slice(0,3));}
 const fset=keep(F),rset=keep(R),pairs=[];
 for(const f of fset)for(const r of rset){const n=r.end0-f.start0,na=n+delta;if(n<p.minProduct||n>p.maxProduct||na<p.minProduct||na>p.maxProduct||Math.abs(f.tm-r.tm)>p.maxTmDiff)continue;
 const cross=complement(f.sequence,r.sequence);if(cross.end>p.maxPairEnd)continue;
 const dF=t-f.end0,dR=r.start0-z;
 const score=merit(f,p)+merit(r,p)+Math.abs(f.tm-r.tm)+0.004*Math.abs(dF-dR)+0.001*Math.abs(n-(p.minProduct+p.maxProduct)/2);
 pairs.push({mode:'sanger',direction:'+',score,productRef:n,productAlt:na,readDistanceF:dF,readDistanceR:dR,crossEnd:cross.end,primers:[mapRef(f,c,'+',seq.length),mapRef(r,c,'+',seq.length)],status:'CANDIDATE_UNVALIDATED'});}
 pairs.sort((a,b)=>a.score-b.score||a.primers[0].start0-b.primers[0].start0||a.primers[1].start0-b.primers[1].start0);
 const chosen=[];for(const p0 of pairs){if(chosen.some(x=>x.primers[0].core===p0.primers[0].core&&x.primers[1].core===p0.primers[1].core))continue;chosen.push(p0);if(chosen.length>=p.numReturn)break;}
 return {candidates:chosen,diagnostics:{leftCandidates:F.length,rightCandidates:R.length,leftAfterBucket:fset.length,rightAfterBucket:rset.length,pairsPassing:pairs.length},warnings:['启发式候选筛选，不是 Primer3 输出；尚未验证 PCR 与读序效果。','仅避让本次提供并通过核对的排除区间；未提供的背景差异和重复注释仍属未知。','杂合 InDel 的混合模板 Sanger 峰图可能难以直接解释；本工具不分析 AB1。']};
}
function makeKasp(c,p){const masks=alleleMasks(c);const tr=trimChange(c),snp=tr.r.length===1&&tr.a.length===1,pure=(tr.r.length===0)!==(tr.a.length===0);
 assert(snp||(pure&&Math.max(tr.r.length,tr.a.length)<=30),'KASP 初筛仅支持 SNP 或去除公共前后缀后 ≤30 bp 的纯 InDel；复杂替换须另行设计。');
 let all=[],diag={alleleRef:0,alleleAlt:0,common:0,tripletsPassing:0};
 for(const dir of p.direction==='both'?['+','-']:[p.direction]){let left=tr.left,r=tr.r,a=tr.a,right=tr.right;if(dir==='-'){left=rc(tr.right);r=rc(tr.r);a=rc(tr.a);right=rc(tr.left);}const rs=left+r+right,as=left+a+right;
 const lists=[];
 for(const [allele,s,change,other] of [['REF',rs,r,as],['ALT',as,a,rs]]){const xs=[],tails=allele===p.famAllele?['FAM',TAILS.FAM]:['HEX',TAILS.HEX];const ends=snp?[left.length+1]:Array.from({length:change.length+4},(_,i)=>left.length+i+1);
  for(const end of ends)for(let len=p.minLen;len<=p.maxLen;len++){const st=end-len;if(st<0||end>s.length)continue;const mask=allele==='REF'?masks.ref:masks.alt;const ms=dir==='+'?st:mask.length-end,me=dir==='+'?end:mask.length-st;if(mask.slice(ms,me).some(Boolean))continue;
   // Allele-specific primer must include the change or straddle the deletion junction.
   if(change.length===0&&!(st<left.length&&end>left.length))continue;
   if(change.length>0&&!(st<left.length+change.length&&end>left.length))continue;
   const core=s.slice(st,end);if(!basic(core,p))continue;
   const oppositeSuffix=suffixMatch(core,other);if(oppositeSuffix>=p.indelSuffix||other.includes(core))continue;
   const o=oligo(core,'AS_'+allele,st,end,'+',p,tails[1],allele);o.channel=tails[0];o.oppositeSuffix=oppositeSuffix;o.orientation=dir;
   if(allele==='ALT'){o.haplotypeStart0=st;o.haplotypeEnd0=end;}
   if(structureOK(o,p))xs.push(o);
  }
  xs.sort((x,y)=>merit(x,p)-merit(y,p)||x.start0-y.start0||x.length-y.length);lists.push(xs.slice(0,16));diag[allele==='REF'?'alleleRef':'alleleAlt']+=xs.length;
 }
 if(!lists[0].length||!lists[1].length)continue;
 const commons=[];
 for(let offset=8;offset<Math.min(p.maxProduct,right.length-p.minLen+1);offset++)for(let len=p.minLen;len<=p.maxLen;len++){if(offset+len>right.length)continue;const qs=left.length+r.length+offset,qe=qs+len;const ms=dir==='+'?qs:rs.length-qe,me=dir==='+'?qe:rs.length-qs;if(masks.ref.slice(ms,me).some(Boolean))continue;const core=rc(right.slice(offset,offset+len));if(!basic(core,p))continue;const o=oligo(core,'COMMON',left.length+r.length+offset,left.length+r.length+offset+len,'-',p);o.rightOffset=offset;if(structureOK(o,p))commons.push(o);}
 commons.sort((x,y)=>merit(x,p)-merit(y,p)||x.start0-y.start0||x.length-y.length);diag.common+=commons.length;
 for(const x of lists[0])for(const y of lists[1]){if(Math.abs(x.tm-y.tm)>p.maxTmDiff)continue;const xy=complement(x.sequence,y.sequence);if(xy.end>p.maxPairEnd)continue;
  for(const z of commons){const nr=z.end0-x.start0,na=left.length+a.length+z.rightOffset+z.length-y.start0;
   if(nr<p.minProduct||nr>p.maxProduct||na<p.minProduct||na>p.maxProduct)continue;
   const diff=Math.max(x.tm,y.tm,z.tm)-Math.min(x.tm,y.tm,z.tm);if(diff>p.maxTmDiff)continue;
   const cross=Math.max(xy.end,complement(x.sequence,z.sequence).end,complement(y.sequence,z.sequence).end);if(cross>p.maxPairEnd)continue;
   const score=merit(x,p)+merit(y,p)+merit(z,p)+diff+0.012*Math.abs((nr+na)/2-(p.minProduct+p.maxProduct)/2);
   const yr={...y,start0:null,end0:null,genomicStart1:null,genomicEnd1:null,strand:dir,coordinateNote:'ALT 构建单倍型坐标；跨 InDel 时不强行映射到连续参考坐标。'};
   all.push({mode:'kasp',direction:dir,score,productRef:nr,productAlt:na,productIncludesTail:false,crossEnd:cross,primers:[mapRef(x,c,dir,rs.length),yr,mapRef(z,c,dir,rs.length)],status:snp?'CANDIDATE_UNVALIDATED':'INDEL_REQUIRES_EXPERT_REVIEW',variantClass:snp?'SNP':'INDEL'});diag.tripletsPassing++;
  }
 }
 }
 all.sort((x,y)=>x.score-y.score||x.direction.localeCompare(y.direction));
 const unique=[];const used=new Set();for(const x of all){const key=x.primers.map(o=>o.sequence).join('|');if(!used.has(key)){unique.push(x);used.add(key);}if(unique.length>=p.numReturn)break;}
 return {candidates:unique,diagnostics:diag,warnings:['KASP 为自定义候选初筛，不是 LGC 专有设计器，也未经厂商或湿实验验证。','REF/ALT 与 FAM/HEX 的绑定由本次参数决定；荧光通道本身不表示野生型、突变型或杂合。','Tm 只按基因组结合区计算；自互补与交叉互补初筛使用包含尾序的完整寡核苷酸。','3′ 后缀筛选不等于等位区分能力。未通过湿实验前不可标记为 validated。',...(snp?[]:['InDel 候选采用断点/插缺序列启发式；本版本不预测错配延伸、滑移或分型聚类，须人工复核。'])]};
}
function design(context,mode,params={}){assert(['sanger','kasp'].includes(mode),'设计模式无效。');const p=validateParams(params,mode),c=validateContext(context);const result=mode==='sanger'?makeSanger(c,p):makeKasp(c,p);result.diagnostics.excludedRegions=c.excludedRegions.length;result.diagnostics.background=c.backgroundAudit||{status:'NOT_PROVIDED'};return {version:VERSION,algorithm:'RPS_ENUMERATION_NN_v1',createdAt:new Date().toISOString(),mode,context:{id:c.id||'target',chrom:c.chrom||null,pos:c.pos||null,ref:c.ref,alt:c.alt,target0:c.target0,start1:c.start1??null,templateLength:c.seq.length,source:c.source||'USER_SEQUENCE'},parameters:p,...result};}
function parseBracket(text,id='sequence_target'){let t=text.split(/\r?\n/).filter(x=>!x.startsWith('>')).join('').replace(/\s/g,'').toUpperCase();const m=t.match(/^([ACGTRYSWKMBDHVN]*)\[([ACGT]+)\/([ACGT]+)\]([ACGTRYSWKMBDHVN]*)$/);assert(m,'请使用一处 [REF/ALT] 标记，例如 ACGT[AC/A]TGCA；仅支持保留 anchor 的序列等位。');return {id,seq:m[1]+m[2]+m[4],target0:m[1].length,ref:m[2],alt:m[3],source:'USER_SEQUENCE',start1:null};}
function parseFai(text){const out=[];for(const line of text.replace(/^\uFEFF/,'').split(/\r?\n/)){if(!line.trim())continue;const f=line.split('\t');assert(f.length>=5,'FAI 需要至少 5 个制表符分隔字段。');const [name,...nums]=f;const [length,offset,lineBases,lineWidth]=nums.map(Number);assert(name&&[length,offset,lineBases,lineWidth].every(Number.isSafeInteger)&&length>0&&offset>=0&&lineBases>0&&lineWidth>=lineBases&&lineWidth-lineBases<=2,'FAI 字段无效；仅支持规则换行、未压缩 FASTA。');out.push({name,length,offset,lineBases,lineWidth});}assert(out.length>0,'FAI 为空。');assert(new Set(out.map(x=>x.name)).size===out.length,'FAI 中存在重复 contig 名。');return out;}
function alias(s){const m=String(s).match(/^(?:chr)?0*(\d+)$/i);return m?String(Number(m[1])):/^chr(?:un|sy)$/i.test(String(s))?String(s).toLowerCase():String(s);}
function resolveContig(entries,name){const exact=entries.find(x=>x.name===name);if(exact)return exact;const matches=entries.filter(x=>alias(x.name)===alias(name));assert(matches.length===1,matches.length?'染色体别名不唯一。':'参考中找不到染色体 '+name);return matches[0];}
async function fetchSeq(file,e,start0,end0){assert(Number.isSafeInteger(start0)&&Number.isSafeInteger(end0)&&start0>=0&&end0>start0&&end0<=e.length,'提取区间越界。');const byte=i=>e.offset+Math.floor(i/e.lineBases)*e.lineWidth+i%e.lineBases;const lo=byte(start0),hi=byte(end0-1)+1;assert(hi<=file.size,'FAI 与 FASTA 大小不匹配。');const s=(await file.slice(lo,hi).text()).replace(/[\r\n]/g,'').toUpperCase();assert(s.length===end0-start0&&/^[ACGTRYSWKMBDHVN]+$/.test(s),'FAI/FASTA 不匹配，或 FASTA 换行不规则。');return s;}
function parseDelimited(text){text=String(text).replace(/^\uFEFF/,'');const rows=[];let row=[],field='',quoted=false;const csv=text.split(/\r?\n/).find(l=>l.trim()&&!l.startsWith('##'))?.includes(',')&&!text.includes('\t');const sep=csv?',':'\t';
 if(!csv&&!text.includes('\t'))text=text.split(/\r?\n/).map(l=>l.trim().replace(/ +/g,'\t')).join('\n');
 for(let i=0;i<=text.length;i++){const ch=text[i]??'\n';if(ch==='"'){if(quoted&&text[i+1]==='"'){field+='"';i++;}else quoted=!quoted;}else if(!quoted&&(ch===sep||ch==='\n')){row.push(field.replace(/\r$/,''));field='';if(ch==='\n'){if(row.some(x=>x.trim()))rows.push(row);row=[];}}else field+=ch;}assert(!quoted,'CSV 引号不闭合。');
 const clean=rows.filter(x=>!x[0].startsWith('##'));assert(clean.length,'输入为空。');let header=clean[0].map(x=>x.replace(/^#/,'').toLowerCase().trim());const hasHeader=header.includes('pos')&&header.includes('ref')&&header.includes('alt');let vcf=header.includes('qual')||header.includes('filter');let data=hasHeader?clean.slice(1):clean;let ix=hasHeader?Object.fromEntries(header.map((h,i)=>[h,i])):{id:0,chrom:1,pos:2,ref:3,alt:4};assert('chrom'in ix&&'pos'in ix&&'ref'in ix&&'alt'in ix,'表头须包含 id、chrom、pos、ref、alt；也可使用标准 VCF。');
 const result=data.filter(x=>!x[0].startsWith('#')).map((x,k)=>{assert(x.length>=5,`第 ${k+1} 行少于 5 列。`);assert((x[ix.chrom]||'').trim(),`第 ${k+1} 行染色体为空。`);const pos=Number(x[ix.pos]),ref=(x[ix.ref]||'').toUpperCase(),alt=(x[ix.alt]||'').toUpperCase();assert(Number.isSafeInteger(pos)&&pos>0,`第 ${k+1} 行 POS 无效。`);assert(/^[ACGT]+$/.test(ref)&&/^[ACGT]+$/.test(alt),`第 ${k+1} 行 REF/ALT 无效；不接受多 ALT、符号等位或缺失 anchor。`);assert(ref!==alt,`第 ${k+1} 行 REF 与 ALT 相同。`);return {id:(ix.id!=null?x[ix.id]:'')&&x[ix.id]!=='.'?x[ix.id]:`${x[ix.chrom]}_${pos}_${ref}_${alt}`,chrom:x[ix.chrom],pos,ref,alt};});assert(result.length>0&&result.length<=100,'每批须为 1–100 个位点。');assert(new Set(result.map(x=>x.id)).size===result.length,'位点 ID 必须唯一。');return result;
}
function demo(kind='snp'){let state=20260910;let seq='';for(let i=0;i<2200;i++){state^=state<<13;state^=state>>>17;state^=state<<5;seq+='ACGT'[(state>>>0)%4];}const t=1100,ref=kind==='del'?seq.slice(t,t+4):seq[t],alt=kind==='del'?ref[0]:'ACGT'[('ACGT'.indexOf(ref)+1)%4];return {id:'SYNTHETIC_'+kind.toUpperCase(),seq,target0:t,ref,alt,source:'SYNTHETIC_DEMO_NOT_MSU',start1:null};}
const api={VERSION,TAILS,DEFAULT,excluded,alleleMasks,rc,tm,gc,poly,complement,hairpin,suffixMatch,validateParams,validateContext,trimChange,design,parseBracket,parseFai,resolveContig,fetchSeq,parseDelimited,demo};
if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.RPS=api;
}
RPSEngineFactory(typeof globalThis!=='undefined'?globalThis:this);

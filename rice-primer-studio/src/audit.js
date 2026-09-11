/* Evidence audits: pure functions shared by the browser and Node regression suite.
 * Intervals are 0-based half-open unless a field explicitly ends in '1'.
 * These checks establish software consistency, not wet-lab performance.
 */
(function(root) {
'use strict';
const R = root.RPS || require('./engine.js');
const V = root.RPSReview || require('./review.js');
const IO = root.RPSIO || require('./io.js');
const RULE_VERSION = 'RPS_AUDIT_1_1_RC4_3';
function requireValue(ok, message) { if (!ok) throw Error(message); }
function normalizeEdit(v) {
  requireValue(typeof v.ref === 'string' && typeof v.alt === 'string', 'REF/ALT missing');
  const ref = v.ref.toUpperCase(), alt = v.alt.toUpperCase();
  requireValue(/^[ACGT]+$/.test(ref) && /^[ACGT]+$/.test(alt) && ref !== alt, 'Invalid sequence alleles');
  const origin0 = Number.isSafeInteger(v.pos) ? v.pos - 1 : v.target0;
  requireValue(Number.isSafeInteger(origin0) && origin0 >= 0, 'Invalid variant position');
  let prefix = 0, suffix = 0;
  while (prefix < ref.length && prefix < alt.length && ref[prefix] === alt[prefix]) prefix++;
  while (suffix < ref.length-prefix && suffix < alt.length-prefix && ref.at(-1-suffix) === alt.at(-1-suffix)) suffix++;
  const changedRef = ref.slice(prefix, ref.length-suffix), changedAlt = alt.slice(prefix, alt.length-suffix);
  const start0 = origin0 + prefix, end0 = start0 + changedRef.length;
  return { recordStart0:origin0, recordEnd0:origin0+ref.length, start0, end0,
    changedRef, changedAlt, prefix, suffix, insertionBreakpoint0:changedRef.length===0?start0:null,
    kind:!changedRef.length?'INSERTION':!changedAlt.length?'DELETION':changedRef.length===1&&changedAlt.length===1?'SNP':'REPLACEMENT',
    normalization:'COMMON_PREFIX_SUFFIX_TRIM_ONLY_NOT_LEFT_NORMALIZATION' };
}
function intervalUnion(regions) {
  const sorted=(regions||[]).map(r=>({start0:r.start0,end0:r.end0}));
  for (const r of sorted) requireValue(Number.isSafeInteger(r.start0)&&Number.isSafeInteger(r.end0)&&r.start0>=0&&r.end0>r.start0,'Invalid interval');
  sorted.sort((a,b)=>a.start0-b.start0||a.end0-b.end0);
  const out=[];
  for(const r of sorted){const last=out.at(-1);if(last&&r.start0<=last.end0)last.end0=Math.max(last.end0,r.end0);else out.push({...r});}
  return out;
}
function overlapBp(start0,end0,regions){
  return intervalUnion(regions).reduce((n,r)=>n+Math.max(0,Math.min(end0,r.end0)-Math.max(start0,r.start0)),0);
}
function sequenceStats(sequence) {
  const seq=sequence.toUpperCase();let known=0,gc=0,run=0,best=0,last='',pairs=0;const counts={};
  for(let i=0;i<seq.length;i++) {
    const b=seq[i],valid='ACGT'.includes(b);
    if(valid){known++;if(b==='G'||b==='C')gc++;run=b===last?run+1:1;best=Math.max(best,run);last=b;}
    else {run=0;last='';}
    if(i&&valid&&'ACGT'.includes(seq[i-1])){const k=seq[i-1]+b;counts[k]=(counts[k]||0)+1;pairs++;}
  }
  let entropy=0;for(const n of Object.values(counts)){const p=n/pairs;entropy-=p*Math.log2(p);}
  return {length:seq.length,knownBases:known,ambiguousBases:seq.length-known,knownFraction:seq.length?known/seq.length:null,
    gc:known?100*gc/known:null,maxHomopolymer:best,dinucleotidePairs:pairs,dinucleotideEntropy:pairs?entropy:null};
}
function gcPeak(seq,window=50) {
  let best=null;
  for(let i=0;i+window<=seq.length;i++) {
    const s=sequenceStats(seq.slice(i,i+window));
    // No GC peak can be called from a window containing N/IUPAC ambiguity.
    if(s.knownBases!==window)continue;
    if(!best||s.gc>best.gc)best={start0:i,end0:i+window,gc:s.gc};
  }
  return best;
}
function ampliconAudit(context,candidate,maxReadSpan=700) {
  const a=V.amplicons(context,candidate),ref=sequenceStats(a.ref),alt=sequenceStats(a.alt);
  const output={ref:{...ref,peak50:gcPeak(a.ref)},alt:{...alt,peak50:gcPeak(a.alt)},includesTail:false};
  if(candidate.mode==='sanger') {
    const f=candidate.primers.find(x=>x.role==='F'),r=candidate.primers.find(x=>x.role==='R');
    const dF=context.target0-f.end0,dR=r.start0-(context.target0+context.ref.length);
    output.readReach={planningLimitBp:maxReadSpan,planningThresholdNotMeasuredReadLength:true,
      forwardFarthestRefBp:dF+context.ref.length,reverseFarthestRefBp:dR+context.ref.length,
      forwardFarthestAltBp:dF+context.alt.length,reverseFarthestAltBp:dR+context.alt.length,
      oneDirectionWithinPlan:Math.min(dF,dR)+Math.max(context.ref.length,context.alt.length)<=maxReadSpan,
      bothDirectionsWithinPlan:Math.max(dF,dR)+Math.max(context.ref.length,context.alt.length)<=maxReadSpan};
    output.readReach.status=output.readReach.oneDirectionWithinPlan?'WITHIN_USER_READ_PLAN':'INTERNAL_SEQUENCING_PRIMER_REVIEW';
  }
  return output;
}
function repeatMetrics(context,oligo,candidate) {
  const regions=intervalUnion(context.repeatRegions||[]),repeatAt=i=>regions.some(r=>r.start0<=i&&i<r.end0);
  if(Number.isSafeInteger(oligo.start0)&&Number.isSafeInteger(oligo.end0)) {
    const lo=oligo.strand==='+'?Math.max(oligo.start0,oligo.end0-8):oligo.start0;
    const hi=oligo.strand==='+'?oligo.end0:Math.min(oligo.end0,oligo.start0+8);
    return {repeatOverlapBp:overlapBp(oligo.start0,oligo.end0,regions),repeatOverlap3prime8Bp:overlapBp(lo,hi,regions),
      repeatProjectionUnknownBp:0,repeatProjection:'REFERENCE_POSITION_ANNOTATION'};
  }
  if(oligo.role!=='AS_ALT'||!Number.isSafeInteger(oligo.haplotypeStart0))return {repeatOverlapBp:null,repeatOverlap3prime8Bp:null,repeatProjectionUnknownBp:oligo.core.length,repeatProjection:'UNKNOWN'};
  // Convert each ALT oligo position to forward ALT-haplotype coordinates, then REF positions.
  const edit=normalizeEdit({ref:context.ref,alt:context.alt,target0:context.target0});
  const delta=context.alt.length-context.ref.length,altLength=context.seq.length+delta;
  let overlap=0,three=0,unknown=0;
  for(let j=0;j<oligo.core.length;j++) {
    const at=candidate.direction==='+'?oligo.haplotypeStart0+j:altLength-1-(oligo.haplotypeStart0+j);
    let refAt=null;
    if(at<edit.start0)refAt=at;
    else if(at>=edit.start0+edit.changedAlt.length)refAt=at-delta;
    else if(edit.changedRef.length===edit.changedAlt.length)refAt=at;
    if(refAt===null){unknown++;continue;}
    if(repeatAt(refAt)){overlap++;if(j>=oligo.core.length-8)three++;}
  }
  return {repeatOverlapBp:overlap,repeatOverlap3prime8Bp:three,repeatProjectionUnknownBp:unknown,
    repeatProjection:unknown?'PARTIAL_REFERENCE_POSITION_PROJECTION':'REFERENCE_POSITION_PROJECTION'};
}
// Recheck saved design constraints independently from candidate generation.
// Tm is recomputed only for our own named NN engine; native Tm is range-checked,
// not relabelled as numerically validated without the native implementation.
function parameterEvidence(result,candidate) {
  const errors=[],rows=[];let tmEvidence='NOT_ASSESSED',pairTmSpan=null,crossEnd=null;
  const need=(ok,msg)=>{if(!ok)errors.push(msg);};
  try {
    const p=R.validateParams(result.parameters,result.mode);
    const names=result.mode==='sanger'?['F','R']:['AS_REF','AS_ALT','COMMON'];
    requireValue(Array.isArray(candidate.primers)&&candidate.primers.length===names.length&&
      new Set(candidate.primers.map(o=>o.role)).size===names.length&&
      names.every(role=>candidate.primers.some(o=>o.role===role)),'Parameter audit requires unique expected roles');
    const oligos=names.map(role=>candidate.primers.find(o=>o.role===role));
    const isEnumeration=['RPS_ENUMERATION_NN_v1','RPS_ENUMERATION_PLUS_PRIMER3_THERMO_REVIEW_v1'].includes(result.algorithm);
    tmEvidence=isEnumeration?'RPS_NN_RECOMPUTED':'REPORTED_TM_RANGE_ONLY_NOT_NATIVE_NUMERICAL_VALIDATION';
    need(Number.isFinite(candidate.score),'CANDIDATE_SCORE_INVALID');
    for(const o of oligos) {
      requireValue(typeof o.core==='string'&&/^[ACGT]{16,32}$/.test(o.core)&&
        typeof o.sequence==='string'&&/^[ACGT]{16,60}$/.test(o.sequence),'Invalid oligo for parameter audit');
      const g=R.gc(o.core),poly=R.poly(o.core),self=R.complement(o.sequence,o.sequence).end,hairpin=R.hairpin(o.sequence);
      need(o.core.length>=p.minLen&&o.core.length<=p.maxLen,'PRIMER_LENGTH_OUT_OF_RANGE_'+o.role);
      need(g>=p.gcMin-1e-7&&g<=p.gcMax+1e-7,'PRIMER_GC_OUT_OF_RANGE_'+o.role);
      need(poly<=p.maxPoly,'PRIMER_HOMOPOLYMER_LIMIT_'+o.role);
      need(Number.isFinite(o.tm)&&o.tm>=p.tmMin-1e-7&&o.tm<=p.tmMax+1e-7,'PRIMER_TM_OUT_OF_RANGE_'+o.role);
      const expectedTm=isEnumeration?R.tm(o.core,p):null;
      if(isEnumeration)need(Number.isFinite(o.tm)&&Math.abs(o.tm-expectedTm)<=1e-6,'PRIMER_TM_RECOMPUTATION_MISMATCH_'+o.role);
      // Some legacy native results omit poly; derive it rather than accepting
      // a missing value as zero. When supplied, it must be correct.
      if(o.poly!=null)need(o.poly===poly,'PRIMER_HOMOPOLYMER_VALUE_MISMATCH_'+o.role);
      need(o.selfEnd===self,'PRIMER_SELF_END_VALUE_MISMATCH_'+o.role);
      need(o.hairpin===hairpin,'PRIMER_HAIRPIN_VALUE_MISMATCH_'+o.role);
      need(self<=p.maxSelfEnd,'PRIMER_SELF_END_LIMIT_'+o.role);
      need(hairpin<=p.maxHairpin,'PRIMER_HAIRPIN_LIMIT_'+o.role);
      rows.push({role:o.role,coreLength:o.core.length,gc:g,homopolymer:poly,selfEnd:self,hairpin,reportedTm:o.tm,expectedTm});
    }
    if(oligos.every(o=>Number.isFinite(o.tm))) {
      pairTmSpan=Math.max(...oligos.map(o=>o.tm))-Math.min(...oligos.map(o=>o.tm));
      need(pairTmSpan<=p.maxTmDiff+1e-7,'PAIR_TM_DIFFERENCE_LIMIT');
    }
    crossEnd=0;
    for(let i=0;i<oligos.length;i++)for(let j=i+1;j<oligos.length;j++)crossEnd=Math.max(crossEnd,R.complement(oligos[i].sequence,oligos[j].sequence).end);
    need(candidate.crossEnd===crossEnd,'PAIR_END_VALUE_MISMATCH');
    need(crossEnd<=p.maxPairEnd,'PAIR_END_LIMIT');
    if(result.mode==='sanger') {
      const ctx=R.validateContext(result.inputContext),[f,r]=oligos;
      need(candidate.readDistanceF===ctx.target0-f.end0&&
           candidate.readDistanceR===r.start0-ctx.target0-ctx.ref.length,'READ_DISTANCE_VALUE_MISMATCH');
    }
  } catch(e) {errors.push('PARAMETER_AUDIT_ERROR:'+e.message);}
  return {status:errors.length?'FAIL':'CONSISTENT',errors:[...new Set(errors)],rows,pairTmSpan,crossEnd,tmEvidence,
    ruleVersion:RULE_VERSION,meaning:'SAVED_DESIGN_CONSTRAINTS_NOT_EXPERIMENTAL_ACCEPTANCE'};
}
function candidateIntegrity(result,candidate) {
  const errors=[];const need=(ok,msg)=>{if(!ok)errors.push(msg);};
  try {
    const c=R.validateContext(result.inputContext),p=R.validateParams(result.parameters,result.mode);
    need(candidate.mode===result.mode,'MODE_MISMATCH');
    const expectedRoles=result.mode==='sanger'?['F','R']:['AS_ALT','AS_REF','COMMON'];
    need(Array.isArray(candidate.primers)&&candidate.primers.map(o=>o.role).sort().join(',')===expectedRoles.join(','),'PRIMER_ROLES_INVALID');
    need(['+','-'].includes(candidate.direction),'DIRECTION_INVALID');
    if(errors.length)return {status:'FAIL',errors};
    const orientedAlt=candidate.direction==='+'?c.seq.slice(0,c.target0)+c.alt+c.seq.slice(c.target0+c.ref.length):R.rc(c.seq.slice(0,c.target0)+c.alt+c.seq.slice(c.target0+c.ref.length));
    for(const o of candidate.primers) {
      need(typeof o.core==='string'&&/^[ACGT]{16,32}$/.test(o.core),'CORE_INVALID_'+o.role);
      need(o.sequence===o.tail+o.core&&o.length===o.core.length&&o.fullLength===o.sequence.length,'SEQUENCE_FIELDS_MISMATCH_'+o.role);
      const strand=result.mode==='sanger'?(o.role==='F'?'+':'-'):(o.role==='COMMON'?(candidate.direction==='+'?'-':'+'):candidate.direction);
      need(o.strand===strand,'ROLE_STRAND_MISMATCH_'+o.role);
      need(Number.isFinite(o.tm)&&Number.isFinite(o.gc)&&Math.abs(o.gc-R.gc(o.core))<1e-7,'PRIMER_STATS_INVALID_'+o.role);
      if(o.role==='AS_ALT') {
        need(Number.isSafeInteger(o.haplotypeStart0)&&Number.isSafeInteger(o.haplotypeEnd0)&&o.haplotypeStart0>=0&&o.haplotypeEnd0<=orientedAlt.length&&o.haplotypeEnd0-o.haplotypeStart0===o.core.length&&orientedAlt.slice(o.haplotypeStart0,o.haplotypeEnd0)===o.core,'ALT_BINDING_MISMATCH');
      } else {
        need(Number.isSafeInteger(o.start0)&&Number.isSafeInteger(o.end0)&&o.start0>=0&&o.end0<=c.seq.length&&o.end0-o.start0===o.core.length,'REF_INTERVAL_INVALID_'+o.role);
        const s=c.seq.slice(o.start0,o.end0);
        need((o.strand==='+'?s:R.rc(s))===o.core,'REF_BINDING_MISMATCH_'+o.role);
        need(!(c.excludedRegions||[]).some(m=>m.start0<o.end0&&m.end0>o.start0),'EXCLUSION_OVERLAP_'+o.role);
        if(c.start1!=null)need(o.genomicStart1===c.start1+o.start0&&o.genomicEnd1===c.start1+o.end0-1,'GENOMIC_COORDINATE_MISMATCH_'+o.role);
      }
      if(o.role==='F'||o.role==='R'||o.role==='COMMON')need(o.tail==='','UNEXPECTED_TAIL_'+o.role);
      else {
        const allele=o.role==='AS_REF'?'REF':'ALT',channel=allele===p.famAllele?'FAM':'HEX';
        need(o.allele===allele&&o.channel===channel&&o.tail===R.TAILS[channel],'ALLELE_CHANNEL_TAIL_MISMATCH_'+o.role);
      }
    }
    if(result.mode==='sanger') {
      const f=candidate.primers.find(o=>o.role==='F'),r=candidate.primers.find(o=>o.role==='R');
      need(f.end0<=c.target0-p.readBuffer&&r.start0>=c.target0+c.ref.length+p.readBuffer,'TARGET_NOT_BRACKETED_WITH_BUFFER');
    } else {
      const tr=R.trimChange(c),x=candidate.primers.find(o=>o.role==='AS_REF'),y=candidate.primers.find(o=>o.role==='AS_ALT');
      if(tr.r.length===1&&tr.a.length===1) {
        const expectedEnd=candidate.direction==='+'?tr.left.length+1:tr.right.length+1;
        const xend=candidate.direction==='+'?x.end0:c.seq.length-x.start0;
        const er=candidate.direction==='+'?tr.r:R.rc(tr.r),ea=candidate.direction==='+'?tr.a:R.rc(tr.a);
        need(xend===expectedEnd&&y.haplotypeEnd0===expectedEnd&&x.core.at(-1)===er&&y.core.at(-1)===ea,'KASP_SNP_TERMINAL_ANCHOR_MISMATCH');
      }
    }
    const a=V.amplicons(c,candidate);
    need(a.ref.length===candidate.productRef&&a.alt.length===candidate.productAlt,'PRODUCT_LENGTH_MISMATCH');
    need([a.ref.length,a.alt.length].every(n=>n>=p.minProduct&&n<=p.maxProduct),'PRODUCT_OUT_OF_RANGE');
    const first=candidate.primers.find(o=>o.role===(result.mode==='sanger'?'F':'AS_REF'));
    const last=candidate.primers.find(o=>o.role===(result.mode==='sanger'?'R':'COMMON'));
    need(a.ref.startsWith(first.core)&&a.ref.endsWith(R.rc(last.core)),'REF_PRODUCT_TERMINI_MISMATCH');
    const altFirst=candidate.primers.find(o=>o.role===(result.mode==='sanger'?'F':'AS_ALT'));
    need(a.alt.startsWith(altFirst.core)&&a.alt.endsWith(R.rc(last.core)),'ALT_PRODUCT_TERMINI_MISMATCH');
  } catch(e) { errors.push('AUDIT_ERROR:'+e.message); }
  const parameters=parameterEvidence(result,candidate);errors.push(...parameters.errors);
  return {status:errors.length?'FAIL':'CONSISTENT',errors:[...new Set(errors)],ruleVersion:RULE_VERSION,parameters};
}
function scanEvidence(r,c) {
  const q=c.approxCheck||c.exactCheck;
  if(!q)return {status:'NOT_RUN',risk:0,truncated:false,detail:[],gaps:['SCAN_NOT_RUN'],errors:[],extraPairProducts:0};
  const errors=[],gaps=[],detail=[],by={};let risk=0,truncated=false,extraPairProducts=0;
  const approximate=!!c.approxCheck, source=approximate?'APPROXIMATE':'EXACT';
  const nonnegative=n=>Number.isSafeInteger(n)&&n>=0;
  const error=s=>errors.push(s);
  const object=x=>x!==null&&typeof x==='object'&&!Array.isArray(x);
  const limit=q.mismatchesAllowed;
  const validLimit=Number.isSafeInteger(limit)&&limit>=0&&limit<=2&&(!approximate?limit===0:true);
  if(!object(q))error('SCAN_PAYLOAD_INVALID');
  if(!validLimit)error('SCAN_MISMATCH_LIMIT_INVALID');
  if(q.gapsAllowed!==false||q.tailIncluded!==false||q.referenceScope!=='UPLOADED_FASTA_ONLY')error('SCAN_MODEL_METADATA_INVALID');
  const expectedMethod=approximate?'DISJOINT_SEED_EXHAUSTIVE_UNGAPPED_HAMMING_v1':'FULL_LENGTH_EXACT_MATCH_BOTH_STRANDS';
  if(q.method!==expectedMethod)error('SCAN_METHOD_UNRECOGNIZED');
  if(!Number.isSafeInteger(q.contigsScanned)||q.contigsScanned<=0||!Number.isSafeInteger(q.basesScanned)||q.basesScanned<=0)error('SCAN_REFERENCE_EXTENT_INVALID');
  if(!nonnegative(q.ambiguousBases)||q.ambiguousBases>q.basesScanned)error('SCAN_AMBIGUITY_COUNT_INVALID');
  if(q.completed!==true)gaps.push('SCAN_COMPLETION_NOT_CONFIRMED');
  if(r.importedReviewOnly)gaps.push('IMPORTED_SCAN_NOT_RECOMPUTED');
  if(!r.reference?.sha256||!q.reference?.sha256)gaps.push('SCAN_REFERENCE_HASH_MISSING');
  else if(r.reference.sha256!==q.reference.sha256)error('SCAN_REFERENCE_HASH_MISMATCH');
  const rows=Array.isArray(q.results)?q.results:[],oligos=Array.isArray(c.primers)?c.primers:[];
  if(!oligos.length||rows.length!==oligos.length||rows.some(x=>!object(x))||new Set(rows.filter(object).map(x=>x.id)).size!==rows.length)error('SCAN_ROLES_INCOMPLETE_OR_DUPLICATED');
  const ctx=r.inputContext||{},hasCoordinates=typeof ctx.chrom==='string'&&ctx.chrom.length>0&&Number.isSafeInteger(ctx.start1)&&ctx.start1>0;
  if(!hasCoordinates)gaps.push('NO_GENOMIC_BINDING_COORDINATES');
  for(const o of oligos) {
    const row=rows.find(x=>object(x)&&x.id===o.role);by[o.role]=row;
    if(!row||row.sequence!==o.core||!nonnegative(row.totalHits)||!Array.isArray(row.hits)) {error('SCAN_QUERY_OR_COUNT_INVALID_'+o.role);continue;}
    if(typeof row.positionsTruncated!=='boolean')error('SCAN_TRUNCATION_FLAG_INVALID_'+o.role);
    if(row.positionsTruncated){truncated=true;gaps.push('POSITIONS_TRUNCATED_'+o.role);}
    if(row.hits.length>row.totalHits||(!row.positionsTruncated&&row.hits.length!==row.totalHits))error('SCAN_POSITION_COUNT_MISMATCH_'+o.role);
    const buckets=row.byMismatches;
    const validBuckets=approximate&&validLimit&&Array.isArray(buckets)&&buckets.length===limit+1&&buckets.every(nonnegative)&&buckets.reduce((a,b)=>a+b,0)===row.totalHits;
    if(approximate&&!validBuckets)error('MISMATCH_COUNT_INVALID_'+o.role);
    const observed=Array(validLimit?limit+1:3).fill(0),seen=new Set();
    for(const h of row.hits) {
      if(!object(h)||typeof h.chrom!=='string'||!h.chrom||!nonnegative(h.start0)||!Number.isSafeInteger(h.end0)||h.end0-h.start0!==o.core.length||!['+','-'].includes(h.strand)) {error('SCAN_HIT_COORDINATE_INVALID_'+o.role);continue;}
      const n=h.mismatches??(approximate?null:0);
      if(!nonnegative(n)||!validLimit||n>limit){error('SCAN_HIT_MISMATCH_INVALID_'+o.role);continue;}
      observed[n]++;
      const key=JSON.stringify([h.chrom,h.start0,h.end0,h.strand]);
      if(seen.has(key))error('DUPLICATE_SCAN_POSITION_'+o.role);seen.add(key);
      if(approximate) {
        const positions=h.mismatchFrom3prime1;
        if(!Array.isArray(positions)||positions.length!==n||new Set(positions).size!==n||positions.some(x=>!Number.isSafeInteger(x)||x<1||x>o.core.length)||h.last5Mismatches!==positions.filter(x=>x<=5).length)error('SCAN_3PRIME_METADATA_INVALID_'+o.role);
      }
      // Recompute positions only where the exact input template supplies the sequence.
      // Outside that window this is an integrity check, not an independent re-alignment.
      if(hasCoordinates&&typeof ctx.seq==='string'&&h.chrom===ctx.chrom&&h.start0>=ctx.start1-1&&h.end0<=ctx.start1-1+ctx.seq.length) {
        const local=ctx.seq.slice(h.start0-ctx.start1+1,h.end0-ctx.start1+1).toUpperCase();
        if(!/^[ACGT]+$/.test(local)){error('SCAN_HIT_IN_AMBIGUOUS_LOCAL_WINDOW_'+o.role);continue;}
        const template=h.strand==='+'?local:R.rc(local),positions=[];
        for(let j=0;j<o.core.length;j++)if(template[j]!==o.core[j])positions.push(o.core.length-j);
        positions.sort((a,b)=>a-b);
        if(positions.length!==n||(approximate&&JSON.stringify(positions)!==JSON.stringify([...(Array.isArray(h.mismatchFrom3prime1)?h.mismatchFrom3prime1:[])].sort((a,b)=>a-b))))error('SCAN_LOCAL_SEQUENCE_MISMATCH_'+o.role);
      }
    }
    if(validBuckets&&observed.some((n,i)=>row.positionsTruncated?n>buckets[i]:n!==buckets[i]))error('SCAN_BUCKET_POSITION_MISMATCH_'+o.role);
    if(!hasCoordinates)continue;
    let expectedStart,expectedEnd,expectedMismatch=0;
    if(o.role!=='AS_ALT') {expectedStart=ctx.start1-1+o.start0;expectedEnd=ctx.start1-1+o.end0;}
    else if(c.variantClass==='SNP') {
      const n=ctx.seq.length;expectedStart=ctx.start1-1+(c.direction==='+'?o.haplotypeStart0:n-o.haplotypeEnd0);
      expectedEnd=expectedStart+o.core.length;expectedMismatch=1;
    } else {gaps.push('KASP_INDEL_ALT_REFERENCE_SCAN_UNRESOLVED');continue;}
    const intended=row.hits.filter(h=>object(h)&&h.chrom===ctx.chrom&&h.start0===expectedStart&&h.end0===expectedEnd&&h.strand===o.strand&&(h.mismatches??0)===expectedMismatch);
    // With >=1 allowed substitution a SNP AS_ALT must have its expected REF-locus
    // one-mismatch hit. Zero-mismatch scans legitimately cannot find that ALT site.
    const requiresIntended=o.role!=='AS_ALT'||(validLimit&&limit>=expectedMismatch);
    if(requiresIntended&&!intended.length) {
      if(row.positionsTruncated)gaps.push('INTENDED_HIT_NOT_RETAINED_'+o.role);
      else error('INTENDED_REFERENCE_HIT_MISSING_'+o.role);
    }
    if(intended.length>1)error('DUPLICATE_INTENDED_HIT_'+o.role);
    const extra=Math.max(0,row.totalHits-intended.length);risk+=extra;
    if(extra)detail.push(o.role+': '+extra+' additional possible binding site(s)');
  }
  const pairs=c.mode==='sanger'?[['F','R']]:[['AS_REF','COMMON'],['AS_ALT','COMMON']];
  // Never pass malformed records to the paired-product routine, and never call
  // skipped product enumeration complete evidence.
  if(hasCoordinates&&!truncated&&!errors.length)for(const [a,b] of pairs) {
    if(!by[a]?.hits||!by[b]?.hits)continue;
    if(by[a].hits.length>1000||by[b].hits.length>1000){gaps.push('PAIR_PRODUCT_ENUMERATION_LIMIT_'+a+'_'+b);continue;}
    const x=oligos.find(o=>o.role===a),y=oligos.find(o=>o.role===b);
    if(!Number.isSafeInteger(x.start0)&&c.variantClass!=='SNP')continue;
    const xs=a==='AS_ALT'?(c.direction==='+'?x.haplotypeStart0:ctx.seq.length-x.haplotypeEnd0):x.start0;
    const xe=xs+x.core.length,start1=ctx.start1+Math.min(xs,y.start0),end1=ctx.start1+Math.max(xe,y.end0)-1;
    const products=IO.exactProducts(by[a],by[b],40,3000);
    if(products.truncated)gaps.push('PAIR_PRODUCT_LIST_TRUNCATED');
    extraPairProducts+=products.products.filter(z=>z.chrom!==ctx.chrom||z.start1!==start1||z.end1!==end1).length;
  }
  if(extraPairProducts)detail.push(extraPairProducts+' unintended inward-facing product(s), bounded ungapped model');
  const status=errors.length?'INVALID':gaps.length?'INCOMPLETE':'COMPLETE_LIMITED';
  return {status,risk,truncated,detail,gaps:[...new Set(gaps)],errors:[...new Set(errors)],extraPairProducts,
    source,mismatchesAllowed:validLimit?limit:null,ambiguousBasesSkipped:q.ambiguousBases??null,
    scope:'UPLOADED_REFERENCE_ONLY_UNGAPPED_SUBSTITUTIONS_NO_PCR_SUCCESS_CLAIM'};
}
function assessReview(r,c,rules={},facts={}) {
  const fail=[];if(!c)return {pass:false,status:'NOT_READY_FOR_MANUAL_REVIEW',fail:['NO_SELECTED_CANDIDATE']};
  const integrity=candidateIntegrity(r,c),scan=scanEvidence(r,c);
  if(r.error)fail.push('RESULT_ERROR');
  if(r.importedReviewOnly)fail.push('IMPORTED_REVIEW_ONLY');
  if(!c.selected)fail.push('EXPLICIT_SELECTION_REQUIRED');
  if((r.context?.source||'').includes('SYNTHETIC')||(r.reference?.filename||'').startsWith('SYNTHETIC_'))fail.push('SYNTHETIC_NOT_FOR_ORDER');
  if(integrity.status!=='CONSISTENT')fail.push(...integrity.errors);
  if(scan.status==='INVALID')fail.push(...scan.errors);
  if(c.mode==='kasp'&&c.variantClass==='INDEL')fail.push('KASP_INDEL_REQUIRES_SEPARATE_EXPERT_SIGNOFF');
  if(rules.requireRefHash&&!/^[a-f0-9]{64}$/.test(r.reference?.sha256||''))fail.push('REFERENCE_SHA256_REQUIRED');
  if(rules.requireRepeat&&r.inputContext?.repeatAudit?.status!=='PROVIDED')fail.push('REPEAT_COORDINATE_SCOPE_NOT_ASSESSED');
  if(rules.requireNear&&(scan.source!=='APPROXIMATE'||scan.status!=='COMPLETE_LIMITED'||scan.mismatchesAllowed<(rules.requiredMismatches??2)))fail.push('MATCHING_COMPLETE_SCAN_REQUIRED');
  if(rules.requireNoConcern&&(facts.concerns||[]).length)fail.push('OBSERVED_CONCERNS_REMAIN');
  if(rules.requireBackground&&(!r.inputContext?.backgroundAudit?.status||r.inputContext.backgroundAudit.status==='NOT_PROVIDED'))fail.push('BACKGROUND_LIST_REQUIRED');
  // A numeric payload is not a completed acceptance review. This release has no accepted native branch.
  if(rules.requireThermo)fail.push('THERMODYNAMIC_ACCEPTANCE_NOT_VALIDATED_IN_THIS_RELEASE');
  return {pass:fail.length===0,status:fail.length?'NOT_READY_FOR_MANUAL_REVIEW':'READY_FOR_MANUAL_REVIEW',
    fail:[...new Set(fail)],integrity,scan,ruleVersion:RULE_VERSION,notOrderAuthorization:true};
}
function motifRelationship(v,start1,end1) {
  const e=normalizeEdit(v),a=start1-1,b=end1;
  const recordOverlaps=e.recordStart0<b&&e.recordEnd0>a;
  const changedReferenceOverlaps=e.start0<e.end0&&e.start0<b&&e.end0>a;
  const insertionInside=e.insertionBreakpoint0!==null&&a<e.insertionBreakpoint0&&e.insertionBreakpoint0<b;
  const insertionAtEdge=e.insertionBreakpoint0!==null&&(a===e.insertionBreakpoint0||b===e.insertionBreakpoint0);
  return {recordOverlaps,changedReferenceOverlaps,insertionInside,insertionAtEdge,
    editDirectlyIntersectsMotif:changedReferenceOverlaps||insertionInside,edit:e};
}
const api={RULE_VERSION,normalizeEdit,intervalUnion,overlapBp,sequenceStats,gcPeak,ampliconAudit,repeatMetrics,parameterEvidence,candidateIntegrity,scanEvidence,assessReview,motifRelationship};
if(typeof module!=='undefined'&&module.exports)module.exports=api;else root.RPSAudit=api;
})(typeof globalThis!=='undefined'?globalThis:this);

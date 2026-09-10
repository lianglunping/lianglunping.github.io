
(function(){
const oldCardHTML=cardHTML;
function trimmedAlleles(r){try{return RPS.trimChange(r.inputContext||r.context)}catch(e){return null}}
function kaspTerminalAudit(r,c){
 if(r.mode!=='kasp')return null;
 const tr=trimmedAlleles(r); const dir=c.direction||c.primers.find(x=>x.role==='AS_REF')?.orientation||'+';
 const ar=c.primers.find(x=>x.role==='AS_REF'), aa=c.primers.find(x=>x.role==='AS_ALT');
 if(!ar||!aa)return {status:'INCOMPLETE'};const consistent=RPSAudit.candidateIntegrity(r,c);if(consistent.status!=='CONSISTENT')return {status:'CANDIDATE_INTEGRITY_FAILURE',errors:consistent.errors};
 const snp=tr&&tr.r?.length===1&&tr.a?.length===1;
 if(!snp)return {status:'INDEL_EXPERT_REVIEW',variantClass:c.variantClass||'INDEL',dir,ref3:ar.core.at(-1),alt3:aa.core.at(-1)};
 const er=dir==='-'?RPS.rc(tr.r):tr.r, ea=dir==='-'?RPS.rc(tr.a):tr.a;
 return {status:(ar.core.at(-1)===er&&aa.core.at(-1)===ea&&er!==ea)?'TERMINAL_ANCHOR_CONFIRMED':'TERMINAL_ANCHOR_MISMATCH',dir,expectedRef3:er,expectedAlt3:ea,ref3:ar.core.at(-1),alt3:aa.core.at(-1),refMinus2:ar.core.at(-2),altMinus2:aa.core.at(-2)};
}
function replaceMinus2(seq,b){if(!b||seq.length<2)return seq;return seq.slice(0,-2)+b+seq.slice(-1)}
function kaspAuditHTML(r,c){
 const a=kaspTerminalAudit(r,c); if(!a)return '';
 const pass=a.status==='TERMINAL_ANCHOR_CONFIRMED';
 let h=`<div class="kasp-audit"><strong>KASP 3′ discrimination audit · <span class="${pass?'audit-pass':'audit-warn'}">${esc(a.status)}</span></strong><div class="hint">方向 ${esc(a.dir||'—')} · AS_REF 3′=${esc(a.ref3||'—')} · AS_ALT 3′=${esc(a.alt3||'—')}${a.expectedRef3?` · expected ${esc(a.expectedRef3)}/${esc(a.expectedAlt3)}`:''}</div>`;
 if($('kaspMismatchPreview')?.checked && c.variantClass==='SNP'){
   const ar=c.primers.find(x=>x.role==='AS_REF'),aa=c.primers.find(x=>x.role==='AS_ALT'); const rb=$('kaspRefMinus2')?.value||'',ab=$('kaspAltMinus2')?.value||'';
   const rp=replaceMinus2(ar.core,rb),ap=replaceMinus2(aa.core,ab);
   h+=`<div class="hint" style="margin-top:7px">人工 −2 预览（不参与候选/订购）：</div><div class="mono">AS_REF ${esc(rp)}${rb?` · −2 ${esc(ar.core.at(-2))}→${esc(rb)}`:' · unchanged'}<br>AS_ALT ${esc(ap)}${ab?` · −2 ${esc(aa.core.at(-2))}→${esc(ab)}`:' · unchanged'}</div>`;
 }
 return h+`<div class="hint" style="margin-top:6px">该模块只确认候选构造与手工预览，不预测错配延伸效率、簇间距或实际分型性能。</div></div>`;
}
cardHTML=function(r,c,ri,ci){return oldCardHTML(r,c,ri,ci).replace('</div>', '</div>') + kaspAuditHTML(r,c)};
function readiness(){
 const total=results.length, withCand=results.filter(r=>r.candidates?.length).length, refHash=results.filter(r=>!!r.reference?.sha256).length, repeat=results.filter(r=>r.inputContext?.repeatAudit?.status==='PROVIDED').length;
 let near=0,thermo=0,kasp=0,kaspAnchor=0;
 for(const r of results){const c=qcCandidateFor(r); if(c?.approxCheck&&RPSAudit.scanEvidence(r,c).status==='COMPLETE_LIMITED')near++; if(c?.thermoCheck)thermo++; if(r.mode==='kasp'&&c){kasp++; if(kaspTerminalAudit(r,c)?.status==='TERMINAL_ANCHOR_CONFIRMED')kaspAnchor++;}}
 const cells=[['候选位点',withCand,total],['reference hash',refHash,total],['repeat 已提供',repeat,total],['near-match 完整有限检查',near,total],['thermo 仅数值未验收',thermo,total],['KASP 3′锚定',kaspAnchor,kasp]];
 const box=$('readinessGrid'); if(box)box.innerHTML=cells.map(([n,a,b])=>`<div class="readiness-item"><strong>${a}/${b||0}</strong><small>${esc(n)}</small></div>`).join('');
}
const oldApply=applyQcFilters; applyQcFilters=function(){oldApply();readiness();};
function exportAuditManifest(){
 const payload={schema:'rps.audit/2',version:'1.1.0-rc3',distribution:globalThis.RPSDeployment?.snapshot()||null,createdAt:new Date().toISOString(),meaning:'review_audit_not_wetlab_validation',qcRules:{requiredMismatches:Number($('scanMismatches').value),requireRefHash:$('qcRequireRefHash')?.checked||false,requireRepeat:$('qcRequireRepeat')?.checked||false,requireNear:$('qcRequireNear')?.checked||false,requireNoConcern:$('qcRequireNoConcern')?.checked||false,requireBackground:$('qcRequireBackground')?.checked||false,requireThermo:$('qcRequireThermo')?.checked||false},reference:runManifest?.reference||null,gff:runManifest?.geneAnnotation||null,currentWorkspaceNotUsedForEvidence:refMeta(),targets:results.map(r=>{const a=adjudicateLocus(r),c=qcCandidateFor(r);return {target_id:r.context.id,mode:r.mode,locus_decision:a.locusTier,preferred_candidate:c?`C${String(r.candidates.indexOf(c)+1).padStart(2,'0')}`:null,selected_candidate:r.candidates.some(x=>x.selected)?`C${String(r.candidates.findIndex(x=>x.selected)+1).padStart(2,'0')}`:null,reference:r.reference||null,annotation:r.context?.annotationAudit||null,candidate_integrity:c?RPSAudit.candidateIntegrity(r,c):null,qc_release:qcReleaseAssessment(r,c),kasp_terminal_audit:c?kaspTerminalAudit(r,c):null,evidence:{repeat:r.inputContext?.repeatAudit?.status||'NOT_PROVIDED',background:r.inputContext?.backgroundAudit?.status||'NOT_PROVIDED',near_match:c?RPSAudit.scanEvidence(r,c).status:'NOT_RUN',thermo:c?.thermoCheck?'NUMERICAL_ONLY_UNACCEPTED':'NOT_RUN',wetlab:'NOT_VALIDATED'}}})};
 save('rice_primer_audit_manifest_v1.1.0-rc3.json',JSON.stringify(payload,null,2),'application/json');
}
$('exportAuditManifest').onclick=exportAuditManifest;
['kaspMismatchPreview','kaspRefMinus2','kaspAltMinus2'].forEach(id=>$(id)?.addEventListener('change',()=>{if(results.some(r=>r.mode==='kasp'))renderResults()}));
const oldExportKasp=exportOrder; exportOrder=function(mode){if(mode==='kasp'&&$('kaspMismatchPreview')?.checked)toast('人工 −2 错配仅为预览，不会写入订购草稿；订购导出仍使用原候选序列。'); return oldExportKasp(mode);};
readiness();
})();

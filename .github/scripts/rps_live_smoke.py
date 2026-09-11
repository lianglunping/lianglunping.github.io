#!/usr/bin/env python3
"""Test the deployed public HTTPS site. No local server or memory fallback.
All entered DNA and variant records are deterministic synthetic fixtures.
"""
from __future__ import annotations
import base64, csv, hashlib, io, json, os, platform, sys, time, traceback
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen
from playwright.sync_api import sync_playwright

SITE = 'https://lianglunping.github.io/rice-primer-studio/'
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'rps-live-evidence'
OUT.mkdir(exist_ok=True)
report = {'schema':'rps.live-test/1', 'harness':'PUBLIC_HTTPS_GITHUB_ACTIONS_CHROMIUM',
          'site':SITE, 'startedAtUtc':datetime.now(timezone.utc).isoformat(),
          'commit':os.getenv('GITHUB_SHA'), 'python':platform.python_version(),
          'playwright':'1.57.0', 'syntheticSeed':20260910, 'checks':[],
          'scientificValidation':'NOT_TESTED', 'status':'RUNNING'}
errors, requests = [], []
page = None

def check(name, actual, expected=True):
    ok = actual == expected
    report['checks'].append({'name':name, 'pass':ok, 'actual':actual, 'expected':expected})
    print(('PASS ' if ok else 'FAIL ') + name, flush=True)
    return ok

def get(path):
    url = SITE + path + ('&' if '?' in path else '?') + 'rps_check=' + str(time.time_ns())
    with urlopen(Request(url,headers={'User-Agent':'Rice-Primer-Studio-Live-Test','Cache-Control':'no-cache'}),timeout=30) as response:
        assert response.status == 200
        assert response.url.startswith(SITE)
        return response.read()

def idle():
    page.wait_for_function('!busy', timeout=120000)

def download(button, name):
    with page.expect_download(timeout=20000) as d:
        page.locator(button).click()
    dest = OUT / name
    d.value.save_as(str(dest))
    return dest.read_text(encoding='utf-8-sig')

def main():
    local = json.loads((ROOT/'rice-primer-studio/release.json').read_text())
    report['version'] = local['version']
    report['sourceSha256'] = local['sourceSha256']
    deadline = time.monotonic()+300
    last_error = ''
    while time.monotonic()<deadline:
        try:
            raw = get('release.json')
            live = json.loads(raw)
            if live == local:
                (OUT/'live_release.json').write_bytes(raw)
                break
            last_error = 'Live release manifest not yet equal to checkout'
        except Exception as exc:
            last_error = str(exc)
        time.sleep(10)
    else:
        raise RuntimeError('Live manifest unavailable: '+last_error)
    check('Live release.json equals checked-out manifest',live,local)
    for name, expected in local['files'].items():
        body = get(name)
        check('HTTPS file length '+name,len(body),expected['bytes'])
        check('HTTPS file SHA256 '+name,hashlib.sha256(body).hexdigest(),expected['sha256'])
    with urlopen('https://lianglunping.github.io/',timeout=30) as response:
        home=response.read()
        check('Existing blog homepage returns HTTP 200',response.status,200)
        check('Existing blog homepage content unchanged',hashlib.sha256(home).hexdigest(),hashlib.sha256((ROOT/'index.html').read_bytes()).hexdigest())
    global page
    with sync_playwright() as p:
        browser=p.chromium.launch(headless=True)
        report['browserVersion']=browser.version
        context=browser.new_context(viewport={'width':1440,'height':1080},accept_downloads=True)
        context.tracing.start(screenshots=True,snapshots=True,sources=True)
        page=context.new_page()
        page.on('pageerror',lambda e:errors.append(str(e)))
        page.on('request',lambda r:requests.append({'url':r.url,'method':r.method,'type':r.resource_type}))
        page.on('console',lambda m: print('BROWSER '+m.type+': '+m.text,flush=True) if m.type=='error' else None)
        try:
            response=page.goto(SITE+'?rps_live='+str(time.time_ns()),wait_until='networkidle',timeout=60000)
            check('Actual Chromium HTTPS navigation',response.status,200)
            page.wait_for_function('typeof RPSDeployment !== "undefined"',timeout=30000)
            check('Actual page remains on public HTTPS',page.url.startswith(SITE))
            check('Runtime version',page.evaluate('RPS.VERSION'),local['version'])
            check('Distribution mode',page.evaluate('RPSDeployment.snapshot().mode'),'PUBLIC_STATIC')
            check('Runtime source digest',page.evaluate('RPSDeployment.snapshot().sourceSha256'),local['sourceSha256'])
            check('All six runtime APIs loaded',page.evaluate('[RPS,RPSIO,RPSReview,RPSGFF,RPSAudit,RPSDeployment].every(x=>x&&typeof x==="object")'))
            check('12 scripts carry SRI',page.locator('script[src][integrity]').count(),12)
            check('Native engine option absent',page.locator('#designEngine option[value=native]').count(),0)
            check('Native backend control hidden',page.locator('#checkBackend').is_visible(),False)
            check('Native origin guard disabled',page.evaluate('nativeOrigin()'),False)
            start_requests=len(requests)
            for mode,kind in [('sanger','snp'),('sanger','del'),('kasp','snp'),('kasp','del')]:
                label=mode+'_'+kind
                page.locator('[data-mode="'+mode+'"]').click()
                page.locator('#demoSnp' if kind=='snp' else '#demoDel').click()
                page.locator('#runDesign').click()
                idle()
                r=page.evaluate('results[0]')
                check(label+' worker produced candidates',not r.get('error') and len(r.get('candidates',[]))>0)
                check(label+' every candidate passes consistency audit',page.evaluate('results[0].candidates.length>0 && results[0].candidates.every(c=>RPSAudit.candidateIntegrity(results[0],c).status==="CONSISTENT")'))
                check(label+' synthetic source clearly labelled',r['context']['source'],'SYNTHETIC_DEMO_NOT_MSU')
                check(label+' Sanger read plan visibility',page.locator('#sangerMaxReadSpan').is_visible(),mode=='sanger')
                page.locator('#autoSelectPreferred').click()
                check(label+' explicit one-group selection',page.evaluate('results[0].candidates.filter(c=>c.selected).length'),1)
                check(label+' synthetic not eligible for order',page.evaluate('!qcReleaseAssessment(results[0],results[0].candidates.find(c=>c.selected)).pass'))
                obj=json.loads(download('#exportJson','SYNTHETIC_'+label+'_run.json'))
                check(label+' actual JSON download schema',obj['schemaVersion'],'rps.run/5')
                check(label+' export records static distribution',obj['manifest']['distribution']['mode'],'PUBLIC_STATIC')
                rows=list(csv.DictReader(io.StringIO(download('#exportSangerOrder' if mode=='sanger' else '#exportKaspOrder','SYNTHETIC_'+label+'_order.csv'))))
                selected=next(c for c in obj['results'][0]['candidates'] if c.get('selected'))
                check(label+' ordering row count',len(rows),2 if mode=='sanger' else 3)
                field='sequence_5to3' if mode=='sanger' else 'full_sequence_5to3'
                check(label+' ordering sequences equal selected candidates',{x['role']:x[field] for x in rows},{x['role']:x['sequence'] for x in selected['primers']})
                check(label+' draft not release',all(x['qc_status']=='DRAFT_NOT_RELEASED' for x in rows))
                if mode=='kasp' and kind=='snp':
                    page.locator('#kaspMismatchPreview').locator('xpath=ancestor::details').locator('summary').click()
                    page.locator('#kaspMismatchPreview').check()
                    page.locator('#kaspRefMinus2').select_option('A')
                    preview_rows=list(csv.DictReader(io.StringIO(download('#exportKaspOrder','SYNTHETIC_mismatch_preview_order.csv'))))
                    check('KASP -2 preview does not rewrite order sequences',{x['role']:x[field] for x in preview_rows},{x['role']:x['sequence'] for x in selected['primers']})
                    page.locator('#kaspMismatchPreview').uncheck()
            check('No HTTP requests during synthetic design and export',[x for x in requests[start_requests:] if x['url'].startswith(('http:','https:'))],[])
            check('No JS runtime errors in all four example runs',errors,[])
            # Load a small synthetic reference through the real file chooser.
            demo=page.evaluate('RPS.demo("snp")')
            seq=demo['seq'];fa=('>Chr01\n'+'\n'.join(seq[i:i+60] for i in range(0,len(seq),60))+'\n').encode()
            page.locator('#toggleRef').click()
            page.locator('#fastaFile').set_input_files({'name':'SYNTHETIC_LIVE_REFERENCE.fa','mimeType':'text/plain','buffer':fa})
            page.wait_for_function('refFile!==null && refFile.name==="SYNTHETIC_LIVE_REFERENCE.fa"')
            page.locator('#buildIndex').click();idle()
            check('FASTA index contig length',page.evaluate('entries[0].length'),len(seq))
            check('FASTA SHA256 agrees with Python',page.evaluate('refSha'),hashlib.sha256(fa).hexdigest())
            page.locator('[data-input="single"]').click()
            page.locator('#targetId').fill('SYNTHETIC_LIVE_SNP')
            page.locator('#chrom').fill('Chr01')
            page.locator('#pos').fill(str(demo['target0']+1))
            page.locator('#ref').fill('')
            page.locator('#refFetchLen').fill('1')
            page.locator('#fetchRefBase').click()
            page.wait_for_function('document.getElementById("ref").value!==""')
            check('REF read from selected FASTA',page.locator('#ref').input_value(),demo['ref'])
            page.locator('#alt').fill(demo['alt'])
            for mode in ['sanger','kasp']:
                page.locator('[data-mode="'+mode+'"]').click()
                page.locator('#runDesign').click();idle()
                check(mode+' local FASTA target produces candidates',page.evaluate('results[0].candidates.length>0&&!results[0].error'))
                page.locator('[data-approx="0"]').click();idle()
                status=page.evaluate('RPSAudit.scanEvidence(results[0],results[0].candidates[0])')
                check(mode+' true Worker 0-2 substitution scan complete',status['status'],'COMPLETE_LIMITED')
                check(mode+' exact expected reference matches accounted for',status['errors'],[])
            page.locator('#autoSelectPreferred').click()
            audit=json.loads(download('#exportAuditManifest','SYNTHETIC_live_audit.json'))
            check('Audit manifest export schema',audit['schema'],'rps.audit/2')
            check('Audit manifest scan state',audit['targets'][0]['evidence']['near_match'],'COMPLETE_LIMITED')
            check('Synthetic FASTA also blocked from review eligibility',audit['targets'][0]['qc_release']['pass'],False)
            page.evaluate("document.documentElement.style.scrollBehavior='auto';window.scrollTo(0,0)")
            page.screenshot(path=str(OUT/'live_desktop.png'))
            page.locator('#results').scroll_into_view_if_needed()
            page.screenshot(path=str(OUT/'live_results.png'))
            check('Desktop no horizontal page overflow',page.evaluate('document.documentElement.scrollWidth<=innerWidth'))
            page.set_viewport_size({'width':390,'height':844})
            page.evaluate('window.scrollTo(0,0)')
            page.screenshot(path=str(OUT/'live_mobile.png'))
            check('Mobile viewport no horizontal page overflow',page.evaluate('document.documentElement.scrollWidth<=innerWidth'))
            page.set_viewport_size({'width':1440,'height':1080})
            wrong=next(b for b in 'ACGT' if b not in [demo['ref'],demo['alt']])
            page.locator('#ref').fill(wrong)
            page.locator('#runDesign').click();idle()
            check('REF mismatch rejects design',page.evaluate('!!results[0].error && results[0].candidates.length===0'))
            check('No uncaught errors after file, scan and export checks',errors,[])
            page.locator('#deploymentHelp').click()
            page.wait_for_url(SITE+'help.html')
            check('Help navigates to real hosted page',page.url,SITE+'help.html')
            page.locator('a[href="index.html"]').click()
            page.wait_for_function('typeof RPSDeployment!=="undefined"')
            check('Return from help reinitializes app',page.evaluate('RPS.VERSION'),local['version'])
        finally:
            report['pageErrors']=errors
            report['requests']=requests
            if page and not page.is_closed():
                page.screenshot(path=str(OUT/'last_page.png'))
            context.tracing.stop(path=str(OUT/'trace.zip'))
            browser.close()

if __name__=='__main__':
    try:
        main()
        report['status']='PASS' if all(x['pass'] for x in report['checks']) else 'FAIL'
    except Exception as exc:
        report.update(status='FAIL',failureType=type(exc).__name__,failure=str(exc))
        traceback.print_exc()
    finally:
        report['completedAtUtc']=datetime.now(timezone.utc).isoformat()
        report['checksTotal']=len(report['checks'])
        report['checksPassed']=sum(x['pass'] for x in report['checks'])
        (OUT/'report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
        # Compact summaries in logs are available even when artifact download is unavailable.
        summary={k:v for k,v in report.items() if k not in ['requests','checks']}
        summary['failedChecks']=[x for x in report['checks'] if not x['pass']]
        print('RPS_LIVE_REPORT '+json.dumps(summary,ensure_ascii=False),flush=True)
        if os.getenv('GITHUB_STEP_SUMMARY'):
            with open(os.environ['GITHUB_STEP_SUMMARY'],'a') as out:
                out.write('## Live HTTPS validation\n\n```json\n'+json.dumps(summary,ensure_ascii=False,indent=2)+'\n```\n')
    sys.exit(0 if report['status']=='PASS' else 1)

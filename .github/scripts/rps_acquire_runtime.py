"""Acquire public, pinned RPS dependencies. No user data, no deployment, no source compilation.
The output is a dependency payload, NOT a complete or scientifically validated application.
"""
from __future__ import annotations
import argparse, base64, csv, email, gzip, hashlib, importlib.metadata, io, json, os
import platform, random, re, shutil, struct, subprocess, sys, tarfile, tempfile
import urllib.parse, urllib.request, zipfile
from pathlib import Path, PurePosixPath
from packaging.tags import cpython_tags, compatible_tags, mac_platforms
from packaging.utils import parse_wheel_filename, canonicalize_name
from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
from packaging.markers import default_environment

PINS={'packaging':'25.0','numpy':'2.3.5','biopython':'1.86','fastapi':'0.128.2',
'uvicorn':'0.48.0','pydantic':'2.13.5','pydantic-core':'2.46.5','python-multipart':'0.0.29',
'starlette':'0.50.0','anyio':'4.13.0','idna':'3.17','typing-extensions':'4.16.0',
 'typing-inspection':'0.4.4','annotated-types':'0.8.0','annotated-doc':'0.0.4',
 'click':'8.1.8','h11':'0.16.0','colorama':'0.4.6'}
PBS='https://github.com/astral-sh/python-build-standalone/releases/download/20261013/'
TARGETS={
 'windows-x64':{'system':'Windows','machine':'amd64','python':'https://www.python.org/ftp/python/3.13.16/python-3.13.16-embed-amd64.zip','python_sha':'97dae5274cc54867065e8d5a3226e48c35017ed332a0fdb0e27d5b5821961297','tracy':'tracy-v0.9.1-windows-amd64.exe','tracy_sha':'a64e19f20044b50fb475657362985ac382cf9aac13386c98bf64f53660cf2484','blast':'x64-win64'},
 'macos-arm64':{'system':'Darwin','machine':'arm64','python':PBS+'cpython-3.13.16+20261013-aarch64-apple-darwin-install_only_stripped.tar.gz','python_sha':'16a5332fb0dcb77f163fc8f8cbccd95be1644cc53e68e5d2a79f50afd6680a0e','tracy':'tracy-v0.9.1-macos-arm64','tracy_sha':'aed8ef2d699a1ed63ec62bfa3697cf73446d5f6da2bd8a30619ec1e9d8d518c3','blast':'aarch64-macosx'},
 'linux-x64':{'system':'Linux','machine':'x86_64','python':PBS+'cpython-3.13.16+20261013-x86_64-unknown-linux-gnu-install_only_stripped.tar.gz','python_sha':'d0080766b6b74aeb81e1949ef0d9a20e6c17f7b160fb295a69b936ccb17660be','tracy':'tracy-v0.9.1-linux-amd64','tracy_sha':'df2dd5016f11aeb16d2df090e06d9fac713d4e74729b44e758165169161bc4d5','blast':'x64-linux'}}
HOSTS={'www.python.org','pypi.org','files.pythonhosted.org','ftp.ncbi.nlm.nih.gov',
 'github.com','raw.githubusercontent.com','release-assets.githubusercontent.com','objects.githubusercontent.com'}
SOURCES=[]

def digest(p,alg='sha256'):
 h=hashlib.new(alg)
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024**2),b''):h.update(b)
 return h.hexdigest()

def checked(url):
 u=urllib.parse.urlsplit(url)
 if u.scheme!='https' or u.hostname not in HOSTS or u.username or u.password or u.port not in (None,443):raise ValueError('Unapproved download host')
 return url

class Redirect(urllib.request.HTTPRedirectHandler):
 def redirect_request(self,req,fp,code,msg,headers,newurl):
  checked(newurl)
  return super().redirect_request(req,fp,code,msg,headers,newurl)

def request(url):
 return urllib.request.build_opener(Redirect()).open(urllib.request.Request(checked(url),headers={'User-Agent':'RPS-portable-runtime/1'}),timeout=120)

def fetch(url):
 with request(url) as r:
  checked(r.geturl());data=r.read(10*1024**2+1)
 if len(data)>10*1024**2:raise ValueError('Metadata too large')
 return data

def download(url,dest,expected,alg='sha256'):
 print('Download',dest.name,flush=True)
 dest.parent.mkdir(parents=True,exist_ok=True)
 with request(url) as r,dest.open('xb') as out:
  checked(r.geturl());total=0
  while b:=r.read(1024**2):
   total+=len(b)
   if total>1024**3:raise ValueError('Asset exceeds 1 GiB')
   out.write(b)
 if digest(dest,alg)!=expected:raise ValueError('Published digest mismatch: '+dest.name)
 SOURCES.append({'url':url,'filename':dest.name,'bytes':dest.stat().st_size,'sha256':digest(dest),'published_algorithm':alg,'published_digest':expected})
 return dest

def safe(name):
 p=PurePosixPath(name)
 if p.is_absolute() or '..' in p.parts or '\\' in name or ':' in name or '\x00' in name:raise ValueError('Unsafe archive member')
 return p

def install_wheel(file,site,runtime):
 with zipfile.ZipFile(file) as z:
  record=[x for x in z.namelist() if x.endswith('.dist-info/RECORD')]
  if len(record)!=1:raise ValueError('Missing wheel RECORD')
  inventory={row[0]:row[1:] for row in csv.reader(io.StringIO(z.read(record[0]).decode()))}
  for i in z.infolist():
   if i.is_dir():continue
   p=safe(i.filename);data=z.read(i)
   if p.suffix.lower() in ('.ttf','.otf','.woff','.woff2','.ttc'):raise ValueError('Fonts not included in this distribution')
   info=inventory.get(i.filename)
   if info is None:raise ValueError('Unlisted wheel member')
   hs,size=info
   if size and len(data)!=int(size):raise ValueError('RECORD length mismatch')
   if hs:
    alg,expected=hs.split('=',1)
    if alg not in ('sha256','sha384','sha512'):raise ValueError('Weak RECORD hash')
    if base64.urlsafe_b64encode(hashlib.new(alg,data).digest()).rstrip(b'=').decode()!=expected:raise ValueError('RECORD digest mismatch')
   elif i.filename!=record[0]:raise ValueError('Missing wheel hash')
   if p.parts[0].endswith('.data'):
    scheme={'purelib':site,'platlib':site,'scripts':runtime/'Scripts','data':runtime,'headers':runtime/'include'}
    dest=scheme[p.parts[1]].joinpath(*p.parts[2:])
   else:dest=site.joinpath(*p.parts)
   dest.parent.mkdir(parents=True,exist_ok=True)
   with dest.open('xb') as f:f.write(data)
  meta=[x for x in z.namelist() if x.endswith('.dist-info/METADATA')]
  m=email.message_from_bytes(z.read(meta[0]))
 return {'name':m['Name'],'version':m['Version'],'requires':m.get_all('Requires-Dist',[]),'requires_python':m.get('Requires-Python',''),'wheel':file.name,'sha256':digest(file)}

def install_blast(archive,dest):
 with tarfile.open(archive,'r:gz') as t:
  for m in t:
   p=safe(m.name);low=p.name.lower()
   if not m.isfile():continue
   binary=('bin' in p.parts and (low in ('blastn','blastp','makeblastdb','blastn.exe','blastp.exe','makeblastdb.exe') or low.endswith(('.dll','.dylib','.so'))))
   library=('lib' in p.parts and ('.so' in low or '.dylib' in low))
   license=('license' in low or 'copying' in low or low.startswith('readme'))
   if not (binary or library or license):continue
   out=dest/('bin' if binary else 'lib' if library else 'licenses')/p.name
   out.parent.mkdir(parents=True,exist_ok=True)
   with t.extractfile(m) as src,out.open('xb') as f:shutil.copyfileobj(src,f)
   if binary or library:out.chmod(0o755)

def abif(seq):
 n=len(seq);positions=[10*i+5 for i in range(n)];channels={b:[0]*(10*n+10) for b in 'ACGT'}
 for b,p in zip(seq,positions):
  for shift,height in [(-2,200),(-1,600),(0,1000),(1,600),(2,200)]:channels[b][p+shift]=height
 tags=[(b'PBAS',2,2,1,n,seq.encode()),(b'PLOC',2,4,2,n,struct.pack('>'+str(n)+'h',*positions)),(b'FWO_',1,2,1,4,b'GATC'),(b'PCON',2,2,1,n,bytes([35]*n))]
 for i,b in enumerate('GATC',9):
  a=channels[b];tags.append((b'DATA',i,4,2,len(a),struct.pack('>'+str(len(a))+'h',*a)))
 offset=34+28*len(tags);dirs=[];payload=[]
 for name,num,code,size,count,data in tags:
  loc=int.from_bytes(data.ljust(4,b'\0'),'big') if len(data)<=4 else offset
  dirs.append(struct.pack('>4sI2H4I',name,num,code,size,count,len(data),loc,0))
  if len(data)>4:payload.append(data);offset+=len(data)
 head=b'ABIF'+struct.pack('>H4sI2H3I',101,b'tdir',1,1023,28,len(tags),28*len(tags),34)+b'\0'*4
 return head+b''.join(dirs)+b''.join(payload)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--target',choices=TARGETS,required=True);a=ap.parse_args();cfg=TARGETS[a.target]
 if platform.system()!=cfg['system'] or platform.machine().lower()!=cfg['machine']:raise RuntimeError('Runner is not target-native')
 root=Path('rps-runtime-work').resolve();root.mkdir(exist_ok=False);downloads=root/'downloads';payload=root/'payload'/a.target;runtime=payload/'runtime/python';runtime.parent.mkdir(parents=True)
 evidence=Path('rps-runtime-output').resolve();evidence.mkdir(exist_ok=True)
 pyname=cfg['python'].rsplit('/',1)[1];archive=download(cfg['python'],downloads/pyname,cfg['python_sha'])
 if a.target=='windows-x64':
  runtime.mkdir()
  with zipfile.ZipFile(archive) as z:
   for n in z.namelist():safe(n)
   z.extractall(runtime)
  (runtime/'python313._pth').write_text('python313.zip\n.\nLib/site-packages\nimport site\n',encoding='ascii')
  executable=runtime/'python.exe';site=runtime/'Lib/site-packages';platforms=['win_amd64']
 else:
  with tarfile.open(archive,'r:gz') as t:t.extractall(runtime.parent,filter='data')
  executable=runtime/'bin/python3';site=runtime/'lib/python3.13/site-packages'
  platforms=list(mac_platforms((14,0),'arm64')) if a.target=='macos-arm64' else ['manylinux_2_'+str(v)+'_x86_64' for v in range(28,16,-1)]+['manylinux2014_x86_64']
 site.mkdir(parents=True,exist_ok=True)
 tags=list(cpython_tags((3,13),abis=['cp313'],platforms=platforms))+list(compatible_tags((3,13),interpreter='cp313',platforms=platforms));rank={tag:i for i,tag in reversed(list(enumerate(tags)))}
 packages=[]
 for name,version in PINS.items():
  metadata=json.loads(fetch('https://pypi.org/pypi/'+name+'/'+version+'/json'));options=[]
  for item in metadata['urls']:
   if item.get('packagetype')!='bdist_wheel' or item.get('yanked'):continue
   _,_,_,wtags=parse_wheel_filename(item['filename']);ranks=[rank[t] for t in wtags if t in rank]
   if ranks:options.append((min(ranks),item['filename'],item))
  if not options:raise RuntimeError('No pinned binary wheel for '+a.target+': '+name+'=='+version)
  item=min(options)[2];whl=download(item['url'],downloads/item['filename'],item['digests']['sha256'])
  packages.append(install_wheel(whl,site,runtime))
 env=default_environment();env.update(python_version='3.13',python_full_version='3.13.16',implementation_version='3.13.16',extra='')
 actual={canonicalize_name(x['name']):x['version'] for x in packages}
 for p in packages:
  if p['requires_python'] and '3.13.16' not in SpecifierSet(p['requires_python']):raise RuntimeError('Python requirement mismatch')
  for text in p['requires']:
   q=Requirement(text)
   if q.marker and not q.marker.evaluate(env):continue
   if q.url or q.extras or canonicalize_name(q.name) not in actual or actual[canonicalize_name(q.name)] not in q.specifier:raise RuntimeError('Unclosed runtime dependency: '+text)
 blasturl='https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-'+cfg['blast']+'.tar.gz'
 md5=re.search(r'\b[0-9a-fA-F]{32}\b',fetch(blasturl+'.md5').decode()).group(0).lower()
 blast=download(blasturl,downloads/blasturl.rsplit('/',1)[1],md5,'md5');install_blast(blast,payload/'tools/blast')
 tracy=payload/'tools/tracy'/('tracy.exe' if a.target=='windows-x64' else 'tracy')
 download('https://github.com/gear-genomics/tracy/releases/download/v0.9.1/'+cfg['tracy'],tracy,cfg['tracy_sha']);tracy.chmod(0o755)
 license_data=fetch('https://raw.githubusercontent.com/gear-genomics/tracy/1d450cee0a832f442042568d23b11cadc76693c0/LICENSE')
 (tracy.parent/'LICENSE').write_bytes(license_data)
 process_env=os.environ.copy();process_env.update(PYTHONNOUSERSITE='1',PYTHONDONTWRITEBYTECODE='1')
 for key in ('PYTHONPATH','PYTHONHOME','VIRTUAL_ENV'):process_env.pop(key,None)
 commands=[]
 def run(argv,cwd):
  p=subprocess.run(list(map(str,argv)),cwd=cwd,env=process_env,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=180)
  commands.append({'argv':[str(v) for v in argv],'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr})
  if p.returncode:raise RuntimeError('Native command failed: '+str(argv[0])+'\n'+p.stdout+'\n'+p.stderr)
  return p.stdout
 with tempfile.TemporaryDirectory(prefix='rps_native_') as tmp:
  tmp=Path(tmp)
  smoke="import sys,json,platform,importlib,importlib.metadata;mods=['numpy','Bio','Bio.Align','Bio.bgzf','fastapi','pydantic','pydantic_core','uvicorn','multipart','sqlite3','ssl','ctypes'];[importlib.import_module(x) for x in mods];from Bio.Seq import Seq;assert str(Seq('ATGAAATAA').translate())=='MK*';pins="+repr(PINS)+";assert platform.python_version()=='3.13.16';assert all(importlib.metadata.version(k)==v for k,v in pins.items());print(json.dumps({'python':platform.python_version(),'modules':mods,'pins':pins,'executable':sys.executable}))"
  run([executable,'-I','-B','-c',smoke],tmp)
  suffix='.exe' if a.target=='windows-x64' else '';bb=payload/'tools/blast/bin'
  for program in ('blastn','blastp','makeblastdb'):
   if '2.17.0' not in run([bb/(program+suffix),'-version'],tmp):raise RuntimeError('BLAST version mismatch')
  rng=random.Random(20261016)
  for typ,seq,program in [('nucl',''.join(rng.choices('ACGT',k=600)),'blastn'),('prot','M'+''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY',k=119)),'blastp')]:
   fa=tmp/(typ+'.fa');fa.write_text('>synthetic_self\n'+seq+'\n',encoding='ascii');gz=tmp/(typ+'.fa.gz');gz.write_bytes(gzip.compress(fa.read_bytes(),mtime=0));db=tmp/(typ+'_db')
   run([bb/('makeblastdb'+suffix),'-in',gz,'-dbtype',typ,'-out',db],tmp)
   cmd=[bb/(program+suffix),'-query',fa,'-db',db,'-outfmt','6 qseqid sseqid pident length']+(['-task','blastn','-dust','no'] if typ=='nucl' else ['-seg','no'])
   rows=[s.split('\t') for s in run(cmd,tmp).splitlines() if s]
   if not any(float(x[2])==100 and int(x[3])==len(seq) for x in rows):raise RuntimeError('Missing exact self hit')
  trace=tmp/'SYNTHETIC_ONLY.ab1';trace.write_bytes(abif(''.join(rng.choices('ACGT',k=180))))
  run([tracy,'--version'],tmp);result=tmp/'trace.json';run([tracy,'basecall','-f','json','-o',result,trace],tmp)
  if not result.is_file() or not json.loads(result.read_text()):raise RuntimeError('Missing Tracy JSON')
 # Relocate the same bytes, then repeat native imports and tool launches.
 relocated=payload.parent/('relocated path '+a.target);payload.rename(relocated);payload=relocated
 executable=payload/'runtime/python'/('python.exe' if a.target=='windows-x64' else 'bin/python3')
 run([executable,'-I','-B','-c',smoke],payload)
 suffix='.exe' if a.target=='windows-x64' else ''
 run([payload/'tools/blast/bin'/('blastn'+suffix),'-version'],payload)
 run([payload/'tools/tracy'/('tracy'+suffix),'--version'],payload)
 receipt={'status':'NATIVE_DEPENDENCY_SMOKE_PASSED','target':a.target,'host':platform.platform(),'python':'3.13.16','sources':SOURCES,'packages':packages,'commands':commands,'seed':20261016,'scope':'Dependencies only; application not included or validated in this job','offline_complete_machine_test':False}
 (payload/'RUNTIME_RECEIPT.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
 (evidence/(a.target+'-receipt.json')).write_text(json.dumps(receipt,indent=2),encoding='utf-8')
 for f in payload.rglob('*'):
  if f.is_file() and f.suffix.lower() in ('.ttf','.otf','.woff','.woff2','.ttc'):raise RuntimeError('Unexpected font in payload')
 target=evidence/('rps-runtime-'+a.target+'.tar.gz')
 with tarfile.open(target,'w:gz',compresslevel=6) as t:t.add(payload,arcname=a.target)
 (evidence/(target.name+'.sha256')).write_text(digest(target)+'  '+target.name+'\n',encoding='ascii')
 print(json.dumps({'status':receipt['status'],'target':a.target,'bytes':target.stat().st_size,'sha256':digest(target)}),flush=True)

if __name__=='__main__':
 try:main()
 except Exception as e:
  out=Path('rps-runtime-output');out.mkdir(exist_ok=True)
  (out/'FAILURE.json').write_text(json.dumps({'error':type(e).__name__+': '+str(e),'sources_acquired':SOURCES},indent=2),encoding='utf-8')
  raise

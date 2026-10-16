"""macOS-specific Tracy source build: no end-user Homebrew dependencies.
Only the dependency build machine uses Homebrew. The application payload retains
exact source/submodule commits, compiler commands and static input hashes.
"""
import hashlib, importlib.util, json, os, platform, shutil, subprocess, sys, traceback
from pathlib import Path
COMMIT='1d450cee0a832f442042568d23b11cadc76693c0'
ROOT=Path('rps-tracy-source-build').resolve()
OUT=Path('rps-runtime-output').resolve()
LOG=[]

def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def run(args,cwd=None,env=None,timeout=1200):
 print('RUN',args,flush=True)
 p=subprocess.run(list(map(str,args)),cwd=cwd,env=env,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=timeout)
 LOG.append({'argv':list(map(str,args)),'cwd':str(cwd),'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr})
 OUT.mkdir(exist_ok=True);(OUT/'macos-tracy-build-progress.json').write_text(json.dumps(LOG,indent=2),encoding='utf-8')
 if p.returncode:raise RuntimeError(p.stdout+'\n'+p.stderr)
 return p.stdout.strip()

def build():
 if platform.system()!='Darwin' or platform.machine()!='arm64':raise RuntimeError('Native macOS arm64 required')
 ROOT.mkdir(exist_ok=False)
 env=os.environ.copy();env.update(HOMEBREW_NO_AUTO_UPDATE='1',HOMEBREW_NO_INSTALL_CLEANUP='1',MACOSX_DEPLOYMENT_TARGET='14.0',CMAKE_POLICY_VERSION_MINIMUM='3.5')
 run(['brew','install','boost','autoconf','automake','cmake','xz','libdeflate'],env=env)
 formulas=json.loads(run(['brew','info','--json=v2','boost','xz','libdeflate','autoconf','automake','cmake'],env=env))
 pref={name:Path(run(['brew','--prefix',name],env=env)) for name in ('boost','xz','libdeflate')}
 source=ROOT/'tracy'
 run(['git','clone','--no-checkout','https://github.com/gear-genomics/tracy.git',source],env=env)
 run(['git','checkout','--detach',COMMIT],cwd=source,env=env)
 if run(['git','rev-parse','HEAD'],cwd=source)!=COMMIT:raise RuntimeError('Tracy source identity mismatch')
 run(['git','submodule','update','--init','--recursive'],cwd=source,env=env)
 submodules=run(['git','submodule','status','--recursive'],cwd=source)
 hts=source/'src/htslib'
 env['CPPFLAGS']=' '.join('-I'+str(pref[n]/'include') for n in ('xz','libdeflate'))
 env['LDFLAGS']=' '.join('-L'+str(pref[n]/'lib') for n in ('xz','libdeflate'))
 env['PKG_CONFIG_PATH']=os.pathsep.join(str(pref[n]/'lib/pkgconfig') for n in ('xz','libdeflate'))
 run(['autoreconf','-i'],cwd=hts,env=env)
 run(['./configure','--disable-s3','--disable-gcs','--disable-libcurl','--disable-plugins'],cwd=hts,env=env)
 run(['make','-j','3','libhts.a'],cwd=hts,env=env)
 sdsl=source/'src/xxsds';sdslbuild=ROOT/'sdsl-build';sdslout=source/'src/sdslLite'
 run(['cmake','-S',sdsl,'-B',sdslbuild,'-DCMAKE_POLICY_VERSION_MINIMUM=3.5','-DCMAKE_INSTALL_PREFIX='+str(sdslout),'-DSDSL_BUILD_TESTS=OFF','-DSDSL_BUILD_EXAMPLES=OFF','-DSDSL_BUILD_TUTORIAL=OFF','-DCMAKE_OSX_DEPLOYMENT_TARGET=14.0'],env=env)
 run(['cmake','--build',sdslbuild,'--target','sdsl','--parallel','3'],env=env)
 run(['cmake','--install',sdslbuild],env=env)
 static=[hts/'libhts.a']
 for name in ('iostreams','filesystem','program_options','date_time','atomic','system'):
  lib=pref['boost']/('lib/libboost_'+name+'.a')
  if lib.exists():static.append(lib)
 for lib in (pref['xz']/'lib/liblzma.a',pref['libdeflate']/'lib/libdeflate.a'):
  if not lib.is_file():raise RuntimeError('Missing static input '+str(lib))
  static.append(lib)
 static=list(sdslout.glob('lib/*.a'))+static
 binary=ROOT/'tracy-macos-arm64'
 args=['clang++','-std=c++17','-O2','-DNDEBUG','-mmacosx-version-min=14.0']
 for inc in (source/'src/jlib',hts,sdslout/'include',pref['boost']/'include',pref['xz']/'include',pref['libdeflate']/'include'):args+=['-isystem',inc]
 args+=[source/'src/tracy.cpp','-o',binary,*static,'-lz','-lbz2','-lpthread']
 run(args,cwd=source,env=env)
 linked=run(['otool','-L',binary])
 for line in linked.splitlines()[1:]:
  lib=line.strip().split(' (',1)[0]
  if not lib.startswith(('/usr/lib/','/System/Library/')):raise RuntimeError('Unexpected non-system dynamic dependency: '+lib)
 # Local ad-hoc signature preserves Mach-O integrity; it is not Developer-ID notarization.
 run(['codesign','--force','--sign','-',binary]);run(['codesign','--verify',binary])
 run([binary,'--version'])
 licenses=ROOT/'licenses';licenses.mkdir()
 for label,p in [('Tracy',source/'LICENSE'),('HTSlib',hts/'LICENSE'),('SDSL',sdsl/'COPYING')]:
  if p.is_file():shutil.copyfile(p,licenses/(label+'-LICENSE.txt'))
 for name in ('boost','xz','libdeflate'):
  paths=[p for p in pref[name].resolve().rglob('*') if p.is_file() and p.name.lower() in ('license','license.txt','copying','copyright','license_1_0.txt','copying.txt')]
  for i,p in enumerate(paths[:15]):shutil.copyfile(p,licenses/(name+'-'+str(i)+'-'+p.name))
 receipt={'status':'SOURCE_BUILT_STATIC_THIRD_PARTY_LIBRARIES','upstream_commit':COMMIT,'submodules':submodules,'binary_sha256':sha(binary),'binary_bytes':binary.stat().st_size,'build_host':platform.platform(),'compiler':run(['clang++','--version']),'brew_formula_metadata':formulas,'static_inputs':[{'name':p.name,'sha256':sha(p),'bytes':p.stat().st_size} for p in static],'dynamic_dependencies':linked,'minimum_macos':'14.0','official_release_binary_used':False,'notarized':False,'commands':LOG}
 (ROOT/'TRACY_SOURCE_BUILD.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
 return binary,licenses,receipt

def main():
 binary,licenses,receipt=build()
 spec=importlib.util.spec_from_file_location('rps_acquire_runtime',Path(__file__).with_name('rps_acquire_runtime.py'))
 m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 original=m.download
 def source_override(url,dest,expected,alg='sha256'):
  if url.endswith('/tracy-v0.9.1-macos-arm64'):
   dest.parent.mkdir(parents=True,exist_ok=True)
   if dest.exists():raise RuntimeError('Refuse to replace Tracy output')
   shutil.copyfile(binary,dest);dest.chmod(0o755)
   shutil.copytree(licenses,dest.parent/'third_party_licenses')
   shutil.copyfile(ROOT/'TRACY_SOURCE_BUILD.json',dest.parent/'TRACY_SOURCE_BUILD.json')
   m.SOURCES.append({'url':'https://github.com/gear-genomics/tracy/tree/'+COMMIT,'filename':dest.name,'bytes':dest.stat().st_size,'sha256':sha(dest),'source_built':True,'published_binary_digest_not_applicable':True,'source_build_receipt':'tools/tracy/TRACY_SOURCE_BUILD.json'})
   return dest
  return original(url,dest,expected,alg)
 m.download=source_override
 m.main()

if __name__=='__main__':
 try:main()
 except Exception as e:
  OUT.mkdir(exist_ok=True);(OUT/'MACOS_BUILD_FAILURE.json').write_text(json.dumps({'error':str(e),'traceback':traceback.format_exc(),'commands':LOG},indent=2),encoding='utf-8');raise

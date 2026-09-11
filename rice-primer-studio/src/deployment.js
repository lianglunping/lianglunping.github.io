/* Distribution metadata only; no primer-design parameters are changed here. */
(function () {
  'use strict';
  const meta = name => document.querySelector(`meta[name="${name}"]`)?.content || '';
  const publicStatic = meta('rps-deployment') === 'PUBLIC_STATIC';
  const metadata = Object.freeze({
    mode: publicStatic ? 'PUBLIC_STATIC' : 'STANDALONE_LOCAL_OPTIONAL',
    appVersion: RPS.VERSION,
    sourceSha256: meta('rps-source-sha256'),
    nativeBackendExposed: false,
    scientificValidation: 'CANDIDATE_REVIEW_ONLY'
  });
  globalThis.RPSDeployment = Object.freeze({
    publicStatic,
    snapshot: () => ({...metadata})
  });
  const description = document.getElementById('deploymentDescription');
  const source = document.getElementById('deploymentSource');
  description.textContent = publicStatic
    ? '公开静态站点模式：选择的序列文件在浏览器处理；本站不提供上传 API，不调用本机 Primer3。候选仍需人工复核。'
    : '本地单文件模式：默认浏览器内计算；Primer3 仅可在配套本机服务中显式启用。真实全基因组与实验效果尚未验收。';
  source.textContent = `v${RPS.VERSION} · source ${metadata.sourceSha256.slice(0, 12)} · 不内置 MSU 全基因组`;
  if (publicStatic) {
    const select = document.getElementById('designEngine');
    select.querySelector('option[value="native"]')?.remove();
    select.value = 'offline';
    document.getElementById('checkBackend').hidden = true;
    document.getElementById('engineStatus').textContent = '本站仅启用离线候选枚举；Primer3 不在 GitHub Pages 上执行，亦不从本站连接本机服务。';
    const help = document.getElementById('deploymentHelp');
    help.href = 'help.html';
    help.textContent = '使用说明与隐私边界';
  }
})();

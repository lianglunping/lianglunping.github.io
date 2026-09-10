# Rice Primer Studio — incomplete publication staging

**NOT A DEPLOYABLE RELEASE. Do not merge this branch into `hexo`.**

This branch stages a subset of the Rice Primer Studio source modules. The upload of the audit module was blocked by the tool platform; publication was stopped rather than bypassing the block or publishing a broken application.

Available source files: `engine.js`, `io.js`, `gff.js`, `review.js`, `deployment.js`, and `audit_ui.js`. These six files were transferred with byte-exact Git tree verification. They are not a complete runnable application. The audit module, application modules, HTML entry point, build files and full tests have not been pushed here.

The existing blog and its `hexo` publishing branch are unchanged. No Rice Primer Studio website has been published by this staging commit. The added CI configuration is a draft and has not passed against a complete remote application; missing-file failures must not be reported as successful validation.

The complete local `1.1.0-rc3` review package is distinct from this partial branch. Tests of that local package do not certify the contents of this branch. Full MSU FASTA end-to-end validation, accepted native Primer3 numerical validation and wet-lab validation remain outstanding.

No sample records, background-variant lists, full reference genome or private research results are included.

## 发布状态

当前只有部分源码进入发布分支，尚未完成代码推送、合并或站点上线。不得将此目录作为完整应用部署。后续必须补齐源码、核对文件摘要、通过完整测试，并另行验证 Pages 的实际部署和页面内容。

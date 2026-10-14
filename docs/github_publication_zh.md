# GitHub 共同作者发布指南

本包对应论文中的十组参数，包含主文、405 补充文件、证书、图目录及精确核验程序。
现有 288/300 研究不属于这篇论文的十项证书，不在本次发布包中追加为已完成成果。

## 1. 选择仓库与协作方式

按最新决定，由 `SelenaJXWANG` 建立个人仓库 `srg-gegenbauer-counting`，再邀请其余作者。
曾考虑的组织名为 `srg-counting-research`，本轮不新建组织。仓库归属与论文作者顺序是不同事项；
引用作者以论文及根目录 `CITATION.cff` 为准。

尚需共同检查时可先建私有仓库；准备公开后再设为 Public。若已决定立即公开，也可
直接创建公开仓库。创建空仓库时不要另选自动生成的 README、许可证或 `.gitignore`，
因为本包已有相应文件。代码采用 MIT，论文与证书采用 CC BY 4.0；适用范围见 `LICENSES.md`。

## 2. 上传完整目录

推荐使用 GitHub Desktop：登录、添加解压后的本地仓库目录，再选择 Publish
repository。上传的是 `srg-gegenbauer-counting` 里面的文件，根目录应直接显示
`README.md`、`paper/`、`code/`、`results/`、`docs/`、`extras/` 和 `CITATION.cff`。

也可在该目录用 Git 上传到已创建的空仓库。下方地址为已经建立的公开仓库。
Git 提交身份应使用作者本人确认的姓名和邮箱；已有 Git 仓库无需再次初始化。

```sh
git init -b main
git add .
git commit -m "Release manuscript and exact SRG certificate snapshot"
git remote add origin https://github.com/SelenaJXWANG/srg-gegenbauer-counting.git
git push -u origin main
```

不要只上传一个 ZIP 到代码区，否则读者无法直接查看论文链接、证书目录和核验脚本。
ZIP 适合作为 Release 附件。环境目录、Python 缓存和生成的实验输出已由 `.gitignore`
排除。文件均小于普通 GitHub 单文件限制，不需要为本包另设 Git LFS。

## 3. 添加共同作者

由仓库所有者在 Settings 中的 Collaborators / Manage access 邀请共同作者的 GitHub
账号。邀请会发送给对方；请由你本人明确选择账号并操作。本文作者名单已保存在
`CITATION.cff`，不需要通过代填贡献比例或虚构 Git 提交历史来显示共同署名。

协作修改使用分支与 Pull Request；可对 `main` 设置要求一名共同作者审阅及精确
核验通过后合并。仓库公开后，GitHub Free 支持公开仓库的受保护分支。

## 4. 验证论文和证书

```sh
python3 -S code/verify.py
```

仅需 Python 3.9 以上标准库，预期最后显示：

```text
PASS: five exact verification stages; fixed inputs unchanged.
```

`.github/workflows/verify.yml` 会在 GitHub push/PR 时执行同一入口。
该工作流的云端运行尚须上传后确认，本机核验不能代替云端执行结果。

## 5. 创建正式版本

共同检查作者、论文版本、许可证及 Actions 结果后，为指定 commit 创建 Release，
建议首次标签为 `v1.0.0`。标题可用：

`Manuscript and exact certificates - v1.0.0`

正文可使用 `docs/release_notes_v1.0.0.md`。附件加入两份 PDF、本包 ZIP 及 SHA-256
清单。核验报告可作为额外附件。修订时创建新标签，保留旧版本供引用者重放。

## 6. 取得 DOI 与引用

如需 DOI，可连接 Zenodo 与 GitHub，在 Zenodo 启用该公开仓库的归档，再发布新的
GitHub Release。也可将固定发布包及元数据手动上传 Zenodo。
论文的预印本 DOI 与代码/证书版本 DOI 可分开维护，应明确指向哪个对象。

取得实际仓库地址、DOI 或 arXiv 编号后，再更新 README、`CITATION.cff` 和文稿中的
公开来源说明。目前没有生成或假定任何 DOI。

## 官方说明

- 建立仓库：https://docs.github.com/en/repositories/creating-and-managing-repositories/quickstart-for-repositories
- 引用文件：https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files
- 邀请合作者：https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/repository-access-and-collaboration/inviting-collaborators-to-a-personal-repository
- 保护分支：https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-protected-branches/about-protected-branches
- Zenodo 归档：https://help.zenodo.org/docs/github/

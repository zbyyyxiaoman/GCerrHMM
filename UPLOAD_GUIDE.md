# GitHub 上传手把手教程

这份教程只针对 `GCerrHMM` 的**干净发布目录**，不要直接把整个开发目录
`project` 推上去。开发目录里可能包含内部脚本、日志、数据路径和本地
Git 历史；发布目录已经过清理和私有路径扫描。

下面的命令以 PowerShell 为例，路径请替换成你自己的解压位置。

## 第 0 步：准备一个干净发布目录

在项目根目录重新生成一次发布包：

```powershell
Set-Location "C:\path\to\project"

python scripts/package_github_release.py `
  --source . `
  --output dist/gcerrhmm_github_20260929 `
  --zip dist/GCerrHMM_github_20260929.zip
```

然后确认发布目录存在：

```powershell
Set-Location "dist\gcerrhmm_github_20260929"
Get-ChildItem
```

你应该看到：

```text
README.md
REPRODUCIBILITY.md
UPLOAD_GUIDE.md
reproduce.sh
src/
scripts/
config/
docs/
tests/
LICENSE
CITATION.cff
environment.yml
pyproject.toml
```

不要在 `project` 根目录执行后面的 `git init`。只在干净发布目录里
创建新仓库。

## 第 1 步：先做发布前检查

在发布目录中运行：

```powershell
python tests/run_tests.py
bash reproduce.sh --static-check
bash reproduce.sh --smoke
$env:PATH_AUDIT_STRICT = "1"
bash reproduce.sh --audit-paths
Remove-Item Env:PATH_AUDIT_STRICT
```

需要看到的信号：

```text
TESTS_OK
STATIC_AUDIT_OK
REPRODUCE_SMOKE_OK
PATH_AUDIT=... findings=0
```

如果这里不是 `findings=0`，不要上传，先把完整错误发给 Codex 处理。

再确认没有超过 GitHub 单文件限制的大文件：

```powershell
Get-ChildItem -Recurse -File |
  Where-Object { $_.Length -gt 100MB } |
  Select-Object FullName, Length
```

正常情况下这条命令不应输出任何文件。

## 第 2 步：创建空的 GitHub 仓库

在浏览器中打开 GitHub：

1. 点右上角 `+`；
2. 选择 `New repository`；
3. Repository name 填 `GCerrHMM`；
4. Description 可填：
   `GC-aware error HMM and reproducible long-read simulation benchmark`
5. 选择 `Public`；
6. **不要勾选** `Add a README file`；
7. **不要勾选** `Add .gitignore`；
8. **不要**选择 licence；
9. 点击 `Create repository`。

这一步之所以不勾选初始化文件，是为了避免首次推送时出现远程和本地
互不相关的历史。

创建后，GitHub 会显示仓库地址，例如：

```text
https://github.com/<你的GitHub用户名>/GCerrHMM.git
```

把它复制下来，下一步要用。

## 第 3 步：配置 Git 身份

在 PowerShell 中执行：

```powershell
git config --global user.name "<你的GitHub显示名>"
git config --global user.email "<你的GitHub邮箱>"
```

如果这台机器已经配置过，可以跳过。

## 第 4 步：在发布目录中初始化本地仓库

确认当前目录是干净发布目录：

```powershell
Set-Location "C:\path\to\gcerrhmm_github_20260929"

git init -b main
git add .
git status
```

检查 `git status`：

- 应看到 `README.md`、`src/`、`scripts/`、`config/`、`docs/`、`tests/`；
- 不应看到 `data/`、`results/`、`logs/`、`dist/`、`.venv/`；
- 不应看到 `.bam`、`.fastq`、`.fastq.gz`。

如果 `git status` 出现数据文件，停止操作，先检查是不是在错误的目录。

## 第 5 步：首次提交

```powershell
git commit -m "Initial release of GCerrHMM v1.0.0"
```

确认提交：

```powershell
git log --oneline -1
```

## 第 6 步：连接远程仓库并推送

把下面的 `<你的GitHub用户名>` 替换成真实用户名：

```powershell
git remote add origin https://github.com/<你的GitHub用户名>/GCerrHMM.git
git remote -v
git push -u origin main
```

第一次推送会要求登录 GitHub：

- 推荐使用 Git Credential Manager / 浏览器登录；
- 如果使用 HTTPS token，只把 token 填进弹出窗口；
- **不要把 token 写进 README、脚本或远程 URL**。

推送成功后，刷新 GitHub 页面，应能看到 `README.md` 和完整目录。

## 第 7 步：打正式版本标签

```powershell
git tag -a v1.0.0 -m "GCerrHMM v1.0.0 reviewer release"
git push origin v1.0.0
```

以后重要的审稿版本可以使用：

```text
v1.0.1
v1.0.2
```

不要删除旧标签；新增标签保留审稿轨迹。

## 第 8 步：检查 GitHub Actions

进入 GitHub 仓库页面：

1. 点 `Actions`；
2. 找到 `reproducibility-smoke`；
3. 确认最近一次运行是绿色；
4. 如果失败，点进失败步骤，把日志发给 Codex。

当前 workflow 会执行：

```text
bash reproduce.sh --static-check
bash reproduce.sh --tests
bash reproduce.sh --smoke
```

## 第 9 步：填写作者和 DOI 信息

首次发布后，编辑：

```text
CITATION.cff
docs/manuscript/author_metadata.template.json
README.md
```

需要确认：

- GitHub 仓库 URL；
- 作者、单位、通讯作者；
- 基金；
- Zenodo DOI；
- 论文 Availability 段落中的链接。

如果还没拿到 Zenodo DOI，可以先写：

```text
DOI: to be assigned
```

拿到后必须替换。

## 第 10 步：以后如何更新仓库

每次修改后，重新生成干净发布目录：

```powershell
Set-Location "C:\path\to\project"

python scripts/package_github_release.py `
  --source . `
  --output dist/gcerrhmm_github_20260929 `
  --zip dist/GCerrHMM_github_20260929.zip
```

进入发布目录，确认差异后提交：

```powershell
Set-Location "dist\gcerrhmm_github_20260929"
git status
git diff
git add .
git commit -m "Update reviewer release"
git push
```

不要直接把开发目录里的 `data/`、`results/` 或 `logs/` 加进 Git。

## 最短命令清单

如果你已经完成环境准备和 GitHub 空仓库创建，最短流程是：

```powershell
Set-Location "C:\path\to\gcerrhmm_github_20260929"

git init -b main
git add .
git commit -m "Initial release of GCerrHMM v1.0.0"
git remote add origin https://github.com/<你的GitHub用户名>/GCerrHMM.git
git push -u origin main
git tag -a v1.0.0 -m "GCerrHMM v1.0.0 reviewer release"
git push origin v1.0.0
```

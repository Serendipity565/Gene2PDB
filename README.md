# Gene2PDB

基于基因名或 PDB ID，快速检索相关蛋白质结构、进行基础物化性质分析，并生成 Markdown 格式的结构分析报告，同时提供 Web 端 3D 可视化界面的小工具。

后端使用 Python + Flask 提供 REST API，前端是一个静态页面（原生 JS + `$3Dmol` + `marked`）用于交互式浏览和展示结果。

---

## 功能简介

- **基因名 → PDB 结构映射**  
  通过 gget + UniProt + PDBe API，将基因名（如 `INS`）映射到对应的 PDB 结构列表。
- **PDB 结构基础信息查询**  
  查询单个 PDB 条目的标题、分辨率、实验方法、来源物种、发布日期等信息。
- **结构物化性质与二级结构分析**  
  使用 BioPython + DSSP，统计链数量、残基数、原子数、二级结构（α-螺旋、β-折叠、线圈）等指标。
- **自动生成 Markdown 分析报告**  
  根据基因名或给定的 PDB ID 列表，生成包含基础信息、物化性质、在线浏览链接等内容的 Markdown 报告，可在网页端直接渲染或下载。
- **一键快速分析**  
  输入基因名或 PDB ID，自动识别类型并给出相应的结构与分析结果。
- **网页端 3D 可视化**  
  使用 `$3Dmol` 在浏览器中渲染 3D 结构，支持 cartoon / stick / sphere / surface 等多种风格与颜色方案。

---

## 仓库结构

- `app.py`：Flask 后端 API 入口。
- `gget_pdb.py`：核心逻辑，包括基因→结构映射、PDB 信息获取、物化性质分析、报告生成、3D 查看等。
- `run_analysis.py`：命令行/Notebook 示例脚本，可用于快速测试后端逻辑。
- `frontend/`
  - `index.html`：前端页面入口。
  - `app.js`：前端交互逻辑，调用后端 API、加载 3D 结构、渲染报告等。
  - `styles.css`：样式文件。
- 若干 `pdbXXXX.ent`：示例/缓存的 PDB 文件，可用于测试本地分析逻辑。

---

## 环境与依赖

### 运行环境

- 操作系统：macOS / Linux / Windows 均可（以下命令以 macOS + zsh 为例）。
- Python：推荐 Python **3.10+**（项目在 3.12 环境下测试通过）。
- 浏览器：Chrome / Edge / Firefox 等现代浏览器。

### Python 依赖

项目根目录下提供了 `requirements.txt`，包含主要依赖，例如：

- `Flask`、`flask-cors` – 提供 Web API 与跨域支持
- `requests` – 访问 RCSB/PDBe/UniProt 等在线接口
- `biopython` – 解析 PDB、进行结构分析
- `pandas`, `numpy` – 数据处理
- `py3Dmol` – Notebook 中 3D 可视化
- `gget` – 基因信息检索

安装方式见下文“安装与启动”。

### 前端依赖

前端为**纯静态页面**，所需库（`$3Dmol`、`marked` 等）通过 CDN 在 `index.html` 中引入，无需额外构建步骤：

- 直接在浏览器中打开 `frontend/index.html` 即可使用（推荐配合简单的 HTTP 静态服务器）。

---

## 安装与启动

以下步骤以项目根目录 `Gene2PDB` 为起点。

### 1. 获取代码

```bash
# 克隆仓库
git clone <你的仓库地址>
cd Gene2PDB
```

也可以直接下载 ZIP 解压后进入目录。

### 2. 创建并激活虚拟环境（推荐）

```bash
# 创建虚拟环境
python3 -m venv .venv

# 激活虚拟环境（macOS / Linux, zsh/bash）
source .venv/bin/activate

# 如果是 Windows（PowerShell）
# .venv\Scripts\Activate.ps1
```

### 3. 安装 Python 依赖

```bash
pip install -r requirements.txt
```

可选：使用下面的命令简单测试核心分析逻辑是否正常工作：

```bash
python run_analysis.py
```

如果能看到类似“找到 INS 的结构”和一段 Markdown 报告输出，说明依赖安装基本正常。

### 4. 启动后端（Flask API）

在项目根目录运行：

```bash
python app.py
```

默认配置：

- 地址：`127.0.0.1`
- 端口：`8080`
- 调试模式：`debug=True`（开发环境适用，正式部署时建议关闭）

启动后终端会打印一段简要 API 文档，例如：

- `GET /api/health`
- `GET /api/gene/structures?gene_name=INS`
- `GET /api/pdb/info/<pdb_id>` 等。

你可以用 curl 或浏览器检查健康状态：

```bash
curl http://localhost:8080/api/health
```

预期返回：

```json
{"status": "ok", "message": "PDB分析服务正常运行"}
```

### 5. 启动前端（可选两种方式）

#### 方式 A：浏览器直接打开静态页面（最简单）

1. 保证后端已在 `8080` 端口运行（`python app.py`）。
2. 使用文件管理器或浏览器，打开：
   - `frontend/index.html`

> 注意：`frontend/app.js` 中将 API 地址写死为 `http://localhost:8080/api`，请确保后端端口与之保持一致，否则需要修改 `API_BASE` 常量。

#### 方式 B：用简单 HTTP 服务器托管前端（推荐）

```bash
cd frontend
python3 -m http.server 8000
```

然后在浏览器访问：

- `http://localhost:8000`

此时前端会通过 `http://localhost:8080/api` 调用后端接口。

页面加载后会自动调用 `/api/health` 检查后端状态，如果未启动，会在页面中提示“无法连接到后端服务，请确保已启动 Flask 服务器 (python app.py)”。

---

## 后端 API 说明（简要）

### 1. 健康检查

- **GET** `/api/health`
- 返回：

```json
{"status": "ok", "message": "PDB分析服务正常运行"}
```

### 2. 基因名 → 结构列表

- **GET** `/api/gene/structures`
- 查询参数：
  - `gene_name`（必填）：基因名，例如 `INS`
  - `species`（可选，默认 `human`）：物种
  - `max_structures`（可选，默认 5）：最多返回的结构数
- 示例：

```bash
curl "http://localhost:8080/api/gene/structures?gene_name=INS&species=human&max_structures=5"
```

- 返回示例（简化）：

```json
{
  "gene_name": "INS",
  "species": "human",
  "structures": ["7s5v", "7s60", ...],
  "count": 2
}
```

### 3. 查询单个 PDB 信息

- **GET** `/api/pdb/info/<pdb_id>`
- 示例：

```bash
curl "http://localhost:8080/api/pdb/info/7s5v"
```

- 返回示例（字段可能略有不同）：

```json
{
  "pdb_id": "7s5v",
  "title": "Some protein structure",
  "resolution": 2.1,
  "method": "X-RAY DIFFRACTION",
  "organism": "Homo sapiens",
  "release_date": "2021-01-01",
  "chains": [],
  "sequence": "...",
  "length": 300
}
```

### 4. 分析 PDB 结构

- **GET** `/api/pdb/analyze/<pdb_id>`
- 示例：

```bash
curl "http://localhost:8080/api/pdb/analyze/7s5v"
```

- 返回示例：

```json
{
  "pdb_id": "7s5v",
  "num_chains": 2,
  "num_residues": 250,
  "num_atoms": 2000,
  "secondary_structure": {
    "helix": 100,
    "beta_sheet": 50,
    "coil": 100
  }
}
```

### 5. 生成分析报告

- **GET** `/api/report`
- 调用方式一：基于基因名

```bash
curl "http://localhost:8080/api/report?gene_name=INS"
```

- 调用方式二：基于指定 PDB ID 列表

```bash
curl "http://localhost:8080/api/report?pdb_ids=7s5v&pdb_ids=7s60"
```

- 返回示例：

```json
{
  "report": "# 🧬 蛋白结构综合分析报告\n...（Markdown 文本）"
}
```

### 6. 一键快速分析

- **GET** `/api/quick`
- 查询参数：
  - `input`：基因名或 PDB ID
- 示例：

```bash
curl "http://localhost:8080/api/quick?input=INS"
```

- 返回示例（根据输入类型可能不同）：

```json
{
  "type": "gene",              // 或 "pdb_id"
  "gene_name": "INS",          // 如果是基因
  "pdb_ids": ["7s5v", ...],
  "info": { ... },              // 第一个结构的基础信息
  "analysis": { ... }           // 物化性质分析
}
```

---

## 前端使用说明

前端逻辑主要在 `frontend/app.js` 中，实现了以下功能：

- 搜索模式切换：
  - **按基因名搜索**：输入基因名，选择物种，列出相关结构，并生成基因级别的分析报告。
  - **按 PDB ID 搜索**：输入 4 位 PDB ID，直接查看该结构的详细信息与分析。
- 结构列表与详情：
  - 左侧显示结构列表（PDB ID、标题、分辨率、实验方法）。
  - 点击某个结构，右侧展示详细信息、二级结构统计、外部链接等。
- 3D 结构查看：
  - 使用 `$3Dmol` 加载 RCSB 上的 PDB 数据。
  - 支持切换显示模式（cartoon / stick / sphere / surface）。
  - 支持切换颜色方案（按谱带、按链、按二级结构）。
  - 提供重置视角按钮。
- 分析报告：
  - 从后端 `/api/report` 获取 Markdown 报告。
  - 使用 `marked` 在页面中渲染为 HTML。
  - 支持一键下载报告为 `.md` 文件。

使用流程示例：

1. 启动后端：`python app.py`。
2. 启动前端：
   - 直接打开 `frontend/index.html`，或
   - `cd frontend && python3 -m http.server 8000`，访问 `http://localhost:8000`。
3. 在页面顶部输入框中输入：
   - 示例基因：`INS`、`TP53` 等。
   - 示例 PDB ID：`7s5v`、`7s60` 等。
4. 浏览结构列表、切换不同结构，查看 3D 结构与详细信息。
5. 在“分析报告”区域查看 Markdown 报告，必要时点击“下载报告”保存到本地。

---

如果在使用或二次开发过程中遇到问题，欢迎在代码中添加注释、完善异常处理，或在 README 中继续补充更多示例与说明。

---
name: chem-structure-recognition
description: 手写实验笔记与非标准化学结构式自纠错工作流：本地高清切片 + RDKit MW/水合物/盐型反向纠偏 + ACS 1996 PNG/SVG 双预览 + 深度联动 Google science 插件（PubChem/NIST/OPSIN/ORD/反应原子守恒）。当用户上传手写化学笔记、反应式、投料表照片，或询问结构识别、化学计量核算、物性/基团效应/实验后处理与安全推演时使用。
---

# 手写化学实验笔记识别与权威自纠错管线 (v2.0)

核心哲学：**“高清多视角读图 + RDKit 确定性反向纠偏 + `science` 插件权威查证 + PNG/SVG 直接预览”**。
不再依赖不可靠的通用 OCR，也不在默认流程弹不兼容的 HTML 窗口；把全部精力放在赋予 Agent **自我发现错误、自我纠正、并用权威数据源托底**的能力上。

---

## 1. 四层自纠错与权威联动架构

1. **第一层：高清预处理与多视角局部读图（`chem_cli.py prep`）**
   - 自动校正 EXIF 与旋转方向，必要时做 CLAHE 对比度增强（抑制背面透印），并将大图切分为重叠的高分辨率局部切片（`top_half`、`bottom_half`、`quad_*` 或自定义 `--crop`）。
   - Agent 必须通过 `read_file` 查看局部高清切片后再确定苯环取代位、杂环原子及表格数字，严禁凭缩略全图盲猜。
2. **第二层：RDKit 确定性校验与 MW 反向纠偏（`chem_cli.py sanitize` / `audit-table`）**
   - `Chem.SanitizeMol()` 严格锁死价态与芳香杂环隐式氢 `[nH]`。
   - **MW 逆向纠错**：将识别出的 SMILES 核算分子量与笔记手写 `M` 自动比对。若偏差 $>0.08$，自动检索常见结晶水/溶剂合物/成盐形态（如 $\cdot\text{H}_2\text{O}$、$\cdot\text{HCl}$、$\cdot\text{TFA}$）及常见拓扑差值（$\pm\text{CH}_2$、$\pm\text{H}_2$、$\pm\text{O}$、氮规则 $\pm 1$ 等），逼迫 Agent 自我排查结构误判或抓出试剂水合差错。
   - 提取 SMARTS 官能团/酸碱中心特征（`functional_motifs`），作为后续性质推演的硬事实依据。
   - `audit-table` 一次性核对整张投料表：限制性试剂基准 `mmol`、各物料按笔记 MW 与按 RDKit MW 的理论值差异、底物浓度与 100% 理论产量。
3. **第三层：深度联动 Google `science` 插件（权威查证，杜绝记忆幻觉）**
   - **反应原子守恒**：调用 `verify_atom_conservation.py` 显式配平反应物、产物与小分子副产物（$\text{H}_2\text{O}$、$\text{Me}_2\text{NH}$、$\text{HCl}$ 等），防止漏算脱水量或氧化还原当量。
   - **IUPAC / 立体化学核验**：调用 `verify_chemical_nomenclature.py`（OPSIN + PubChem）校验系统命名与 $(E)/(Z)$、$(R)/(S)$ 构型。
   - **实测物性与安全**：调用 `pubchem_api.py`（CID、GHS 安全警示、理化性质）与 `nist_webbook_fetcher.py`（NIST 实测 $T_{\text{boil}}$ / $T_{\text{fus}}$ / 热化学参数）。
   - **反应先例与文献**：调用 `precedent_lookup.py`（ORD/USPTO/Google Patents）或 `pubmed-database` / `literature-search-openalex` 查同类反应的真实当量、温度与后处理先例。
   - **新分子/基团性质推演（无现成数据时）**：由 RDKit 拆解母核与官能团，去 `science` 工具查最接近母核的实测常数（如 $\text{p}K_a$、极性），再结合电子效应/氢键效应自主推演，**严禁拍脑袋编造**。
4. **第四层：Markdown + 内联预览图交付（PNG + SVG）**
   - 默认不生成定制 HTML 交互弹窗；所有分子与反应式均生成**白底 ACS 1996 PNG**（供 Markdown `![名称](file:///...)` 直接内联预览）与 **SVG** 矢量文件（供下载/排版），并在存在结构歧义时（如 `3-F` vs `4-F`）直接并用 `draw_molecules.py --align_core` 生成并排对比预览图供用户一眼确认。
   - 仅当用户明确要求亲自手动筛选/微调且多轮对话仍未调准时，才考虑使用 `rdkit-general/scripts/interactive_grid_viewer.py` 作为人工兜底。

---

## 2. 工具链速查

### A. 本插件前置纠错 CLI（Python + RDKit 本地执行）
脚本路径（相对本插件根目录）：`skills/chem-structure-recognition/scripts/chem_cli.py`
（若作为 Antigravity 插件安装，通常位于 `~/.gemini/config/plugins/chem-structure-recognition/skills/chem-structure-recognition/scripts/chem_cli.py` 或 Windows 下的 `%USERPROFILE%\.gemini\config\plugins\chem-structure-recognition\skills\chem-structure-recognition\scripts\chem_cli.py`）

1. **图片转正与高清切片**：
   ```bash
   python "<PLUGIN_DIR>/skills/chem-structure-recognition/scripts/chem_cli.py" prep "<IMG_PATH>" --out-dir "<SCRATCH_DIR>/prep" --rotate 0 --enhance
   ```
   （若图片横置/倒置，改 `--rotate 90|180|270`；需局部特写可加 `--crop x1,y1,x2,y2`）

2. **结构校验、MW 反向纠偏与 PNG+SVG 双导出**：
   ```bash
   python "<PLUGIN_DIR>/skills/chem-structure-recognition/scripts/chem_cli.py" sanitize "<SMILES>" --note-mw <NOTE_MW> --out-dir "<ARTIFACT_DIR>" --name "S1"
   ```

3. **反应式 PNG+SVG 渲染与元素差值速查**：
   ```bash
   python "<PLUGIN_DIR>/skills/chem-structure-recognition/scripts/chem_cli.py" render-reaction --reactants "<S1>" "<S2>" --products "<P>" --out-dir "<ARTIFACT_DIR>" --name "reaction"
   ```

4. **批量投料表审计**：
   ```bash
   python "<PLUGIN_DIR>/skills/chem-structure-recognition/scripts/chem_cli.py" audit-table --json-file "<SCRATCH_DIR>/table.json"
   ```

### B. 联动的 Google `science` 插件核心工具（可选增强）
（若已安装 `science` 插件，其脚本通常位于 `<SCIENCE_PLUGIN_DIR>/skills`，例如 `~/.gemini/config/plugins/science/skills`）

- **反应原子守恒核查（`python` 直接运行；若副产物是 SMILES 如二甲胺请加 `smiles:` 前缀并可加 `--strict-hydrogen`）**：
  ```bash
  python "<SCIENCE_PLUGIN_DIR>/skills/cheminformatics_reactions/scripts/verify_atom_conservation.py" --precursor "<S1>" "<S2>" --product "<P>" --byproducts "2*H2O" "smiles:CN(C)" --strict-hydrogen --json
  ```
- **并排对齐对比/高亮图（`python` 直接运行，适合歧义异构体对比）**：
  ```bash
  python "<SCIENCE_PLUGIN_DIR>/skills/rdkit_general/scripts/draw_molecules.py" --smiles "<SMI_A>" "<SMI_B>" --legends "A: 3-F" "B: 4-F" --align_core "<CORE_SMILES>" --output "<ARTIFACT_DIR>/compare.png"
  ```
- **描述符与结构滤网（`python` 直接运行）**：
  ```bash
  python "<SCIENCE_PLUGIN_DIR>/skills/cheminformatics_core/scripts/compute_descriptors.py" --smiles "<SMILES>" --json
  ```
- **PubChem 权威身份/性质/GHS 安全（`uv run`）**：
  ```bash
  uv run "<SCIENCE_PLUGIN_DIR>/skills/pubchem_database/scripts/pubchem_api.py" resolve --name "<NAME>" --output "<SCRATCH_DIR>/pubchem.json"
  uv run "<SCIENCE_PLUGIN_DIR>/skills/pubchem_database/scripts/pubchem_api.py" safety --cid <CID> --output "<SCRATCH_DIR>/safety.json"
  ```
- **NIST WebBook 实测热力学与沸点/熔点（`uv run`）**：
  ```bash
  uv run "<SCIENCE_PLUGIN_DIR>/skills/cheminformatics_core/scripts/nist_webbook_fetcher.py" --smiles "<SMILES>" --output-json "<SCRATCH_DIR>/nist.json"
  ```
- **IUPAC 命名与立体化学校验（`uv run`）**：
  ```bash
  uv run "<SCIENCE_PLUGIN_DIR>/skills/rdkit_general/scripts/verify_chemical_nomenclature.py" --smiles "<SMILES>" --claimed_name "<IUPAC_NAME>" --json
  ```
- **反应先例检索（ORD / USPTO / Google Patents，`uv run`）**：
  ```bash
  uv run "<SCIENCE_PLUGIN_DIR>/skills/cheminformatics_doc2mol/scripts/precedent_lookup.py" --reactants "<S1>" "<S2>" --product "<P>" --api all --output-json "<SCRATCH_DIR>/precedents.json"
  ```

---

## 3. 标准执行 SOP（4 步完成全局掌控）

1. **Step 1 · 切图看清**：收到实验笔记照片后，先用 `chem_cli.py prep` 转正并生成局部切片，用 `read_file` 逐块看清结构骨架、取代基位置与投料表数字，标出所有低置信度歧义点。
2. **Step 2 · 算账纠偏**：对每个分子运行 `chem_cli.py sanitize --note-mw ... --out-dir ...` 生成 PNG+SVG 并核对分子量/水合物；有投料表则跑 `audit-table`；有反应式则用 `verify_atom_conservation.py` 配平原子并生成 `render-reaction` 反应预览图。若发现 MW 或原子不守恒，**必须先自查结构是否正确，再下结论**。
3. **Step 3 · 权威查证与自主推演**：按问题需要调用 `pubchem_api.py`、`nist_webbook_fetcher.py`、`verify_chemical_nomenclature.py`、`precedent_lookup.py` 或文献检索工具获取实测/先例依据；涉及后处理、极性、酸碱萃取或安全注意事项时，基于 RDKit 识别出的基团与权威查证结果自主推演。
4. **Step 4 · 带图交付与证据标注**：
   - 在 Markdown 回复中用 `![标签](file:///.../name.png)` **直接内联渲染结构/反应式预览图**，并附带 `.svg` 链接；
   - 若有异构体或断笔歧义，一并给出 `draw_molecules.py` 对齐对比预览图；
   - 关键结论与性质参数须明确标出来源层级（`[RDKit 实算]`、`[PubChem/NIST/OPSIN 实测或校验]`、`[ORD/文献先例]`、`[基于基团的推演]`）。

# Chemical Structure Recognition & Self-Correction Suite (v2.0)

[![Version](https://img.shields.io/badge/version-2.0.0-blue.svg)](https://github.com/enenhengrusty/chem-structure-recognition)
[![RDKit](https://img.shields.io/badge/RDKit-2023.9%2B-green.svg)](https://www.rdkit.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

面向手写实验笔记与化学结构分析的**自纠错工作流插件（跨平台支持 Windows / macOS / Linux）**：聚焦**高清多视角读图 + RDKit 确定性反向纠偏 + ACS 1996 PNG/SVG 直接预览 + 可选联动 Google `science` 插件**。

---

## 🌟 v2.0 重铸设计理念

1. **丢弃无效 OCR 与硬编码 HTML 弹窗**
   - 实践证明整页通用 OCR 对手写小数点/分子式误识率高，真正保证结构与数字准确率的是 **VLM 局部高清切片读图 + RDKit 分子量反向碰撞**；
   - 默认用 **Markdown + 白底 ACS 1996 PNG 内联预览图 + SVG 矢量文件** 交付，结构直观可见且兼容性最好。
2. **赋予 Agent 自我纠错能力**
   - 通过 [chem_cli.py](skills/chem-structure-recognition/scripts/chem_cli.py) 对 `SMILES` 进行价态消毒、官能团 SMARTS 拆解、`note_mw` 偏差排查（自动匹配水合物/常见盐/多碳少氢拓扑差）以及全投料表批量核账，发现矛盾立刻自查修正。
3. **可选联动 Google `science` 插件，杜绝化学幻觉**
   - `chem_cli.py` 本身可**完全独立运行**；若宿主环境同时安装了 `science` 插件，可进一步联动：
     - 反应原子守恒 → `verify_atom_conservation.py`
     - IUPAC/手性/双键命名 → `verify_chemical_nomenclature.py` (OPSIN + PubChem)
     - 物性与 GHS 安全 → `pubchem_api.py` & `nist_webbook_fetcher.py`
     - 反应先例与条件 → `precedent_lookup.py` (ORD/USPTO/Google Patents) & 文献检索
     - 对齐/高亮对比预览 → `draw_molecules.py`

---

## 📦 安装与使用

### 1. 环境要求与依赖安装（Python 3.10+）

```bash
git clone https://github.com/enenhengrusty/chem-structure-recognition.git
cd chem-structure-recognition
pip install -r requirements.txt
# 或使用 uv：
# uv pip install -r requirements.txt
```

### 2. 作为 Antigravity 插件安装

Windows PowerShell：
```powershell
git clone https://github.com/enenhengrusty/chem-structure-recognition.git "$env:USERPROFILE\.gemini\config\plugins\chem-structure-recognition"
```

macOS / Linux：
```bash
git clone https://github.com/enenhengrusty/chem-structure-recognition.git "$HOME/.gemini/config/plugins/chem-structure-recognition"
```

（若你使用其他支持 `SKILL.md` 的 Agent，可将 `skills/chem-structure-recognition` 复制到对应 `skills/` 目录即可。）

---

## 🛠️ 核心 CLI 工具 (`skills/chem-structure-recognition/scripts/chem_cli.py`)

支持 5 个子命令，输出标准 JSON：

- **`prep`**：EXIF 自动转正、CLAHE 对比度增强（抑制背面透印）及重叠高清局部切片（`top_half`、`bottom_half`、`quad_*` 或 `--crop`）：
  ```bash
  python skills/chem-structure-recognition/scripts/chem_cli.py prep note.jpg --out-dir ./out/prep --enhance
  ```
- **`sanitize`**：RDKit 价态/芳香氢校验、SMARTS 官能团识别、`--note-mw` 反向纠偏（自动匹配结晶水/成盐/常见拓扑差）、白底 ACS 1996 **PNG + SVG** 双导出：
  ```bash
  python skills/chem-structure-recognition/scripts/chem_cli.py sanitize "CC(=O)Oc1ccccc1C(=O)O" --note-mw 180.16 --out-dir ./out --name aspirin
  ```
- **`render-reaction`**：反应式 PNG + SVG 渲染及元素差值核查：
  ```bash
  python skills/chem-structure-recognition/scripts/chem_cli.py render-reaction --reactants "Oc1ccccc1C(=O)O" "CC(=O)OC(C)=O" --products "CC(=O)Oc1ccccc1C(=O)O" --out-dir ./out --name rxn
  ```
- **`audit-table`**：JSON 批量投料表核账、限制性试剂基准计算、浓度与理论产量计算：
  ```bash
  python skills/chem-structure-recognition/scripts/chem_cli.py audit-table --json-file table.json
  ```
- **`stoichiometry`**：单步理论产量与 $>100\%$ 残余溶剂/无机盐警示：
  ```bash
  python skills/chem-structure-recognition/scripts/chem_cli.py stoichiometry --mw 180.16 --mmol 10 --actual-g 1.62
  ```

---

## 🧪 运行单元测试

```bash
python -m unittest discover -s tests -p "test_*.py"
```

---

## 📁 目录结构

```text
chem-structure-recognition/
├── plugin.json                                      # Antigravity 插件清单
├── README.md                                        # 项目说明与使用指南
├── LICENSE                                          # MIT License
├── requirements.txt                                 # Python 依赖清单
├── pyproject.toml                                   # 项目元数据
├── assets/
│   └── logo.svg                                     # 插件图标
├── rules/
│   └── AGENTS.md                                    # 零幻觉自纠错协议规则
├── skills/
│   └── chem-structure-recognition/
│       ├── SKILL.md                                 # 4 步自纠错 SOP 技能文档
│       └── scripts/
│           └── chem_cli.py                          # 核心 RDKit + OpenCV CLI
└── tests/
    └── test_chem_cli.py                             # 单元测试套件
```

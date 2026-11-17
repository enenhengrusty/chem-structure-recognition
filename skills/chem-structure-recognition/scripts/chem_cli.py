"""
chem_cli.py - Antigravity Chemistry Self-Correction & Preview CLI (v2.0)
========================================================================
Deterministic local Python 3.13 + RDKit + OpenCV/PIL helper providing:
  1. prep:            EXIF-aware upright rotation, contrast enhancement & high-res tile crops
  2. sanitize:        RDKit valence/Huckel validation, functional motif SMARTS,
                      MW reverse-correction (hydrates/salts/deltas), and ACS 1996 PNG+SVG previews
  3. render-reaction: Reaction scheme PNG+SVG rendering & quick elemental delta check
  4. audit-table:     Batch stoichiometry table & limiting-reagent mass-balance auditor
  5. stoichiometry:   Quick theoretical yield & >100% solvent/salt anomaly detector
"""

import argparse
from collections import Counter
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import cv2
import numpy as np
from PIL import Image, ImageOps
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

# Common hydrate / solvate / salt mass offsets (g/mol)
COMMON_ADDUCTS: List[Tuple[str, float]] = [
    ("·0.5H2O", 9.0076),
    ("·H2O", 18.0153),
    ("·1.5H2O", 27.0229),
    ("·2H2O", 36.0306),
    ("·3H2O", 54.0459),
    ("·HCl", 36.4609),
    ("·2HCl", 72.9219),
    ("·HBr", 80.9119),
    ("·TFA (C2HF3O2)", 114.0233),
    ("·AcOH", 60.0520),
    ("·0.5EtOAc", 44.0526),
    ("·MeOH", 32.0419),
    ("·EtOH", 46.0684),
    ("Na+ salt (-H +Na)", 21.9818),
    ("K+ salt (-H +K)", 38.0904),
]

# Common topological misrecognition mass deltas (g/mol)
STRUCTURAL_DELTAS: List[Tuple[str, float]] = [
    ("+CH2 (少画一个亚甲基)", 14.0266),
    ("-CH2 (多画一个亚甲基)", -14.0266),
    ("+H2 (漏画不饱和度/双键)", 2.0159),
    ("-H2 (多画不饱和度/双键)", -2.0159),
    ("+O (漏画一个氧原子)", 15.9994),
    ("-O (多画一个氧原子)", -15.9994),
    ("+N vs CH (氮规则: 多一个N少一个CH)", 1.0005),
    ("-N vs CH (氮规则: 少一个N多一个CH)", -1.0005),
    ("+F vs H (漏画一个氟)", 17.9905),
    ("-F vs H (多画一个氟)", -17.9905),
    ("+Cl vs H (漏画一个氯)", 34.4451),
    ("-Cl vs H (多画一个氯)", -34.4451),
]

# Deterministic SMARTS patterns for acid-base & reactivity motifs
FUNCTIONAL_MOTIF_SMARTS: Dict[str, str] = {
    "primary_aliphatic_amine": "[NX3;H2;!$(NC=O);!$(Nc)]",
    "secondary_aliphatic_amine": "[NX3;H1;!$(NC=O);!$(Nc)]([#6])[#6]",
    "tertiary_aliphatic_amine": "[NX3;H0;!$(NC=O);!$(Nc)]([#6])([#6])[#6]",
    "aniline_nh2": "[NX3;H2]c1ccccc1",
    "aromatic_aza_n_pyridine_like": "[nX2;r5,r6]",
    "aromatic_nh_pyrrole_like": "[nH;r5,r6]",
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "sulfonic_acid": "[SX4](=O)(=O)[OX2H1]",
    "phenol": "[OX2H1]c1ccccc1",
    "aliphatic_alcohol": "[OX2H1][CX4]",
    "aldehyde": "[CX3H1](=O)[#6]",
    "ketone": "[#6][CX3](=O)[#6]",
    "ester": "[#6][CX3](=O)[OX2][#6]",
    "amide": "[NX3][CX3](=O)[#6]",
    "enaminone": "[NX3][CX3]=[CX3][CX3](=O)",
    "aryl_halide": "[F,Cl,Br,I]c1ccccc1",
    "boronic_acid_or_ester": "[BX3]([OX2])[OX2]",
    "azide": "[NX1]~[NX2]~[NX2]",
    "peroxide_bond": "[OX2][OX2]",
    "nitro": "[NX3+](=O)[O-]",
    "nitrile": "[CX2]#[NX1]",
}


def _to_posix(path: Path) -> str:
    return path.resolve().as_posix()


def detect_functional_motifs(mol: Chem.Mol) -> List[str]:
    found: List[str] = []
    for name, smarts in FUNCTIONAL_MOTIF_SMARTS.items():
        pat = Chem.MolFromSmarts(smarts)
        if pat is not None and mol.HasSubstructMatch(pat):
            found.append(name)
    return found


def audit_mw_against_note(calc_mw: float, note_mw: float) -> dict:
    diff = round(note_mw - calc_mw, 4)
    abs_diff = abs(diff)
    result: Dict[str, Any] = {
        "calc_mw": round(calc_mw, 3),
        "note_mw": note_mw,
        "delta_note_minus_calc": diff,
    }
    if abs_diff <= 0.08:
        result["status"] = "MATCH"
        result["message"] = "手写分子量与 RDKit 计算值高度一致。"
        return result

    adduct_hits = []
    for label, offset in COMMON_ADDUCTS:
        if abs((calc_mw + offset) - note_mw) <= 0.15:
            adduct_hits.append({
                "adduct": label,
                "expected_mw": round(calc_mw + offset, 3),
                "residual": round((calc_mw + offset) - note_mw, 3),
            })

    struct_hits = []
    for label, offset in STRUCTURAL_DELTAS:
        if abs((calc_mw + offset) - note_mw) <= 0.25:
            struct_hits.append({
                "hypothesis": label,
                "expected_mw": round(calc_mw + offset, 3),
                "residual": round((calc_mw + offset) - note_mw, 3),
            })

    if adduct_hits:
        result["status"] = "ADDUCT_OR_SALT_MATCH"
        result["message"] = "手写 MW 与游离态不符，但高度匹配结晶水/成盐形态，请核对试剂瓶签。"
        result["adduct_matches"] = adduct_hits
    elif abs_diff <= 0.5:
        result["status"] = "MINOR_ROUNDING_OR_ATOM_WEIGHT_DIFF"
        result["message"] = f"存在轻微偏差 ({diff:+.3f})，多为四舍五入或原子量表精度差异。"
    else:
        result["status"] = "MISMATCH_REVIEW_REQUIRED"
        result["message"] = f"手写 MW 与结构偏差 {diff:+.3f}，请优先排查结构识别错误、结晶水或成盐状态。"

    if struct_hits:
        result["structural_hypotheses"] = struct_hits
    return result


def _render_mol_previews(
    mol: Chem.Mol, out_dir: Path, name: str, width: int = 420, height: int = 280
) -> Dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    AllChem.Compute2DCoords(mol)
    mean_bl = rdMolDraw2D.MeanBondLength(mol) if mol.GetNumBonds() > 0 else 1.5

    # White-background ACS 1996 PNG (universal inline preview in dark/light UI)
    png_path = out_dir / f"{name}.png"
    d_png = rdMolDraw2D.MolDraw2DCairo(width, height)
    opts_png = d_png.drawOptions()
    rdMolDraw2D.SetACS1996Mode(opts_png, mean_bl)
    opts_png.fixedBondLength = -1.0
    opts_png.clearBackground = True
    opts_png.padding = 0.08
    d_png.DrawMolecule(mol)
    d_png.FinishDrawing()
    png_path.write_bytes(d_png.GetDrawingText())

    # White-background ACS 1996 SVG
    svg_path = out_dir / f"{name}.svg"
    d_svg = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts_svg = d_svg.drawOptions()
    rdMolDraw2D.SetACS1996Mode(opts_svg, mean_bl)
    opts_svg.fixedBondLength = -1.0
    opts_svg.clearBackground = True
    opts_svg.padding = 0.08
    d_svg.DrawMolecule(mol)
    d_svg.FinishDrawing()
    svg_text = d_svg.GetDrawingText()
    idx = svg_text.find("<svg")
    if idx != -1:
        svg_text = svg_text[idx:]
    svg_path.write_text(svg_text, encoding="utf-8")

    return {
        "preview_png": _to_posix(png_path),
        "preview_svg": _to_posix(svg_path),
    }


def prep_image(
    img_path: str,
    out_dir: str,
    rotate_deg: int = 0,
    enhance: bool = False,
    custom_crop: Optional[str] = None,
) -> dict:
    src = Path(img_path)
    if not src.is_file():
        return {"error": f"Image file not found: {img_path}"}

    dst_dir = Path(out_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)

    img = ImageOps.exif_transpose(Image.open(src)).convert("RGB")
    if rotate_deg in (90, 180, 270):
        img = img.rotate(rotate_deg, expand=True)

    w, h = img.size
    upright_path = dst_dir / "upright.png"
    img.save(upright_path)

    outputs: Dict[str, Any] = {
        "upright": _to_posix(upright_path),
        "width": w,
        "height": h,
    }

    # Contrast stats & optional CLAHE enhancement for faint pencil strokes
    gray = np.array(img.convert("L"), dtype=np.float32)
    bg_median = float(np.median(gray))
    stroke_dark = float(np.percentile(gray, 3))
    outputs["bg_median_gray"] = round(bg_median, 1)
    outputs["stroke_contrast"] = round(bg_median - stroke_dark, 1)

    if enhance:
        bgr = cv2.cvtColor(np.array(img), cv2.COLOR_RGB2BGR)
        lab = cv2.cvtColor(bgr, cv2.COLOR_BGR2LAB)
        l_chan, a_chan, b_chan = cv2.split(lab)
        clahe = cv2.createCLAHE(clipLimit=2.2, tileGridSize=(8, 8))
        l_enh = clahe.apply(l_chan)
        enh_rgb = cv2.cvtColor(cv2.merge((l_enh, a_chan, b_chan)), cv2.COLOR_LAB2RGB)
        enh_path = dst_dir / "enhanced.png"
        Image.fromarray(enh_rgb).save(enh_path)
        outputs["enhanced"] = _to_posix(enh_path)

    if custom_crop:
        parts = [float(x.strip()) for x in custom_crop.split(",")]
        if len(parts) != 4:
            return {"error": "custom_crop must be x1,y1,x2,y2 (normalized 0..1)"}
        x1, y1, x2, y2 = parts
        box = (int(w * x1), int(h * y1), int(w * x2), int(h * y2))
        crop_path = dst_dir / "custom_crop.png"
        img.crop(box).save(crop_path)
        outputs["custom_crop"] = _to_posix(crop_path)
        return outputs

    # Generate overlapping high-resolution crops so VLM inspects structures & tables cleanly
    regions = {
        "top_half": (0.0, 0.0, 1.0, 0.55),
        "bottom_half": (0.0, 0.45, 1.0, 1.0),
        "quad_tl": (0.0, 0.0, 0.55, 0.55),
        "quad_tr": (0.45, 0.0, 1.0, 0.55),
        "quad_bl": (0.0, 0.45, 0.55, 1.0),
        "quad_br": (0.45, 0.45, 1.0, 1.0),
    }
    crops: Dict[str, str] = {}
    for key, (x1, y1, x2, y2) in regions.items():
        box = (int(w * x1), int(h * y1), int(w * x2), int(h * y2))
        p = dst_dir / f"{key}.png"
        img.crop(box).save(p)
        crops[key] = _to_posix(p)
    outputs["crops"] = crops
    return outputs


def sanitize_structure(
    smiles: str,
    note_mw: Optional[float] = None,
    out_dir: Optional[str] = None,
    name: str = "mol",
    width: int = 420,
    height: int = 280,
) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "valid": False,
            "error": f"Chem.MolFromSmiles failed to parse SMILES: {smiles}",
        }

    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        return {
            "valid": False,
            "error": f"Sanitization / valence error: {exc}",
        }

    mw_raw = Descriptors.MolWt(mol)
    chiral = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

    res: Dict[str, Any] = {
        "valid": True,
        "canon_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
        "inchikey": Chem.MolToInchiKey(mol),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(mw_raw, 3),
        "exact_mw": round(Descriptors.ExactMolWt(mol), 4),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "rot_bonds": Descriptors.NumRotatableBonds(mol),
        "rings": rdMolDescriptors.CalcNumRings(mol),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
        "chiral_centers": [{"atom_idx": i, "descriptor": d} for i, d in chiral],
        "functional_motifs": detect_functional_motifs(mol),
    }

    if note_mw is not None:
        res["mw_audit"] = audit_mw_against_note(mw_raw, note_mw)

    if out_dir:
        res.update(_render_mol_previews(mol, Path(out_dir), name, width=width, height=height))

    return res


def _sum_elements(smiles_list: List[str]) -> Counter:
    total: Counter = Counter()
    for smi in smiles_list:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise ValueError(f"Invalid SMILES in reaction: {smi}")
        m_h = Chem.AddHs(m)
        for atom in m_h.GetAtoms():
            total[atom.GetSymbol()] += 1
    return total


def render_reaction(
    reactants: List[str],
    products: List[str],
    out_dir: str,
    agents: Optional[List[str]] = None,
    name: str = "reaction",
    width: int = 1000,
    height: int = 300,
) -> dict:
    agents = agents or []
    rxn_smarts = f"{'.'.join(reactants)}>{'.'.join(agents)}>{'.'.join(products)}"
    rxn = AllChem.ReactionFromSmarts(rxn_smarts, useSmiles=True)

    r_elems = _sum_elements(reactants)
    p_elems = _sum_elements(products)
    delta = {}
    for el in sorted(set(r_elems) | set(p_elems)):
        diff = p_elems[el] - r_elems[el]
        if diff != 0:
            delta[el] = diff

    dst_dir = Path(out_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)

    png_path = dst_dir / f"{name}.png"
    d_png = rdMolDraw2D.MolDraw2DCairo(width, height)
    d_png.drawOptions().clearBackground = True
    d_png.drawOptions().padding = 0.08
    d_png.DrawReaction(rxn)
    d_png.FinishDrawing()
    png_path.write_bytes(d_png.GetDrawingText())

    svg_path = dst_dir / f"{name}.svg"
    d_svg = rdMolDraw2D.MolDraw2DSVG(width, height)
    d_svg.drawOptions().clearBackground = True
    d_svg.drawOptions().padding = 0.08
    d_svg.DrawReaction(rxn)
    d_svg.FinishDrawing()
    svg_text = d_svg.GetDrawingText()
    idx = svg_text.find("<svg")
    if idx != -1:
        svg_text = svg_text[idx:]
    svg_path.write_text(svg_text, encoding="utf-8")

    return {
        "valid": True,
        "reaction_smiles": rxn_smarts,
        "preview_png": _to_posix(png_path),
        "preview_svg": _to_posix(svg_path),
        "element_delta_product_minus_reactant": delta,
        "note": "If delta is non-zero (e.g. -H2O, -Me2NH), verify full balance via science verify_atom_conservation.py.",
    }


def audit_table(payload: dict) -> dict:
    """
    Expected JSON payload:
    {
      "ref_name": "S2",
      "solvent_vol_ml": 15.0,
      "actual_product_g": null,
      "rows": [
        {"role": "reactant", "name": "S1", "smiles": "...", "note_mw": 123.16, "eq": 1.0, "mass_mg": 560.89},
        {"role": "reactant", "name": "S2", "smiles": "...", "note_mw": 193.23, "eq": 1.0, "mass_mg": 880.0},
        {"role": "product",  "name": "P",  "smiles": "...", "eq": 1.0}
      ]
    }
    """
    rows = payload.get("rows", [])
    ref_name = payload.get("ref_name")
    if not rows:
        return {"error": "Payload must contain non-empty 'rows' list."}

    enriched: List[Dict[str, Any]] = []
    ref_row = None

    for row in rows:
        item = dict(row)
        smi = item.get("smiles")
        if smi:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                return {"error": f"Invalid SMILES in row {item.get('name')}: {smi}"}
            Chem.SanitizeMol(m)
            rd_mw = Descriptors.MolWt(m)
            item["rdkit_mw"] = round(rd_mw, 3)
            item["formula"] = rdMolDescriptors.CalcMolFormula(m)
            if item.get("note_mw") is not None:
                item["mw_audit"] = audit_mw_against_note(rd_mw, float(item["note_mw"]))
        else:
            item["rdkit_mw"] = None

        # Effective mass in mg (supports mass_mg, mass_g, or vol_ml * density)
        mass_mg = item.get("mass_mg")
        if mass_mg is None and item.get("mass_g") is not None:
            mass_mg = float(item["mass_g"]) * 1000.0
        elif mass_mg is None and item.get("vol_ml") is not None and item.get("density") is not None:
            mass_mg = float(item["vol_ml"]) * float(item["density"]) * 1000.0
        purity = float(item.get("purity", 1.0))
        if mass_mg is not None:
            mass_mg = float(mass_mg) * purity
        item["effective_mass_mg"] = round(mass_mg, 4) if mass_mg is not None else None

        if item.get("name") == ref_name:
            ref_row = item
        enriched.append(item)

    if ref_row is None:
        return {"error": f"Reference row '{ref_name}' not found in rows."}
    if ref_row["effective_mass_mg"] is None:
        return {"error": f"Reference row '{ref_name}' has no mass/volume+density."}

    ref_eq = float(ref_row.get("eq", 1.0))
    ref_note_mw = float(ref_row["note_mw"]) if ref_row.get("note_mw") else ref_row["rdkit_mw"]
    ref_rd_mw = ref_row["rdkit_mw"] or ref_note_mw

    n_base_note = (ref_row["effective_mass_mg"] / ref_note_mw) / ref_eq
    n_base_rdkit = (ref_row["effective_mass_mg"] / ref_rd_mw) / ref_eq

    audited_rows = []
    theoretical_products = []
    for item in enriched:
        eq = float(item.get("eq", 1.0))
        mw_note = float(item["note_mw"]) if item.get("note_mw") else item["rdkit_mw"]
        mw_rd = item["rdkit_mw"] or mw_note

        if mw_note is not None:
            exp_mmol_note = n_base_note * eq
            exp_mass_note = exp_mmol_note * mw_note
            item["expected_mmol_by_note_mw"] = round(exp_mmol_note, 4)
            item["expected_mass_mg_by_note_mw"] = round(exp_mass_note, 2)
        if mw_rd is not None:
            exp_mmol_rd = n_base_rdkit * eq
            exp_mass_rd = exp_mmol_rd * mw_rd
            item["expected_mmol_by_rdkit_mw"] = round(exp_mmol_rd, 4)
            item["expected_mass_mg_by_rdkit_mw"] = round(exp_mass_rd, 2)

        if item["effective_mass_mg"] is not None and mw_note is not None:
            delta_mg = item["effective_mass_mg"] - item["expected_mass_mg_by_note_mw"]
            rel_pct = (delta_mg / item["expected_mass_mg_by_note_mw"] * 100.0) if item["expected_mass_mg_by_note_mw"] else 0.0
            item["written_vs_expected_note_mg_delta"] = round(delta_mg, 3)
            item["written_vs_expected_note_rel_pct"] = round(rel_pct, 3)

        if item.get("role") == "product" and mw_rd is not None:
            theo_mg = n_base_rdkit * eq * mw_rd
            theo_g = theo_mg / 1000.0
            prod_summary = {
                "name": item.get("name"),
                "theoretical_100pct_mg": round(theo_mg, 2),
                "theoretical_100pct_g": round(theo_g, 4),
            }
            actual_g = payload.get("actual_product_g")
            if actual_g is not None and theo_g > 0:
                yp = round((float(actual_g) / theo_g) * 100.0, 2)
                prod_summary["actual_g"] = float(actual_g)
                prod_summary["yield_pct"] = yp
                if yp > 100.0:
                    prod_summary["anomaly_warning"] = (
                        f"产率 {yp}% (>100%)：优先排查高沸点溶剂残留（如 DMSO/DMF/甲苯）或未洗净无机盐。"
                    )
            theoretical_products.append(prod_summary)
        audited_rows.append(item)

    solvent_ml = payload.get("solvent_vol_ml")
    conc_m = round(n_base_rdkit / float(solvent_ml), 4) if solvent_ml else None

    return {
        "ref_name": ref_name,
        "base_1eq_mmol_by_note_mw": round(n_base_note, 6),
        "base_1eq_mmol_by_rdkit_mw": round(n_base_rdkit, 6),
        "solvent_vol_ml": solvent_ml,
        "substrate_molarity_M": conc_m,
        "rows": audited_rows,
        "theoretical_products": theoretical_products,
    }


def calculate_stoichiometry(mw: float, mmol: float, actual_g: Optional[float] = None) -> dict:
    theory_g = round(mmol * 1e-3 * mw, 4)
    res: Dict[str, Any] = {
        "mw": mw,
        "mmol": mmol,
        "theoretical_yield_g": theory_g,
    }
    if actual_g is not None:
        yield_pct = round((actual_g / theory_g) * 100.0, 2) if theory_g > 0 else 0.0
        res["actual_yield_g"] = actual_g
        res["yield_percentage"] = yield_pct
        if yield_pct > 105.0:
            res["anomaly_warning"] = (
                f"产率达到 {yield_pct}% (>100%)：极可能是高沸点溶剂残留（如甲苯、DMF、DMSO）或无机盐未洗净。"
            )
        elif yield_pct < 30.0:
            res["notice"] = f"产率较低 ({yield_pct}%)，建议排查水相未萃尽、反应未完全或副反应。"
        else:
            res["status"] = f"产率处于合理合成区间 ({yield_pct}%)。"
    return res


def main() -> None:
    parser = argparse.ArgumentParser(description="Antigravity Chemistry Self-Correction & Preview CLI (v2.0)")
    sub = parser.add_subparsers(dest="command")

    # 1. prep
    p_prep = sub.add_parser("prep", help="Rotate upright, enhance contrast, and slice high-res crops")
    p_prep.add_argument("image_path", help="Source image path")
    p_prep.add_argument("--out-dir", required=True, help="Output directory (e.g. scratch/prep)")
    p_prep.add_argument("--rotate", type=int, default=0, choices=[0, 90, 180, 270], help="CCW rotation degrees")
    p_prep.add_argument("--enhance", action="store_true", help="Also write CLAHE contrast-enhanced image")
    p_prep.add_argument("--crop", help="Optional normalized crop box x1,y1,x2,y2 (0..1)")

    # 2. sanitize
    p_san = sub.add_parser("sanitize", help="Validate SMILES, audit MW vs note, detect motifs, and render PNG+SVG")
    p_san.add_argument("smiles", help="Input SMILES string")
    p_san.add_argument("--note-mw", type=float, help="Handwritten MW to audit against RDKit MW")
    p_san.add_argument("--out-dir", help="Directory to save PNG and SVG previews")
    p_san.add_argument("--name", default="mol", help="Base filename for PNG/SVG previews")
    p_san.add_argument("--width", type=int, default=420)
    p_san.add_argument("--height", type=int, default=280)

    # 3. render-reaction
    p_rxn = sub.add_parser("render-reaction", help="Render reaction scheme PNG+SVG and check elemental delta")
    p_rxn.add_argument("--reactants", nargs="+", required=True, help="Reactant SMILES")
    p_rxn.add_argument("--products", nargs="+", required=True, help="Product SMILES")
    p_rxn.add_argument("--agents", nargs="*", default=[], help="Optional reagent/catalyst SMILES")
    p_rxn.add_argument("--out-dir", required=True, help="Directory to save reaction PNG and SVG")
    p_rxn.add_argument("--name", default="reaction", help="Base filename")
    p_rxn.add_argument("--width", type=int, default=1000)
    p_rxn.add_argument("--height", type=int, default=300)

    # 4. audit-table
    p_tbl = sub.add_parser("audit-table", help="Audit full stoichiometry table from JSON")
    p_tbl.add_argument("--json-file", help="Path to JSON payload")
    p_tbl.add_argument("--json-str", help="Inline JSON string")

    # 5. stoichiometry
    p_st = sub.add_parser("stoichiometry", help="Quick theoretical yield and >100% anomaly check")
    p_st.add_argument("--mw", type=float, required=True)
    p_st.add_argument("--mmol", type=float, required=True)
    p_st.add_argument("--actual-g", type=float)

    args = parser.parse_args()

    if args.command == "prep":
        res = prep_image(args.image_path, args.out_dir, rotate_deg=args.rotate, enhance=args.enhance, custom_crop=args.crop)
    elif args.command == "sanitize":
        res = sanitize_structure(args.smiles, note_mw=args.note_mw, out_dir=args.out_dir, name=args.name, width=args.width, height=args.height)
    elif args.command == "render-reaction":
        res = render_reaction(args.reactants, args.products, args.out_dir, agents=args.agents, name=args.name, width=args.width, height=args.height)
    elif args.command == "audit-table":
        if args.json_file:
            payload = json.loads(Path(args.json_file).read_text(encoding="utf-8-sig"))
        elif args.json_str:
            payload = json.loads(args.json_str)
        else:
            print(json.dumps({"error": "Provide --json-file or --json-str"}))
            sys.exit(1)
        res = audit_table(payload)
    elif args.command == "stoichiometry":
        res = calculate_stoichiometry(args.mw, args.mmol, args.actual_g)
    else:
        parser.print_help()
        return

    print(json.dumps(res, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

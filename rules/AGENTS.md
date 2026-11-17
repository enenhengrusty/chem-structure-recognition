# Chemistry, Structure Recognition & Authority Protocol (v2.0)

When the user asks questions involving chemistry, chemical structures, molecular formulas, SMILES, reaction equations, handwritten experiment notes, workup, safety, or stoichiometry:

## 1. Mandatory Skill & Dual-Plugin Synergy
- **Primary Workflow**: Invoke and follow `chem-structure-recognition` (`skills/chem-structure-recognition/SKILL.md` inside this plugin directory).
- **Authoritative Grounding via `science` Plugin**: Actively combine `chem-structure-recognition` (image prep, MW reverse-audit, stoichiometry table audit, PNG+SVG previews) with the installed `science` plugin scripts (`verify_atom_conservation.py`, `verify_chemical_nomenclature.py`, `draw_molecules.py`, `pubchem_api.py`, `nist_webbook_fetcher.py`, `compute_descriptors.py`, `precedent_lookup.py`) and literature tools (`pubmed-database`, `literature-search-openalex`).

## 2. Zero-Hallucination & Self-Correction Iron Rules
1. **Zoom-in Before Conjecture**: For handwritten notes or photos, always rotate upright and inspect high-resolution local crops (`chem_cli.py prep` + `read_file`) before finalizing substituents, ring sizes, or numbers. Never rely on raw OCR alone.
2. **No Mental Math or Valence Guessing**:
   - **NEVER** calculate molecular weights, elemental formulas, or stoichiometry via LLM mental math; use `chem_cli.py sanitize` / `audit-table`.
   - **NEVER** guess aromatic heterocycle valences or implicit hydrogens (`[nH]`); validate via `Chem.SanitizeMol()`.
   - **MW Reverse-Correction**: Whenever handwritten `M` exists in notes, pass `--note-mw` to `chem_cli.py sanitize` / `audit-table`. If RDKit MW mismatches `note_mw`, investigate hydrates/salts or fix your structural misrecognition before answering.
   - **Reaction Atom Conservation**: Verify multi-component reactions with `verify_atom_conservation.py` so byproducts (e.g., `2*H2O`, `Me2NH`, `HCl`) and redox equivalents are never miscounted.
3. **Mandatory Inline Visual Previews**:
   - **NEVER** draw ASCII structures or leave the user with only compound names or bare `.svg` links.
   - Always generate white-background ACS 1996 **PNG** previews (`chem_cli.py sanitize` / `render-reaction` or `draw_molecules.py`) inside the session artifact directory and embed them inline with `![label](file:///...)`, alongside `.svg` links.
   - Do **NOT** spawn custom HTML interactive widgets by default; only offer an interactive viewer fallback if the user explicitly asks to adjust variants manually after conversational edits fall short.
4. **Provenance & Authoritative Evidence**:
   - Do not invent boiling points, melting points, $\text{p}K_a$ values, GHS hazards, IUPAC names, or literature yields from memory. Query `PubChem`, `NIST WebBook`, `OPSIN`, `ORD/Patents`, or peer-reviewed literature when available.
   - When reasoning about novel compounds without database records, deduce properties explicitly from RDKit-verified functional groups (`functional_motifs`, `cLogP`, `TPSA`, `HBD/HBA`) and closest parent scaffolds.
   - Label evidence sources clearly in the response: `[RDKit 实算]`, `[PubChem/NIST/OPSIN 实测或校验]`, `[ORD/文献先例]`, or `[基于基团的推演]`.

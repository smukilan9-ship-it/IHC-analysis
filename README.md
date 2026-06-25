# IHC Analyzer — *Spatial biology for the price of a stain*

Recover cell-resolution spatial immune information from ordinary **brightfield
immunohistochemistry** — the ~$10 stain every pathology lab already runs — instead of a
six-figure multiplex machine, with an instrument designed so it **cannot report a
conclusion it cannot defend**.

The analysis core is deterministic (QuPath + InstanSeg for segmentation; classical spatial
statistics) and every sensitive conclusion is gated behind an explicit validity check:
uncertified registration, uncertified cell correspondence, and shared-tissue-preference
artifacts are **blocked or flagged, never silently reported**.

## The three pipelines

1. **Quantification** — per-image DAB-positive cell counting from brightfield IHC. The
   validated, deterministic core.
2. **Serial-section Spatial Association** — for two markers on *adjacent* sections:
   landmark-certified registration + intensity-reweighted cross-type Ripley's K against a
   calibrated null + a global DCLF significance test. A **population-level** statistic — it
   does **not** assert single-cell co-expression (serial sections are different Z-planes).
3. **Same-section Restained Co-expression** — for stains imaged on the *same* physical
   section (stain → image → strip → restain → image): segment once, reuse the exact cell
   coordinates across captures, and call genuine **single-cell** co-expression with a 2×2
   contingency test, gated behind manual coordinate-correspondence certification.

## Architecture

```text
Raw brightfield images
  -> QuPath headless + InstanSeg (brightfield_nuclei)
  -> CSV / GeoJSON / JSON exports
  -> Python: overlays, dashboard, Excel, spatial association, restained co-expression
  -> pywebview desktop UI (Quantification / Spatial / Restained tabs)
```

## Requirements

- Python 3.10+
- QuPath 0.7.x with the InstanSeg extension
- InstanSeg `brightfield_nuclei-0.1.1` model downloaded locally
- macOS is the currently targeted desktop environment

## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cp .env.example .env
cp config.example.yaml config.yaml   # then edit local paths
```

Edit `config.yaml` with local paths for `input_dir`, `output_dir`, `dashboard_dir`,
`qupath_binary`, and `instanseg_model`. The desktop app stores user setup in
`~/.ihc_analyzer/`. An optional in-app results-chat assistant can use an LLM provider key
(`GEMINI_API_KEY` or `ANTHROPIC_API_KEY`); it is unrelated to the deterministic analysis core.

## Usage

```bash
python app.py                                   # desktop UI (all three tabs)
python run_pipeline.py --config config.yaml     # quantification
python run_pipeline.py --config config.yaml --mode spatial   # serial-section spatial association
```

## Validation

The statistics and certification machinery are validated separately from the biological
datasets. Harnesses live in `validation/` (cross-K vs. brute force; DCLF calibration; the
reweighted-null 500-run size/power calibration; cross-validation against R `spatstat`;
registration certification on ANHIR/CIMA expert landmarks; real-data CODEX controls;
restained segmentation against expert masks). The consolidated results table is in
**`ihc.md` → Appendix A**.

**Honest scope:** the *method* is validated; a finished CD8/TIM-3 biological result is not —
no serial-section pair has yet passed registration certification on the target cohort, which
is a data-acquisition limit, not a software one. See the documents below.

## Documentation

- **`moonshot_paper.md`** — the problem, the first-principles insight, the technical
  foundations, the validation, and the future.
- **`ihc.md`** — full technical reference + consolidated validation ledger (Appendix A) +
  corrections log (Appendix B).
- **`learn.md`** — plain-language companion explaining the whole project with no assumed
  background in pathology, microscopy, or statistics.

## Resources

- Vision deck (~50 slides): https://canva.link/agmetzji7s05ify
- Supporting materials — outputs, overlays, segmentation validation (Google Drive):
  https://drive.google.com/drive/folders/1N3wTEH9Won0i12BUm7qTa2qOzFoo7s2J?usp=share_link

## Credits

Built on QuPath, the InstanSeg `brightfield_nuclei` model, and the scientific Python stack
(NumPy, SciPy, scikit-image, OpenCV, Shapely, SimpleITK, matplotlib, pywebview). Validated
against, and crediting, the Schürch et al. 2020 colorectal CODEX dataset, the ANHIR/CIMA
expert-landmark benchmark, the HNSCC mIF/mIHC comparison dataset, and the spatial-statistics
literature (Ripley's K; Baddeley–Møller–Waagepetersen inhomogeneous K; the
Diggle–Cressie–Loosmore–Ford envelope test). Cross-validated against R `spatstat`.

## License

MIT — see [LICENSE](LICENSE).


# Quantum Bioenergetics Mapping (Fixed)

Drop‑in repo you can push to GitHub and run locally or on Streamlit.

## Quick Start

```bash
pip install -r requirements.txt
python QB_6_main_pipeline_fixed.py
```

Outputs appear in `outputs/`.

## Streamlit (optional)

```bash
streamlit run streamlit_app.py
```

This shows the ENAQT curve and lets you tweak parameters.

## Notes

- Real TCGA/GEO loading stubs are present; the pipeline uses synthetic data unless you wire in GDC/GEO.
- The ENAQT physics uses a Lindblad master equation with sink flux for ETE.

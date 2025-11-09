
# Quantum Bioenergetics Mapping (ENAQT-validated)

## Quick Start
```
pip install -r requirements.txt
python QB_6_main_pipeline_fixed.py
```

## Streamlit
```
streamlit run streamlit_app.py
streamlit run streamlit_upload.py   # upload your own omics
```

## CLI
```
python qemd.py simulate --expr ./data/templates/sample_omics.csv --graph ./data/templates/graph.json --mapping ./data/templates/mapping.csv --out ./outputs/metrics.json
python qemd.py cohort --cohort ./data/templates/cohort.csv --graph ./data/templates/graph.json --mapping ./data/templates/mapping.csv --out ./outputs/cohort_metrics.parquet
python analyze_cohort.py --parquet ./outputs/cohort_metrics.parquet
```

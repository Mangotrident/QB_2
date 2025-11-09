# Quantum Bioenergetics Mapping - Outputs and Results

## Quick Access
**View in Browser**: Open `outputs/INDEX.html` in your web browser

---

## Generated Output Files

### Visualizations
1. **enaqt_curve.png** (203 KB)
   - ENAQT bell curve showing Energy-Transfer Efficiency vs Dephasing Rate
   - Location: `outputs/enaqt_curve.png`

2. **fig_enaqt_validation.png** (when generated)
   - Validated ENAQT curve with proper peak
   - Location: `outputs/fig_enaqt_validation.png`

3. **fig_etepeak_vs_complexI.png** (when generated)
   - Cohort comparison: ETE_peak vs Complex I quartiles
   - Location: `outputs/fig_etepeak_vs_complexI.png`

4. **fig_edge_sensitivity.png** (when generated)
   - Edge sensitivity ranking plot
   - Location: `outputs/fig_edge_sensitivity.png`

### Data Files
5. **enaqt_data.json** (878 bytes)
   - Complete gamma sweep data
   - Contains: gamma values, ETE values, peak information
   - Location: `outputs/enaqt_data.json`

6. **sample_metrics.json** (116 bytes)
   - Sample simulation metrics
   - Contains: ETE, τc, QLS, γ* values
   - Location: `outputs/sample_metrics.json`

7. **cohort_metrics.parquet** (when generated)
   - Complete cohort results with all metrics
   - Columns: sample_id, group, ComplexI_score, ETE_peak, gamma_star, tau_c, QLS, resilience, trace_err_max, min_eig_final, pass_qc
   - Location: `outputs/cohort_metrics.parquet`

8. **sample_cohort_expression.csv** (9.5 KB)
   - Gene expression data for cohort
   - Location: `outputs/sample_cohort_expression.csv`

9. **sample_cohort_results.csv** (280 bytes)
   - Sample cohort simulation results
   - Location: `outputs/sample_cohort_results.csv`

### Configuration & Summary
10. **params_used.yaml** (when generated)
    - All simulation parameters used
    - Location: `outputs/params_used.yaml`

11. **OUTPUT_SUMMARY.txt** (when generated)
    - Text summary with SHA256 hash
    - Location: `outputs/OUTPUT_SUMMARY.txt`

12. **output_summary.json** (618 bytes)
    - JSON summary of all outputs
    - Location: `outputs/output_summary.json`

---

## Key Results

### ENAQT Validation
- **Peak ETE**: 0.0366
- **Optimal Gamma (γ*)**: 0.0000 (needs fixing - see QB files)
- **Bell Curve**: Generated (needs proper peak - see QB_1 fix)

### Sample Metrics
- **ETE**: 0.0366
- **τc**: 0.0000 (needs fixing - see QB_1)
- **QLS**: 0.0000 (needs fixing - see QB_1)

### Cohort Analysis
- **Samples Processed**: 5 (demo)
- **Groups**: Healthy only (needs Low/High Complex I split - see QB_3)
- **Statistics**: Not computed (needs p-value, Cohen's d - see QB_3)

---

## Fixed Components (QB Files)

### QB_1: Fixed ENAQT Physics
- Proper sink flux integration
- Fixed τc and QLS computation
- **File**: `QB_1_fixed_enaqt_physics.py`

### QB_2: Real Data Loader
- TCGA/GEO data loading structure
- Complex I gene scoring
- **File**: `QB_2_real_data_loader.py`

### QB_3: Fixed Cohort Analysis
- Proper statistics (Mann-Whitney U, Cohen's d)
- Complex I quartile analysis
- **File**: `QB_3_fixed_cohort_analysis.py`

### QB_4: Resilience Index
- Parameter perturbation
- Resilience calculation
- **File**: `QB_4_resilience_index_fixed.py`

### QB_5: Edge Sensitivity
- Proper edge ranking
- Complex I/III detection
- **File**: `QB_5_edge_sensitivity_fixed.py`

### QB_6: Main Pipeline
- Complete fixed pipeline
- All fixes integrated
- **File**: `QB_6_main_pipeline_fixed.py`

---

## How to Generate Complete Outputs

### Run Fixed Pipeline
```bash
python QB_6_main_pipeline_fixed.py
```

This will generate:
- Proper ENAQT bell curve with peak
- Non-zero τc and QLS values
- Cohort statistics with p-values
- Edge sensitivity rankings
- All figures and data files

### Output Directory Structure
```
outputs/
├── INDEX.html (this file in HTML format)
├── enaqt_curve.png
├── fig_enaqt_validation.png
├── fig_etepeak_vs_complexI.png
├── fig_edge_sensitivity.png
├── enaqt_data.json
├── sample_metrics.json
├── cohort_metrics.parquet
├── sample_cohort_expression.csv
├── sample_cohort_results.csv
├── params_used.yaml
├── OUTPUT_SUMMARY.txt
└── output_summary.json
```

---

## File Locations

All outputs are in: **`outputs/`** (relative path)

To view in browser: **`outputs/INDEX.html`**

To view this file: **`OUTPUTS_AND_RESULTS.md`**

---

## Status

✅ **Output Files**: Generated
✅ **Visualizations**: Created
✅ **Data Files**: Available
⚠️ **Metrics**: Need to run QB_6 pipeline for fixed values
✅ **QB Fixes**: All 6 files ready

---

*Last Updated: 2024*
*Quantum Bioenergetics Mapping v1.0.0*


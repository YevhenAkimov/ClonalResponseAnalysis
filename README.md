# ClonalResponseAnalysis

> Growth-rate analysis for clonal/barcode tracing experiments.

---

## What it does

The pipeline is designed for the multi-readout bracode-tracing experiments. Typically, it involves the lentiviral DNA barcode delivery and pool expansion. Expanded pool is devided into treatment/control samples.
<img width="665" height="309" alt="image" src="https://github.com/user-attachments/assets/1ef05e09-123b-4e7f-9292-08830bf6a84a" />


Input is a raw barcode counts table for a set of sample. The algorithm uses a weighted mixed-effects model to estimates **per-lineage growth rates** for each treatment. The cGR score is per-lineage perturbation/treatment effect estimated via constrained-lasso step. It residualizes each treatment which isolates drug-specific effects from baseline drift.

Per lineage × condition you get:

- **`growth_rate`** — absolute growth rate under treatment
- **`centered_growth_rate`** — growth rate of each lineage minus population growth rate
- **`growth_difference`** — treatment vs vehicle growth rate difference (negative = sensitive, positive = resistant)
- **`cGR`** — vehicle-residualized centered growth rate (the cleanest per-lineage treatment effect)
- **`FC`**, **`fraction`**, **`fraction_control`**, **`fraction_t0`**, **`p_value`**

---

## Quick start

```r
source("ResponseAnalysis.R")

# counts:   matrix — rows = lineages (barcodes), cols = samples
# metadata: data.frame — see schema below

obj <- ClonalAnalysis$new(metadata, counts)
obj$fit_model()

res <- obj$results()
head(res$result)       # per-lineage × condition metrics (long format)

res$analysis_info      # which samples were used as baseline / vehicle per rep
res$last_meta          # the exact metadata slice used for this fit
```

For a wide/matrix layout (one column per condition):

```r
res <- obj$results(format = "list")
res$result$cGR         # lineages × conditions matrix for cGR
```

To combine multiple fits (e.g. different drug panels):

```r
merge_clonal_results(res_A, res_B, id = c("panel_A", "panel_B"))
```

---

## Metadata table

To streamline the analysis we use an automatic sample matching via the metadata table.

### Required columns

| Column | What | Details |
| :--- | :--- | :--- |
| `sample_name` | Unique sample ID | Must **exactly match** a single column name in the counts matrix |
| `treatment` | `control0`, `control1`, or a drug/condition label | `control0` = baseline, `control1` = vehicle |
| `time` |  Sampling time (numeric) | Consistent units across the table |
| `rep` | Replicate ID | Groups samples from the same experimental run |
| `fold_expansion` | Bulk cell-count fold change vs seeding | `> 0`, finite. Skip if you supply `growth_rate` directly |


### Rules

- Each `rep` needs **≥1 `control0`** and **≥1 `control1`**.
- Non-`control0` samples must be sampled **after** the replicate's baseline.
- `fold_expansion` must be `> 0` for every non-`control0` row.

---

##  Baseline (t0) and vehicle (t1) selection

| Sample type | How it's chosen | Why |
| :--- | :--- | :--- |
| **Baseline** (`control0`) | Earliest `control0` per replicate | All lineage log-fold-changes are measured against it |
| **Vehicle** (`control1`) | One per replicate, closest in time to the average treatment timepoint | Keeps the vehicle a fair temporal comparator; extras are dropped |

---

##  Minimal example

```text
sample_name    treatment   time   rep   fold_expansion
S1_baseline    control0    0      R1    NA
S1_dmso_14     control1    14     R1    12.0
S1_drugA_14    drugA       14     R1    3.5
S2_baseline    control0    0      R2    NA
S2_dmso_14     control1    14     R2    11.5
S2_drugA_14    drugA       14     R2    3.8
```

One baseline + one vehicle + one treatment per replicate. Anything more elaborate is just more rows obeying the same rules.

---

<details>
<summary><b> Narrowing or reshaping a fit without touching metadata</b></summary>

<br>

`fit_model()` takes three optional arguments that subset or relabel samples **for a single run**. Your underlying metadata is never modified.

- **`control0 = c(...)`** — restrict the baseline to specific samples. Useful when metadata carries multiple baselines (batches, passages) and you want this fit pinned to one.

- **`control1 = c(...)`** — restrict the vehicle pool the closest-in-time picker can choose from. Useful when vehicles span many timepoints or batches.

- **`treatments = list(new_label = c(...))`** — selects *and* relabels treatment samples. The list names become the condition labels the model sees. Typical uses:
  - **Pool** batches under one label — `list(drugA = c("drugA_b1", "drugA_b2"))`
  - **Split** a label into subgroups — e.g. early vs late samples of the same drug
  - **Filter** — anything not listed is dropped from this fit

</details>

# ClonalResponseAnalysis

> Lineage-level growth-rate analysis for clonal/barcode tracing experiments.

---

## 📋 Metadata table

The analysis is driven by a single metadata table that pairs each sample with its experimental context. Everything downstream — baseline, vehicle matching, growth rates — is derived from it.

### Required columns

| Column | What it holds | Notes |
| :--- | :--- | :--- |
| `sample_name` | Unique sample ID | Must **exactly match** a column in the counts matrix |
| `treatment` | `control0`, `control1`, or a drug/condition label | `control0` = baseline, `control1` = vehicle |
| `time` | Numeric sampling time | Consistent units across the table |
| `rep` | Replicate ID | Groups samples from the same experimental run |
| `fold_expansion` | Bulk cell-count fold change vs seeding | `> 0`, finite. Skip if you supply `growth_rate` directly |

> [!NOTE]
> Derived automatically: `condition`, `delta_t`, `growth_rate`, baseline time.

### Rules

- Each `rep` needs **≥1 `control0`** and **≥1 `control1`**.
- Non-`control0` samples must be sampled **after** the replicate's baseline.
- `fold_expansion` must be `> 0` for every non-`control0` row.

---

## 🎯 Baseline (t0) and vehicle (t1) selection

| Sample type | How it's chosen | Why |
| :--- | :--- | :--- |
| **Baseline** (`control0`) | Earliest `control0` per replicate | All lineage log-fold-changes are measured against it |
| **Vehicle** (`control1`) | One per replicate, closest in time to the average treatment timepoint | Keeps the vehicle a fair temporal comparator; extras are dropped |

---

## 🧪 Minimal example

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
<summary><b>🔧 Narrowing or reshaping a fit without touching metadata</b></summary>

<br>

`fit_model()` takes three optional arguments that subset or relabel samples **for a single run**. Your underlying metadata is never modified.

- **`control0 = c(...)`** — restrict the baseline to specific samples. Useful when metadata carries multiple baselines (batches, passages) and you want this fit pinned to one.

- **`control1 = c(...)`** — restrict the vehicle pool the closest-in-time picker can choose from. Useful when vehicles span many timepoints or batches.

- **`treatments = list(new_label = c(...))`** — selects *and* relabels treatment samples. The list names become the condition labels the model sees. Typical uses:
  - **Pool** batches under one label — `list(drugA = c("drugA_b1", "drugA_b2"))`
  - **Split** a label into subgroups — e.g. early vs late samples of the same drug
  - **Filter** — anything not listed is dropped from this fit

</details>

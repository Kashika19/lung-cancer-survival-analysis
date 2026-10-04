# Dataset guidance

The source dataset is intentionally not included in this public repository because its provenance and redistribution licence have not been confirmed.

Place an authorised local copy at:

```text
data/raw/Lung Cancer.csv
```

The analysis normalises column names to lowercase snake case. It requires a binary `survived` target with values such as `Yes`/`No`, `1`/`0`, or `TRUE`/`FALSE`.

The academic script also references these optional fields:

| Field | Intended use |
|---|---|
| `id` | Identifier; excluded from modelling |
| `age` | Numeric predictor and age band |
| `gender` | Categorical predictor |
| `smoking_status` | Categorical predictor |
| `bmi` | Numeric predictor and BMI band |
| `cholesterol_level` | Numeric predictor and cholesterol band |
| `cancer_stage` | Categorical predictor |
| `treatment_type` | Categorical predictor |
| `country` | Categorical predictor |
| `family_history` | Categorical predictor |
| `diagnosis_date` | Used with treatment end date, then excluded |
| `end_treatment_date` | Used to derive treatment duration, then excluded |
| `asthma` | Optional comorbidity indicator |
| `cirrhosis` | Optional comorbidity indicator |
| `other_cancer` | Optional comorbidity indicator |
| `survived` | Required binary target |

Do not commit patient-identifiable, confidential, or unlicensed data. The repository `.gitignore` excludes `data/raw/` and `data/private/`.


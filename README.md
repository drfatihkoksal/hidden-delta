# HIDDEN-DELTA

Code for a preregistered matched case-control test of whether one 12-lead ECG without pre-excitation identifies patients who have ventricular pre-excitation on other ECGs (intermittent pre-excitation).

- **Preregistration:** https://doi.org/10.17605/OSF.IO/GCVJT. The frozen analysis plan (version 1.2) is in `frozen/`.
- **Model freeze:** `FROZEN_MANIFEST.sha256` lists the plan, the code, the final models with their calibration and threshold, and the matched test sets, as frozen before any test-set ECG was scored. The SHA-256 of the manifest itself is `a16e66cfcc36b878c2c572454baeb0fc8d607237316f11ab06ba23b51948b264` and was added to the registration as an update. Code entries can be checked with `sha256sum -c FROZEN_MANIFEST.sha256 2>/dev/null | grep -v '^/'`; model and data entries point to local files that are not shared.

## Data

The data are not included and cannot be shared by the author. They are available to credentialed users under the providers' data use agreements:

| Source | Version | Role | Access |
|---|---|---|---|
| Harvard-Emory ECG Database (HEEDB), Massachusetts General Hospital (I0001) | 5.0 (https://doi.org/10.60508/765d-f248) | Development | Brain Data Science Platform |
| HEEDB, Emory (I0006) | 5.0 | Second test | Brain Data Science Platform |
| MIMIC-IV-ECG | 1.0 (https://doi.org/10.13026/4nqg-sb35) | Primary test | PhysioNet |
| MIMIC-IV | 3.1 (https://doi.org/10.13026/kpb9-mt58) | Diagnosis and procedure codes | PhysioNet |
| MIMIC-IV-Note | 2.2 (https://doi.org/10.13026/1n74-ne17) | AVNRT subgroup (fixed text patterns; notes are never quoted) | PhysioNet |

This repository holds code and aggregate results only. Paths to the local copies of the data and to intermediate files are set as constants at the top of the scripts (for example `/home/fatih/data/hidden-delta`); edit them for another machine.

## Pipeline

| Folder | What it does | Main outputs (aggregate) |
|---|---|---|
| `stage0/` | Counts of candidate patients and ECG statements; diagnosis and procedure code flags | `stage0/out/` |
| `stage1/` | Labels: pre-excitation in at least 2 of 3 HEEDB sources (12SL version 24, acquisition statement, physician overread); definite machine statements in MIMIC-IV-ECG; pre-excitation-free ECGs; cohorts | `stage1/out/labels_summary.*` |
| `stage2/` | Signal processing (`hdsig.py`): median-filter baseline removal, beat detection and clustering, median beats, QRS and P boundaries, features, the discrete pre-QRS wave; device pre-ringing from pacemaker spikes; checks against machine intervals | `stage2/out/` |
| `stage3/` | Matching of cases to controls (exact strata and calipers); removal of patients read in earlier studies from the test sets | `stage3/out/match_*` |
| `stage4/` | Development in MGH only: feature model (M0), self-supervised 1D ResNet on median beats (M1, primary), M1 embedding plus features (M2), demographic model, raw 10-s model; 5-fold patient-level cross-validation; final models | `stage4/out/` |
| `stage5/` | Test preparation, freeze manifest, single scoring run, registered analyses, Holm correction | `stage5/out/results.md`, `results_holm.json` |

Code comments are in Turkish.

### Order of runs

```
python3 stage0/s0_mimic.py && python3 stage0/s0_heedb.py
python3 stage1/s1_labels.py
python3 stage2/s2_extract.py <site>            # I0001, I0006, MIMIC
python3 stage3/s3_match.py <site>
python3 stage4/s4_data.py
python3 stage4/s4_m0.py && python3 stage4/s4_m1.py && python3 stage4/s4_m2.py
python3 stage4/s4_m1raw.py build && python3 stage4/s4_m1raw.py train
python3 stage5/s5_prepare.py MIMIC I0006
python3 stage5/s5_prepare.py agg MIMIC I0006 && python3 stage5/s5_prepare.py raw MIMIC I0006
python3 stage5/s5_freeze.py                    # writes the manifest; refuses to overwrite it
python3 stage5/s5_score.py MIMIC && python3 stage5/s5_score.py I0006   # checks the manifest; runs once
python3 stage5/s5_analysis.py site MIMIC && python3 stage5/s5_analysis.py site I0006
python3 stage5/s5_analysis.py combine
```

`s5_score.py` and `s5_analysis.py` also have a rehearsal mode (`--rehearsal`) that uses randomly initialised models and, in the test sets, replaces the case in each matched set with a random member.

On one RTX 5090 GPU, training M1 took about 9 minutes and the raw 10-s model about 15 minutes; scoring and the registered analyses take seconds.

### Analyses written after the registered run (post hoc)

| Script | What it adds |
|---|---|
| `stage5/s5_supplement.py` | Descriptive outputs named in the plan but missing from the registered analysis code: attribution share in the P wave, discrete-wave polarity, lead profile and frequency by recording filter, AUROC difference by filter setting, concealed-pathway counts |
| `stage5/s5_self_similarity.py` | The initial-QRS agreement measure between two pre-excitation-free ECGs of the same person, as a reference for the registered own-pathway test |
| `stage5/s5_before_after.py` | AUROC by timing of the index ECG relative to the first pre-excited ECG |
| `stage4/s4_label_strength.py` | AUROC by the number of label sources that stated pre-excitation, and in ECGs with strict PR and QRS limits |
| `stage5/s5_by_sex.py` | AUROC in women and men (matched sets are exact on sex) |

## Environment

Python 3.13 with the packages in `requirements.txt`; PyTorch with CUDA for stages 4 and 5.

## Licence and citation

MIT licence (`LICENSE`). Citation details are in `CITATION.cff`.

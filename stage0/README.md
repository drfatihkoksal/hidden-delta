# Aşama 0 — fizibilite sayımı

concept.md §10'daki sayımlar. Pipeline'dan önce G1 kararına veri sağlar.

| Script | Kaynak | Çıktı |
|---|---|---|
| `s0_mimic.py` | MIMIC-IV-ECG 1.0, MIMIC-IV 3.1 hosp, MIMIC-IV-Note | `out/mimic_stage0.{md,json}` |
| `s0_heedb.py` | HEEDB 5.0 (I0001 MGH, I0006 Emory) | `out/heedb_stage0.{md,json}` |

Ortak yollar, ifade kalıpları ve 12SL/ICD kod kümeleri `common.py`'de.

```
python3 stage0/s0_mimic.py   # ~1 dk
python3 stage0/s0_heedb.py   # MGH metadata 3 GB; birkaç dakika
```

**DUA:** `out/` yalnızca agregat sayı içerir. Hasta düzeyindeki ara listeler (`*_case_candidates`,
`*_icd_only_pool`, `*_candidate_ecgs`, `mimic_pathway_ablation`, `mimic_k2`) proje dışında,
`~/data/hidden-delta/stage0/` altında tutulur ve makineden çıkmaz. MIMIC-IV-Note metni yalnızca
kalıp eşleştirmede kullanılır; hiçbir çıktıya girmez.

## Tanımlar (makine raporu / kod düzeyi; uzman doğrulaması öncesi)

- **Preeksitasyon EKG'si.** MIMIC'te rapor metninde `WOLFF|WPW|PRE-?EXCIT|DELTA WAVE` geçiyorsa. HEEDB'de 12SL kodu 300/302/303/304 varsa; v24 yeniden işlenmiş kodlar, kayıt anı yazılım kodları ya da hekim overread'i, herhangi biri yeterli.
- **Kısa PR.** Ayrı tutulur (`SHORT PR` / kod 102); çoğu WPW değildir.
- **T0 adayı EKG.** Preeksitasyon ya da kısa PR ifadesi yok, sinüs ritmi, AF/flutter/pace yok. QRS <120 ms olmalı (yalnızca MIMIC'te makine ölçümüyle); her iki veri setinde de tam dal bloğu ve IVCD dışlanır. Yaş ≥18.
- **Katı intermittan aday.** ≥1 preeksitasyon EKG'si + ≥1 T0 adayı EKG'si.
- **Yalnızca-ICD havuzu.** I45.6/426.7 var, hiç preeksitasyon ifadesi yok, ≥1 T0 adayı.
- **Yolak ablasyonu vekili (yalnızca MIMIC).** Ablasyon kodu (02583ZZ/37.34) ile birlikte I45.6/426.7 ya da ablasyon öncesi bir preeksitasyon EKG'si. İlk ablasyon tarihi esas alınır; ablasyonla aynı günkü EKG'ler sayılmaz.
- **K2.** SVT kodu (I47.1x, I47.11 hariç / 427.0); preeksitasyon ifadesi ve I45.6/426.7 yok; ≥1 T0 adayı. MIMIC'te alt tip taburcu özeti kalıplarıyla kabaca ayrılır; olumsuzlama ayıklanmaz.

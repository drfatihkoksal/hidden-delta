# HIDDEN-DELTA — kayıtlı test sonuçları

ANALYSIS_PLAN v1.2; OSF https://osf.io/kvn28. Yalnızca agregat.

## H1 (birincil, MIMIC)

AUROC 0.721 (%95 GA 0.660–0.779); 103 vaka, 461 kontrol. M1 − demografik model +0.201 (+0.143 – +0.259). **H1 karşılandı.**

## H2 ailesi (Holm, tek yönlü α 0.05)

| Site | Test | Tahmin / ortalama r | Ham p | Holm p |
|---|---|---|---|---|
| MIMIC | 5.1_onset | +0.0166 | 0.1449 | 1.0000 |
| MIMIC | 5.1_stt | -0.0567 | 0.9785 | 1.0000 |
| MIMIC | 5.2_onset | -0.0322 | 1.0000 | 1.0000 |
| MIMIC | 5.2_stt | -0.0056 | 0.7946 | 1.0000 |
| MIMIC | 5.3_delta | +0.2283 | 0.0002 | 0.0022 |
| MIMIC | 5.4_tmemory | -0.1797 | 0.9617 | 1.0000 |
| I0006 | 5.1_onset | +0.0262 | 0.0075 | 0.0750 |
| I0006 | 5.1_stt | -0.0519 | 0.9920 | 1.0000 |
| I0006 | 5.2_onset | -0.0420 | 1.0000 | 1.0000 |
| I0006 | 5.2_stt | -0.0058 | 0.8666 | 1.0000 |
| I0006 | 5.3_delta | +0.3210 | 0.0001 | 0.0012 |
| I0006 | 5.4_tmemory | -0.0373 | 0.8453 | 1.0000 |

## Ayrık dalga
- MIMIC: OR 1.138 (0.689–1.880); sıklık vaka 0.650, kontrol 0.653
- I0006: OR 1.614 (1.144–2.278); sıklık vaka 0.538, kontrol 0.436
- Doğrulama (iki sitede de OR > 1, GA 1'i dışlıyor): hayır

## H4 (MIMIC, AVNRT − K1, M1 olasılığı)

Fark +0.0014 (-0.0157 – +0.0207); destekleniyor.

Ayrıntılar: site başına `analysis.json` (ikincil sonlanımlar, denetimler, duyarlılık analizleri).

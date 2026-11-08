# Aşama 3 — eşleştirme (MGH)

concept.md §5.3 "Eşleştirme uygulaması"; ANALYSIS_PLAN §2–3. Yalnızca agregat.

Birincil vaka hastası (uygun indeks T0'ı olan): 1,013; ilk T2'den önce T0'ı olan: 442.

## K1 (1:5)

Havuz 19,920; eşleşen 852 / 1,013 (%84.1); kontrol 3,198; tam (1:5) küme 467.

**Denge** (SMD; hedef |SMD| < 0.1):

| Değişken | Eşleşmemiş tüm | Havuz | SMD önce | Eşleşen | Kontrol (eşleşen) | SMD sonra |
|---|---|---|---|---|---|---|
| yaş | 57.8 | 60.66 | -0.157 | 58.7 | 60.11 | -0.08 |
| KH | 78.16 | 78.33 | -0.008 | 77.74 | 76.5 | 0.063 |
| EKG sayısı | 39.12 | 20.0 | 0.496 | 31.6 | 27.1 | 0.152 |
| dönem (klasör yılı) | 2014.09 | 2015.24 | -0.217 | 2014.68 | 2015.34 | -0.133 |
| kaydırılmış yıl | 2013.78 | 2015.16 | -0.247 | 2014.39 | 2014.99 | -0.112 |
| log(EKG sayısı) | 3.12 | 2.44 | 0.635 | 3.0 | 2.9 | 0.102 |
| erkek (%) | 59.7 | 50.0 | 0.197 | 60.6 | 60.7 | -0.002 |

**Eşleşme sonrası betimsel karşılaştırma** (vaka − küme kontrol ortalaması; küme bootstrap %95 GA):

| Özellik | Vaka | Kontrol | Fark | %95 GA | SMD |
|---|---|---|---|---|---|
| Septal q yok (I, aVL, V5, V6) | 8.6% | 4.7% | +3.7 puan | [+1.6 puan, +5.9 puan] | 0.158 |
| V6'da q var | 71.7% | 77.5% | -6.1 puan | [-9.3 puan, -2.6 puan] | -0.133 |
| Yavaş başlangıç (ms) | 12.2 | 10.9 | +1.22 | [+0.311, +2.06] | 0.114 |
| QRS süresi (ms) | 94.5 | 88.4 | +5.88 | [+4.35, +7.36] | 0.351 |
| PR (ms; +42 ms sapmalı) | 172 | 172 | -0.354 | [-4.05, +3.46] | -0.004 |
| V1 R/S | 0.467 | 0.328 | +0.133 | [+0.00679, +0.278] | 0.089 |
| QRS-T açısı (°) | 85.4 | 76.7 | +8.27 | [+4.9, +11.7] | 0.185 |
| QRS öncesi genlik, düşük eşik (mV) | 0.0276 | 0.0205 | +0.00728 | [+0.00497, +0.00973] | 0.219 |
| QRS öncesi genlik, görünür başlangıç (mV) | 0.0862 | 0.0613 | +0.0235 | [+0.0169, +0.0309] | 0.31 |

## K2 (1:2)

Havuz 2,995; eşleşen 456 / 1,013 (%45.0); kontrol 705; tam (1:2) küme 249.

**Denge** (SMD; hedef |SMD| < 0.1):

| Değişken | Eşleşmemiş tüm | Havuz | SMD önce | Eşleşen | Kontrol (eşleşen) | SMD sonra |
|---|---|---|---|---|---|---|
| yaş | 57.8 | 63.51 | -0.323 | 62.74 | 63.88 | -0.072 |
| KH | 78.16 | 77.23 | 0.043 | 76.69 | 75.03 | 0.088 |
| EKG sayısı | 39.12 | 31.24 | 0.184 | 33.18 | 32.74 | 0.015 |
| dönem (klasör yılı) | 2014.09 | 2016.57 | -0.497 | 2015.84 | 2016.19 | -0.085 |
| kaydırılmış yıl | 2013.78 | 2017.37 | -0.707 | 2015.37 | 2016.3 | -0.204 |
| log(EKG sayısı) | 3.12 | 2.9 | 0.196 | 3.1 | 3.12 | -0.013 |
| erkek (%) | 59.7 | 38.3 | 0.438 | 53.1 | 51.9 | 0.023 |

**Eşleşme sonrası betimsel karşılaştırma** (vaka − küme kontrol ortalaması; küme bootstrap %95 GA):

| Özellik | Vaka | Kontrol | Fark | %95 GA | SMD |
|---|---|---|---|---|---|
| Septal q yok (I, aVL, V5, V6) | 9.2% | 4.4% | +4.8 puan | [+1.8 puan, +7.7 puan] | 0.192 |
| V6'da q var | 68.6% | 75.2% | -6.6 puan | [-11.6 puan, -1.8 puan] | -0.146 |
| Yavaş başlangıç (ms) | 12.2 | 11.8 | +0.456 | [-0.941, +1.65] | 0.031 |
| QRS süresi (ms) | 96.1 | 89.9 | +6.42 | [+4.24, +8.48] | 0.329 |
| PR (ms; +42 ms sapmalı) | 171 | 172 | -1.28 | [-6.54, +4.1] | -0.036 |
| V1 R/S | 0.555 | 0.316 | +0.233 | [+0.0226, +0.492] | 0.13 |
| QRS-T açısı (°) | 88.8 | 75.8 | +13.4 | [+7.97, +18.5] | 0.279 |
| QRS öncesi genlik, düşük eşik (mV) | 0.033 | 0.0201 | +0.0128 | [+0.00518, +0.0249] | 0.157 |
| QRS öncesi genlik, görünür başlangıç (mV) | 0.0844 | 0.0566 | +0.0265 | [+0.0162, +0.0371] | 0.342 |

**Not.** Betimseldir; karar ve hipotez testi Aşama 4'te. HEEDB'de cart/filtre bilgisi yok; dönem, klasör yılından vekil olarak alındı.

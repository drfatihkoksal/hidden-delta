# Aşama 3 — eşleştirme (MIMIC)

concept.md §5.3 "Eşleştirme uygulaması"; ANALYSIS_PLAN §2–3. Yalnızca agregat.

Birincil vaka hastası (uygun indeks T0'ı olan): 109; ilk T2'den önce T0'ı olan: 83.

Test seti: yalnızca denge raporlanır; vaka-kontrol özellik karşılaştırması dondurmadan sonra.

Önceki IAB çalışmalarında EKG'si okunmuş 149 hasta eşleştirmeden önce çıkarıldı (ANALYSIS_PLAN v1.2).

## K1 (1:5)

Havuz 19,681; eşleşen 103 / 109 (%94.5); kontrol 461; tam (1:5) küme 83.

**Denge** (SMD; hedef |SMD| < 0.1):

| Değişken | Eşleşmemiş tüm | Havuz | SMD önce | Eşleşen | Kontrol (eşleşen) | SMD sonra |
|---|---|---|---|---|---|---|
| yaş | 61.88 | 61.38 | 0.028 | 62.01 | 62.97 | -0.053 |
| KH | 80.52 | 78.62 | 0.088 | 79.86 | 77.96 | 0.084 |
| EKG sayısı | 16.99 | 10.65 | 0.4 | 14.56 | 14.41 | 0.014 |
| kaydırılmış yıl | 2158.92 | 2155.03 | 0.157 | 2157.8 | 2157.21 | 0.024 |
| log(EKG sayısı) | 2.45 | 1.95 | 0.563 | 2.38 | 2.39 | -0.012 |
| erkek (%) | 54.1 | 47.8 | 0.126 | 53.4 | 53.8 | -0.008 |

## K2 (1:2)

Havuz 2,149; eşleşen 72 / 109 (%66.1); kontrol 125; tam (1:2) küme 53.

**Denge** (SMD; hedef |SMD| < 0.1):

| Değişken | Eşleşmemiş tüm | Havuz | SMD önce | Eşleşen | Kontrol (eşleşen) | SMD sonra |
|---|---|---|---|---|---|---|
| yaş | 61.88 | 67.05 | -0.299 | 66.6 | 67.68 | -0.077 |
| KH | 80.52 | 82.75 | -0.099 | 77.61 | 77.48 | 0.005 |
| EKG sayısı | 16.99 | 13.45 | 0.199 | 14.88 | 14.7 | 0.016 |
| kaydırılmış yıl | 2158.92 | 2153.62 | 0.215 | 2159.46 | 2156.21 | 0.127 |
| log(EKG sayısı) | 2.45 | 2.03 | 0.415 | 2.46 | 2.43 | 0.037 |
| erkek (%) | 54.1 | 46.0 | 0.163 | 55.6 | 54.4 | 0.023 |

## H4: AVNRT (K2a) → K1 (1:5)

Havuz 19,681; eşleşen 120 / 133 (%90.2); kontrol 543; tam (1:5) küme 98.

**Denge** (SMD; hedef |SMD| < 0.1):

| Değişken | Eşleşmemiş tüm | Havuz | SMD önce | Eşleşen | Kontrol (eşleşen) | SMD sonra |
|---|---|---|---|---|---|---|
| yaş | 68.58 | 61.38 | 0.44 | 67.91 | 67.89 | 0.001 |
| KH | 84.02 | 78.62 | 0.24 | 83.33 | 82.77 | 0.024 |
| EKG sayısı | 17.25 | 10.65 | 0.355 | 14.31 | 12.7 | 0.121 |
| kaydırılmış yıl | 2156.26 | 2155.03 | 0.051 | 2155.22 | 2157.05 | -0.077 |
| log(EKG sayısı) | 2.33 | 1.95 | 0.409 | 2.28 | 2.2 | 0.091 |
| erkek (%) | 35.3 | 47.8 | -0.255 | 37.5 | 37.9 | -0.009 |

**Not.** MIMIC'te dönem tabakası yok (kaydırılmış tarih); cihaz filtre ayarı kesin tabakadır. İndeks T0 ablasyon öncesidir.

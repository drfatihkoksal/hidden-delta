# Aşama 3 — eşleştirme (Emory)

concept.md §5.3 "Eşleştirme uygulaması"; ANALYSIS_PLAN §2–3. Yalnızca agregat.

Birincil vaka hastası (uygun indeks T0'ı olan): 185; ilk T2'den önce T0'ı olan: 93.

Test seti: yalnızca denge raporlanır; vaka-kontrol özellik karşılaştırması dondurmadan sonra.

Önceki IAB çalışmalarında EKG'si okunmuş 43 hasta eşleştirmeden önce çıkarıldı (ANALYSIS_PLAN v1.2).

## K1 (1:5)

Havuz 19,934; eşleşen 173 / 185 (%93.5); kontrol 755; tam (1:5) küme 131.

**Denge** (SMD; hedef |SMD| < 0.1):

| Değişken | Eşleşmemiş tüm | Havuz | SMD önce | Eşleşen | Kontrol (eşleşen) | SMD sonra |
|---|---|---|---|---|---|---|
| yaş | 49.94 | 60.91 | -0.65 | 50.95 | 53.06 | -0.116 |
| KH | 76.95 | 74.12 | 0.153 | 77.5 | 76.42 | 0.058 |
| EKG sayısı | 9.24 | 5.92 | 0.382 | 7.7 | 6.76 | 0.158 |
| dönem (klasör yılı) | 2014.76 | 2014.68 | 0.036 | 2014.76 | 2014.72 | 0.018 |
| kaydırılmış yıl | 2014.74 | 2014.69 | 0.018 | 2014.72 | 2014.78 | -0.022 |
| log(EKG sayısı) | 1.84 | 1.52 | 0.417 | 1.76 | 1.66 | 0.136 |
| erkek (%) | 54.1 | 50.5 | 0.07 | 54.3 | 55.5 | -0.023 |

## K2 (1:2)

Havuz 2,992; eşleşen 127 / 185 (%68.6); kontrol 223; tam (1:2) küme 96.

**Denge** (SMD; hedef |SMD| < 0.1):

| Değişken | Eşleşmemiş tüm | Havuz | SMD önce | Eşleşen | Kontrol (eşleşen) | SMD sonra |
|---|---|---|---|---|---|---|
| yaş | 49.94 | 59.29 | -0.532 | 54.96 | 56.16 | -0.07 |
| KH | 76.95 | 82.44 | -0.259 | 77.44 | 77.57 | -0.007 |
| EKG sayısı | 9.24 | 10.53 | -0.123 | 7.48 | 7.69 | -0.037 |
| dönem (klasör yılı) | 2014.76 | 2015.28 | -0.227 | 2015.03 | 2015.01 | 0.01 |
| kaydırılmış yıl | 2014.74 | 2015.31 | -0.243 | 2015.02 | 2015.02 | 0.001 |
| log(EKG sayısı) | 1.84 | 1.97 | -0.158 | 1.77 | 1.8 | -0.039 |
| erkek (%) | 54.1 | 45.0 | 0.183 | 48.8 | 46.2 | 0.053 |

**Not.** HEEDB'de cart/filtre bilgisi yok; dönem, klasör yılından vekil olarak alındı.

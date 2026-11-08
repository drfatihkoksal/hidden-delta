# Aşama 0 — HEEDB sayımları

HEEDB 5.0. Yalnızca agregat. Preeksitasyon kodları 12SL 300/302/303/304; kaynaklar v24 yeniden işlenmiş 12SL, kayıt anı yazılım ve hekim overread. T0 adayı: hiçbir kaynakta preeksitasyon/kısa PR kodu yok, sinüs kodu var, AF/flutter/pace/dal bloğu/IVCD kodu yok, yaş ≥18. HEEDB'de prosedür kodu ve WFDB başlığında filtre bilgisi yok.

## MGH (I0001)

- EKG: 10,608,417; hasta: 1,818,247; tarih aralığı (kaydırılmış): 1961-06-06 – 2105-09-01
- 12SL v24 kodu olan EKG: 10,471,531; kayıt anı kodu olan: 10,608,415; hekim overread'i olan: 8,608,018; yaşı eksik: 2,405,753; zamanı eksik: 14,743
- T0 adayı EKG: 5,782,623 (1,181,915 hasta)

| Tanım | EKG | Hasta |
|---|---|---|
| Preeksitasyon kodu (300/302/303/304), herhangi kaynak | 20,535 | 9,984 |
|   v24 yeniden işlenmiş 12SL | 7,549 | 3,731 |
|   kayıt anı yazılım | 12,845 | 5,367 |
|   kayıt anı hekim overread | 10,235 | 5,964 |
| ALTWPW (303: fusion/intermittan preeksitasyon) | 6,490 | 4,278 |
| Yalnızca kısa PR (102), preeksitasyon kodu yok | 262,826 | 143,274 |

| Grup | Hasta |
|---|---|
| ≥1 preeksitasyon kodlu EKG'si olan hasta | 9,984 |
|   bunlardan ≥2 EKG'si olan | 9,041 |
| **Aday (katı): ≥1 preeksitasyon EKG + ≥1 T0 adayı EKG** | 4,724 |
|   adayların T0 adayı EKG sayısı | 61,215 |
|   adayların preeksitasyon EKG sayısı | 10,675 |
|   ilk preeksitasyon EKG'sinden önce T0'ı olan (prospektif varyant) | 3,281 |
|   ≥2 T0 adayı EKG'si olan | 3,930 |
|   hekim overread'inde preeksitasyon olan | 2,960 |
|   ALTWPW (303) kodlu EKG'si olan | 1,842 |
|   I45.6/426.7 kodu da olan | 312 |
| I45.6/426.7 kodlu hasta | 1,092 |
|   bunlardan EKG'si olan | 1,092 |
|   bunlardan ≥1 preeksitasyon EKG'si olan | 576 |
|   **yalnızca ICD (hiç preeksitasyon kodu yok) + ≥1 T0 adayı** | 413 |
| Uzman okuma havuzu: (katı aday ∪ yalnızca-ICD) | 5,137 |
| Yalnızca kısa PR + ≥1 T0 adayı (geniş havuz) | 74,132 |
| SVT kodlu hasta (I47.1 / 427.0) | 16,911 |
| K2 adayı: SVT kodu, preeksitasyon kodu ve I45.6/426.7 yok, ≥1 T0 adayı | 15,024 |
|   K2 adaylarının T0 adayı EKG sayısı | 187,662 |
|   kısa PR kodlu EKG'si olan | 2,299 |

Katı adayların ilk T0 EKG'sindeki medyan yaş: 60.6.

Katı adayların preeksitasyon + T0 EKG'lerinde örnekleme hızı (WFDB başlığı):

| fs (Hz) | EKG |
|---|---|
| 500 | 39,715 |
| 250 | 32,175 |

Süre: 10 s.

## Emory (I0006)

- EKG: 998,844; hasta: 349,547; tarih aralığı (kaydırılmış): 2009-01-05 – 2020-11-26
- 12SL v24 kodu olan EKG: 974,172; kayıt anı kodu olan: 998,844; hekim overread'i olan: 998,168; yaşı eksik: 3,030; zamanı eksik: 0
- T0 adayı EKG: 737,106 (301,578 hasta)

| Tanım | EKG | Hasta |
|---|---|---|
| Preeksitasyon kodu (300/302/303/304), herhangi kaynak | 1,568 | 936 |
|   v24 yeniden işlenmiş 12SL | 1,049 | 711 |
|   kayıt anı yazılım | 993 | 560 |
|   kayıt anı hekim overread | 923 | 524 |
| ALTWPW (303: fusion/intermittan preeksitasyon) | 36 | 34 |
| Yalnızca kısa PR (102), preeksitasyon kodu yok | 14,566 | 10,029 |

| Grup | Hasta |
|---|---|
| ≥1 preeksitasyon kodlu EKG'si olan hasta | 937 |
|   bunlardan ≥2 EKG'si olan | 707 |
| **Aday (katı): ≥1 preeksitasyon EKG + ≥1 T0 adayı EKG** | 334 |
|   adayların T0 adayı EKG sayısı | 22,256 |
|   adayların preeksitasyon EKG sayısı | 687 |
|   ilk preeksitasyon EKG'sinden önce T0'ı olan (prospektif varyant) | 205 |
|   ≥2 T0 adayı EKG'si olan | 208 |
|   hekim overread'inde preeksitasyon olan | 197 |
|   ALTWPW (303) kodlu EKG'si olan | 19 |
|   I45.6/426.7 kodu da olan | 110 |
| I45.6/426.7 kodlu hasta | 1,057 |
|   bunlardan EKG'si olan | 1,057 |
|   bunlardan ≥1 preeksitasyon EKG'si olan | 307 |
|   **yalnızca ICD (hiç preeksitasyon kodu yok) + ≥1 T0 adayı** | 592 |
| Uzman okuma havuzu: (katı aday ∪ yalnızca-ICD) | 926 |
| Yalnızca kısa PR + ≥1 T0 adayı (geniş havuz) | 5,974 |
| SVT kodlu hasta (I47.1 / 427.0) | 20,177 |
| K2 adayı: SVT kodu, preeksitasyon kodu ve I45.6/426.7 yok, ≥1 T0 adayı | 16,766 |
|   K2 adaylarının T0 adayı EKG sayısı | 74,706 |
|   kısa PR kodlu EKG'si olan | 1,241 |

Katı adayların ilk T0 EKG'sindeki medyan yaş: 55.0.

Katı adayların preeksitasyon + T0 EKG'lerinde örnekleme hızı (WFDB başlığı):

| fs (Hz) | EKG |
|---|---|
| 500 | 1,746 |
| 250 | 124 |

Süre: 2 s.

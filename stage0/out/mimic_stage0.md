# Aşama 0 — MIMIC sayımları

MIMIC-IV-ECG 1.0 + MIMIC-IV 3.1 hosp + MIMIC-IV-Note. Yalnızca agregat. Çalışma süresi 2 s.

- EKG: 800,035; hasta: 161,352; EKG tarih aralığı (kaydırılmış): 2087-05-08 – 2224-04-14
- Yaşı hesaplanamayan EKG: 54
- T0 adayı EKG (preeksitasyon/kısa PR ifadesi yok, sinüs, AF/flutter/pace yok, QRS <120 ms ve dal bloğu/IVCD yok, yaş ≥18): 555,605 EKG, 141,638 hasta

## 1. Rapor regex'i

| Tanım | EKG | Hasta |
|---|---|---|
| Concept regex (preeksitasyon ifadesi VEYA kısa PR) | 11,984 | 7,446 |
| Preeksitasyon ifadesi (katı: WPW/Wolff/pre-excit/delta) | 713 | 506 |
| Yalnızca kısa PR (preeksitasyon ifadesi yok) | 11,271 | 7,077 |

Yalnızca kısa PR ifadesi olup hiç preeksitasyon ifadesi olmayan hasta: 6,940. MIMIC makine raporlarında preeksitasyonla ilgili seri karşılaştırma ifadesi ("no longer evident" vb.) yok.

## 2. İntermittan aday havuzları (makine raporu düzeyinde; uzman doğrulaması öncesi)

| Grup | Hasta |
|---|---|
| ≥1 preeksitasyon ifadeli EKG olan hasta | 506 |
|   bunlardan ≥2 EKG'si olan | 473 |
| **Aday (katı): ≥1 preeksitasyon EKG + ≥1 T0 adayı EKG** | 352 |
|   adayların T0 adayı EKG sayısı | 3,562 |
|   adayların preeksitasyon EKG sayısı | 487 |
|   ilk preeksitasyon EKG'sinden önce T0'ı olan (prospektif varyant) | 293 |
|   ≥2 T0 adayı EKG'si olan | 298 |
|   I45.6/426.7 kodu da olan | 22 |
| I45.6/426.7 kodlu hasta | 178 |
|   bunlardan EKG'si olan | 139 |
|   bunlardan ≥1 preeksitasyon EKG'si olan | 45 |
|   **yalnızca ICD (hiç preeksitasyon ifadesi yok) + ≥1 T0 adayı** | 80 |
| Uzman okuma havuzu: (katı aday ∪ yalnızca-ICD) | 432 |
| Yalnızca kısa PR + ≥1 T0 adayı (geniş havuz, çoğu WPW değil) | 5,345 |

Katı adayların ilk T0 EKG'sindeki yaş: medyan 61.0, aralık 18–92.

## 3. Ablasyon (H3)

Ablasyon kodu (02583ZZ / 37.34; yolağa özgü değil) olan hasta: 1,695; EKG'si olan: 1,478.

| Grup | Hasta |
|---|---|
| **Yolak ablasyonu vekili** (ablasyon + [I45.6/426.7 ya da ablasyon öncesi preeksitasyon EKG'si]) | 39 |
|   ablasyon öncesi preeksitasyon EKG'si olan | 18 |
|   yalnızca ICD ile tanımlanan | 21 |
|   ablasyon öncesi ve sonrası EKG'si olan | 23 |
|   öncesi preeksitasyon EKG + sonrası T0 adayı | 9 |
|     sonrası T0 0–90 gün (bellek olası) | 6 |
|     sonrası T0 >90 gün | 5 |
|   ablasyon sonrası preeksitasyon EKG'si olan (başarısız/nüks?) | 2 |
|   aynı zamanda katı intermittan aday olan | 14 |

## 4. Filtre / cart envanteri

Farklı cart_id: tümünde 156, katı adayların preeksitasyon + T0 EKG'lerinde 128.

Tüm EKG'ler:

| filtering | bandwidth | EKG |
|---|---|---|
| 60 Hz notch Baseline filter | 0.005-150 Hz | 622,477 |
| <not specified> | 0.0005-150 Hz | 107,398 |
| 60 Hz notch Baseline filter | 0.05-150 Hz | 51,006 |
| Baseline filter | 0.005-150 Hz | 12,505 |
| 60 Hz notch Baseline filter | 0.0005-150 Hz | 3,932 |
| 50 Hz notch Baseline filter | 0.005-150 Hz | 2,215 |
| 50 Hz notch Baseline filter | 0.05-150 Hz | 502 |

Katı adayların preeksitasyon + T0 EKG'leri:

| filtering | bandwidth | EKG |
|---|---|---|
| 60 Hz notch Baseline filter | 0.005-150 Hz | 3,241 |
| <not specified> | 0.0005-150 Hz | 411 |
| 60 Hz notch Baseline filter | 0.05-150 Hz | 289 |
| Baseline filter | 0.005-150 Hz | 67 |
| 60 Hz notch Baseline filter | 0.0005-150 Hz | 20 |
| 50 Hz notch Baseline filter | 0.005-150 Hz | 15 |
| 50 Hz notch Baseline filter | 0.05-150 Hz | 6 |

## 5. K2 — SVT kontrolü

SVT kodlu hasta (I47.1x, I47.11 hariç / 427.0): 4,040.

| Grup | Hasta |
|---|---|
| K2 adayı: SVT kodu, preeksitasyon ifadesi ve I45.6/426.7 yok, ≥1 T0 adayı | 2,392 |
|   T0 adayı EKG sayısı | 22,833 |
|   kısa PR ifadeli EKG'si olan | 268 |
|   ablasyon kodu olan | 214 |
| Taburcu özeti olan K2 hastası | 2,198 |
|   yalnızca AVNRT kalıbı | 178 |
|   yalnızca aksesuar yol / AVRT kalıbı | 14 |
|   ikisi birden | 34 |
|   hiçbiri | 1,972 |
|   notta WPW/preeksitasyon sözü geçen (kontaminasyon adayı) | 10 |
|   notta EPS/ablasyon sözü geçen | 434 |
|   **K2a adayı: yalnızca AVNRT, WPW sözü yok** | 178 |
|   **K2b adayı: yalnızca aksesuar yol/AVRT, WPW sözü yok** | 12 |

Not kalıpları olumsuzlamayı ayıklamaz ("no evidence of AVNRT" da sayılır); alt tip ancak okuma ile kesinleşir.

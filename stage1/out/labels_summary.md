# Aşama 1 — etiketler, kesin G1 ve kaynaklar arası uyum

concept.md §5.2 tanımları; manuel okuma yok. Yalnızca agregat.

## G1 özeti

| Tanım | Vaka hastası | T2'den önce T0'ı olan | ≥1 sıkı T0'ı olan |
|---|---|---|---|
| MIMIC (dış test) — kesin T2 | 117 | 89 | 88 |
| MIMIC — herhangi T2 (duyarlılık) | 352 | 293 | 265 |
| MGH — ≥2 kaynak uzlaşısı (birincil), 500 Hz T0 | 1,015 | 526 | — |
| MGH — hekim overread (duyarlılık), 500 Hz T0 | 1,992 | 1,073 | — |
| MGH — v24 12SL (duyarlılık), 500 Hz T0 | 1,260 | 843 | — |
| Emory — ≥2 kaynak uzlaşısı (birincil), 500 Hz T0 | 186 | 96 | — |
| Emory — hekim overread (duyarlılık), 500 Hz T0 | 193 | 92 | — |
| Emory — v24 12SL (duyarlılık), 500 Hz T0 | 240 | 148 | — |

G1 eşikleri: geliştirme ≥100 (MGH, uzlaşı: 1,015), dış test ≥40 (MIMIC kesin: 117).

## MIMIC ayrıntı

| Grup | Hasta |
|---|---|
| Kesin T2 vaka hastası | 117 |
|   T0 EKG sayısı | 1,085 |
|   ilk T2'den önce T0'ı olan (prospektif) | 89 |
|   ≥1 sıkı T0 (PR ≥140 ms, QRS <110 ms) | 88 |
|   ilk T2'den önce sıkı T0'ı olan | 63 |
| Herhangi T2 vaka hastası (duyarlılık) | 352 |
|   ≥1 sıkı T0 | 265 |
| H3: kesin T2 ablasyon öncesi + ablasyon sonrası T0 | 5 |
|   sonrası T0 0–90 gün | 5 |
|   sonrası T0 >90 gün | 1 |
| H3, herhangi T2 ile | 9 |
| K2 hastası (≥1 T0) | 2,392 |
|   ≥1 sıkı T0 | 2,083 |
|   K2a (yalnızca AVNRT kalıbı) | 178 |
|     ≥1 sıkı T0 | 151 |
|   K2b (yalnızca aksesuar yol kalıbı) | 12 |

T0 EKG'lerinde makine PR ölçümü olan: 553,958 / 555,605.

Kesin T2 vakalarının T2 + T0 EKG'lerinde filtre ayarı:

| filtering | bandwidth | EKG |
|---|---|---|
| 60 Hz notch Baseline filter | 0.005-150 Hz | 992 |
| <not specified> | 0.0005-150 Hz | 146 |
| 60 Hz notch Baseline filter | 0.05-150 Hz | 88 |
| Baseline filter | 0.005-150 Hz | 20 |
| 60 Hz notch Baseline filter | 0.0005-150 Hz | 6 |
| 50 Hz notch Baseline filter | 0.005-150 Hz | 4 |

## MGH (I0001)

| Grup | Hasta |
|---|---|
| **T2 = ≥2 kaynak uzlaşısı (birincil) + ≥1 500 Hz T0** | 1,015 |
|   (fs kısıtı öncesi) | 1,401 |
|   500 Hz T0 EKG sayısı | 7,734 |
|   ilk T2'den önce 500 Hz T0'ı olan (prospektif) | 526 |
|   ≥1 500 Hz T2'si de olan (H2 delta vektörü için) | 707 |
| T2 = hekim overread (duyarlılık) + ≥1 500 Hz T0 | 1,992 |
|   (fs kısıtı öncesi) | 2,792 |
|   500 Hz T0 EKG sayısı | 15,199 |
|   ilk T2'den önce 500 Hz T0'ı olan (prospektif) | 1,073 |
|   ≥1 500 Hz T2'si de olan (H2 delta vektörü için) | 1,134 |
| T2 = v24 12SL (duyarlılık) + ≥1 500 Hz T0 | 1,260 |
|   (fs kısıtı öncesi) | 1,771 |
|   500 Hz T0 EKG sayısı | 9,955 |
|   ilk T2'den önce 500 Hz T0'ı olan (prospektif) | 843 |
|   ≥1 500 Hz T2'si de olan (H2 delta vektörü için) | 903 |
| Uzlaşı ve overread vakası birlikte | 798 |
| Yalnızca overread'e göre vaka (uzlaşı değil) | 1,194 |
| Yalnızca v24'e göre vaka (uzlaşı değil) | 535 |
| K2 hastası (≥1 overread'li T0; fs kontrol edilmedi) | 14,246 |

Başlığı okunamayan / beklenmeyen fs: 0 EKG.

**Kaynaklar arası uyum** (üç kaynağın da olduğu 8,471,177 EKG, preeksitasyon kodu):

| Çift | Cohen κ | Pozitif uyum |
|---|---|---|
| v24–yazılım | 0.396 | 0.396 |
| v24–overread | 0.306 | 0.307 |
| yazılım–overread | 0.547 | 0.547 |

Fleiss κ (üç kaynak): 0.428. Desen sayıları (v24, yazılım, overread): {'000': 8453844, '001': 4359, '010': 3953, '011': 3286, '100': 2434, '101': 95, '110': 878, '111': 2328}.

## Emory (I0006)

| Grup | Hasta |
|---|---|
| **T2 = ≥2 kaynak uzlaşısı (birincil) + ≥1 500 Hz T0** | 186 |
|   (fs kısıtı öncesi) | 190 |
|   500 Hz T0 EKG sayısı | 661 |
|   ilk T2'den önce 500 Hz T0'ı olan (prospektif) | 96 |
|   ≥1 500 Hz T2'si de olan (H2 delta vektörü için) | 179 |
| T2 = hekim overread (duyarlılık) + ≥1 500 Hz T0 | 193 |
|   (fs kısıtı öncesi) | 196 |
|   500 Hz T0 EKG sayısı | 695 |
|   ilk T2'den önce 500 Hz T0'ı olan (prospektif) | 92 |
|   ≥1 500 Hz T2'si de olan (H2 delta vektörü için) | 185 |
| T2 = v24 12SL (duyarlılık) + ≥1 500 Hz T0 | 240 |
|   (fs kısıtı öncesi) | 244 |
|   500 Hz T0 EKG sayısı | 852 |
|   ilk T2'den önce 500 Hz T0'ı olan (prospektif) | 148 |
|   ≥1 500 Hz T2'si de olan (H2 delta vektörü için) | 227 |
| Uzlaşı ve overread vakası birlikte | 161 |
| Yalnızca overread'e göre vaka (uzlaşı değil) | 32 |
| Yalnızca v24'e göre vaka (uzlaşı değil) | 86 |
| K2 hastası (≥1 overread'li T0; fs kontrol edilmedi) | 16,766 |

Başlığı okunamayan / beklenmeyen fs: 0 EKG.

**Kaynaklar arası uyum** (üç kaynağın da olduğu 973,517 EKG, preeksitasyon kodu):

| Çift | Cohen κ | Pozitif uyum |
|---|---|---|
| v24–yazılım | 0.625 | 0.625 |
| v24–overread | 0.550 | 0.550 |
| yazılım–overread | 0.786 | 0.786 |

Fleiss κ (üç kaynak): 0.651. Desen sayıları (v24, yazılım, overread): {'000': 971981, '001': 152, '010': 122, '011': 216, '100': 400, '101': 17, '110': 109, '111': 520}.

# Keşifsel: vaka etiketinin gücü, sıkı T0 ve AUROC (kayıt sonrası)

| Site | Vaka | 3/3 kaynak | ≥2 konsensüs T2 |
|---|---|---|---|
| MGH | 852 | 39.8% | 53.5% |
| Emory | 173 | 59.5% | 39.9% |

MGH OOF AUROC: birleştirilmiş 0.665, kat içi ortalama 0.666.

| MGH tabakası | Vaka | AUROC (%95 GA, küme bootstrap) |
|---|---|---|
| 3/3 kaynak: evet | 339 | 0.709 (0.675–0.742) |
| 3/3 kaynak: hayır | 513 | 0.636 (0.607–0.663) |
| ≥2 konsensüs T2: evet | 456 | 0.670 (0.642–0.698) |
| ≥2 konsensüs T2: hayır | 396 | 0.659 (0.625–0.689) |

### MGH — M1 OOF (geliştirme)

| Tabaka | Vaka | Kontrol | AUROC (%95 GA, küme bootstrap) |
|---|---|---|---|
| Tüm vakalar | 852 | 3198 | 0.665 (0.643–0.686) |
| Tüm vakalar, sıkı T0 | 144 | 248 | 0.605 (0.540–0.661) |
| 3/3 kaynak | 339 | 1280 | 0.709 (0.675–0.742) |
| 3/3 kaynak, sıkı T0 | 48 | 86 | 0.654 (0.545–0.761) |
| 2/3 kaynak | 513 | 1918 | 0.636 (0.607–0.663) |
| 2/3 kaynak, sıkı T0 | 96 | 162 | 0.581 (0.504–0.656) |

### Emory — dondurulmuş M1 (test seti; kayıt sonrası keşifsel)

| Tabaka | Vaka | Kontrol | AUROC (%95 GA, küme bootstrap) |
|---|---|---|---|
| Tüm vakalar | 173 | 755 | 0.736 (0.693–0.779) |
| Tüm vakalar, sıkı T0 | 28 | 52 | 0.720 (0.609–0.823) |
| 3/3 kaynak | 103 | 445 | 0.739 (0.686–0.791) |
| 3/3 kaynak, sıkı T0 | 17 | 32 | 0.711 (0.565–0.849) |
| 2/3 kaynak | 70 | 310 | 0.734 (0.662–0.798) |
| 2/3 kaynak, sıkı T0 | 11 | 20 | 0.736 (0.582–0.909) |

Sıkı T0 oranı (vaka / kontrol): MGH 24.3% / 31.4%; Emory 23.7% / 32.6%.
MGH vakalarında sıkı T0 oranı: 3/3 19.5%, 2/3 27.5%.
Emory vakalarında sıkı T0 oranı: 3/3 26.2%, 2/3 20.0%.

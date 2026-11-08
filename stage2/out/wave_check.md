# QRS öncesi ayrık dalga: ekran filtresi testi

Gözlenen dalga: ~20 ms, ~10 µV, QRS'ten ayrı, Q'dan hemen önce (IAB etiketleme çalışması; medyan vurular 0.5–40 Hz 2. derece Butterworth ile sıfır fazlı filtrelenmiş sinyalden).

## 1. Filtrenin keskin kenara yanıtı

| girdi                       |   QRS öncesi en büyük sapma (µV) |   zamanı (ms) |   yarı yükseklik genişliği (ms) | ilk QRS sapmasına göre   |   −4 ms'deki değer (µV) |
|:----------------------------|---------------------------------:|--------------:|--------------------------------:|:-------------------------|------------------------:|
| R 1 mV (üçgen, 40 ms)       |                            -18.4 |            -8 |                               8 | ters                     |                    -1.7 |
| qR: q −0.1 mV + R 1 mV      |                            -18.9 |            -4 |                               4 | ters                     |                   -18.9 |
| rS: r 0.2 mV + S −1 mV      |                              6.9 |            -4 |                              10 | aynı                     |                     6.9 |
| geniş R 1 mV (üçgen, 80 ms) |                             -9.4 |            -8 |                              10 | ters                     |                    -1   |

## 2. Aynı kayıtta filtresiz ve IAB filtreli medyan vuru: filtrenin QRS öncesine eklediği sapma

Derivasyon başına (8 bağımsız derivasyon), düşük eşikli QRS başlangıcından önceki −30…−4 ms'de, başlangıçtan 80–50 ms önceki düzeye göre en büyük sapma. P kuyruğu iki zincirde aynı olduğundan fark filtreye aittir. Baskın QRS: başlangıç–bitiş arasında en büyük mutlak değerli örnek. Kontroller preeksitasyonu olmayan hastalardır; etiket kullanılmaz.

| veri                     |   derivasyon |   |sapma| medyanı, filtresiz (µV) |   |sapma| medyanı, IAB filtreli (µV) |   ≥10 µV, filtresiz (%) |   ≥10 µV, IAB filtreli (%) |   filtrenin eklediği, baskın QRS'e ters (%) |   aynısı, baskın QRS ≥0.2 mV olan derivasyonlarda (%) |   eklenen sapma / baskın QRS (medyan, %) |   Spearman r (|eklenen|, |baskın QRS|) |
|:-------------------------|-------------:|----------------------------------:|-------------------------------------:|------------------------:|---------------------------:|--------------------------------------------:|------------------------------------------------------:|-----------------------------------------:|---------------------------------------:|
| sentetik normal EKG (20) |          160 |                               4.6 |                                  4.5 |                       0 |                        0   |                                        68.8 |                                                  66.7 |                                     0.15 |                                 -0.089 |
| MGH K1 kontrol (1000)    |         8000 |                              25   |                                 25   |                      82 |                       79.9 |                                        49.9 |                                                  49.9 |                                     0.5  |                                  0.11  |

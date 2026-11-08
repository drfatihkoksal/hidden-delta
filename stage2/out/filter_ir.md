# Pace spike'larından cihaz dürtü yanıtı (yalnızca spike öncesi)

Sıkı spike ölçütleri: ≥0.3 mV; spike'tan 24–8 ms önce düz taban çizgisi (SD ≤ genliğin %5'i); ±6 ms'de genliğin yarısının altına iniş; ≥3 derivasyonda görünürlük. h: spike tepesine (t = 0, ms) göre normalize medyan yanıt (önce kayıt içi, sonra kayıtlar arası medyan). Spike sonrası kısım ardından gelen P/QRS ve polarizasyon artığıyla kirlenebilir; pre-ringing yalnızca t < 0'dan okunur.

| veri   | grup                                        |   kayıt |   spike'lı kayıt |   spike |   genişlik (örnek) |   h(−8) |   h(−6) |   h(−4) |   h(−2) |   h(+2) |   max|h| −20…−10 |   kayıtlar arası SD (−8…−2) |
|:-------|:--------------------------------------------|--------:|-----------------:|--------:|-------------------:|--------:|--------:|--------:|--------:|--------:|-----------------:|----------------------------:|
| Emory  | Emory                                       |    1200 |              752 |    8765 |               2.23 |   0     |  -0.003 |   0     |   0.041 |   0.027 |            0     |                       0.183 |
| MGH    | MGH                                         |    1200 |              428 |    4763 |               1.82 |   0     |  -0.001 |  -0.002 |   0.008 |   0.002 |            0     |                       0.147 |
| MIMIC  | 50 Hz notch Baseline filter | 0.005-150 Hz  |     109 |               57 |     488 |               2.8  |   0.037 |   0.035 |  -0.145 |   0.313 |   0.135 |            0     |                       0.136 |
| MIMIC  | 50 Hz notch Baseline filter | 0.05-150 Hz   |      14 |                7 |      39 |               3.33 |   0.039 |   0.042 |  -0.109 |   0.565 |   0.332 |            0.009 |                       0.176 |
| MIMIC  | 60 Hz notch Baseline filter | 0.0005-150 Hz |     198 |               95 |     763 |               3.16 |   0.033 |   0.024 |  -0.153 |   0.291 |  -0.027 |            0.001 |                       0.179 |
| MIMIC  | 60 Hz notch Baseline filter | 0.005-150 Hz  |     300 |              176 |    1526 |               2.83 |   0.033 |   0.03  |  -0.137 |   0.271 |   0.136 |            0     |                       0.146 |
| MIMIC  | 60 Hz notch Baseline filter | 0.05-150 Hz   |     300 |              179 |    1672 |               2.69 |   0.033 |   0.032 |  -0.152 |   0.299 |   0.03  |            0     |                       0.131 |
| MIMIC  | <not specified> | 0.0005-150 Hz             |     300 |              190 |    1871 |               2.47 |   0     |   0     |  -0.065 |   0.137 |   0     |            0     |                       0.111 |
| MIMIC  | Baseline filter | 0.005-150 Hz              |     300 |              157 |    1397 |               2.81 |   0.032 |   0.031 |  -0.148 |   0.296 |   0.087 |            0     |                       0.134 |

**Yorum.** MIMIC'te t = −4 ms'de spike'ın tersi işaretli ~%15'lik bir lob ve −6…−8 ms'de ~%3'lük bir lob var: doğrusal fazlı bir zincirin pre-ringing'i. Lob filtre ayarından bağımsız; kaynağı örnekleme/anti-aliasing zinciri. "<not specified> | 0.0005-150 Hz" grubunda lob daha küçük (farklı cart). HEEDB'de spike'lar tek örnek genişliğinde ve öncesi tam sıfır: büyük olasılıkla cihazın eklediği yapay pace işaretleri. HEEDB'nin dürtü yanıtı bu yolla tahmin edilemez.

Eğriler: `stage2/out/filter_ir_curves.csv`.

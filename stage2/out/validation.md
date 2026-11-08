# Aşama 2 — pipeline doğrulaması

## 1. Sentetik testler

| vaka   |   slow_init_ms |   qrs_ms |   pr_ms |   q_V6 |   pqhigh_mag |   pqhigh_fq_corr |   pqhigh_opp_frac |   pqlow_fq_corr |   n_beats |   n_clusters |
|:-------|---------------:|---------:|--------:|-------:|-------------:|-----------------:|------------------:|----------------:|----------:|-------------:|
| normal |              6 |       80 |     192 |  0.113 |        0.004 |           -0.375 |               nan |          -0.039 |        11 |            1 |
| delta  |             42 |      128 |     142 |  0     |        0.207 |            0.948 |                 0 |           0.834 |        11 |            1 |
| gibbs  |             34 |      110 |     178 |  0.12  |        0.044 |           -0.996 |                 1 |           0.709 |        11 |            1 |

Beklenen: delta'da `slow_init_ms` uzun ve `pqhigh_fq_corr` > 0; gibbs'te `pqhigh_fq_corr` < 0.

## 2. MIMIC makine ölçümleriyle uyum (rastgele 1000 T0 EKG, preeksitasyon/kısa PR hastaları hariç)

| ölçüm                          |   n |   fark medyanı (ms) |   mutlak fark medyanı (ms) |   |fark| ≤20 ms (%) |   Spearman r |   bizim medyan |   makine medyanı |
|:-------------------------------|----:|--------------------:|---------------------------:|--------------------:|-------------:|---------------:|-----------------:|
| QRS süresi (düşük eşik)        | 990 |                  -2 |                          8 |                83.4 |        0.528 |             88 |               90 |
| QRS süresi (görünür başlangıç) | 990 |                 -14 |                         14 |                80.4 |        0.678 |             76 |               90 |
| PR aralığı                     | 952 |                  44 |                         46 |                19.9 |        0.476 |            200 |              162 |

**Kalite:**

```
{
 "n": 1000,
 "başarılı": 990,
 "hatalar": {
  "SignalError": 10
 },
 "baskın küme oranı medyanı": 1.0,
 "küme sayısı dağılımı": {
  "1.0": 868,
  "2.0": 74,
  "3.0": 17,
  "4.0": 9,
  "5.0": 5,
  "6.0": 5,
  "7.0": 1,
  "8.0": 5,
  "9.0": 1,
  "10.0": 1,
  "11.0": 1,
  "13.0": 1,
  "15.0": 1,
  "16.0": 1
 },
 "başlangıç sınıra dayanan (%)": 6.57,
 "PR-sonu taban penceresi eğimi medyanı (mV)": 0.0171,
 "yavaş başlangıç medyanı (ms)": 10.0
}
```

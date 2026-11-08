# HIDDEN-DELTA — kayıt sonrası betimsel ekler

Resmî sonuçlardan sonra yazıldı (`stage5/s5_supplement.py`). Planın betimsel olarak istediği ama analiz kodunda eksik kalan çıktılar. Test ve karar kuralı yok. Yeni skorlama yok. Yalnızca agregat.

## MIMIC

- **5.2 P dalgası atıf payı:** vaka 0.099, kontrol 0.119; küme farkı -0.0163 (-0.0263 – -0.0062).
- **Ayrık dalga, ayrıntı** (dalgalı EKG: 67 vaka, 301 kontrol; yeniden hesap resmî bayrakla aynı):
  - İlk QRS sapmasıyla aynı polaritedeki derivasyon oranı: vaka 0.407, kontrol 0.464. Çoğunlukla ters polarite: vaka 0.493, kontrol 0.458.
  - Medyan genlik 17.5 / 15.0 µV; genişlik 9.0 / 9.0 ms; tepe -32.0 / -32.5 ms (vaka / kontrol; QRS başlangıcına göre).
  - Kendi delta yönüyle uyum (dalgalı vakalar, n 67): ortalama r -0.157, permütasyon -0.027; tek yönlü p 0.9935.
- **Dalga sıklığı, filtre grubuna göre:**
  - 60 Hz notch Baseline filter | 0.005-150 Hz: 85 vaka, 412 kontrol; genel 0.694, vaka 0.706, kontrol 0.692
  - 60 Hz notch Baseline filter | 0.05-150 Hz: 6 vaka, 15 kontrol; genel 0.714, vaka 0.833, kontrol 0.667
  - <not specified> | 0.0005-150 Hz: 12 vaka, 34 kontrol; genel 0.174, vaka 0.167, kontrol 0.176
- **7.2 ikinci koşul:** AUROC farkı (<not specified> | 0.0005-150 Hz − 60 Hz notch Baseline filter | 0.005-150 Hz) -0.085 (-0.266 – +0.071).
- **7.4 K2b (gizli yol):** notlarda 12 hasta, uygun T0'ı olan 8, eşleşmiş kümelerde 0. Betimsel skor raporlanamaz.

## Emory

- **5.2 P dalgası atıf payı:** vaka 0.096, kontrol 0.106; küme farkı -0.0085 (-0.0148 – -0.0025).
- **Ayrık dalga, ayrıntı** (dalgalı EKG: 93 vaka, 329 kontrol; yeniden hesap resmî bayrakla aynı):
  - İlk QRS sapmasıyla aynı polaritedeki derivasyon oranı: vaka 0.426, kontrol 0.417. Çoğunlukla ters polarite: vaka 0.538, kontrol 0.526.
  - Medyan genlik 14.8 / 13.4 µV; genişlik 12.0 / 10.0 ms; tepe -30.5 / -34.0 ms (vaka / kontrol; QRS başlangıcına göre).
  - Kendi delta yönüyle uyum (dalgalı vakalar, n 90): ortalama r -0.047, permütasyon 0.004; tek yönlü p 0.8579.


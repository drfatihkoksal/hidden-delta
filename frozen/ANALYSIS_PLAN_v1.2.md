# HIDDEN-DELTA — Analiz planı (Aşama 4–5)

**Sürüm:** v1.2 — 2026-11-08
**Durum:** Sonuçlar görülmeden önce yazıldı ve kullanıcı onayıyla donduruldu (2026-11-08; dosya özeti `FROZEN.sha256`). Sonraki değişiklikler gerekçe ve tarihle §13'e eklenir; dondurulmuş sürümün yerine geçmez.
**Bağlı belgeler:** `concept.md` (v0.2; §3 hipotezler, §5 etiketler ve gruplar, §8 denetimler), `stage1–3/README.md`.

> **Bu plandan önce görülmüş veriler:** MGH'de G2 pilotu ve eşleşme sonrası betimsel farklar (`stage2/out/g2_pilot.md`, `stage3/out/match_I0001.md`). Bu nedenle MGH, hipotez testi için değil, model geliştirme ve keşif içindir. Hipotezler **Emory** (iç-dış) ve **MIMIC** (dış) üzerinde test edilir; ikisinde de etiketli analiz yapılmadı.

---

## 1. Hipotezler ve birincil sonlanım

| | Hipotez | Test edildiği veri |
|---|---|---|
| **H1 (birincil)** | Preeksitasyonsuz (T0) EKG'den üretilen skor, vakaları eşleştirilmiş K1 kontrollerinden ayırır | MIMIC (birincil), Emory (ikincil) |
| H2a | Skoru sürükleyen sinyal QRS başlangıcındadır ve hastanın kendi delta vektörüyle uyumludur | MIMIC + Emory |
| H2b | Sinyalin bir kısmı repolarizasyondadır (kardiyak bellek); T vektörü hastanın kendi preeksite QRS vektörüyle uyumludur | MIMIC + Emory |
| H3 | Ablasyondan sonra skor düşer | MIMIC; ~5–9 hasta, yalnızca tanımlayıcı olgu serisi |
| H4 | AVNRT'de skor yükselmez | MIMIC (K2a, taburcu özetinden) |
| Keşif | Kullanıcının gözlediği QRS öncesi ayrık dalga vakaya özgüdür (§6) | MGH (keşif) → Emory/MIMIC (doğrulama) |

**Birincil sonlanım:** MIMIC dış testinde, birincil modelin (M1; §4) hasta düzeyinde AUROC'u (%95 GA, hasta düzeyinde bootstrap).
- H1 desteklendi sayılır: AUROC'un %95 GA alt sınırı > 0.60 ve demografik kısa yol modelinden (§7.6) anlamlı biçimde yüksek (eşleştirilmiş bootstrap farkı, alt sınır > 0).
- %95 özgüllükte duyarlılık, AUPRC ve kalibrasyon ikincil sonlanımlardır.
- Birincil sonlanımın AUROC olması kullanıcı tarafından onaylandı (2026-11-08; concept §13 soru 6).

**Güç (yaklaşık):** MIMIC'te 117 vaka ve 1:5 kontrolle, AUROC 0.70 için SE ≈ 0.029 (Hanley–McNeil). %95 GA yarı genişliği ≈ ±0.056. Alt sınır 0.60'ı geçmek için gerçek AUROC'un ≈ 0.66'nın üstünde olması gerekir.

---

## 2. Veri, bölmeler ve dondurma

| Veri | Rol | Erişim kuralı |
|---|---|---|
| MGH (HEEDB I0001) | Geliştirme: hasta düzeyinde 5 katlı çapraz doğrulama ile model seçimi; final model tüm MGH ile eğitilir | Serbest |
| Emory (HEEDB I0006) | İç-dış test | Dondurmaya kadar yalnızca etiketsiz sinyal çıkarımı ve eşleştirme |
| MIMIC-IV-ECG | Dış test (birincil) | Dondurmaya kadar yalnızca etiketsiz sinyal çıkarımı ve eşleştirme; skor hesaplanmaz |

- **Dondurma:** Final modeller, eşik (%95 özgüllük, MGH çapraz doğrulama tahminlerinden), ön işleme kodu ve bu plan özetlenip `FROZEN.sha256`'ya yazılır. Ondan sonra Emory ve MIMIC skorlanır. Dondurmadan sonra model değişirse yeni sonuçlar "dondurma sonrası, keşifsel" diye ayrı raporlanır.
- **Ön eğitim:** Self-supervised ön eğitim yalnızca MGH eğitim katlarının etiketsiz EKG'leriyle yapılır. HEEDB'de etiketle eğitilmiş hazır modeller (ör. ECGFounder) kullanılmaz (concept §7, sızıntı).
- **Eşleştirme:** Emory ve MIMIC'te Aşama 3'ün kuralları aynen uygulanır. MIMIC'te ek kesin tabaka olarak cihaz filtre ayarı (`filtering × bandwidth`) kullanılır. MIMIC'te dönem yerine kaydırılmış yıl yoktur; dönem tabakası atlanır.

---

## 3. Etiketler, indeks EKG ve birleştirme

- **Vaka:** concept §5.2 birincil tanımı. HEEDB'de üç kaynaktan ≥2'si; MIMIC'te kesin ifade. ≥1 T0 + ≥1 T2.
- **Kontrol:** K1 (1:5, eşleştirilmiş). H4 için ayrıca K2.
- **Birincil birim — tek indeks EKG:** Vakada rastgele bir T0, kontrolde tek T0 (Aşama 3 ile aynı). Bu, iki grupta EKG sayısından kaynaklanan "fırsat" farkını ortadan kaldırır.
- **İkincil birim — hasta düzeyinde ortalama ve maksimum:** Vakada tüm T0'lar; kontrolde, vakanın T0 sayısına eşit sayıda rastgele T0. Bunun için Aşama 4'ün ilk adımında kontrollerden ek T0 çıkarımı yapılır.
- **Hariç tutulanlar:** Pipeline başarısız (%0.5), `onset_capped`, yaşı bilinmeyen. Sayılar akış diyagramında raporlanır.
- **IAB okumaları (v1.2):** Yazarın önceki IAB çalışmalarında (AUTO-IAB, LEAD-IAB, P-SPLIT Emory işaret seti) herhangi bir EKG'sini okuduğu hastalar Emory ve MIMIC test setlerinden eşleştirmeden önce çıkarılır (`stage3/s3_iab_overlap.py`). Gerekçe: Gözlenen dalga bu okumalarda görüldü.

---

## 4. Modeller

Tümü MGH'de hasta düzeyinde 5 katlı çapraz doğrulama ile seçilir. Hiperparametre arama alanı aşağıdadır ve genişletilmez.

| Model | Girdi | Ayrıntı |
|---|---|---|
| **M1 (birincil)** | Baskın küme medyan vurusu, 12 × 500 örnek (500 Hz, −400…+600 ms; Aşama 2) | 1D ResNet (4–8 blok). MGH eğitim katlarında SSL ön eğitimi (maskelenmiş yeniden yapılandırma) → ince ayar ya da linear probe; ikisi arasında seçim çapraz doğrulama AUROC'una göre. Filtre augmentasyonu (concept §7). Platt kalibrasyonu |
| M0 | §4.1'deki önceden belirlenmiş özellikler + yaş + cinsiyet | Doğal kübik spline'lı (3 düğüm) L2 lojistik regresyon |
| M2 | M1 gömmesi + M0 özellikleri | Lojistik regresyon |
| M-demo | Yaş, cinsiyet, KH | Demografik kısa yol kontrolü (§7.6) |

- **Neden M1 medyan vuru üzerinde, ham 10 s üzerinde değil:** Mekanizma denetimleri (oklüzyon, atıf) medyan vuru üzerinde tanımlıdır. Ham 10 s'lik bir M1 duyarlılık analizidir.
- **M1 ile M0 arasındaki birincil karşılaştırma tektir:** MIMIC'te AUROC farkı, eşleştirilmiş bootstrap ile.

### 4.1 M0 özellikleri (sabit liste)

- **Başlangıç:**
  - lateral derivasyonlarda septal q: I, aVL, V5 ve V6'da genlik, ayrıca dördünün hiçbirinde q olmaması;
  - V1'de R/S ve rSr';
  - QRS süresi;
  - `slow_init_ms`, `svr10`;
  - başlangıç vektörleri (`v20_*`, `v40_*`, `i20_*`).
- **QRS öncesi:** `pqlow_mag`, `pqlow_fq_corr` ve §6'daki ayrık dalga bayrağı.
- **Repolarizasyon:** QRS-T açısı, T integral vektörü.
- **Diğer:** PR (Aşama 2 ölçümü; +42 ms sapmalı, bütün gruplarda aynı), KH.

---

## 5. Mekanizma analizleri (H2a, H2b)

Bütün testler MIMIC ve Emory'de, dondurmadan sonra, birincil M1 ile yapılır. MGH'deki karşılıkları yalnızca betimseldir.

| Analiz | Ölçü | Gerçek sinyal lehine ölçüt |
|---|---|---|
| 5.1 Oklüzyon | AUROC düşüşü: QRS başlangıcı ±40 ms maskelenince; ST-T (QRS bitişi +40…+400 ms) maskelenince; aynı genişlikte rastgele bir pencere maskelenince (kontrol) | H2a: başlangıç maskesi AUROC'u, rastgele pencereden anlamlı biçimde daha çok düşürür. H2b: ST-T maskesi için aynı |
| 5.2 Atıf zamanlaması | Integrated gradients (taban: sıfır sinyal, 64 adım). EKG başına mutlak atfın QRS başlangıcı ±40 ms'ye düşen payı ve ST-T'ye (bitiş +40…+400 ms) düşen payı | Bu pay vaka EKG'lerinde, eşleştirilmiş kontrollerden yüksek: küme düzeyinde ortalama farkın bootstrap %95 GA alt sınırı > 0 (H2a için başlangıç payı, H2b için ST-T payı). P dalgasında yoğunluk betimsel olarak raporlanır |
| 5.3 Kendi delta vektörüyle uyum | Vakanın T0'ındaki (a) başlangıç vektörü `i20`, (b) atıf ağırlıklı başlangıç profili ile aynı hastanın T2'sindeki ilk 20 ms eğim profilinin korelasyonu (8 derivasyon) | Ortalama korelasyon, vakalar arası permütasyon dağılımının üstünde (tek yönlü p < 0.05) |
| 5.4 Kardiyak bellek | T integral vektörünün hastanın kendi T2 QRS integral vektörüyle korelasyonu; skorun son T2'ye uzaklıkla ilişkisi (0–30 gün, 30–365 gün, >365 gün) | Permütasyonun üstünde uyum; yakın T0'larda daha yüksek skor (H2b) |
| 5.5 Polarite kuralı | `pqlow_fq_corr` dağılımı, vaka ve kontrol | Vaka EKG'lerinde Gibbs benzeri oran kontrolden yüksek değil |

**Çoklu karşılaştırma:** 5.1–5.4'ün iki veri setindeki testleri tek bir aile; Holm düzeltmesi uygulanır.

---

## 6. Kullanıcının gözlediği QRS öncesi ayrık dalga (keşif → doğrulama)

**Tanım** (medyan vuru, 8 bağımsız derivasyon):
- **Pencere:** Düşük eşikli QRS başlangıcından önceki 60 ms. Taban çizgisi, derivasyon başına pencerenin medyanı.
- **Bir derivasyonda dalga:** Tepesi en geç başlangıçtan 4 ms önce olan, |genlik| ≥ 8 µV'lik yerel bir ekstremum. Ayrıca:
  - Tepeden önce ve tepeden sonra (QRS başlamadan) sapma, tepenin yarısının altına iner.
  - Yarı yükseklikteki genişlik 6–30 ms arasıdır.
  - Tepenin iki yanında türevin işareti değişir; P'nin tekdüze inişi elenir.
- **EKG'de dalga:** ≥2 bitişik derivasyonda (I–aVL, II–III–aVF, V1…V6 komşulukları), tepe zamanları birbirinden en fazla 6 ms farklı dalga var.
- **Kaydedilenler:** Zamanlama, genlik, genişlik, ilk QRS sapmasına göre polarite, 12 derivasyon profili.

**Analiz:**
- **MGH (keşif, eşleştirilmiş kümeler):** Sıklık; vaka ile kontrol arasında koşullu lojistik regresyon. Tabakalamalar: P süresi (IAB vekili: Aşama 2 P süresi ≥120 ms), PR ve cinsiyet.
- **Doğrulama:** Emory ve MIMIC'te aynı tanımla, dondurmadan sonra. Doğrulama ölçütü: eşleştirilmiş OR > 1 ve %95 GA 1'i dışlıyor.
- **Ek kontroller:** Polarite (ilk QRS sapmasının tersi mi), kendi delta vektörüyle uyum, MIMIC'te cart grubuna göre sıklık (ölçülen pre-ringing farklı olan "<not specified>" grubu).
- **Sınırlılık:** ~10 µV, kuantizasyon adımının (≈5 µV) iki katı. Medyan vuru alt-LSB çözünürlük sağlar, ama tek derivasyonlu dalgalar gürültüye duyarlıdır. Bu yüzden ≥2 derivasyon şartı konuldu.

---

## 7. Artefakt ve kısa yol denetimleri (concept §8)

| # | Denetim | Veri | Geçme ölçütü (önceden) |
|---|---|---|---|
| 7.1 | **T0 sıkılaştırma (etiket gürültüsü)** | MIMIC: makine PR ≥140 ms ve QRS <110 ms. HEEDB: Aşama 2 QRS <110 ms ve düzeltilmiş PR (pr_ms − 42) ≥140 ms | Sıkı T0'da AUROC, tüm T0'daki AUROC'un %95 GA'sı içinde ya da ≥0.60 |
| 7.2 | Filtre/cart tabakalama | MIMIC filtre gruplarına göre AUROC | Gruplar arası heterojenlik testi p > 0.05; "<not specified>" grubunda AUROC düşmüyor |
| 7.3 | Aynı kayıt içi (vuru-vuru intermittan) | ALTWPW (12SL 303) kayıtlarında normal ve preeksite küme | Normal kümenin QRS başlangıcında görünür iz yok (negatif kontrol); varsa §5.4 yönünde |
| 7.4 | K2 / AVNRT (H4) | MIMIC K2a (yalnızca AVNRT kalıbı, WPW sözü yok) | K2a skoru K1'den yüksek değil (fark %95 GA'sı 0'ı içeriyor ya da negatif); K2b ayrıca betimsel |
| 7.5 | Ablasyon (H3) | MIMIC, kesin T2 ile 5 hasta (herhangi T2 ile 9) | Yalnızca tanımlayıcı: ablasyon öncesi, 0–90 gün ve >90 gün skorları; test yok |
| 7.6 | Demografik kısa yol | M-demo | M1'in AUROC'u M-demo'dan yüksek (eşleştirilmiş bootstrap farkının alt sınırı > 0) |
| 7.7 | Eşlik eden hastalık | MIMIC'te ICD: kardiyomiyopati, HCM, Ebstein, konjenital kalp hastalığı, dal bloğu öyküsü | Bu hastalar çıkarıldığında AUROC %95 GA içinde kalıyor |

---

## 8. Duyarlılık analizleri

1. Etiket varyantları: HEEDB'de yalnızca overread ve yalnızca v24 12SL; MIMIC'te herhangi T2 (352 vaka).
2. Prospektif varyant: İlk T2'den önceki T0'lar (MGH'de 442 vaka uygun).
3. Son T2'ye uzaklık tabakaları (H2b ile ortak).
4. Birleştirme: Hasta düzeyinde ortalama ve maksimum (§3).
5. Eşleştirme sonrası kalan dengesizlik: dönem, kaydırılmış yıl ve log(EKG sayısı) kovaryatlı lojistik model.
6. Yaş <40.
7. Ham 10 s girişli M1.

---

## 9. İstatistik

- **Belirsizlik:** Hasta düzeyinde (eşleştirilmiş kümelerde küme düzeyinde) bootstrap, 2,000 tekrar, sabit tohum.
- **AUROC:** %95 GA bootstrap ile; DeLong yalnızca karşılaştırma için.
- **Kalibrasyon:** Kesişim, eğim ve ECE. Dış veride yeniden kalibrasyon yapılmaz; ayrıca yeniden kalibre edilmiş sonuç ikincil olarak verilir.
- **PPV:** %0.1, %1 ve %10 prevalans senaryoları (concept §9).
- **Eşleştirilmiş betimsel farklar:** Koşullu lojistik regresyon.
- **Çoklu karşılaştırma:** Tek birincil sonlanım (H1, MIMIC AUROC). H2 ailesi Holm ile düzeltilir. Diğerleri GA ile, düzeltmesiz ve "ikincil/keşifsel" etiketiyle raporlanır.
- **Raporlama:** TRIPOD+AI.

---

## 10. Çıktılar

- `stage4/` kodu ve `stage4/out/` (agregat).
- Model ağırlıkları `~/data/hidden-delta/stage4/` (paylaşılmaz).
- `FROZEN.sha256` (plan, ön işleme ve model dosyalarının özetleri).
- Sonuç raporu: birincil sonlanım, mekanizma tablosu, denetim tablosu, duyarlılık tablosu ve akış diyagramı.

---

## 11. Sırayla iş adımları

1. Emory ve MIMIC'te etiketsiz sinyal çıkarımı (Aşama 2 kodu).
2. Emory ve MIMIC'te eşleştirme (MIMIC'te filtre tabakasıyla); MGH kontrollerinde ek T0 çıkarımı (§3 ikincil birim).
3. MGH'de M0/M1/M2 geliştirme (5 katlı çapraz doğrulama), eşik seçimi.
4. Dondurma (`FROZEN.sha256`).
5. Emory ve MIMIC skorlama; §1, §5, §6, §7 ve §8 analizleri.
6. Rapor.

---

## 12. Kararlar (kullanıcı onayı, 2026-11-08)

1. Birincil sonlanım: MIMIC'te hasta düzeyinde AUROC.
2. Birincil model: M1 (medyan vuru üzerinde derin öğrenme). M0 önemli ikincil modeldir.
3. Plan OSF'ye ön kayıt olarak yüklenecek. Form metinleri `osf/` altında hazırlanır; yükleme kullanıcı hesabıyla yapılır. Emory ve MIMIC'te etiketli herhangi bir analizden önce tamamlanır.

---

## 13. Değişiklik günlüğü

- v1.0 taslak (2026-11-08): ilk yazım.
- v1.0 (2026-11-08): §12'deki üç karar kullanıcı tarafından onaylandı; plan donduruldu.
- v1.1 (2026-11-08, sonuç görülmeden, OSF metni hazırlanırken): §5.2'nin ölçütü test edilebilir hale getirildi (atıf payı, vaka − eşleştirilmiş kontrol, bootstrap alt sınır > 0). Başka değişiklik yok. Dondurulmuş kopyalar `frozen/` altında.
- v1.2 (2026-11-08, sonuç görülmeden): §3'e IAB okuma dışlaması eklendi. Gözlenen dalganın görüldüğü IAB çalışmalarının Emory ve MIMIC'ten olduğu anlaşıldı. O okumalarda EKG'si bulunan hastalar test setlerinden çıkarıldı (MIMIC'te 149, Emory'de 43 hasta; ikisinde de en fazla 4 vaka). Eşleştirme yeniden yapıldı: MIMIC 103 vaka ve 461 kontrol, Emory 173 vaka ve 755 kontrol.

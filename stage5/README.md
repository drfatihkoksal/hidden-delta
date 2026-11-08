# Aşama 5 — kayıtlı test (MIMIC, Emory)

ANALYSIS_PLAN v1.2 (`frozen/ANALYSIS_PLAN_v1.2.md`); OSF https://osf.io/kvn28 (ambargolu).

| Script | İş | Not |
|---|---|---|
| `common5.py` | Veri yükleme, türetilmiş özellikler, dondurulmuş modeller, oklüzyon, integrated gradients, koşullu lojistik, Holm, DeLong, kalibrasyon | — |
| `s5_prepare.py` | Test girdileri: kümeler (birincil, K2, AVNRT→K1, prospektif, etiket varyantları), indeks EKG'leri, T2 profilleri, sıkı T0, yapısal kalp hastalığı kodları, son T2'ye uzaklık, ablasyon serisi, ALTWPW | Skor yok, vaka-kontrol karşılaştırması yok |
| `s5_freeze.py` | `FROZEN_MANIFEST.sha256`: plan, kod, final modeller, test girdileri | Manifest özeti OSF'ye güncelleme olarak eklenir |
| `s5_score.py` | M1/M0/M2/M-demo skorları; oklüzyon (başlangıç ±40 ms, ST-T, aynı genişlikte rastgele); integrated gradients payları; ablasyon serisi skorları | Manifest doğrulanmadan çalışmaz; bir kez |
| `s5_supplement.py` | **Kayıt sonrası, resmî sonuçlardan sonra yazıldı.** Planın betimsel istediği ama analizde eksik kalanlar: 5.2 P atıf payı, §6 dalga ek kontrolleri (polarite, kendi delta yönü, MIMIC filtre grubu), 7.2 ikinci koşul, K2b sayısı → `out/supplement.md` | Test yok; yeni skorlama yok |
| `s5_analysis.py` | `site`: H1 + ikincil, H2 (5.1–5.4, ham p), 5.5, H4, ayrık dalga (koşullu lojistik), 7.1/7.2/7.3/7.5/7.7, duyarlılık. `combine`: Holm ve karar kuralları → `out/results.md` | Bir kez |

## Sıra

```
python3 stage5/s5_prepare.py MIMIC I0006          # yapıldı (2026-11-08)
python3 stage5/s5_freeze.py                       # manifest; özeti OSF'ye güncelleme olarak eklenir
# — OSF güncellemesi yapıldıktan sonra —
python3 stage5/s5_score.py MIMIC && python3 stage5/s5_score.py I0006
python3 stage5/s5_analysis.py site MIMIC && python3 stage5/s5_analysis.py site I0006
python3 stage5/s5_analysis.py combine
python3 stage5/s5_supplement.py                   # kayıt sonrası betimsel ekler
```

## Prova (`--rehearsal`)

Çıktı `~/data/hidden-delta/stage5/rehearsal/`.
- **MGH:** gerçek modelle ve gerçek etiketlerle (geliştirme verisi; final modeller MGH ile eğitildiği için örneklem içi).
- **MIMIC ve Emory:** rastgele ağırlıklı modelle, ayrıca her kümede "vaka" rastgele bir üyeyle değiştirilerek. Dondurmadan önce hiçbir gerçek test skoru üretilmez.
- **Yan etki:** Prova, ayrık dalganın test setlerindeki genel sıklığını (vaka-kontrol ayrımı olmadan) gösterdi: MIMIC ~%60–66, Emory ~%37–47. OSF güncellemesinde açıklanır.

## Uygulama notları (plan metninin yorumlandığı yerler)

- **5.1:** Maskeleme, medyan vuruda pencerenin sıfıra (taban çizgisi) çekilmesidir. Rastgele pencere başlangıcı tekdüze ve sabit tohumludur. Test istatistiği: AUROC(rastgele maskeli) − AUROC(hedef maskeli). p, küme bootstrap dağılımının ≤ 0 oranıdır.
- **5.2:** Mutlak IG'nin pencere payı. Küme başına "vaka − kontrol ortalaması"; küme bootstrap'ı.
- **5.3:** Birincil test (a): ilk 20 ms derivasyon profili ile hastanın T2 eğim profili; 10,000 permütasyon. Atıf ağırlıklı profil (b) uygulanmadı.
- **5.4:** İndeks EKG'nin ST-T integral profili (bitiş +40…+400 ms, 8 derivasyon) ile hastanın T2 QRS integral profili; permütasyon. Zaman tabakaları betimsel.
- **7.7:** Plan listesi uygulandı: kardiyomiyopati, HCM, Ebstein, konjenital kalp hastalığı ve dal bloğu kodları. OSF metni dal bloğunu saymıyor; plan esas alındı.

## Dondurmadan önce eklenen duyarlılık analizleri (plan §8)

- **§8.4 hasta düzeyinde ortalama/maksimum** (`s5_prepare.py agg`):
  - Vakanın bütün T0'ları kullanılır; her kontrolden vakanın T0 sayısı kadar rastgele T0 alınır.
  - Eksik kontrol EKG'leri etiketsiz olarak çıkarılır: MIMIC'te 2,582, Emory'de 1,603.
  - M1 logitinin hasta düzeyinde ortalaması ve maksimumu üzerinden AUROC.
- **§8.7 ham 10 s girişli M1** (`stage4/s4_m1raw.py`; `s5_prepare.py raw`):
  - M1'in yapılandırması korunur (8 blok, linear probe); girdi 12 × 5000.
  - MGH'de kat başına SSL ve final model.

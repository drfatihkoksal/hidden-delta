# Aşama 2 — sinyal pipeline'ı

concept.md §6 ve §6.1. Etiketten bağımsızdır.

| Dosya | İş |
|---|---|
| `hdsig.py` | Kütüphane: yükleme, taban çizgisi, vuru tespiti, hizalama, morfoloji kümeleme, medyan vuru, global QRS/P sınırları, Kors XYZ, §6.1 özellikleri |
| `s2_validate.py` | Sentetik testler (normal / delta / Gibbs) ve MIMIC makine aralıklarıyla uyum → `out/validation.md` |
| `s2_extract.py` | Toplu çıkarım (HEEDB sitesi): vaka EKG'leri, K1 havuzu, K2 örneklemi → `~/data/hidden-delta/stage2/<site>/` |
| `s2_filter_ir.py` | Pace spike'larından cihaz zincirinin dürtü yanıtı ve pre-ringing'i → `out/filter_ir.md` |
| `s2_g2_pilot.py` | G2 pilot polarite denetimi (concept.md §10; ölçütler önceden yazıldı) → `out/g2_pilot.md` |

```
python3 stage2/s2_validate.py --n 1000
python3 stage2/s2_extract.py I0001        # MGH; kaldığı parçadan sürer
python3 stage2/s2_filter_ir.py
python3 stage2/s2_g2_pilot.py
```

## Zincir

1. **Yükleme:** mV birimi, standart derivasyon sırası (MIMIC'in I, II, III, aVR, aVF, aVL sırası düzeltilir). Yalnızca 500 Hz.
2. **Taban çizgisi:** Medyan filtreyle (0.2 s + 0.6 s) kestirilip çıkarılır. Medyan filtre doğrusal değildir ve keskin kenarlarda Gibbs salınımı üretmez. Hiçbir doğrusal ya da sıfır fazlı alçak geçiren filtre uygulanmaz.
3. **Vuru tespiti ve medyan vuru:**
   - Tespit, 8 bağımsız derivasyonun (I, II, V1–V6) uzaysal hızından yapılır; 40 ms ortalama, 250 ms en kısa aralık.
   - Vurular, uzaysal hız şablonuna ±20 ms içinde çapraz korelasyonla hizalanır.
   - Morfoloji kümeleme: QRS penceresinde (−80…+100 ms) korelasyon ≥0.90, ortalama bağlantılı hiyerarşik kümeleme. Baskın küme ve ≥3 vurulu ikinci küme için medyan vuru alınır (−400…+600 ms). İkinci küme, vuru-vuru intermittan kayıtlarda preeksite ve normal vuruları ayırmak içindir.
4. **Global QRS sınırları** (medyan vurunun uzaysal hızından):
   - `on_low`: gürültü tabanı + 4 SD ya da tepenin %4'ü; ilk elektriksel aktivite.
   - `on_high`: tepenin %15'i; görünür başlangıç.
   - `slow_init_ms = on_high − on_low`: yavaş başlangıç bileşeni.
   - Bitiş: son anlamlı QRS tepesinden sonra eşiğin altına ilk düşüş.
   - P sınırları: P tepe hızının %20'si.
5. **Referans düzeyler:** İki ayrı referans kullanılır.
   - Morfoloji ve vektör özellikleri (q, V1, başlangıç vektörleri, QRS/T integralleri): kayıt düzeyindeki medyan filtre taban çizgisi.
   - QRS öncesi lob özellikleri: yerel PR-sonu düzeyi, yani düşük eşikli başlangıçtan 40–22 ms önce; ölçüm penceresi −20…−2 ms. Kısa PR'de bu pencere P'ye taşabilir. Penceredeki eğim `pq_base_slope` olarak saklanır.
   - İlk sürümdeki TP'ye geçiş kuralı kayıtların %73–85'inde tetikleniyordu; Ta dalgasını sahte pre-QRS lobu olarak gösterme riski nedeniyle kaldırıldı.
6. **Kors XYZ** (Kors 1990; X+ sol, Y+ aşağı, Z+ arka). Katsayılar Vondrak ve Penhakert 2022 (IEEE JTEHM, Tablo 1) ile birebir doğrulandı.

## Özellikler (`features_*.parquet`)

Zamanlar fiducial'a göre ms, genlikler mV.

- **QRS sınırları:** `on_low_ms`, `on_high_ms`, `off_ms`, `slow_init_ms`, `qrs_ms`, `qrs_high_ms`, `onset_capped`
- **P ve PR:** `p_on_ms`, `p_end_ms`, `pr_ms`, `pr_seg_ms`; QRS öncesi referans penceresinin kalitesi: `pq_base_slope`, `pq_base_level`
- **Başlangıç vektörleri:** `v{10,20,40}_{x,y,z,mag}` (anlık), `i{20,40}_{x,y,z}` (ortalama), `slope20_<derivasyon>`, `svr{5,10,20}`
- **Septal q ve V1:** `q_{I,aVL,V5,V6,V1}`; `v1_r`, `v1_s`, `v1_rs`, `v1_rsr`
- **QRS öncesi defleksiyon** (`pqlow_*` düşük eşikli, `pqhigh_*` görünür başlangıca göre):
  - derivasyon başına değer, XYZ ve büyüklük;
  - `*_fq_corr`: QRS öncesi lobun 12 derivasyonluk profiliyle o başlangıçtan sonraki ilk 15 ms'lik QRS sapmasının korelasyonu. Filtre salınımında negatif, delta'da pozitiftir.
  - `*_opp_frac`: zıt polariteli derivasyonların oranı.
- **Repolarizasyon:** QRS ve T integral vektörleri, T tepe vektörü, `qrst_angle`, `vg_mag`
- **Kayıt düzeyi:** `n_beats`, `n_clusters`, `dom_frac`, `within_corr`, `hr`, `rr_sd`, `n_flat`. İkinci küme özellikleri `c2_` önekiyle.

## Doğrulama (`out/validation.md`)

- **Sentetik testler:** Delta'da `pqhigh_fq_corr` +0.95, zıt polarite oranı 0 ve yavaş başlangıç 42 ms. Gibbs'te `pqhigh_fq_corr` −1.00 ve zıt polarite oranı 1.00. Normal vurularda QRS öncesi lob yok (0.004 mV). Delta'da sahte septal q yok.
- **MIMIC (1,000 rastgele T0, etiket kullanılmadan):**
  - QRS süresi makineyle uyumlu: medyan fark −2 ms, mutlak fark medyanı 8 ms, %83'ü ±20 ms içinde.
  - PR'yi makineden sistematik olarak +42–44 ms uzun ölçüyoruz (Spearman ~0.5). P eşiğini 0.2–0.5 arasında değiştirmek uyumu anlamlı biçimde iyileştirmedi, bu yüzden 0.20 korundu. Sentetik testte doğru ölçüldüğü için fark büyük olasılıkla P başlangıcı tanımından geliyor.
  - Başlangıç %6.6 kayıtta 140 ms sınırına dayanıyor (`onset_capped` ile işaretli). MGH'de bu oran %0.8.

## MGH çıkarımı (2026-11-08)

`~/data/hidden-delta/stage2/I0001/` (2.2 GB): 47,301 EKG işlendi, %99.6'sı hatasız.

| Rol | EKG | Başarı |
|---|---|---|
| Vaka T0 | 21,243 | %99.6 |
| Vaka T2 | 3,058 | %99.1 |
| K1 havuzu | 20,000 | %99.6 |
| K2 | 3,000 | %99.8 |

- Hatalar: az vuru ve düz derivasyon.
- Kayıtların %13.2'sinde birden çok morfoloji kümesi var; 1,676 kayıtta ikinci küme medyan vurusu çıkarıldı.
- QRS süresi medyanı 88 ms; yavaş başlangıç medyanı 10 ms; başlangıcın 140 ms sınırına dayandığı kayıt %0.8.

## Cihaz dürtü yanıtı (`out/filter_ir.md`)

- **MIMIC:** Pace spike'larında t = −4 ms'de spike'ın tersi işaretli ~%15'lik bir lob ve −6…−8 ms'de ~%3'lük bir lob var. Bu doğrusal fazlı bir pre-ringing.
  - Lob notch/baseline/bant ayarından bağımsız; kaynağı örnekleme zinciri.
  - "<not specified> | 0.0005–150 Hz" grubunda lob daha küçük (−0.065). Farklı bir cart modeli olması muhtemel; §8'deki filtre tabakalaması için doğal bir karşılaştırma.
- **HEEDB:** Spike'lar tek örnek genişliğinde ve öncesi tam sıfır. Büyük olasılıkla cihazın eklediği yapay pace işaretleri; dürtü yanıtı bu yolla tahmin edilemiyor.
- **İlk sürüm geçersizdi:** Gevşek dedektör R tepeleri gibi başka keskin yapıları da topladığı için tahmin geçersizdi. Sıkı ölçütlerle düzeltildi; ölçütler raporun başında.

## Bilinen sınırlılıklar

- PR ölçümü yukarıdaki nedenle yalnızca göreli karşılaştırmalarda kullanılmalı. MIMIC'teki sıkı T0 tanımı makine PR'sini kullanır.
- HEEDB'de filtre ayarı bilgisi yok ve pace spike'ları yapay işaret gibi görünüyor. Bu yüzden HEEDB'de beklenen filtre salınımı modellenemez; filtre artefaktı denetimi orada polarite kuralına (`pqhigh_fq_corr`, `pqhigh_opp_frac`) dayanır.

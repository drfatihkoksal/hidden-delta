# Aşama 4 — model geliştirme (yalnızca MGH)

ANALYSIS_PLAN v1.2 §3–4 (dondurulmuş: `frozen/ANALYSIS_PLAN_v1.2.md`). OSF: https://osf.io/kvn28.
Emory ve MIMIC bu aşamada skorlanmaz; modeller ve test kodu dondurulup manifest özeti OSF'ye eklendikten sonra skorlanır.

| Script | İş | Çıktı |
|---|---|---|
| `s4_data.py` | Eşleşmiş MGH kümeleri (vaka + K1, 1:5), özellikler, medyan vurular, ayrık dalga bayrağı, 5 kat (kümeye göre); SSL havuzu (T2 hariç tüm başarılı EKG'ler) | `~/data/hidden-delta/stage4/I0001_{dev,ssl}.parquet`, `*_mb.npy` |
| `s4_m0.py` | M0 (sabit özellikler + yaş + cinsiyet; 3 düğümlü kısıtlı kübik spline, L2 lojistik) ve M-demo (yaş, cinsiyet, KH) | `out/m0_mgh.md`, `I0001_oof_m0.parquet` |
| `s4_m1.py` | M1: 1D ResNet, derinlik {4, 8} × {ince ayar, linear probe}. Katın test hastaları hariç SSL; filtre augmentasyonu; bf16; Platt; %95 özgüllük eşiği; final model tüm MGH ile | `out/m1_mgh.md`, `~/data/hidden-delta/stage4/m1/` |
| `s4_m2.py` | M2: M1 SSL kodlayıcı gömmesi (128) + M0 özellikleri, L2 lojistik | `out/m2_mgh.md`, `final/M2.json` |
| `s4_m1raw.py` | §8.7 duyarlılık: ham 10 s girişli M1 (M1 yapılandırması; 12 × 5000) | `out/m1raw_mgh.md`, `~/data/hidden-delta/stage4/m1raw/` |
| `m0lib.py`, `statlib.py` | Ortak tasarım (spline, ölçek, lojistik; parametreler JSON'da) ve hızlı küme bootstrap'ı / AUROC | — |

```
python3 stage4/s4_data.py
python3 stage4/s4_m0.py
python3 stage4/s4_m1.py          # GPU; ~9 dk
python3 stage4/s4_m2.py          # ~15 s
```

## Uygulama notları

- **Ayrık dalga (plan §6):** Plan hem "8 bağımsız derivasyon" hem "I–aVL, II–III–aVF komşulukları" diyor. Komşuluk şartı türetilmiş derivasyonlar olmadan anlamsız olduğundan dalga 12 derivasyonda aranır (`hdsig.discrete_wave`, `WAVE_PAIRS`).
- **V1 R/S:** S < 0.01 mV ise 0.01 mV ile bölünür (`v1_rs_c`). Bu bir tanım; impütasyon değil.
- **M0'da eksik özellik:** Eksik değeri olan hastalar M0'dan çıkar (MGH'de 131; çoğu P bulunamadığı için PR eksik). İmpütasyon yok.
- **T2 EKG'leri** hiçbir aşamada, SSL dahil, kullanılmaz (concept §7).

## MGH geliştirme sonuçları (hipotez testi değil)

- **Ayrık dalga sıklığı** (eşleşmiş kümeler): vaka %59.6, kontrol %57.6.
- **M0** OOF AUROC 0.649 (0.628–0.670); **M-demo** 0.537 (0.528–0.547); fark +0.113 (0.090–0.138). Ayrıntı: `out/m0_mgh.md`.
- **M1:** seçilen yapılandırma 8 blok, linear probe (SSL kodlayıcısı + L2 lojistik). OOF AUROC 0.665 (katlarda 0.645–0.687). %95 özgüllük eşiği 0.362; bu eşikte duyarlılık 0.185.
- **M2:** OOF AUROC 0.676 (0.655–0.697); M2 − M1 +0.011 (−0.001 – +0.023).

| Model (aynı 3,919 hasta) | OOF AUROC | %95 GA |
|---|---|---|
| M2 | 0.676 | 0.655–0.697 |
| M1 (birincil) | 0.665 | 0.644–0.686 |
| M0 | 0.649 | 0.627–0.671 |
| M-demo (4,050 hasta) | 0.537 | 0.528–0.547 |
| M1 ham 10 s (§8.7; 4,050 hasta) | 0.658 | 0.638–0.679 |

## Hız notları

- **bf16 karma hassasiyet:** M1'de ölçülebilir hızlanma sağlamadı. Model küçük ve nvidia-smi'deki %99 değeri SM doluluğunu değil, çekirdek çalıştığını gösteriyor.
- **Bootstrap:** pandas birleştirmesi yerine numpy indeksleri ve vektörel sıra (`statlib`); 2,000 tekrar 0.6 s (önce dakikalar).
- **Lojistik CV:** bütün çekirdeklerde (`n_jobs=-1`). M2 13.5 s'de bitti; önceki sürüm 12 dakikada bitmemişti.

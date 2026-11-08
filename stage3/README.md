# Aşama 3 — vaka-kontrol eşleştirmesi

concept.md §5.3 "Eşleştirme uygulaması" (kurallar sonuçlardan önce yazıldı; sapmalar orada).

```
python3 stage3/s3_match.py I0001    # MGH (betimsel karşılaştırma dahil)
python3 stage3/s3_match.py I0006    # Emory (yalnızca denge)
python3 stage3/s3_match.py MIMIC    # MIMIC (yalnızca denge; filtre ayarı tabakası; ablasyon öncesi indeks; H4 eşleştirmesi)
```
Deterministik: sabit tohum, sıralı girdi. Genelleştirmeden sonra MGH eşleşme dosyasının md5 özeti değişmedi.

Test setlerinde (Emory, MIMIC) vaka-kontrol özellik karşılaştırması model dondurulana kadar yapılmaz (ANALYSIS_PLAN §2).

- **Vaka:** Birincil uzlaşı tanımı. İndeks EKG, hasta başına bir rastgele başarılı 500 Hz T0 (yaşı bilinen).
- **Kontrol:** K1 (1:5) ve K2 (1:2); Aşama 2'de çıkarılmış tek T0.
- **Kesin tabakalar:** Cinsiyet; KH bandı; EKG sayısı bandı (2–3, 4–6, 7–12, 13–24, 25–49, ≥50). HEEDB'de ayrıca dönem (WFDB klasör yılı, 5 yıllık bant); MIMIC'te ayrıca cihaz filtre ayarı (`filtering | bandwidth`).
- **MIMIC'e özgü:** Vakanın indeks T0'ı ablasyon öncesidir. H4 için AVNRT (K2a) hastaları K1'e 1:5 eşleştirilir (`pool = k2a_vs_k1`).
- **Kaliperler:** Yaş ±3 yıl, log(EKG sayısı) ±0.25. Açgözlü, yerine koymadan.
- **Çıktılar:**
  - `stage3/out/match_<site>.{md,json}`: denge tablosu (agregat). MGH'de eşleşme sonrası betimsel farklar da var.
  - `~/data/hidden-delta/stage3/<site>_matched.parquet`: `set_id`, `pid`, `k`, `grp` (case/k1/k2/k2a), `prospective`, `pool`. Hasta düzeyi; paylaşılmaz.

Varsayım: HEEDB tarihleri kaydırılmış. WFDB klasör yılı gerçek takvim yılı gibi görünüyor (1987–2022), ama kaydırılmış tarihle farkı hasta içinde sabit değil. Bu yüzden yalnızca dönem vekili olarak kullanılıyor.

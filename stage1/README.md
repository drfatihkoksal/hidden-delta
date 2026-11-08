# Aşama 1 — etiketler (manuel okuma yok)

concept.md §5.2: Etiketler makine yorumlarından ve HEEDB'deki kayıt anı hekim overread'inden tanımlanır.

```
python3 stage1/s1_labels.py        # Aşama 0 PRIVATE çıktılarını kullanır; birkaç dakika
```

| Veri | T2 (birincil) | T2 (duyarlılık) | T0 (indeks) |
|---|---|---|---|
| HEEDB (MGH, Emory) | Üç kaynaktan ≥2'sinde preeksitasyon kodu (300/302/303/304): v24 12SL, kayıt anı yazılım, hekim overread | Yalnızca hekim overread; yalnızca v24 12SL | Overread var; hiçbir kaynakta preeksitasyon ya da kısa PR kodu yok; sinüs; AF/flutter, pace ve dal bloğu yok; yaş ≥18; **500 Hz** |
| MIMIC | Kesin preeksitasyon ifadesi ("Possible / Consider / Borderline" içermeyen) | Herhangi preeksitasyon ifadesi | Aşama 0 T0 tanımı (MIMIC'in tamamı 500 Hz) |

- **Sıkı T0** (§8 etiket gürültüsü denetimi; yalnızca MIMIC): T0 + makine PR ≥140 ms + QRS <110 ms. HEEDB'de aralık ölçümü yok; Aşama 2 pipeline'ı hesaplayacak.
- **Vaka hastası:** ≥1 T2 + ≥1 T0.
- **Kaynaklar arası uyum** (HEEDB): üç kaynağın da olduğu EKG'lerde preeksitasyon kodu için ikili Cohen κ, pozitif uyum ve Fleiss κ.

Çıktılar:
- `stage1/out/labels_summary.{md,json}`: yalnızca agregat.
- `~/data/hidden-delta/stage1/`: hasta düzeyi, paylaşılmaz.
  - `mimic_cohort.parquet`, `mimic_ecg_labels.parquet`
  - `heedb_{I0001,I0006}_cohort.parquet`, `heedb_{I0001,I0006}_ecg_labels.parquet`

## Kullanılmayan okuma aracı

`s1_build_queue.py` ve `reader/` bir körleştirilmiş okuma arayüzüdür. 2026-11-08'de manuel okuma
kaldırıldığı için kullanılmıyor. Kuyruklar `~/data/hidden-delta/reading/` altında ve hiç etiket
içermiyor. Silinebilir; ileride bir alt örneklem denetimi gerekirse diye duruyor.

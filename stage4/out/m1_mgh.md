# M1 — MGH 5 katlı çapraz doğrulama (geliştirme; hipotez testi değil)

|   derinlik | mod      |   OOF AUROC | kat AUROC'ları                    | en iyi epoch'lar   |
|-----------:|:---------|------------:|:----------------------------------|:-------------------|
|          4 | finetune |      0.6434 | 0.662, 0.693, 0.657, 0.655, 0.613 | 5, 4, 3, 12, 11    |
|          4 | probe    |      0.6601 | 0.664, 0.689, 0.674, 0.658, 0.612 | —                  |
|          8 | finetune |      0.6625 | 0.689, 0.708, 0.689, 0.646, 0.595 | 4, 9, 2, 12, 3     |
|          8 | probe    |      0.6649 | 0.678, 0.687, 0.653, 0.665, 0.645 | —                  |

Seçilen yapılandırma: derinlik 8, probe (OOF AUROC en yüksek).
Platt kalibrasyonu OOF skorlarla; %95 özgüllük eşiği (OOF kontrollerin 95. persentili) 0.3622; bu eşikte OOF duyarlılık 0.185.
Final model tüm MGH ile eğitildi. Süre 512 s.

Not: geliştirme sayılarıdır; birincil test MIMIC'te dondurulmuş modelle (ANALYSIS_PLAN §1).

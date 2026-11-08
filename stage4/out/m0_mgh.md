# M0 ve M-demo — MGH 5 katlı çapraz doğrulama (geliştirme; hipotez testi değil)

| model   |   hasta |   vaka |   çıkan (eksik özellik) |   AUROC | %95 GA      |
|:--------|--------:|-------:|------------------------:|--------:|:------------|
| M0      |    3919 |    820 |                     131 |   0.649 | 0.628–0.670 |
| M-demo  |    4050 |    852 |                       0 |   0.537 | 0.528–0.547 |

M0 − M-demo AUROC farkı (ortak 3,919 hasta): +0.113, küme bootstrap %95 GA +0.090 – +0.138.

Final modeller tüm MGH ile eğitildi (`~/data/hidden-delta/stage4/final/`).

Not: bu sayılar geliştirme setindendir; birincil test MIMIC'te dondurulmuş modelle yapılır (ANALYSIS_PLAN §1).

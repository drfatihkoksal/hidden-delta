# M2 — MGH 5 katlı çapraz doğrulama (geliştirme; hipotez testi değil)

M1 gömmesi: derinlik 8 SSL kodlayıcısı (128 boyut) + M0 özellikleri; ortak hastalarda karşılaştırma.

| model   |   hasta (ortak) |   OOF AUROC | %95 GA      |
|:--------|----------------:|------------:|:------------|
| M2      |            3919 |       0.676 | 0.655–0.697 |
| M1      |            3919 |       0.665 | 0.644–0.686 |
| M0      |            3919 |       0.649 | 0.627–0.671 |

M2 − M1 AUROC farkı: +0.011 (küme bootstrap %95 GA -0.001 – +0.023).

Final M2 tüm MGH ile eğitildi (`final/M2.json`; kodlayıcı: M1 probe ise `m1/final_encoder.pt`).

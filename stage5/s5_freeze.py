"""HIDDEN-DELTA — model dondurma manifesti (ANALYSIS_PLAN §2; OSF https://osf.io/kvn28 güncellemesi için).

Plan, kod, final modeller, eşik/kalibrasyon ve hazırlanmış test girdileri SHA-256 ile `FROZEN_MANIFEST.sha256`'ya yazılır
(`sha256sum -c` biçimi). Manifestin kendi özeti OSF'ye güncelleme olarak eklenir; `s5_score.py` resmî çalıştırmadan önce
manifesti doğrular. Var olan manifestin üzerine yazmaz.

Çalıştırma: python3 stage5/s5_freeze.py
"""
import sys

import common5 as C

CODE = ["frozen/ANALYSIS_PLAN_v1.2.md",
        "stage0/common.py", "stage0/s0_mimic.py", "stage0/s0_heedb.py", "stage1/s1_labels.py",
        "stage2/hdsig.py", "stage2/s2_extract.py", "stage3/s3_match.py", "stage3/s3_iab_overlap.py",
        "stage4/s4_data.py", "stage4/m0lib.py", "stage4/statlib.py", "stage4/s4_m0.py", "stage4/s4_m1.py", "stage4/s4_m2.py", "stage4/s4_m1raw.py",
        "stage5/common5.py", "stage5/s5_prepare.py", "stage5/s5_score.py", "stage5/s5_analysis.py", "stage5/s5_freeze.py"]
MODELS = ["stage4/m1/final_encoder.pt", "stage4/m1/final_probe.npz", "stage4/m1/final_meta.json",
          "stage4/final/M0.json", "stage4/final/M-demo.json", "stage4/final/M2.json",
          "stage4/m1raw/final_encoder.pt", "stage4/m1raw/final_probe.npz", "stage4/m1raw/final_meta.json"]
TEST = [f"stage5/{s}/{f}" for s in ("MIMIC", "I0006") for f in ("sets.parquet", "index.parquet", "index_mb.npy", "t2.parquet",
                                                                 "agg.parquet", "agg_mb.npy", "index_raw.npy", "index_raw_ok.npy")] + \
       ["stage5/MIMIC/h3.parquet", "stage5/MIMIC/h3_mb.npy", "stage5/I0006/alt.parquet",
        "stage3/MIMIC_matched.parquet", "stage3/I0006_matched.parquet", "stage3/iab_read_patients.parquet"]


def main():
    out = C.ROOT / "FROZEN_MANIFEST.sha256"
    if out.exists():
        sys.exit("FROZEN_MANIFEST.sha256 zaten var; üzerine yazılmaz.")
    lines = ["# HIDDEN-DELTA model dondurma manifesti (sha256sum -c FROZEN_MANIFEST.sha256 ile doğrulanır)",
             "# plan ve kod (proje içi göreli yol)"]
    lines += [f"{C.sha256(C.ROOT / f)}  {f}" for f in CODE]
    lines.append("# final modeller, kalibrasyon, eşik (mutlak yol)")
    lines += [f"{C.sha256(C.D / f)}  {C.D / f}" for f in MODELS]
    lines.append("# test girdileri: eşleşmiş kümeler, indeks EKG'leri, T2 profilleri (mutlak yol)")
    lines += [f"{C.sha256(C.D / f)}  {C.D / f}" for f in TEST]
    out.write_text("\n".join(lines) + "\n")
    out.chmod(0o444)
    print(f"{out}\nmanifest SHA-256: {C.sha256(out)}")


if __name__ == "__main__":
    main()

# balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288

## Source Data
- archived results dir: `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/archived_results/balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288_20260505_021306`
- params csv: `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/archived_results/balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288_20260505_021306/params/run_params.csv`
- filter used by `peak_pr_stats.py`: `(gc_exp == 0 and acc_exp == 0) OR (gc_exp > 0 and acc_exp > 0)`
- included per-run rows: `288`

## Fixed Parameters
- `acc_exp`: `0.5`
- `acc_key`: `bedA`
- `aligner`: `bowtie2`
- `fragment_length`: `150`
- `gc_exp`: `0.6`
- `gc_key`: `gcA`
- `genome`: `ce11_1pct`
- `macs2_mode`: `narrow`
- `map_coverage_pct`: `3`
- `map_enrich`: `10`
- `map_exp`: `1.0`
- `map_sigma`: `1.5`
- `nb_k`: `1000000`
- `peakcaller`: `macs2`
- `read_length`: `38`
- `tf_exp`: `1.0`
- `tf_peak_count_ctrl`: `0`
- `tf_peak_count_treat`: `5`
- `tf_sigma`: `5`
- `use_control`: `True`

## Swept Parameters
- `coverage_ctrl`: `0.5, 1.0, 12.0, 16.0, 2.0, 24.0, 4.0, 8.0`
- `coverage_treat`: `10, 20, 5`
- `id_ctrl`: `0001_con, 0002_con, 0003_con, 0004_con, 0005_con, 0006_con, 0007_con, 0008_con, 0009_con, 0010_con, 0011_con, 0012_con, 0013_con, 0014_con, 0015_con, 0016_con, 0017_con, 0018_con, 0019_con, 0020_con, 0021_con, 0022_con, 0023_con, 0024_con, 0025_con, 0026_con, 0027_con, 0028_con, 0029_con, 0030_con, 0031_con, 0032_con, 0033_con, 0034_con, 0035_con, 0036_con, 0037_con, 0038_con, 0039_con, 0040_con, 0041_con, 0042_con, 0043_con, 0044_con, 0045_con, 0046_con, 0047_con, 0048_con, 0049_con, 0050_con, 0051_con, 0052_con, 0053_con, 0054_con, 0055_con, 0056_con, 0057_con, 0058_con, 0059_con, 0060_con, 0061_con, 0062_con, 0063_con, 0064_con, 0065_con, 0066_con, 0067_con, 0068_con, 0069_con, 0070_con, 0071_con, 0072_con, 0073_con, 0074_con, 0075_con, 0076_con, 0077_con, 0078_con, 0079_con, 0080_con, 0081_con, 0082_con, 0083_con, 0084_con, 0085_con, 0086_con, 0087_con, 0088_con, 0089_con, 0090_con, 0091_con, 0092_con, 0093_con, 0094_con, 0095_con, 0096_con, 0097_con, 0098_con, 0099_con, 0100_con, 0101_con, 0102_con, 0103_con, 0104_con, 0105_con, 0106_con, 0107_con, 0108_con, 0109_con, 0110_con, 0111_con, 0112_con, 0113_con, 0114_con, 0115_con, 0116_con, 0117_con, 0118_con, 0119_con, 0120_con, 0121_con, 0122_con, 0123_con, 0124_con, 0125_con, 0126_con, 0127_con, 0128_con, 0129_con, 0130_con, 0131_con, 0132_con, 0133_con, 0134_con, 0135_con, 0136_con, 0137_con, 0138_con, 0139_con, 0140_con, 0141_con, 0142_con, 0143_con, 0144_con, 0145_con, 0146_con, 0147_con, 0148_con, 0149_con, 0150_con, 0151_con, 0152_con, 0153_con, 0154_con, 0155_con, 0156_con, 0157_con, 0158_con, 0159_con, 0160_con, 0161_con, 0162_con, 0163_con, 0164_con, 0165_con, 0166_con, 0167_con, 0168_con, 0169_con, 0170_con, 0171_con, 0172_con, 0173_con, 0174_con, 0175_con, 0176_con, 0177_con, 0178_con, 0179_con, 0180_con, 0181_con, 0182_con, 0183_con, 0184_con, 0185_con, 0186_con, 0187_con, 0188_con, 0189_con, 0190_con, 0191_con, 0192_con, 0193_con, 0194_con, 0195_con, 0196_con, 0197_con, 0198_con, 0199_con, 0200_con, 0201_con, 0202_con, 0203_con, 0204_con, 0205_con, 0206_con, 0207_con, 0208_con, 0209_con, 0210_con, 0211_con, 0212_con, 0213_con, 0214_con, 0215_con, 0216_con, 0217_con, 0218_con, 0219_con, 0220_con, 0221_con, 0222_con, 0223_con, 0224_con, 0225_con, 0226_con, 0227_con, 0228_con, 0229_con, 0230_con, 0231_con, 0232_con, 0233_con, 0234_con, 0235_con, 0236_con, 0237_con, 0238_con, 0239_con, 0240_con, 0241_con, 0242_con, 0243_con, 0244_con, 0245_con, 0246_con, 0247_con, 0248_con, 0249_con, 0250_con, 0251_con, 0252_con, 0253_con, 0254_con, 0255_con, 0256_con, 0257_con, 0258_con, 0259_con, 0260_con, 0261_con, 0262_con, 0263_con, 0264_con, 0265_con, 0266_con, 0267_con, 0268_con, 0269_con, 0270_con, 0271_con, 0272_con, 0273_con, 0274_con, 0275_con, 0276_con, 0277_con, 0278_con, 0279_con, 0280_con, 0281_con, 0282_con, 0283_con, 0284_con, 0285_con, 0286_con, 0287_con, 0288_con`
- `id_treat`: `0001_treat, 0002_treat, 0003_treat, 0004_treat, 0005_treat, 0006_treat, 0007_treat, 0008_treat, 0009_treat, 0010_treat, 0011_treat, 0012_treat, 0013_treat, 0014_treat, 0015_treat, 0016_treat, 0017_treat, 0018_treat, 0019_treat, 0020_treat, 0021_treat, 0022_treat, 0023_treat, 0024_treat, 0025_treat, 0026_treat, 0027_treat, 0028_treat, 0029_treat, 0030_treat, 0031_treat, 0032_treat, 0033_treat, 0034_treat, 0035_treat, 0036_treat, 0037_treat, 0038_treat, 0039_treat, 0040_treat, 0041_treat, 0042_treat, 0043_treat, 0044_treat, 0045_treat, 0046_treat, 0047_treat, 0048_treat, 0049_treat, 0050_treat, 0051_treat, 0052_treat, 0053_treat, 0054_treat, 0055_treat, 0056_treat, 0057_treat, 0058_treat, 0059_treat, 0060_treat, 0061_treat, 0062_treat, 0063_treat, 0064_treat, 0065_treat, 0066_treat, 0067_treat, 0068_treat, 0069_treat, 0070_treat, 0071_treat, 0072_treat, 0073_treat, 0074_treat, 0075_treat, 0076_treat, 0077_treat, 0078_treat, 0079_treat, 0080_treat, 0081_treat, 0082_treat, 0083_treat, 0084_treat, 0085_treat, 0086_treat, 0087_treat, 0088_treat, 0089_treat, 0090_treat, 0091_treat, 0092_treat, 0093_treat, 0094_treat, 0095_treat, 0096_treat, 0097_treat, 0098_treat, 0099_treat, 0100_treat, 0101_treat, 0102_treat, 0103_treat, 0104_treat, 0105_treat, 0106_treat, 0107_treat, 0108_treat, 0109_treat, 0110_treat, 0111_treat, 0112_treat, 0113_treat, 0114_treat, 0115_treat, 0116_treat, 0117_treat, 0118_treat, 0119_treat, 0120_treat, 0121_treat, 0122_treat, 0123_treat, 0124_treat, 0125_treat, 0126_treat, 0127_treat, 0128_treat, 0129_treat, 0130_treat, 0131_treat, 0132_treat, 0133_treat, 0134_treat, 0135_treat, 0136_treat, 0137_treat, 0138_treat, 0139_treat, 0140_treat, 0141_treat, 0142_treat, 0143_treat, 0144_treat, 0145_treat, 0146_treat, 0147_treat, 0148_treat, 0149_treat, 0150_treat, 0151_treat, 0152_treat, 0153_treat, 0154_treat, 0155_treat, 0156_treat, 0157_treat, 0158_treat, 0159_treat, 0160_treat, 0161_treat, 0162_treat, 0163_treat, 0164_treat, 0165_treat, 0166_treat, 0167_treat, 0168_treat, 0169_treat, 0170_treat, 0171_treat, 0172_treat, 0173_treat, 0174_treat, 0175_treat, 0176_treat, 0177_treat, 0178_treat, 0179_treat, 0180_treat, 0181_treat, 0182_treat, 0183_treat, 0184_treat, 0185_treat, 0186_treat, 0187_treat, 0188_treat, 0189_treat, 0190_treat, 0191_treat, 0192_treat, 0193_treat, 0194_treat, 0195_treat, 0196_treat, 0197_treat, 0198_treat, 0199_treat, 0200_treat, 0201_treat, 0202_treat, 0203_treat, 0204_treat, 0205_treat, 0206_treat, 0207_treat, 0208_treat, 0209_treat, 0210_treat, 0211_treat, 0212_treat, 0213_treat, 0214_treat, 0215_treat, 0216_treat, 0217_treat, 0218_treat, 0219_treat, 0220_treat, 0221_treat, 0222_treat, 0223_treat, 0224_treat, 0225_treat, 0226_treat, 0227_treat, 0228_treat, 0229_treat, 0230_treat, 0231_treat, 0232_treat, 0233_treat, 0234_treat, 0235_treat, 0236_treat, 0237_treat, 0238_treat, 0239_treat, 0240_treat, 0241_treat, 0242_treat, 0243_treat, 0244_treat, 0245_treat, 0246_treat, 0247_treat, 0248_treat, 0249_treat, 0250_treat, 0251_treat, 0252_treat, 0253_treat, 0254_treat, 0255_treat, 0256_treat, 0257_treat, 0258_treat, 0259_treat, 0260_treat, 0261_treat, 0262_treat, 0263_treat, 0264_treat, 0265_treat, 0266_treat, 0267_treat, 0268_treat, 0269_treat, 0270_treat, 0271_treat, 0272_treat, 0273_treat, 0274_treat, 0275_treat, 0276_treat, 0277_treat, 0278_treat, 0279_treat, 0280_treat, 0281_treat, 0282_treat, 0283_treat, 0284_treat, 0285_treat, 0286_treat, 0287_treat, 0288_treat`
- `map_seed`: `11, 23, 37, 53, 71, 89`
- `seed`: `11, 23, 37, 53, 71, 89`
- `tf_enrich`: `1500, 2500`
- `tf_seed`: `11, 23, 37, 53, 71, 89`

## Plot Point Counts

| coverage_treat | coverage_ctrl | n_runs | total_called | total_planted | precision | recall | f1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 5 | 0.5 | 12 | 266 | 360 | 0.9211 | 0.6944 | 0.7919 |
| 5 | 1 | 12 | 344 | 360 | 0.9128 | 0.8889 | 0.9007 |
| 5 | 2 | 12 | 595 | 360 | 0.5765 | 0.9694 | 0.7230 |
| 5 | 4 | 12 | 1685 | 360 | 0.2059 | 0.9806 | 0.3404 |
| 5 | 8 | 12 | 1685 | 360 | 0.2059 | 0.9806 | 0.3404 |
| 5 | 12 | 12 | 1685 | 360 | 0.2059 | 0.9806 | 0.3404 |
| 5 | 16 | 12 | 1685 | 360 | 0.2059 | 0.9806 | 0.3404 |
| 5 | 24 | 12 | 1685 | 360 | 0.2059 | 0.9806 | 0.3404 |
| 10 | 0.5 | 12 | 267 | 360 | 0.9288 | 0.7000 | 0.7983 |
| 10 | 1 | 12 | 347 | 360 | 0.9251 | 0.9083 | 0.9166 |
| 10 | 2 | 12 | 420 | 360 | 0.8048 | 0.9556 | 0.8737 |
| 10 | 4 | 12 | 1538 | 360 | 0.2282 | 0.9917 | 0.3710 |
| 10 | 8 | 12 | 2448 | 360 | 0.1434 | 0.9944 | 0.2506 |
| 10 | 12 | 12 | 2448 | 360 | 0.1434 | 0.9944 | 0.2506 |
| 10 | 16 | 12 | 2448 | 360 | 0.1434 | 0.9944 | 0.2506 |
| 10 | 24 | 12 | 2448 | 360 | 0.1434 | 0.9944 | 0.2506 |
| 20 | 0.5 | 12 | 269 | 360 | 0.9294 | 0.7056 | 0.8021 |
| 20 | 1 | 12 | 342 | 360 | 0.9357 | 0.9056 | 0.9204 |
| 20 | 2 | 12 | 387 | 360 | 0.8786 | 0.9611 | 0.9180 |
| 20 | 4 | 12 | 1011 | 360 | 0.3442 | 0.9833 | 0.5099 |
| 20 | 8 | 12 | 2262 | 360 | 0.1556 | 0.9944 | 0.2691 |
| 20 | 12 | 12 | 2946 | 360 | 0.1195 | 0.9944 | 0.2133 |
| 20 | 16 | 12 | 2946 | 360 | 0.1195 | 0.9944 | 0.2133 |
| 20 | 24 | 12 | 2946 | 360 | 0.1195 | 0.9944 | 0.2133 |

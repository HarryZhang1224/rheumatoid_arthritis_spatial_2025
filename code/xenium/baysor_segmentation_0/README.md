# Xenium In Situ data segmentation with Baysor (v0.7.0)

---

## Installation
Refer to [Baysor GitHub](https://github.com/kharchenkolab/Baysor) and [10x Genomics Xenium Ranger (3.0.1)](https://www.10xgenomics.com/support/software/xenium-ranger/downloads/previous-versions).

---

## Pipeline Steps

1. **Filter Transcripts** (`filter_transcripts_0`)
   - This step filters out negative control probes and transcripts with a Q-score below 20 prior to running Baysor.

2. **Write Baysor parameters and run Baysor** (`run_baysor_1`)
   - This step writes a `config.toml` Baysor parameter configuration file and runs Baysor for each sample.

3. **Process Baysor outputs** (`process_baysor_outputs_2`)
   - This step process Baysor outputs to make sure they are compatible with Xenium Ranger input file format.

4. **Import Baysor output using Xenium Ranger** (`run_xenium_ranger_3`)
    - This step imports Baysor outputs using Xenium Ranger for downstream analysis.
<!-- Workflow Management -->
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0.0-brightgreen.svg?style=flat&logo=snakemake)](https://snakemake.readthedocs.io)
[![Workflow](https://img.shields.io/badge/workflow-reproducible-blue.svg?style=flat)](https://snakemake.readthedocs.io)

<!-- License and Documentation -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-passing-success.svg?style=flat)](./docs/)

<!-- Bioinformatics Tools -->
[![BWA](https://img.shields.io/badge/aligner-BWA--MEM-orange.svg?style=flat)](http://bio-bwa.sourceforge.net/)
[![Bowtie2](https://img.shields.io/badge/aligner-Bowtie2-orange.svg?style=flat)](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
[![Novoalign](https://img.shields.io/badge/aligner-Novoalign-orange.svg?style=flat)](http://www.novocraft.com/products/novoalign/)

<!-- Variant Callers -->
[![GATK](https://img.shields.io/badge/caller-GATK%20HaplotypeCaller-red.svg?style=flat)](https://gatk.broadinstitute.org/)
[![DeepVariant](https://img.shields.io/badge/caller-DeepVariant-red.svg?style=flat&logo=tensorflow)](https://github.com/google/deepvariant)
[![FreeBayes](https://img.shields.io/badge/caller-FreeBayes-red.svg?style=flat)](https://github.com/freebayes/freebayes)

<!-- Reference Data -->
[![GIAB](https://img.shields.io/badge/benchmark-GIAB-purple.svg?style=flat)](https://www.nist.gov/programs-projects/genome-bottle)
[![Reference](https://img.shields.io/badge/reference-GRCh38-purple.svg?style=flat)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/)

<!-- ML/AI -->
[![DeepLearning](https://img.shields.io/badge/ML-CNN-blueviolet.svg?style=flat&logo=tensorflow)](https://www.tensorflow.org/)
[![Python](https://img.shields.io/badge/python-3.8+-blue.svg?style=flat&logo=python)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-≥4.0-blue.svg?style=flat&logo=r)](https://www.r-project.org/)

<!-- Environment & Dependencies -->
[![Conda](https://img.shields.io/badge/package-conda-green.svg?style=flat&logo=anaconda)](https://docs.conda.io/)
[![Bioconda](https://img.shields.io/badge/install-bioconda-green.svg?style=flat)](https://bioconda.github.io/)
[![Container](https://img.shields.io/badge/container-ready-2496ED.svg?style=flat&logo=docker)](https://www.docker.com/)

<!-- Code Quality -->
[![Code Style](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat)](https://github.com/psf/black)
[![Linting](https://img.shields.io/badge/linting-pylint-yellowgreen.svg?style=flat)](https://www.pylint.org/)

<!-- Project Status -->
[![Status](https://img.shields.io/badge/status-active-success.svg?style=flat)]()
[![Maintained](https://img.shields.io/badge/maintained-yes-green.svg?style=flat)]()
[![Course](https://img.shields.io/badge/course-BINF6310-informational.svg?style=flat)]()

---

# Systematic Benchmark of Variant Calling Pipelines (HG002, GRCh38)

Reproducibility project for **BINF6310 (Northeastern University)** based on:

> Barbitoff, Y. A., Abasov, R., Tvorogova, V. E., Glotov, A. S., & Predeus, A. V. (2022).  
> *Systematic benchmark of state-of-the-art variant calling pipelines identifies major factors affecting accuracy of coding sequence variant discovery.* BMC Genomics, 23(1), 155. https://doi.org/10.1186/s12864-022-08365-3

Our goal is to reproduce (on a smaller scale) the WGS/WES variant-calling benchmark using **HG002** and two pipelines:

- **GATK HaplotypeCaller**
- **DeepVariant**

and then compare them against **GIAB truth sets** using **hap.py**.

---

## Dataset Strategy: Synthetic Normal/Tumor Data

### Background
The original plan for this project was to benchmark variant calling using the Genome in a Bottle (GIAB) samples HG001 (NA12878) and HG002 (NA24385), which are considered the gold-standard datasets for evaluating variant calling pipelines.

These samples provide high-confidence truth sets for SNPs and indels. However, the raw FASTQ datasets are multiple hundreds of gigabytes in size. The HPC environment used in this project restricts external downloads, bandwidth, and storage allocation.

### Synthetic Dataset Approach
To overcome these limitations while keeping the benchmarking scientifically valid, we created a synthetic dataset that includes:

- 2× Normal FASTQ samples  
- 2× Tumor FASTQ samples  
- Known somatic mutations injected into the tumor reads  
- Realistic Illumina-like sequencing noise  
- 150 bp paired-end reads  

Reads were aligned to the full GRCh38 reference genome using BWA-MEM and processed with samtools. The known mutations act as the ground truth for evaluating the accuracy of the variant calling pipeline.

### Advantages of Synthetic Data
- Exact control over truth variant set  
- Full reproducibility  
- Works within HPC constraints  
- Enables clean benchmarking (TP, FP, FN)  
- Small dataset footprint (<1 GB)  


# System Architecture Diagram 
```
  flowchart TD

subgraph Synthetic Data Layer
A[GRCh38 Reference] --> B[Extract Region for Simulation]
B --> C[Generate Normal Reads (2 samples)]
B --> D[Inject Somatic Mutations]
D --> E[Generate Tumor Reads (2 samples)]
C --> F[FASTQ Files]
E --> F
end

subgraph Processing Layer (HPC)
F --> G[bwa mem (SLURM)]
G --> H[samtools view]
H --> I[samtools sort + index]
I --> J[Sorted BAMs + QC]
end

subgraph Analysis Layer
J --> K[Merge Replicates Normal & Tumor]
K --> L[bcftools mpileup]
L --> M[bcftools call]
M --> N[VCF Output]
N --> O[Compare to Truth VCF]
O --> P[Precision/Recall/F1]
end

Q[Logs + Reports] --> P

```

### Future Work
If GIAB FASTQ datasets become available via shared storage or instructor provision, the pipeline can be applied directly to HG001/HG002 datasets with minimal changes.



## 1. Project Structure (HPC)

All paths are relative to:

```bash
BASE=/scratch/islam.mdtar/vbench_group/variant-calling-benchmark
```

## 📋 Group Project Overview

This project reproduces and extends the benchmarking study from:

> **Barbitoff et al. (2022)** - *Systematic benchmark of state-of-the-art variant calling pipelines identifies major factors affecting accuracy of coding sequence variant discovery*  
> DOI: [10.1186/s12864-022-08365-3](https://doi.org/10.1186/s12864-022-08365-3)

### 🎯 Key Objectives

- ✅ Reproduce variant calling benchmark on GIAB samples (HG001, HG002)
- ✅ Compare performance of **3 aligners** and **3 variant callers** (18 combinations)
- ✅ Evaluate precision, recall, and F1-score using built-in ML tools (DeepVariant, GATK CNN)
- ✅ Provide fully reproducible Snakemake workflow
- ✅ Generate comprehensive performance visualizations

---

## 🏗️ Project Structure

```
variant-calling-benchmark/
├── README.md
├── docs/
│   ├── rationale_synthetic_data.pdf
│   ├── architecture_diagram.png
│   ├── workflow_diagram.png
│   └── presentation_slides.pptx
├── data/
│   ├── reference/
│   │   ├── GRCh38.fa
│   │   ├── GRCh38.fa.fai
│   │   ├── GRCh38.dict
│   │   └── bwa_index/  (pac, bwt, ann, amb, sa)
│   ├── fastq/
│   │   ├── normal1_R1.fastq.gz
│   │   ├── normal1_R2.fastq.gz
│   │   ├── normal2_R1.fastq.gz
│   │   ├── normal2_R2.fastq.gz
│   │   ├── tumor1_R1.fastq.gz
│   │   ├── tumor1_R2.fastq.gz
│   │   ├── tumor2_R1.fastq.gz
│   │   ├── tumor2_R2.fastq.gz
│   └── truth/
│       ├── truth_variants.vcf   (positions you injected)
│       └── truth_bed.bed        (optional)
├── bam/
│   ├── normal1.sorted.bam
│   ├── normal1.sorted.bam.bai
│   ├── normal2.sorted.bam
│   ├── normal2.sorted.bam.bai
│   ├── tumor1.sorted.bam
│   ├── tumor1.sorted.bam.bai
│   ├── tumor2.sorted.bam
│   ├── tumor2.sorted.bam.bai
│   ├── merged_normal.bam
│   ├── merged_tumor.bam
│   └── flagstat/
│       ├── normal1.flagstat
│       ├── normal2.flagstat
│       ├── tumor1.flagstat
│       └── tumor2.flagstat
├── results/
│   ├── raw.vcf
│   ├── raw.vcf.gz
│   ├── raw.vcf.gz.csi
│   ├── filtered.vcf
│   ├── comparison/
│   │   ├── true_positive.txt
│   │   ├── false_positive.txt
│   │   └── false_negative.txt
│   └── metrics/
│       ├── precision_recall.csv
│       └── plots/
│           ├── ROC.png
│           └── PR_curve.png
├── scripts/
│   ├── generate_synthetic_reads.py
│   ├── alignment.slurm
│   ├── merge_bams.slurm
│   ├── call_variants.slurm
│   ├── compare_truth.ipynb
│   └── utils/
│       ├── generate_truth_vcf.py
│       └── qc_summary.py
├── env/
│   ├── environment.yml
│   └── requirements.txt
└── logs/
    ├── alignment_*.out
    ├── alignment_*.err
    ├── call_*.out
    ├── call_*.err
    └── qc_reports/

```

---

## 🚀 Quick Start

### Prerequisites

- **Conda** or **Mamba** (recommended for faster package resolution)
- **Snakemake** ≥7.0.0
- **Minimum 100GB** free disk space
- **16GB RAM** minimum (32GB recommended)
- **8+ CPU cores** for efficient parallel processing

### Installation

```bash
# 1. Clone repository
git clone git@github.com:mtariqi/variant-calling-benchmark.git
cd variant-calling-benchmark

# 2. Create conda environment
conda env create -f environment.yml
conda activate variant-benchmark

# 3. Set up project and download test data
bash scripts/setup_project.sh

# 4. Verify installation
snakemake --version
```

### Running the Pipeline

```bash
# Full workflow with 8 cores
snakemake --cores 8 --use-conda

# Dry run to preview steps
snakemake -n

# Generate workflow diagram
snakemake --dag | dot -Tpdf > workflow.pdf

# Run specific stage
snakemake --cores 8 --use-conda results/alignments/all

# Continue after interruption
snakemake --cores 8 --use-conda --rerun-incomplete
```

### Quick Test Run

```bash
# Test with reduced dataset
snakemake --cores 4 --use-conda --config subset=chr21
```

---

## 👥 Team Members & Contributions

| Member | Role | Primary Responsibilities | GitHub |
|--------|------|--------------------------|--------|
| **Md Tariqul Islam** | Data Management Lead | • Data acquisition & preprocessing<br>• Alignment pipeline development<br>• QC implementation | [@atra](https://github.com/) |
| **Raghad Al-Ampudi** | Pipeline Integration Lead | • Variant calling workflows<br>• Snakemake integration<br>• Environment management | [@raghad](https://github.com/) |
| **Atra Alimoradian** | Analysis & ML Lead | • Benchmarking analysis<br>• ML insights & visualization<br>• Final reporting | [@mtariqi](https://github.com/mtariqi) |

---

## 🔧 Pipeline Components

### Tools Evaluated

#### 🧬 Sequence Aligners
| Tool | Version | Algorithm | Speed | Accuracy |
|------|---------|-----------|-------|----------|
| **BWA-MEM** | 0.7.17 | Burrows-Wheeler Transform | Fast | High |
| **Bowtie2** | 2.4.5 | FM-index (local mode) | Very Fast | Moderate-High |
| **Novoalign** | 4.03.07 | Full dynamic programming | Moderate | Very High |

#### 🔍 Variant Callers
| Tool | Version | Method | ML Integration | Best For |
|------|---------|--------|----------------|----------|
| **GATK HaplotypeCaller** | 4.3.0 | Local de novo assembly | CNN filtering | SNPs & Indels |
| **DeepVariant** | 1.4.0 | Deep CNN | Native | High accuracy SNPs |
| **FreeBayes** | 1.3.6 | Bayesian haplotype-based | None | Population studies |

#### 📊 Benchmarking Tools
- **hap.py** (v0.3.14) - GIAB standardized benchmarking
- **bcftools** (v1.15) - VCF manipulation and stats
- **R/Bioconductor** - Statistical analysis and visualization

---
## Figure — Mutation Injection Strategy
ASCII Diagram (simple)
```
Reference Genome (GRCh38)
        |
        |  Extract Target Region (e.g., chr20:100000-100500)
        v
   Synthetic Normal Reads ──> FASTQ (R1/R2)
            |
            |  Inject Somatic Mutations
            v
   Synthetic Tumor Reads ───> FASTQ (R1/R2)

Mutations introduced:
    • SNPs at known positions
    • Small indels
    • Allele frequency modulation
    • Replicate variance
```

```
Reference Genome (GRCh38)
        |
        |  Extract Target Region (e.g., chr20:100000-100500)
        v
   Synthetic Normal Reads ──> FASTQ (R1/R2)
            |
            |  Inject Somatic Mutations
            v
   Synthetic Tumor Reads ───> FASTQ (R1/R2)

Mutations introduced:
    • SNPs at known positions
    • Small indels
    • Allele frequency modulation
    • Replicate variance
```


## 📈 Workflow Architecture

```
┌──────────────┐
│  Raw FASTQ   │
│  (GIAB Data) │
└──────┬───────┘
       │
       ├─────────────────────────────────┐
       │                                 │
       ▼                                 ▼
┌──────────────┐                  ┌──────────────┐
│  Quality QC  │                  │  Reference   │
│   (FastQC)   │                  │   Genome     │
└──────┬───────┘                  └──────┬───────┘
       │                                 │
       └────────────┬────────────────────┘
                    │
        ┌───────────┴───────────┐
        │   Parallel Alignment  │
        ├───────────────────────┤
        │  • BWA-MEM           │
        │  • Bowtie2           │
        │  • Novoalign         │
        └───────────┬───────────┘
                    │
        ┌───────────┴───────────┐
        │   BAM Processing      │
        │  • Sort               │
        │  • Mark Duplicates    │
        │  • Base Recalibration │
        └───────────┬───────────┘
                    │
        ┌───────────┴───────────┐
        │  Parallel Variant     │
        │  Calling              │
        ├───────────────────────┤
        │  • GATK HC           │
        │  • DeepVariant       │
        │  • FreeBayes         │
        └───────────┬───────────┘
                    │
        ┌───────────┴───────────┐
        │  Variant Filtering    │
        │  • GATK CNN          │
        │  • Hard Filters       │
        └───────────┬───────────┘
                    │
        ┌───────────┴───────────┐
        │  Benchmarking         │
        │  (hap.py vs GIAB)    │
        └───────────┬───────────┘
                    │
        ┌───────────┴───────────┐
        │  Statistical Analysis │
        │  & Visualization      │
        │  (R/ggplot2)         │
        └───────────────────────┘
```

---

## 📊 Expected Outcomes

### 1. Performance Metrics
- ✅ **Precision, Recall, F1-score** for each pipeline combination
- ✅ **Runtime and resource usage** benchmarks
- ✅ **Error profile analysis** (SNP vs Indel performance)

### 2. Comparative Analysis
- 📊 Heatmaps of pipeline performance across metrics
- 📈 ROC curves for ML-enabled callers
- 📉 Error distribution plots by variant type

### 3. Deliverables
- 📄 **Comprehensive report** (R Markdown → PDF/HTML)
- 🎨 **Publication-ready figures** (vector graphics)
- 💾 **Reproducible workflow** (Snakemake + Conda)
- 📦 **Archived results** on Zenodo (optional)

---

## 🗂️ File Structure with Ownership

```
variant-calling-benchmark/
├── 📁 workflows/
│   ├── main.smk                    👤 Member 2 (Integrator)
│   ├── alignment.smk               👤 Member 1 (Owner)
│   ├── variant_calling.smk         👤 Member 2 (Owner) 
│   ├── benchmarking.smk            👤 Member 3 (Owner)
│   └── qc.smk                      👤 Member 1 (Owner)
│
├── 📁 scripts/
│   ├── download_data.sh            👤 Member 1 (Owner)
│   ├── setup_project.sh            👤 Member 2 (Owner)
│   ├── alignment/                  👤 Member 1 (Owner)
│   ├── variant_calling/            👤 Member 2 (Owner)
│   └── benchmarking/               👤 Member 3 (Owner)
│
├── 📁 config/
│   ├── config.yaml                 👤 Member 2 (Integrator)
│   ├── alignment_config.yaml       👤 Member 1 (Owner)
│   ├── calling_config.yaml         👤 Member 2 (Owner)
│   └── benchmark_config.yaml       👤 Member 3 (Owner)
│
├── 📁 analysis/
│   ├── notebooks/
│   │   ├── exploratory_analysis.Rmd   👤 Member 3 (Owner)
│   │   ├── performance_analysis.Rmd   👤 Member 3 (Owner)
│   │   └── ml_analysis.Rmd            👤 Member 3 (Owner)
│   └── plots/                         👤 Member 3 (Owner)
│
├── 📁 docs/
│   ├── project_plan.md             👥 All (Collaborative)
│   ├── task_allocation.md          👥 All (Collaborative)
│   └── final_report.Rmd            👤 Member 3 (Lead)
│
└── environment.yml                  👤 Member 2 (Owner)
```

---

## 🧪 Testing & Validation

### Unit Tests
```bash
# Test individual workflow rules
snakemake --cores 1 --use-conda test_alignment
snakemake --cores 1 --use-conda test_variant_calling
```

### Continuous Integration
- Automated testing via GitHub Actions
- Workflow validation on push/PR
- Container builds for reproducibility

---

## 📚 Documentation

| Document | Description |
|----------|-------------|
| [Project Plan](docs/project_plan.md) | Detailed project roadmap and milestones |
| [Task Allocation](docs/task_allocation.md) | Responsibility matrix and timeline |
| [Configuration Guide](docs/configuration.md) | How to customize pipeline parameters |
| [Troubleshooting](docs/troubleshooting.md) | Common issues and solutions |

---

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Workflow
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## 📖 Citation

If you use this workflow in your research, please cite:

```bibtex
@software{variant_calling_benchmark2024,
  title = {Variant Calling Pipeline Benchmark},
  author = {Alimoradian, Atra and Al-Ampudi, Raghad and Islam, Md Tariqul},
  year = {2024},
  url = {https://github.com/mtariqi/variant-calling-benchmark},
  note = {BINF6300 Group Project}
}
```

**Original Study:**
```bibtex
@article{barbitoff2022systematic,
  title={Systematic benchmark of state-of-the-art variant calling pipelines identifies major factors affecting accuracy of coding sequence variant discovery},
  author={Barbitoff, Yury A and others},
  journal={BMC Genomics},
  volume={23},
  number={1},
  pages={155},
  year={2022},
  publisher={Springer}
}
```

---

## 📝 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- **GIAB Consortium** for providing gold-standard benchmarking datasets
- **Broad Institute** for GATK
- **Google Health** for DeepVariant
- **Snakemake** community for workflow management tools
- **Bioconda** for simplified tool installation
- **Course Instructor & TAs** for guidance and feedback

---

## 📞 Contact

- **Project Repository**: [github.com/mtariqi/variant-calling-benchmark](https://github.com/mtariqi/variant-calling-benchmark)
- **Issues**: [Report bugs or request features](https://github.com/mtariqi/variant-calling-benchmark/issues)
- **Discussions**: [Join the conversation](https://github.com/mtariqi/variant-calling-benchmark/discussions)

---

## 🔗 Useful Links

- [Snakemake Documentation](https://snakemake.readthedocs.io/)
- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
- [DeepVariant GitHub](https://github.com/google/deepvariant)
- [GIAB FTP Server](ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/)
- [hap.py Documentation](https://github.com/Illumina/hap.py)

---

<div align="center">

**⭐ Star this repo if you find it useful! ⭐**

Made with ❤️ by the BINF6310 Team

</div>

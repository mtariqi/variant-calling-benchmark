# Variant Calling Pipeline Benchmark

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A reproducible benchmark analysis of state-of-the-art variant calling pipelines for coding sequence variant discovery.

## 📋 Project Overview

This project reproduces and extends the benchmarking study from:
**Barbitoff et al. (2022)** - *Systematic benchmark of state-of-the-art variant calling pipelines identifies major factors affecting accuracy of coding sequence variant discovery*

**Key Objectives:**
- Reproduce variant calling benchmark on GIAB samples (HG001, HG002)
- Compare performance of 3 aligners and 3 variant callers
- Evaluate precision, recall, and F1-score using built-in ML tools (DeepVariant, GATK CNN)
- Provide fully reproducible Snakemake workflow

## 🏗️ Project Structure
```
variant-calling-benchmark/
├── workflows/ # Snakemake workflow definitions
├── scripts/ # Helper scripts for each analysis step
├── config/ # Configuration files
├── data/ # Data directory (external data goes here)
├── analysis/ # R notebooks and analysis code
├── docs/ # Documentation
├── results/ # Generated results (not in version control)
└── environment.yml # Conda environment specification
```
text

## 🚀 Quick Start

### Prerequisites
- Conda or Mamba
- Snakemake ≥7.0.0

### Installation
```bash
# Clone repository
git clone git@github.com:mtariqi/variant-calling-benchmark.git
cd variant-calling-benchmark

# Set up environment and download test data
bash scripts/setup_project.sh

# Run workflow
snakemake --cores 4 --use-conda
👥 Team Members
Member 1: Atra Alimoradian - Data management & alignment pipelines

Member 2: Raghad Al-Ampudi - Variant calling & workflow development

Member 3: Md Tariqul Islam - Benchmarking analysis & ML insights

🔧 Pipeline Components
Tools Evaluated
Aligners: BWA-MEM, Bowtie2 (local), Novoalign

Variant Callers: GATK HaplotypeCaller, DeepVariant, FreeBayes

ML Integration: DeepVariant (CNN), GATK CNN filtering

Benchmarking: hap.py, custom R analysis

Workflow
text
Raw Reads → Alignment → Variant Calling → Filtering → Benchmarking → Analysis
    ↓           ↓           ↓           ↓         ↓         ↓
  FASTQ       BAM         VCF        Filtered  hap.py   R/ML Analysis
📊 Expected Outcomes
Comparative performance ranking of 18 pipeline combinations

ML-powered variant calling vs traditional methods comparison

Reproducible Snakemake workflow for variant calling benchmarks

Comprehensive performance report with visualizations

📝 License
MIT License - see LICENSE file for details.# variant-calling-benchmark
Reproducible benchmark of variant calling pipelines using GIAB data - BINF6300 Group Project

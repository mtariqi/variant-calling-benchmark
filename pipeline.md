# Variant Calling Pipeline Workflow (Mermaid Version)

```mermaid
graph TD
    subgraph Input["📥 INPUT DATA"]
        A1[GIAB Reference Samples<br/>HG001-HG007]
        A2[Raw WGS/WES Reads<br/>14 samples]
        A3[Reference Genome<br/>GRCh37/GRCh38]
        A4[GIAB Gold Standards<br/>High-confidence variants]
    end

    subgraph Alignment["🧭 ALIGNMENT & PREPROCESSING"]
        B1[Quality Control<br/>FastQC, MultiQC]
        B2[Alignment Tools]
        B2A[BWA-MEM]
        B2B[Bowtie2 Local]
        B2C[Bowtie2 E2E]
        B2D[Novoalign]
        B2E[Isaac4]
        B3[Post-processing<br/>Sort, MarkDuplicates]
    end

    subgraph Calling["🔍 VARIANT CALLING & FILTERING"]
        C1[Variant Callers]
        C1A[DeepVariant]
        C1B[GATK HaplotypeCaller]
        C1C[FreeBayes]
        C1D[Sentieon]
        C1E[Octopus]
        C1F[Clair3]
        C2[Filtering Methods]
        C2A[GATK CNN 1D]
        C2B[GATK CNN 2D]
        C2C[GATK Hard Filter]
        C2D[Octopus RF]
        C2E[Octopus Std]
    end

    subgraph Benchmark["📊 BENCHMARKING & EVALUATION"]
        D1[Performance Metrics<br/>Precision, Recall, F1]
        D2[Benchmarking Strata<br/>SNPs vs INDELs]
        D3[Variant Comparison<br/>hap.py, vcfeval]
        D4[Stratified Analysis]
    end

    subgraph Output["📈 OUTPUT & ANALYSIS"]
        E1[Comparative Performance<br/>Pipeline Rankings]
        E2[Factor Impact Analysis]
        E3[Reproducible Workflow<br/>Snakemake Pipeline]
        E4[Final Report]
    end

    Input --> Alignment
    Alignment --> B2A & B2B & B2C & B2D & B2E
    B2A & B2B & B2C & B2D & B2E --> B3
    B3 --> Calling
    Calling --> C1A & C1B & C1C & C1D & C1E & C1F
    C1A & C1B & C1C & C1D & C1E & C1F --> C2A & C2B & C2C & C2D & C2E
    C2A & C2B & C2C & C2D & C2E --> Benchmark
    Benchmark --> D1 --> D2 --> D3 --> D4
    D4 --> Output
    Output --> E1 & E2 & E3 & E4

    classDef inputStyle fill:#89CFF0,stroke:#5DADE2,stroke-width:3px,color:#000
    classDef alignStyle fill:#DDA0DD,stroke:#BA68C8,stroke-width:3px,color:#000
    classDef callStyle fill:#98D8C8,stroke:#6BCB77,stroke-width:3px,color:#000
    classDef benchStyle fill:#FFD93D,stroke:#F4A460,stroke-width:3px,color:#000
    classDef outStyle fill:#FFB6C1,stroke:#FF69B4,stroke-width:3px,color:#000

    class Input inputStyle
    class Alignment alignStyle
    class Calling callStyle
    class Benchmark benchStyle
    class Output outStyle
```



```
┌─────────────────────────────────────────────────────────────┐
│                    📥 INPUT DATA                            │
├─────────────────────────────────────────────────────────────┤
│  • GIAB Samples (HG001-HG007)                              │
│  • Raw WGS/WES Reads (14 samples)                          │
│  • Reference Genome (GRCh37/GRCh38)                        │
│  • GIAB Gold Standards                                      │
└────────────────────────┬────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│           🧭 ALIGNMENT & PREPROCESSING                      │
├─────────────────────────────────────────────────────────────┤
│  Aligners: BWA-MEM | Bowtie2 (Local/E2E) | Novoalign      │
│  QC: FastQC, MultiQC                                       │
│  Post-processing: Sort, MarkDuplicates                     │
└────────────────────────┬────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│         🔍 VARIANT CALLING & FILTERING                      │
├─────────────────────────────────────────────────────────────┤
│  Callers: DeepVariant | GATK HC | FreeBayes | Octopus     │
│  Filters: CNN 1D/2D | Hard Filters | Random Forest        │
└────────────────────────┬────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│         📊 BENCHMARKING & EVALUATION                        │
├─────────────────────────────────────────────────────────────┤
│  Metrics: Precision, Recall, F1-Score                     │
│  Tools: hap.py, vcfeval                                    │
│  Strata: SNPs/INDELs, CDS, Coverage levels                │
└────────────────────────┬────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────┐
│             📈 OUTPUT & ANALYSIS                            │
├─────────────────────────────────────────────────────────────┤
│  • Comparative Performance Rankings                        │
│  • Factor Impact Analysis                                  │
│  • Reproducible Snakemake Workflow                         │
│  • Final Benchmarking Report                               │
└─────────────────────────────────────────────────────────────┘
```

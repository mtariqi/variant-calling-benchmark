---

# 🧬 **Project Report Summary — Job 1: Reference Preparation, FASTQ Retrieval, and Alignment (HG002)**

## **1. Overview**

Job 1 focused on preparing all foundational inputs required for the downstream variant calling workflow. This included organizing a reproducible shared working environment, downloading raw sequencing reads for the HG002 sample, retrieving and indexing the human reference genome, and performing short-read alignment using BWA-MEM. The final output of Job 1 is a high-quality, sorted, and indexed BAM file ready for variant calling.

This stage establishes the computational groundwork for the entire pipeline and ensures full reproducibility for all group members.

---

# **2. Shared Project Environment Setup**

A shared working directory was created on the Northeastern HPC under:

```
/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/
```

A setup script applied ACL permissions so all team members could:

* Read/write FASTQs
* Access reference files
* Submit SLURM jobs
* Work collaboratively without permission errors

Core subdirectories created:

```
raw_fastq/
reference/
results/
logs/
slurm_scripts/
```

This standardized layout ensures reproducibility and team-based workflow.

---

# **3. Downloading Raw FASTQ Files (HG002)**

Because NCBI SRA prefetch downloads were blocked by the cluster proxy, an HTTPS-based FASTQ download workflow was implemented using European Nucleotide Archive (ENA), which hosts mirrored versions of the same data.

The following paired-end reads were successfully downloaded:

```
HG002_R1.fastq.gz  
HG002_R2.fastq.gz
```

Located at:

```
/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/raw_fastq/HG002/
```

A dedicated SLURM script (`job1_fastq_download_https.slurm`) automated:

* Creating the output directory
* Downloading both FASTQs
* Logging all progress

The download completed successfully, and the files were verified with `ls -lh`.

---

# **4. Reference Genome Preparation**

The reference genome used is:

```
Homo_sapiens_assembly38.fasta  (GRCh38-no-alt)
```

Files prepared include:

* `.fasta`
* `.fasta.fai` (samtools index)
* `.fasta.amb/.ann/.bwt/.pac/.sa` (BWA index)
* `.dict` (Picard sequence dictionary)

These files are stored in:

```
/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/reference/
```

The indexing was successful with no errors, ensuring compatibility with BWA-MEM, GATK, and downstream tools.

---

# **5. Read Alignment With BWA-MEM and BAM Processing**

A dedicated SLURM script (`job1_alignment.slurm`) performed all alignment and BAM processing steps:

### **Steps executed:**

1. **BWA-MEM alignment**

```
bwa mem -t 8 REF.fa R1.fastq.gz R2.fastq.gz > HG002.raw.sam
```

2. **Convert SAM → BAM**

```
samtools view -@ 8 -bS HG002.raw.sam > HG002.raw.bam
```

3. **Sort BAM**

```
samtools sort -@ 8 HG002.raw.bam -o HG002.sorted.bam
```

4. **Index BAM**

```
samtools index HG002.sorted.bam
```

5. **Cleanup**

Temporary `.sam` files removed.

### **Outputs generated:**

Located under:

```
/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/results/bwa/
```

Final alignment files:

* `HG002.sorted.bam`
* `HG002.sorted.bam.bai`

These are the required inputs for Job 2 (variant calling).

---

# **6. Verification of Successful Execution**

The SLURM log (`job1_alignment.out`) reported:

```
==== JOB 1: ALIGNMENT STARTED ====
[INFO] Running BWA MEM...
[INFO] Converting SAM --> BAM...
[INFO] Sorting BAM...
[INFO] Indexing BAM...
[INFO] Cleaning up...
==== JOB 1: ALIGNMENT COMPLETE ====
```

No `.err` file was generated, confirming:

* No errors occurred
* All modules loaded correctly
* All alignment and processing steps completed successfully

---

# **7. Status Summary (as of now)**

| Stage                     | Status      | Notes                              |
| ------------------------- | ----------- | ---------------------------------- |
| Shared directory setup    | ✔ Completed | ACL permissions verified           |
| Reference genome download | ✔ Completed | All indexes (FAI, DICT, BWA) ready |
| FASTQ download (HTTPS)    | ✔ Completed | R1/R2 verified                     |
| Read alignment            | ✔ Completed | Sorted + indexed BAM generated     |
| Log files                 | ✔ Completed | Located under `/logs/`             |

**Job 1 is fully completed and validated.**

---

# **8. What Comes Next (For Job 2)**

Job 2 will perform:

* Variant calling (GATK HaplotypeCaller OR DeepVariant)
* Generate raw VCF or gVCF
* Store outputs in `/results/`

Your teammate responsible for Job 2 can now begin safely, using:

```
HG002.sorted.bam
Homo_sapiens_assembly38.fasta
```

Both are ready.

---


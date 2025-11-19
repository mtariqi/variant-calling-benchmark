```
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Variant Calling Pipeline Workflow</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 40px 20px;
            min-height: 100vh;
        }

        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 20px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            overflow: hidden;
        }

        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }

        .header h1 {
            font-size: 32px;
            font-weight: 700;
            margin-bottom: 10px;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.2);
        }

        .header p {
            font-size: 16px;
            opacity: 0.95;
        }

        .workflow {
            padding: 40px;
        }

        .stage {
            margin-bottom: 30px;
            position: relative;
        }

        .stage-title {
            font-size: 20px;
            font-weight: 700;
            margin-bottom: 15px;
            padding: 12px 20px;
            border-radius: 10px;
            text-align: center;
            color: white;
            text-transform: uppercase;
            letter-spacing: 1px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
        }

        .stage-subtitle {
            font-size: 14px;
            font-weight: 600;
            margin-bottom: 12px;
            padding: 8px 15px;
            border-radius: 8px;
            text-align: center;
            background: #f8f9fa;
            color: #495057;
            border-left: 4px solid;
        }

        /* Stage-specific colors */
        .input-stage .stage-title { background: linear-gradient(135deg, #89CFF0, #5DADE2); }
        .alignment-stage .stage-title { background: linear-gradient(135deg, #DDA0DD, #BA68C8); }
        .calling-stage .stage-title { background: linear-gradient(135deg, #98D8C8, #6BCB77); }
        .benchmark-stage .stage-title { background: linear-gradient(135deg, #FFD93D, #F4A460); }
        .output-stage .stage-title { background: linear-gradient(135deg, #FFB6C1, #FF69B4); }

        .input-stage .stage-subtitle { border-left-color: #5DADE2; }
        .alignment-stage .stage-subtitle { border-left-color: #BA68C8; }
        .calling-stage .stage-subtitle { border-left-color: #6BCB77; }
        .benchmark-stage .stage-subtitle { border-left-color: #F4A460; }

        .tools-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 15px;
            margin-bottom: 20px;
        }

        .tool-box {
            background: linear-gradient(135deg, #ffffff 0%, #f8f9fa 100%);
            border: 2px solid #e9ecef;
            border-radius: 10px;
            padding: 15px;
            text-align: center;
            transition: all 0.3s ease;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08);
        }

        .tool-box:hover {
            transform: translateY(-5px);
            box-shadow: 0 8px 20px rgba(0,0,0,0.15);
            border-color: #667eea;
        }

        .tool-name {
            font-weight: 700;
            font-size: 15px;
            color: #2c3e50;
            margin-bottom: 5px;
        }

        .tool-desc {
            font-size: 12px;
            color: #6c757d;
            line-height: 1.4;
        }

        .data-box {
            background: linear-gradient(135deg, #E3F2FD, #BBDEFB);
            border: 2px solid #90CAF9;
            border-radius: 12px;
            padding: 20px;
            margin-bottom: 15px;
            text-align: center;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        }

        .data-title {
            font-weight: 700;
            font-size: 16px;
            color: #1565C0;
            margin-bottom: 8px;
        }

        .data-desc {
            font-size: 13px;
            color: #455A64;
            line-height: 1.5;
        }

        .arrow {
            text-align: center;
            margin: 20px 0;
            position: relative;
        }

        .arrow::before {
            content: '▼';
            font-size: 30px;
            color: #667eea;
            display: block;
            animation: bounce 2s infinite;
        }

        @keyframes bounce {
            0%, 100% { transform: translateY(0); }
            50% { transform: translateY(10px); }
        }

        .process-box {
            background: #f8f9fa;
            border-left: 4px solid #667eea;
            border-radius: 8px;
            padding: 15px 20px;
            margin-bottom: 15px;
        }

        .process-title {
            font-weight: 700;
            color: #2c3e50;
            font-size: 15px;
            margin-bottom: 5px;
        }

        .process-desc {
            font-size: 13px;
            color: #6c757d;
        }

        .metrics-grid {
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 15px;
            margin-bottom: 20px;
        }

        .metric-box {
            background: linear-gradient(135deg, #FFF9E6, #FFF3CD);
            border: 2px solid #FFE69C;
            border-radius: 10px;
            padding: 15px;
            text-align: center;
        }

        .metric-title {
            font-weight: 700;
            font-size: 14px;
            color: #856404;
            margin-bottom: 5px;
        }

        .metric-desc {
            font-size: 12px;
            color: #856404;
        }

        .analysis-grid {
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 15px;
        }

        .analysis-box {
            background: linear-gradient(135deg, #FFE4E8, #FFC1CC);
            border: 2px solid #FFB3C1;
            border-radius: 10px;
            padding: 20px;
            text-align: center;
        }

        .analysis-title {
            font-weight: 700;
            font-size: 15px;
            color: #C41E3A;
            margin-bottom: 8px;
        }

        .analysis-desc {
            font-size: 13px;
            color: #8B0000;
            line-height: 1.4;
        }

        .footer {
            background: #2c3e50;
            color: white;
            padding: 20px;
            text-align: center;
            font-size: 13px;
        }

        .badge {
            display: inline-block;
            padding: 4px 10px;
            border-radius: 12px;
            font-size: 11px;
            font-weight: 600;
            margin: 3px;
            background: #667eea;
            color: white;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Variant Calling Pipeline Workflow</h1>
            <p>Comprehensive Benchmarking of State-of-the-Art Variant Calling Methods</p>
            <div style="margin-top: 15px;">
                <span class="badge">GIAB Benchmark</span>
                <span class="badge">18 Pipeline Combinations</span>
                <span class="badge">ML-Enhanced</span>
                <span class="badge">Reproducible</span>
            </div>
        </div>

        <div class="workflow">
            <!-- INPUT DATA STAGE -->
            <div class="stage input-stage">
                <div class="stage-title">📥 Input Data</div>
                
                <div class="data-box">
                    <div class="data-title">GIAB Reference Samples</div>
                    <div class="data-desc">HG001-HG007</div>
                </div>

                <div class="data-box">
                    <div class="data-title">Raw WGS/WES Reads</div>
                    <div class="data-desc">14 samples</div>
                </div>

                <div class="data-box">
                    <div class="data-title">Reference Genome</div>
                    <div class="data-desc">GRCh37/GRCh38</div>
                </div>

                <div class="data-box">
                    <div class="data-title">GIAB Gold Standards</div>
                    <div class="data-desc">High-confidence variants & regions</div>
                </div>
            </div>

            <div class="arrow"></div>

            <!-- ALIGNMENT STAGE -->
            <div class="stage alignment-stage">
                <div class="stage-title">🧭 Alignment & Preprocessing</div>
                
                <div class="stage-subtitle">Alignment Tools</div>
                <div class="tools-grid">
                    <div class="tool-box">
                        <div class="tool-name">BWA-MEM</div>
                        <div class="tool-desc">BW</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">Bowtie2 Local</div>
                        <div class="tool-desc">BT-LOC</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">Bowtie2 End-to-end</div>
                        <div class="tool-desc">BT-E2E</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">Novoalign</div>
                        <div class="tool-desc">NO</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">Isaac4</div>
                        <div class="tool-desc">IS</div>
                    </div>
                </div>

                <div class="process-box">
                    <div class="process-title">Quality Control</div>
                    <div class="process-desc">FastQC, MultiQC</div>
                </div>

                <div class="process-box">
                    <div class="process-title">Read Alignment</div>
                    <div class="process-desc">Map reads to reference genome</div>
                </div>

                <div class="process-box">
                    <div class="process-title">Post-processing</div>
                    <div class="process-desc">Sorting, MarkDuplicates</div>
                </div>
            </div>

            <div class="arrow"></div>

            <!-- VARIANT CALLING STAGE -->
            <div class="stage calling-stage">
                <div class="stage-title">🔍 Variant Calling & Filtering</div>
                
                <div class="stage-subtitle">Filtering Methods</div>
                <div class="tools-grid">
                    <div class="tool-box">
                        <div class="tool-name">GATK CNN 1D</div>
                        <div class="tool-desc">G1</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">GATK CNN 2D</div>
                        <div class="tool-desc">G2</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">GATK Hard Filter</div>
                        <div class="tool-desc">GH</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">Octopus Random Forest</div>
                        <div class="tool-desc">OF</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">Octopus Standard</div>
                        <div class="tool-desc">OS</div>
                    </div>
                </div>

                <div class="stage-subtitle">Calling Tools</div>
                <div class="tools-grid">
                    <div class="tool-box">
                        <div class="tool-name">DeepVariant</div>
                        <div class="tool-desc">DV</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">GATK HaplotypeCaller</div>
                        <div class="tool-desc">GATK</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">FreeBayes</div>
                        <div class="tool-desc">FB</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">Sentieon</div>
                        <div class="tool-desc">ST</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">Octopus</div>
                        <div class="tool-desc">OF/OS</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">Clair3</div>
                        <div class="tool-desc">CL</div>
                    </div>
                </div>

                <div class="process-box">
                    <div class="process-title">Variant Calling</div>
                    <div class="process-desc">Identify SNPs and Indels from aligned reads</div>
                </div>

                <div class="process-box">
                    <div class="process-title">Variant Filtering</div>
                    <div class="process-desc">Apply quality filters and ML-based classification</div>
                </div>
            </div>

            <div class="arrow"></div>

            <!-- BENCHMARKING STAGE -->
            <div class="stage benchmark-stage">
                <div class="stage-title">📊 Benchmarking & Evaluation</div>
                
                <div class="stage-subtitle">Performance Metrics</div>
                <div class="metrics-grid">
                    <div class="metric-box">
                        <div class="metric-title">Precision</div>
                        <div class="metric-desc">TP/(TP+FP)</div>
                    </div>
                    <div class="metric-box">
                        <div class="metric-title">Recall</div>
                        <div class="metric-desc">TP/(TP+FN)</div>
                    </div>
                    <div class="metric-box">
                        <div class="metric-title">F1-Score</div>
                        <div class="metric-desc">2*Precision*Recall/(Precision+Recall)</div>
                    </div>
                </div>

                <div class="stage-subtitle">Benchmarking Strata</div>
                <div class="tools-grid">
                    <div class="tool-box">
                        <div class="tool-name">Variant Type</div>
                        <div class="tool-desc">SNPs vs INDELs</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">Genomic Context</div>
                        <div class="tool-desc">CDS, promoters</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">Functional Impact</div>
                        <div class="tool-desc">synonymous, missense</div>
                    </div>
                    <div class="tool-box">
                        <div class="tool-name">Coverage Levels</div>
                        <div class="tool-desc">low vs high coverage</div>
                    </div>
                </div>

                <div class="process-box">
                    <div class="process-title">Variant Comparison</div>
                    <div class="process-desc">hap.py, vcfeval</div>
                </div>

                <div class="process-box">
                    <div class="process-title">Stratified Analysis</div>
                    <div class="process-desc">Performance across different genomic contexts</div>
                </div>

                <div class="process-box">
                    <div class="process-title">Metric Calculation</div>
                    <div class="process-desc">Compute precision, recall, F1-score for each pipeline</div>
                </div>
            </div>

            <div class="arrow"></div>

            <!-- OUTPUT & ANALYSIS STAGE -->
            <div class="stage output-stage">
                <div class="stage-title">📈 Output & Analysis</div>
                
                <div class="analysis-grid">
                    <div class="analysis-box">
                        <div class="analysis-title">Comparative Performance</div>
                        <div class="analysis-desc">Pipeline rankings</div>
                    </div>
                    <div class="analysis-box">
                        <div class="analysis-title">Factor Impact Analysis</div>
                        <div class="analysis-desc">Key accuracy drivers</div>
                    </div>
                    <div class="analysis-box">
                        <div class="analysis-title">Reproducible Workflow</div>
                        <div class="analysis-desc">Snakemake pipeline</div>
                    </div>
                    <div class="analysis-box">
                        <div class="analysis-title">Final Report</div>
                        <div class="analysis-desc">Benchmarking conclusions</div>
                    </div>
                </div>
            </div>
        </div>

        <div class="footer">
            <p><strong>Variant Calling Pipeline Benchmark</strong> | BINF6300 Group Project</p>
            <p style="margin-top: 10px; opacity: 0.8;">Reproducing Barbitoff et al. (2022) - Systematic benchmark of variant calling pipelines</p>
        </div>
    </div>
</body>
</html>
```

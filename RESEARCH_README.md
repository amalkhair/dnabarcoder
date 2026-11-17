# Fungal Biomarker Research

## Computational Analysis and Validation

**Researcher:** Amal Khairunnisa
**Institution:** Fontys University of Applied Sciences
**Research collaboration:** Westerdijk Fungal Biodiversity Institute
**Research branch:** `publication-scale-analysis`

## 1. Research Overview

This project investigates alternative DNA biomarkers for fungal species where the Internal Transcribed Spacer (ITS) region provides insufficient species-level resolution.

**Main research question:**

Which alternative biomarkers are promising for fungal species where ITS does not provide sufficient species-level resolution?

The research evaluates marker performance using DNABarcoder, supported by taxonomic coverage, species ambiguity, species-matched comparisons, phylogenetic analysis, and barcode-gap calculations.

The computational analyses cover 12 genetic loci, with more detailed evaluation of selected fungal genera.

## 2. Software Contributions

The research builds on the original DNABarcoder software developed by Duong Vu and collaborators.

The additional research scripts support:

- Preparation of standardized FASTA and classification files
- Automated marker characterization and result processing
- Genus-level biomarker performance comparisons
- Identification and quantification of species ambiguity
- Preparation of species-matched and balanced datasets
- Phylogenetic monophyly evaluation
- Genetic-distance and barcode-gap calculations
- Computational quality and consistency checks

The research scripts support the analysis; they do not constitute a new DNA barcoding algorithm or a replacement for DNABarcoder.

## 3. Implementation Overview

### Data Preparation

`scripts/prepare_dnabarcoder_inputs.py`

Generates standardized FASTA and corresponding classification files from previously processed marker datasets.

### Automated Analysis

`scripts/run_marker_characterization.py`

Automates marker characterization, checks required input files, executes DNABarcoder clustering, and processes the resulting outputs.

`run_all_markers.sh`

Runs DNABarcoder similarity-threshold prediction across the selected genetic markers.

### Species Ambiguity

`scripts/plot_candidate_species_ambiguity.py`

Analyzes species involved in 100%-identity multi-species groups and calculates species ambiguity percentages.

### Species-Matched Datasets

`scripts/build_fusarium_its_cam_matched.py`

Generates balanced Fusarium ITS and CaM datasets using shared species and reproducible sequence sampling.

`scripts/prepare_phylo_files_for_duong.py`

Prepares genus-specific ITS and alternative-marker FASTA datasets for further phylogenetic analysis.

### Phylogenetic Validation

`scripts/analyze_matched_tree.py`

Evaluates species monophyly using phylogenetic tree bipartitions.

### Barcode-Gap Analysis

`scripts/analyze_fusarium_barcode_gap.py`

Calculates maximum intraspecific and nearest interspecific genetic distances to evaluate species-level barcode gaps.

## 4. Research Outputs

The analyses generated:

- Biomarker performance and coverage summaries
- Genus-level F-measure comparisons
- Species ambiguity results
- Matched FASTA datasets
- Preliminary phylogenetic validation results
- Fusarium ITS-CaM barcode-gap results
- Scientific figures and summary tables

Selected lightweight outputs are maintained under:

- `results/final_summary/`
- `results/figures/`

Large and intermediate datasets may be excluded from version control.

## 5. Computational Workflow

The research workflow consists of:

1. Preparing marker-specific sequence inputs
2. Running DNABarcoder analyses
3. Processing genus-level performance and ambiguity results
4. Selecting candidate markers and shared species
5. Generating matched datasets
6. Performing phylogenetic and barcode-gap analyses
7. Reviewing results and computational quality checks

Processing decisions and analysis parameters are documented in the corresponding scripts.

## 6. Reproducibility and Limitations

Selected workflows use fixed random seeds, input-file checks, sequence-count validation, and structured output files.

The analyses depend on the availability of appropriate sequence datasets and external scientific software, including DNABarcoder, BLAST, MAFFT, and IQ-TREE.

Not every workflow has been packaged as a single end-to-end executable pipeline. Some analyses require intermediate datasets or separately generated outputs.

Phylogenetic results and biomarker recommendations remain subject to further scientific review.

## 7. Research Documentation

The accompanying Research Progress Report documents:

- Research methodology
- Computational results
- Selected implementation excerpts
- Methodological decisions
- Current limitations and next steps

The scientific figures included in the computational report represent analysis outputs. Research concerning the selection, design, and evaluation of visualization techniques is documented separately.

## 8. Attribution

This repository builds on DNABarcoder developed by Duong Vu and collaborators.

Original software:
https://github.com/vuthuyduong/dnabarcoder

Research branch:
https://github.com/amalkhair/dnabarcoder/tree/publication-scale-analysis

Reference:

Vu, D., Nilsson, R. H., & Verkley, G. J. M. (2022).
Dnabarcoder: An open-source software package for analysing
and predicting DNA sequence similarity cutoffs for fungal
sequence identification. Molecular Ecology Resources.
https://doi.org/10.1111/1755-0998.13651

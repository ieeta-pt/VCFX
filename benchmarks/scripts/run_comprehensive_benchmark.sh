#!/bin/bash
# Comprehensive VCFX Performance Benchmark
# Tests all tools on chr21 1000 Genomes Phase 3 data (4GB, 427K variants, 2504 samples)

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="$ROOT_DIR/build/src"
DATA_DIR="$ROOT_DIR/benchmarks/data"
RESULTS_DIR="$ROOT_DIR/benchmarks/results"
DATA_FILE="$DATA_DIR/chr21.1kg.phase3.v5a.vcf"
REF_FASTA="$DATA_DIR/chr21.fa"
POP_FREQ_FILE="$DATA_DIR/pop_frequencies.tsv"

# Create results directory
mkdir -p "$RESULTS_DIR"

# =============================================================================
# Generate helper files for tools that need additional data
# =============================================================================
echo "=== Preparing helper files ==="

# Download chr21 reference FASTA if not present
if [ ! -f "$REF_FASTA" ]; then
    echo "Downloading chr21 reference FASTA..."
    curl -s "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr21.fa.gz" | gunzip > "$REF_FASTA" 2>/dev/null || {
        echo "Warning: Could not download reference FASTA, creating minimal version"
        echo ">21" > "$REF_FASTA"
        # Generate minimal reference from VCF positions (for benchmark purposes)
        awk -F'\t' '!/^#/ {print $2, $4}' "$DATA_FILE" | head -1000 | while read pos ref; do
            printf "%s" "$ref"
        done >> "$REF_FASTA"
        echo "" >> "$REF_FASTA"
    }
fi

# Generate population frequencies file from the VCF (extract AF from first 10K variants)
if [ ! -f "$POP_FREQ_FILE" ]; then
    echo "Generating population frequencies file..."
    awk -F'\t' 'BEGIN {OFS="\t"}
    !/^#/ && NR <= 10000 {
        # Parse AF from INFO field
        split($8, info, ";")
        af = "0.5"
        for (i in info) {
            if (info[i] ~ /^AF=/) {
                split(info[i], afval, "=")
                af = afval[2]
                gsub(/,.*/, "", af)  # Take first AF value for multiallelic
            }
        }
        # Assign population based on AF range (simulated)
        if (af < 0.1) pop = "EUR"
        else if (af < 0.2) pop = "AFR"
        else if (af < 0.3) pop = "EAS"
        else if (af < 0.4) pop = "SAS"
        else pop = "AMR"
        print $1, $2, $4, $5, pop, af
    }' "$DATA_FILE" > "$POP_FREQ_FILE"
fi
echo ""

# Results file
RESULTS_CSV="$RESULTS_DIR/comprehensive_benchmark_$(date +%Y%m%d_%H%M%S).csv"
echo "tool,category,time_seconds,status,description" > "$RESULTS_CSV"

echo "=============================================="
echo "   VCFX Comprehensive Performance Benchmark"
echo "=============================================="
echo ""
echo "Test file: $DATA_FILE"
echo "Size: $(ls -lh "$DATA_FILE" | awk '{print $5}')"
echo "Date: $(date)"
echo "Results: $RESULTS_CSV"
echo ""

# Function to benchmark a tool (default 5 min timeout)
benchmark() {
    local name="$1"
    local category="$2"
    local cmd="$3"
    local desc="$4"
    local timeout_secs="${5:-300}"  # Default 5 minutes

    echo -n "Testing $name... "

    # Run with timeout
    local start=$(python3 -c "import time; print(time.time())")
    if timeout $timeout_secs bash -c "$cmd" > /dev/null 2>&1; then
        local end=$(python3 -c "import time; print(time.time())")
        local elapsed=$(python3 -c "print(f'{$end - $start:.3f}')")
        echo "${elapsed}s"
        echo "$name,$category,$elapsed,success,$desc" >> "$RESULTS_CSV"
    else
        echo "FAILED/TIMEOUT"
        echo "$name,$category,0,failed,$desc" >> "$RESULTS_CSV"
    fi
}

# Warm-up run
echo "=== Warm-up (priming file cache) ==="
cat "$DATA_FILE" > /dev/null
echo ""

echo "=== Category 1: Basic I/O & Validation ==="
benchmark "variant_counter" "basic_io" \
    "$BUILD_DIR/VCFX_variant_counter/VCFX_variant_counter < $DATA_FILE" \
    "Count variants"

benchmark "validator_stdin" "basic_io" \
    "$BUILD_DIR/VCFX_validator/VCFX_validator < $DATA_FILE" \
    "Validate VCF (stdin)"

benchmark "validator_mmap" "basic_io" \
    "$BUILD_DIR/VCFX_validator/VCFX_validator -i $DATA_FILE" \
    "Validate VCF (mmap)"

benchmark "header_parser" "basic_io" \
    "$BUILD_DIR/VCFX_header_parser/VCFX_header_parser < $DATA_FILE" \
    "Parse header"

benchmark "reformatter" "basic_io" \
    "$BUILD_DIR/VCFX_reformatter/VCFX_reformatter < $DATA_FILE" \
    "Reformat VCF"

echo ""
echo "=== Category 2: Filtering Tools ==="
benchmark "phred_filter" "filtering" \
    "$BUILD_DIR/VCFX_phred_filter/VCFX_phred_filter -p 30 < $DATA_FILE" \
    "Filter QUAL>=30"

benchmark "record_filter" "filtering" \
    "$BUILD_DIR/VCFX_record_filter/VCFX_record_filter --filter 'FILTER==PASS' < $DATA_FILE" \
    "Filter PASS only"

benchmark "nonref_filter" "filtering" \
    "$BUILD_DIR/VCFX_nonref_filter/VCFX_nonref_filter < $DATA_FILE" \
    "Non-reference filter"

benchmark "gl_filter" "filtering" \
    "$BUILD_DIR/VCFX_gl_filter/VCFX_gl_filter --filter 'GQ>20' < $DATA_FILE" \
    "GQ filter (GQ>20)"

benchmark "allele_balance_filter" "filtering" \
    "$BUILD_DIR/VCFX_allele_balance_filter/VCFX_allele_balance_filter --min-ab 0.2 --max-ab 0.8 < $DATA_FILE" \
    "Allele balance filter"

benchmark "probability_filter" "filtering" \
    "$BUILD_DIR/VCFX_probability_filter/VCFX_probability_filter --min-gp 0.9 < $DATA_FILE" \
    "GP filter"

benchmark "phase_quality_filter" "filtering" \
    "$BUILD_DIR/VCFX_phase_quality_filter/VCFX_phase_quality_filter --min-pq 30 < $DATA_FILE" \
    "Phase quality filter"

benchmark "impact_filter" "filtering" \
    "$BUILD_DIR/VCFX_impact_filter/VCFX_impact_filter --filter-impact HIGH < $DATA_FILE" \
    "Impact filter (HIGH)"

echo ""
echo "=== Category 3: Analysis & Calculation Tools ==="
benchmark "allele_freq_calc" "analysis" \
    "$BUILD_DIR/VCFX_allele_freq_calc/VCFX_allele_freq_calc < $DATA_FILE" \
    "Allele frequencies"

benchmark "variant_classifier" "analysis" \
    "$BUILD_DIR/VCFX_variant_classifier/VCFX_variant_classifier < $DATA_FILE" \
    "Classify variants"

# NOTE: These O(variants × samples) tools need longer timeout (60 min = 3600s)
benchmark "hwe_tester" "analysis" \
    "$BUILD_DIR/VCFX_hwe_tester/VCFX_hwe_tester < $DATA_FILE" \
    "HWE test" 3600

benchmark "allele_counter" "analysis" \
    "$BUILD_DIR/VCFX_allele_counter/VCFX_allele_counter < $DATA_FILE" \
    "Count alleles" 3600

benchmark "allele_balance_calc" "analysis" \
    "$BUILD_DIR/VCFX_allele_balance_calc/VCFX_allele_balance_calc < $DATA_FILE" \
    "Allele balance calc" 3600

benchmark "inbreeding_calculator" "analysis" \
    "$BUILD_DIR/VCFX_inbreeding_calculator/VCFX_inbreeding_calculator < $DATA_FILE" \
    "Inbreeding coef" 3600

benchmark "dosage_calculator" "analysis" \
    "$BUILD_DIR/VCFX_dosage_calculator/VCFX_dosage_calculator < $DATA_FILE" \
    "Dosage calc" 3600

benchmark "distance_calculator" "analysis" \
    "$BUILD_DIR/VCFX_distance_calculator/VCFX_distance_calculator < $DATA_FILE" \
    "Genetic distance"

echo ""
echo "=== Category 4: Quality Control Tools ==="
benchmark "missing_detector" "qc" \
    "$BUILD_DIR/VCFX_missing_detector/VCFX_missing_detector < $DATA_FILE" \
    "Detect missing"

# Phase checker may take time on large files
benchmark "phase_checker" "qc" \
    "$BUILD_DIR/VCFX_phase_checker/VCFX_phase_checker < $DATA_FILE" \
    "Check phasing" 600

benchmark "outlier_detector" "qc" \
    "$BUILD_DIR/VCFX_outlier_detector/VCFX_outlier_detector < $DATA_FILE" \
    "Detect outliers"

# Cross-sample concordance with specific samples
benchmark "cross_sample_concordance" "qc" \
    "$BUILD_DIR/VCFX_cross_sample_concordance/VCFX_cross_sample_concordance --samples HG00096,HG00097 < $DATA_FILE" \
    "Cross-sample concordance"

# Alignment checker needs reference FASTA
benchmark "alignment_checker" "qc" \
    "$BUILD_DIR/VCFX_alignment_checker/VCFX_alignment_checker --alignment-discrepancy $DATA_FILE $REF_FASTA" \
    "Alignment check"

echo ""
echo "=== Category 5: Transformation Tools ==="
benchmark "multiallelic_splitter" "transformation" \
    "$BUILD_DIR/VCFX_multiallelic_splitter/VCFX_multiallelic_splitter < $DATA_FILE" \
    "Split multiallelic"

# Indel normalizer may need extra time
benchmark "indel_normalizer" "transformation" \
    "$BUILD_DIR/VCFX_indel_normalizer/VCFX_indel_normalizer < $DATA_FILE" \
    "Normalize indels" 900

benchmark "duplicate_remover" "transformation" \
    "$BUILD_DIR/VCFX_duplicate_remover/VCFX_duplicate_remover < $DATA_FILE" \
    "Remove duplicates"

# Sorter needs longer timeout for large files
benchmark "sorter" "transformation" \
    "$BUILD_DIR/VCFX_sorter/VCFX_sorter < $DATA_FILE" \
    "Sort variants" 1800

benchmark "quality_adjuster" "transformation" \
    "$BUILD_DIR/VCFX_quality_adjuster/VCFX_quality_adjuster < $DATA_FILE" \
    "Adjust quality"

# Missing data handler with fill option
benchmark "missing_data_handler" "transformation" \
    "$BUILD_DIR/VCFX_missing_data_handler/VCFX_missing_data_handler --fill-missing < $DATA_FILE" \
    "Handle missing"

benchmark "sv_handler" "transformation" \
    "$BUILD_DIR/VCFX_sv_handler/VCFX_sv_handler < $DATA_FILE" \
    "Handle SVs"

echo ""
echo "=== Category 6: Extraction & Subsetting Tools ==="
benchmark "sample_extractor" "extraction" \
    "$BUILD_DIR/VCFX_sample_extractor/VCFX_sample_extractor --samples HG00096,HG00097,HG00099 < $DATA_FILE" \
    "Extract samples"

benchmark "field_extractor" "extraction" \
    "$BUILD_DIR/VCFX_field_extractor/VCFX_field_extractor --fields CHROM,POS,REF,ALT < $DATA_FILE" \
    "Extract fields"

benchmark "genotype_query" "extraction" \
    "$BUILD_DIR/VCFX_genotype_query/VCFX_genotype_query --genotype-query '0/1' < $DATA_FILE" \
    "Query genotypes"

benchmark "position_subsetter" "extraction" \
    "$BUILD_DIR/VCFX_position_subsetter/VCFX_position_subsetter --region 21:1-10000000 < $DATA_FILE" \
    "Subset by region"

benchmark "af_subsetter" "extraction" \
    "$BUILD_DIR/VCFX_af_subsetter/VCFX_af_subsetter --af-filter 0.01-0.99 < $DATA_FILE" \
    "Subset by AF"

benchmark "subsampler" "extraction" \
    "$BUILD_DIR/VCFX_subsampler/VCFX_subsampler --fraction 0.1 < $DATA_FILE" \
    "Random subsample"

echo ""
echo "=== Category 7: Annotation & INFO Tools ==="
benchmark "info_parser" "annotation" \
    "$BUILD_DIR/VCFX_info_parser/VCFX_info_parser --info 'AC,AN,AF' < $DATA_FILE" \
    "Parse INFO"

benchmark "info_summarizer" "annotation" \
    "$BUILD_DIR/VCFX_info_summarizer/VCFX_info_summarizer --info 'AC,AN,AF' < $DATA_FILE" \
    "Summarize INFO"

# Info aggregator needs --aggregate-info argument
benchmark "info_aggregator" "annotation" \
    "$BUILD_DIR/VCFX_info_aggregator/VCFX_info_aggregator --aggregate-info 'AC,AN,AF' < $DATA_FILE" \
    "Aggregate INFO"

benchmark "metadata_summarizer" "annotation" \
    "$BUILD_DIR/VCFX_metadata_summarizer/VCFX_metadata_summarizer < $DATA_FILE" \
    "Summarize metadata"

benchmark "annotation_extractor" "annotation" \
    "$BUILD_DIR/VCFX_annotation_extractor/VCFX_annotation_extractor --annotation-extract 'AF' < $DATA_FILE" \
    "Extract annotations"

echo ""
echo "=== Category 8: Haplotype Tools ==="
benchmark "haplotype_phaser" "haplotype" \
    "$BUILD_DIR/VCFX_haplotype_phaser/VCFX_haplotype_phaser --streaming < $DATA_FILE" \
    "Phase haplotypes (streaming)"

benchmark "haplotype_extractor" "haplotype" \
    "$BUILD_DIR/VCFX_haplotype_extractor/VCFX_haplotype_extractor --streaming < $DATA_FILE" \
    "Extract haplotypes (streaming)"

echo ""
echo "=== Category 9: Ancestry Tools ==="
# Ancestry inferrer needs population frequencies file
benchmark "ancestry_inferrer" "ancestry" \
    "$BUILD_DIR/VCFX_ancestry_inferrer/VCFX_ancestry_inferrer --frequency $POP_FREQ_FILE < $DATA_FILE" \
    "Infer ancestry"

# Ancestry assigner needs population frequencies file
benchmark "ancestry_assigner" "ancestry" \
    "$BUILD_DIR/VCFX_ancestry_assigner/VCFX_ancestry_assigner --pop-freqs $DATA_DIR/pop_frequencies.txt < $DATA_FILE" \
    "Assign ancestry"

echo ""
echo "=== Category 9b: Population Tools ==="
# Population filter needs a population map file
POP_MAP_FILE="$DATA_DIR/pop_map.tsv"
if [ ! -f "$POP_MAP_FILE" ]; then
    echo "Generating population map file from VCF samples..."
    # Extract sample names and assign populations based on sample ID patterns
    head -1000 "$DATA_FILE" | grep "^#CHROM" | cut -f10- | tr '\t' '\n' | awk '{
        # Assign populations based on sample patterns (simulated for 1000G samples)
        n = NR % 5
        if (n == 0) pop = "EUR"
        else if (n == 1) pop = "AFR"
        else if (n == 2) pop = "EAS"
        else if (n == 3) pop = "SAS"
        else pop = "AMR"
        print $0 "\t" pop
    }' > "$POP_MAP_FILE"
fi
benchmark "population_filter" "population" \
    "$BUILD_DIR/VCFX_population_filter/VCFX_population_filter --population EUR --pop-map $POP_MAP_FILE < $DATA_FILE" \
    "Filter by population"

echo ""
echo "=== Category 10: File Management Tools ==="
benchmark "indexer" "file_management" \
    "$BUILD_DIR/VCFX_indexer/VCFX_indexer $DATA_FILE" \
    "Create index (mmap)"

benchmark "indexer_stdin" "file_management" \
    "$BUILD_DIR/VCFX_indexer/VCFX_indexer < $DATA_FILE" \
    "Create index (stdin)"

# Compressor needs --compress flag
benchmark "compressor" "file_management" \
    "$BUILD_DIR/VCFX_compressor/VCFX_compressor --compress < $DATA_FILE" \
    "Compress VCF"

# File splitter
benchmark "file_splitter" "file_management" \
    "$BUILD_DIR/VCFX_file_splitter/VCFX_file_splitter --chunks 10 < $DATA_FILE" \
    "Split into chunks"

# Merger (merge with itself as a test)
benchmark "merger" "file_management" \
    "$BUILD_DIR/VCFX_merger/VCFX_merger $DATA_FILE $DATA_FILE" \
    "Merge VCF files"

# Region subsampler - needs BED file with regions
REGION_BED="$DATA_DIR/regions.bed"
if [ ! -f "$REGION_BED" ]; then
    echo "Creating regions BED file..."
    echo -e "21\t0\t10000000\n21\t20000000\t30000000" > "$REGION_BED"
fi
benchmark "region_subsampler" "file_management" \
    "$BUILD_DIR/VCFX_region_subsampler/VCFX_region_subsampler --region-bed $REGION_BED < $DATA_FILE" \
    "Subsample by region"

echo ""
echo "=== Category 11: Format Conversion ==="
benchmark "format_converter_bed" "conversion" \
    "$BUILD_DIR/VCFX_format_converter/VCFX_format_converter --to-bed < $DATA_FILE" \
    "Convert to BED"

# FASTA converter is O(variants × samples) - needs longer timeout
benchmark "fasta_converter" "conversion" \
    "$BUILD_DIR/VCFX_fasta_converter/VCFX_fasta_converter < $DATA_FILE" \
    "Convert to FASTA" 3600

echo ""
echo "=== Category 12: Comparison & Diff Tools ==="
# Diff tool needs two VCF files with --file1 and --file2 flags
benchmark "diff_tool" "comparison" \
    "$BUILD_DIR/VCFX_diff_tool/VCFX_diff_tool --file1 $DATA_FILE --file2 $DATA_FILE --assume-sorted" \
    "Diff VCF files"

# Concordance checker - reads from stdin, compares two samples within the VCF
# Extract first two sample names from VCF header
SAMPLE1=$(head -1000 "$DATA_FILE" | grep "^#CHROM" | cut -f10)
SAMPLE2=$(head -1000 "$DATA_FILE" | grep "^#CHROM" | cut -f11)
benchmark "concordance_checker" "comparison" \
    "$BUILD_DIR/VCFX_concordance_checker/VCFX_concordance_checker --samples \"$SAMPLE1 $SAMPLE2\" < $DATA_FILE" \
    "Check concordance"

# Reference comparator needs reference FASTA
benchmark "ref_comparator" "comparison" \
    "$BUILD_DIR/VCFX_ref_comparator/VCFX_ref_comparator --reference $REF_FASTA < $DATA_FILE" \
    "Compare to reference"

echo ""
echo "=== Category 13: Advanced Analysis ==="
# LD calculator - computationally intensive
benchmark "ld_calculator" "advanced_analysis" \
    "$BUILD_DIR/VCFX_ld_calculator/VCFX_ld_calculator --window 1000 < $DATA_FILE" \
    "Calculate LD" 3600

# Custom annotator
benchmark "custom_annotator" "advanced_analysis" \
    "$BUILD_DIR/VCFX_custom_annotator/VCFX_custom_annotator --annotation 'TEST=1' < $DATA_FILE" \
    "Custom annotation"

echo ""
echo "=== Baseline Comparison: bcftools ==="
if command -v bcftools &> /dev/null; then
    benchmark "bcftools_view" "baseline" \
        "bcftools view -H $DATA_FILE | wc -l" \
        "bcftools count"

    benchmark "bcftools_query" "baseline" \
        "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $DATA_FILE" \
        "bcftools query"

    benchmark "bcftools_filter" "baseline" \
        "bcftools view -i 'QUAL>=30' $DATA_FILE" \
        "bcftools filter"
else
    echo "bcftools not installed - skipping baseline"
fi

echo ""
echo "=== Baseline Comparison: vcftools ==="
if command -v vcftools &> /dev/null; then
    benchmark "vcftools_freq" "baseline" \
        "vcftools --vcf $DATA_FILE --freq --stdout 2>/dev/null" \
        "vcftools freq"

    benchmark "vcftools_missing" "baseline" \
        "vcftools --vcf $DATA_FILE --missing-indv --stdout 2>/dev/null" \
        "vcftools missing"
else
    echo "vcftools not installed - skipping baseline"
fi

echo ""
echo "=============================================="
echo "   Benchmark Complete!"
echo "=============================================="
echo ""
echo "Results saved to: $RESULTS_CSV"
echo ""

# Print summary
echo "=== Summary by Category ==="
echo ""
python3 << PYSCRIPT
import csv
from collections import defaultdict

results = defaultdict(list)
with open("$RESULTS_CSV") as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['status'] == 'success':
            results[row['category']].append(float(row['time_seconds']))

print(f"{'Category':<25} {'Tools':<8} {'Avg (s)':<10} {'Min (s)':<10} {'Max (s)':<10}")
print("-" * 63)
for cat in sorted(results.keys()):
    times = results[cat]
    print(f"{cat:<25} {len(times):<8} {sum(times)/len(times):<10.3f} {min(times):<10.3f} {max(times):<10.3f}")

total_time = sum(t for times in results.values() for t in times)
total_tools = sum(len(times) for times in results.values())
print("-" * 63)
print(f"{'TOTAL':<25} {total_tools:<8} {total_time:<10.3f}")
PYSCRIPT

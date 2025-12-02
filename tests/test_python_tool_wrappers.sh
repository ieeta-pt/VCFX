#!/usr/bin/env bash
set -e
set -o pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"
BUILD_DIR="${ROOT_DIR}/build/python_tool_wrappers"

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake -DPYTHON_BINDINGS=ON ../..
make -j

# Extract the Python executable from CMake cache to ensure we use the same Python
# that the extension was built against
PYTHON_CONFIG=$(grep '_Python3_CONFIG:INTERNAL' CMakeCache.txt 2>/dev/null | cut -d= -f2)
if [ -n "$PYTHON_CONFIG" ] && [ -x "$PYTHON_CONFIG" ]; then
    # python3.13-config -> python3.13
    PYTHON_EXECUTABLE="${PYTHON_CONFIG%-config}"
fi

if [ -z "$PYTHON_EXECUTABLE" ] || [ ! -x "$PYTHON_EXECUTABLE" ]; then
    # Fallback: try to find Python from the include dir
    PYTHON_INCLUDE=$(grep '_Python3_INCLUDE_DIR:INTERNAL' CMakeCache.txt 2>/dev/null | cut -d= -f2)
    if [[ "$PYTHON_INCLUDE" =~ python([0-9]+\.[0-9]+) ]]; then
        PYTHON_VERSION="${BASH_REMATCH[1]}"
        # Try common locations
        for candidate in /usr/local/Cellar/python@${PYTHON_VERSION}/*/bin/python${PYTHON_VERSION} \
                         /usr/local/bin/python${PYTHON_VERSION} \
                         /opt/homebrew/Cellar/python@${PYTHON_VERSION}/*/bin/python${PYTHON_VERSION} \
                         /opt/homebrew/bin/python${PYTHON_VERSION}; do
            # Use eval to expand the glob
            for real_candidate in $candidate; do
                if [ -x "$real_candidate" ]; then
                    PYTHON_EXECUTABLE="$real_candidate"
                    break 2
                fi
            done
        done
    fi
fi

# Final fallback to python3 in PATH
if [ -z "$PYTHON_EXECUTABLE" ] || [ ! -x "$PYTHON_EXECUTABLE" ]; then
    PYTHON_EXECUTABLE="python3"
fi

echo "Using Python: $PYTHON_EXECUTABLE"

# ensure tools on PATH for the wrappers
for exe in "${BUILD_DIR}"/src/VCFX_*/*; do
    if [ -x "$exe" ]; then
        PATH="$(dirname "$exe"):${PATH}"
    fi
done
export PATH

cd "$SCRIPT_DIR"

PYTHONPATH="${BUILD_DIR}/python" "$PYTHON_EXECUTABLE" - <<'PY'
import vcfx

rows = vcfx.alignment_checker("data/align_Y.vcf", "data/align_refY.fa")
assert rows and rows[0].Discrepancy_Type == "ALT_MISMATCH"

counts = vcfx.allele_counter("data/allele_counter_A.vcf")
assert counts[0].Sample == "S1"
assert isinstance(counts[0].Ref_Count, int)

n = vcfx.variant_counter("data/variant_counter_normal.vcf")
assert n == 5

freqs = vcfx.allele_freq_calc("data/allele_freq_calc/simple.vcf")
assert abs(freqs[0].Allele_Frequency - 0.5) < 1e-6

annotated = vcfx.info_aggregator("data/aggregator/basic.vcf", ["DP"])
assert "#AGGREGATION_SUMMARY" in annotated

parsed = vcfx.info_parser("data/info_parser/basic.vcf", ["DP"])
assert parsed[0]["DP"] == "10"

summary = vcfx.info_summarizer("data/info_summarizer/basic.vcf", ["DP"])
assert abs(summary[0].Mean - 20.0) < 1e-6
assert isinstance(summary[0].Median, float)

fasta = vcfx.fasta_converter("data/fasta_converter/basic.vcf")
assert fasta.startswith(">")

balance = vcfx.allele_balance_calc("data/allele_balance_calc_A.vcf")
assert abs(balance[0].Allele_Balance - 1.0) < 1e-6

conc = vcfx.concordance_checker("data/concordance_input.vcf", "SAMPLE1", "SAMPLE2")
assert conc[0].Concordance == "Concordant"

filtered = vcfx.genotype_query("data/genotype_query/sample.vcf", "0/1")
assert filtered.startswith("##")

dedup = vcfx.duplicate_remover("data/allele_balance_calc_A.vcf")
assert dedup.startswith("##")

subset = vcfx.af_subsetter("data/af_subsetter_A.vcf", "0.01-0.1")
assert subset.startswith("##")

flagged = vcfx.missing_detector("data/concordance_missing_data.vcf")
assert "MISSING_GENOTYPES=1" in flagged

hwe_rows = vcfx.hwe_tester("data/hwe_tester/basic_hwe.vcf")
assert isinstance(hwe_rows[0].HWE_pvalue, float)

coeff = vcfx.inbreeding_calculator(
    "data/inbreeding_calculator/single_sample_excludeSample_false.vcf",
    freq_mode="excludeSample",
)
assert isinstance(coeff[0].InbreedingCoefficient, float)

classes = vcfx.variant_classifier("data/classifier_mixed.vcf")
assert classes[0].Classification

xconc = vcfx.cross_sample_concordance("data/concordance_some_mismatch.vcf")
assert xconc[0].Concordance_Status
assert isinstance(xconc[0].Num_Samples, int)

fields = vcfx.field_extractor("data/field_extractor_input.vcf", ["CHROM", "POS"])
assert fields[0]["CHROM"] == "chr1"

assign = vcfx.ancestry_assigner(
    "data/ancestry_assigner/input.vcf",
    "data/ancestry_assigner/freq.tsv",
)
assert assign[0].Assigned_Population == "EUR"

dos = vcfx.dosage_calculator("data/dosage_calculator/basic.vcf")
assert dos[0].Dosages == "0,1,2"

# Additional wrappers
inf = vcfx.ancestry_inferrer(
    "data/ancestry_inferrer/eur_samples.vcf",
    "data/ancestry_inferrer/population_freqs.txt",
)
assert inf[0].Inferred_Population == "EUR"

dist = vcfx.distance_calculator("data/variant_counter_normal.vcf")
assert isinstance(dist[1].DISTANCE, int)

index_rows = vcfx.indexer("data/indexer/basic.vcf")
assert isinstance(index_rows[0].FILE_OFFSET, int)

norm_text = vcfx.indel_normalizer("data/basic_indel.vcf")
assert norm_text.startswith("##")

val_msg = vcfx.validator("data/variant_counter_normal.vcf")
assert "VCF" in val_msg
print("Python tool wrappers OK")
PY

#include "VCFX_format_converter.h"
#include "vcfx_core.h"
#include "vcfx_io.h"
#include <algorithm>
#include <cctype>
#include <sstream>

// -----------------------------------------------------------------------
// printHelp
// -----------------------------------------------------------------------
void printHelp() {
    std::cout << "VCFX_format_converter\n"
              << "Usage: VCFX_format_converter [OPTIONS]\n\n"
              << "Options:\n"
              << "  --to-bed             Convert VCF to BED format.\n"
              << "  --to-csv             Convert VCF to CSV format.\n"
              << "  --help, -h           Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Converts VCF files to specified formats (BED or CSV).\n\n"
              << "Example:\n"
              << "  ./VCFX_format_converter --to-bed < input.vcf > output.bed\n"
              << "  ./VCFX_format_converter --to-csv < input.vcf > output.csv\n";
}

// -----------------------------------------------------------------------
// parseArguments
// -----------------------------------------------------------------------
bool parseArguments(int argc, char *argv[], OutputFormat &format) {
    format = OutputFormat::UNKNOWN;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--to-bed") {
            format = OutputFormat::BED;
        } else if (arg == "--to-csv") {
            format = OutputFormat::CSV;
        } else if (arg == "--help" || arg == "-h") {
            // We'll handle help outside
        }
    }
    return (format != OutputFormat::UNKNOWN);
}

// -----------------------------------------------------------------------
// convertVCFtoBED
//   - For each variant line, output a single BED line:
//     chrom, start=(pos-1 clamped to >=0), end=(pos-1 + ref.size()), name=id
// -----------------------------------------------------------------------
void convertVCFtoBED(std::istream &in, std::ostream &out) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            // Skip header or empty
            continue;
        }

        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string f;
            while (std::getline(ss, f, '\t')) {
                fields.push_back(f);
            }
        }

        if (fields.size() < 5) {
            // Not enough columns to parse properly
            continue;
        }

        const std::string &chrom = fields[0];
        int pos = 0;
        try {
            pos = std::stoi(fields[1]);
        } catch (...) {
            // invalid pos => skip
            continue;
        }
        const std::string &id = fields[2];
        const std::string &ref = fields[3];
        // alt = fields[4] if needed, but not used in basic bed

        // Start coordinate
        // If pos==1 => start=0
        int start = std::max(0, pos - 1);  // clamp to >=0
        int end = start + (int)ref.size(); // simplistic approach if ref is multi-base

        // Format: chrom, start, end, name
        out << chrom << "\t" << start << "\t" << end << "\t" << id << "\n";
    }
}

// -----------------------------------------------------------------------
// Helper function: CSV-escape a single field.
//   - If field has comma or quote, enclose in quotes and double the quotes.
// -----------------------------------------------------------------------
static std::string csvEscape(const std::string &field) {
    bool needQuotes = false;
    // check if it contains a comma or a quote
    if (field.find(',') != std::string::npos || field.find('"') != std::string::npos) {
        needQuotes = true;
    }
    // if it has any control chars or whitespace, you might also consider quoting
    // but let's just do commas and quotes

    if (!needQuotes) {
        return field;
    }
    // we do need quotes -> double all existing quotes
    std::string tmp;
    tmp.reserve(field.size() + 2);
    tmp.push_back('"');
    for (char c : field) {
        if (c == '"') {
            // double it by writing two quotes
            tmp.push_back('"');
            tmp.push_back('"');
        } else {
            tmp.push_back(c);
        }
    }
    tmp.push_back('"');
    return tmp;
}

// -----------------------------------------------------------------------
// convertVCFtoCSV
//   - For each variant line, produce a CSV line of the same columns
//   - We do minimal CSV: if a field has comma or quote, we enclose in quotes
//     and double any quotes inside.
// -----------------------------------------------------------------------
void convertVCFtoCSV(std::istream &in, std::ostream &out) {
    std::string line;
    bool wroteHeader = false;
    std::vector<std::string> headerCols; // from #CHROM line if we want to produce a CSV header

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '#') {
            // Optionally parse #CHROM line to produce a CSV header.
            // If you'd like a CSV header row, you can do something like:
            if (!wroteHeader && line.rfind("#CHROM", 0) == 0) {
                // parse the columns
                std::stringstream ss(line.substr(1)); // drop '#'
                std::vector<std::string> hdrTokens;
                {
                    std::string t;
                    while (std::getline(ss, t, '\t')) {
                        hdrTokens.push_back(t);
                    }
                }
                // output as CSV header
                for (size_t i = 0; i < hdrTokens.size(); i++) {
                    out << csvEscape(hdrTokens[i]);
                    if (i + 1 < hdrTokens.size())
                        out << ",";
                }
                out << "\n";
                wroteHeader = true;
            }
            // Skip the actual variant lines in the header
            continue;
        }

        // This is a data line
        std::stringstream ss(line);
        std::vector<std::string> fields;
        {
            std::string f;
            while (std::getline(ss, f, '\t')) {
                fields.push_back(f);
            }
        }

        // Join them as CSV, escaping if needed
        for (size_t i = 0; i < fields.size(); i++) {
            out << csvEscape(fields[i]);
            if (i + 1 < fields.size())
                out << ",";
        }
        out << "\n";
    }
}

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_format_converter", show_help))
        return 0;
    OutputFormat format;
    bool valid = parseArguments(argc, argv, format);

    // Check if user asked for help
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--help" || a == "-h") {
            printHelp();
            return 0;
        }
    }

    if (!valid || format == OutputFormat::UNKNOWN) {
        std::cerr << "No valid output format specified (--to-bed or --to-csv).\n";
        printHelp();
        return 1;
    }

    switch (format) {
    case OutputFormat::BED:
        convertVCFtoBED(std::cin, std::cout);
        break;
    case OutputFormat::CSV:
        convertVCFtoCSV(std::cin, std::cout);
        break;
    default:
        std::cerr << "Unsupported output format.\n";
        return 1;
    }

    return 0;
}

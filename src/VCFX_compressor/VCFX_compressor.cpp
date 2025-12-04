#include "vcfx_core.h" 
#include "vcfx_io.h"
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <string>
#include <zlib.h>

// ---------------------------------------------------------------------------
// Show help
// ---------------------------------------------------------------------------
static void printHelp() {
    std::cout << "VCFX_compressor\n"
              << "Usage: VCFX_compressor [OPTIONS]\n\n"
              << "Options:\n"
              << "  --compress, -c         Compress the input VCF file (to stdout).\n"
              << "  --decompress, -d       Decompress the input VCF file (from stdin).\n"
              << "  --help, -h             Display this help message and exit.\n\n"
              << "Description:\n"
              << "  Compresses or decompresses data using zlib's raw DEFLATE (similar to gzip).\n"
              << "  Note that for .vcf.gz indexing via tabix, one typically needs BGZF blocks,\n"
              << "  which is not implemented here.\n\n"
              << "Examples:\n"
              << "  Compress:\n"
              << "    ./VCFX_compressor --compress < input.vcf > output.vcf.gz\n\n"
              << "  Decompress:\n"
              << "    ./VCFX_compressor --decompress < input.vcf.gz > output.vcf\n";
}

// ---------------------------------------------------------------------------
// compressDecompressVCF
//   compress = true  => read from 'in', produce gzip to 'out'
//   compress = false => read gzip from 'in', produce plain text to 'out'
// ---------------------------------------------------------------------------
static bool compressDecompressVCF(std::istream &in, std::ostream &out, bool compress) {
    constexpr int CHUNK = 16384;
    char inBuffer[CHUNK];
    char outBuffer[CHUNK];

    if (compress) {
        // Initialize for deflate
        z_stream strm;
        std::memset(&strm, 0, sizeof(strm));
        if (deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
            std::cerr << "Error: deflateInit2 failed.\n";
            return false;
        }

        // read -> deflate -> out
        int flush = Z_NO_FLUSH;
        do {
            in.read(inBuffer, CHUNK);
            std::streamsize bytesRead = in.gcount();
            flush = in.eof() ? Z_FINISH : Z_NO_FLUSH;

            strm.avail_in = static_cast<uInt>(bytesRead);
            strm.next_in = reinterpret_cast<Bytef *>(inBuffer);

            // compress until input is used up
            do {
                strm.avail_out = CHUNK;
                strm.next_out = reinterpret_cast<Bytef *>(outBuffer);

                int ret = deflate(&strm, flush);
                if (ret == Z_STREAM_ERROR) {
                    std::cerr << "Error: deflate failed.\n";
                    deflateEnd(&strm);
                    return false;
                }
                // # of bytes written
                size_t have = CHUNK - strm.avail_out;
                if (have > 0) {
                    out.write(outBuffer, have);
                    if (!out.good()) {
                        std::cerr << "Error: write to output stream failed.\n";
                        deflateEnd(&strm);
                        return false;
                    }
                }
            } while (strm.avail_out == 0);

        } while (flush != Z_FINISH);

        deflateEnd(&strm);
        return true;

    } else {
        // Initialize for inflate
        z_stream strm;
        std::memset(&strm, 0, sizeof(strm));
        // 15+32 to allow auto-detect of gzip/zlib
        if (inflateInit2(&strm, 15 + 32) != Z_OK) {
            std::cerr << "Error: inflateInit2 failed.\n";
            return false;
        }

        int ret = Z_OK;
        do {
            in.read(inBuffer, CHUNK);
            std::streamsize bytesRead = in.gcount();
            if (bytesRead == 0 && !in.eof()) {
                // Possibly an I/O error
                if (in.bad()) {
                    std::cerr << "Error: reading input stream.\n";
                    inflateEnd(&strm);
                    return false;
                }
            }
            strm.avail_in = static_cast<uInt>(bytesRead);
            strm.next_in = reinterpret_cast<Bytef *>(inBuffer);

            // decompress until output buffer not needed
            do {
                strm.avail_out = CHUNK;
                strm.next_out = reinterpret_cast<Bytef *>(outBuffer);

                ret = inflate(&strm, Z_NO_FLUSH);
                if (ret == Z_STREAM_ERROR || ret == Z_NEED_DICT || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR) {
                    std::cerr << "Error: inflate failed with code " << ret << "\n";
                    inflateEnd(&strm);
                    return false;
                }
                size_t have = CHUNK - strm.avail_out;
                if (have > 0) {
                    out.write(outBuffer, have);
                    if (!out.good()) {
                        std::cerr << "Error: write to output stream failed.\n";
                        inflateEnd(&strm);
                        return false;
                    }
                }
            } while (strm.avail_out == 0);

        } while (ret != Z_STREAM_END && !in.eof());

        // Clean up
        inflateEnd(&strm);
        // If we didn't get Z_STREAM_END, the compressed data might be incomplete
        if (ret != Z_STREAM_END) {
            // Possibly truncated
            std::cerr << "Warning: stream ended prematurely or was truncated.\n";
        }
        return true;
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
static void show_help() { printHelp(); }

int main(int argc, char *argv[]) {
    vcfx::init_io();  // Performance: disable sync_with_stdio
    if (vcfx::handle_common_flags(argc, argv, "VCFX_compressor", show_help))
        return 0;
    bool compress = false;
    bool decompress = false;

    static struct option long_options[] = {{"compress", no_argument, 0, 'c'},
                                           {"decompress", no_argument, 0, 'd'},
                                           {"help", no_argument, 0, 'h'},
                                           {0, 0, 0, 0}};

    while (true) {
        int optIndex = 0;
        int c = getopt_long(argc, argv, "cdh", long_options, &optIndex);
        if (c == -1)
            break;

        switch (c) {
        case 'c':
            compress = true;
            break;
        case 'd':
            decompress = true;
            break;
        case 'h':
            printHelp();
            return 0;
        default:
            printHelp();
            return 1;
        }
    }

    if ((compress && decompress) || (!compress && !decompress)) {
        std::cerr << "Error: must specify exactly one of --compress or --decompress.\n";
        return 1;
    }

    if (!compressDecompressVCF(std::cin, std::cout, compress)) {
        std::cerr << "Error: Compression/Decompression failed.\n";
        return 1;
    }
    return 0;
}

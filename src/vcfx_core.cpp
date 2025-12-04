#include "vcfx_core.h"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <zlib.h>

namespace vcfx {

std::string trim(const std::string &str) {
    auto first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) {
        return "";
    }
    auto last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}

std::vector<std::string> split(const std::string &str, char delimiter) {
    std::vector<std::string> result;
    std::istringstream iss(str);
    std::string item;
    while (std::getline(iss, item, delimiter)) {
        result.push_back(item);
    }
    return result;
}

bool flag_present(int argc, char *argv[], const char *long_flag, const char *short_flag) {
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], long_flag) == 0 || (short_flag && std::strcmp(argv[i], short_flag) == 0)) {
            return true;
        }
    }
    return false;
}

void print_error(const std::string &msg, std::ostream &os) { os << "Error: " << msg << '\n'; }

void print_version(const std::string &tool, const std::string &version, std::ostream &os) {
    os << tool << " version " << version << '\n';
}

// ------------------------------------------------------------
// Internal helper: decompress gzip/BGZF data from 'in' into 'out'
// ------------------------------------------------------------
static bool decompress_gzip_stream(std::istream &in, std::string &out) {
    constexpr int CHUNK = 16384;
    char inBuf[CHUNK];
    char outBuf[CHUNK];

    z_stream strm;
    std::memset(&strm, 0, sizeof(strm));
    if (inflateInit2(&strm, 15 + 32) != Z_OK) {
        return false;
    }

    int ret = Z_OK;
    do {
        in.read(inBuf, CHUNK);
        strm.avail_in = static_cast<uInt>(in.gcount());
        if (strm.avail_in == 0 && in.eof()) {
            break;
        }
        strm.next_in = reinterpret_cast<Bytef *>(inBuf);

        do {
            strm.avail_out = CHUNK;
            strm.next_out = reinterpret_cast<Bytef *>(outBuf);
            ret = inflate(&strm, Z_NO_FLUSH);
            if (ret == Z_STREAM_ERROR || ret == Z_NEED_DICT || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR) {
                inflateEnd(&strm);
                return false;
            }
            size_t have = CHUNK - strm.avail_out;
            if (have > 0) {
                out.append(outBuf, have);
            }
        } while (strm.avail_out == 0);
    } while (ret != Z_STREAM_END);

    inflateEnd(&strm);
    return ret == Z_STREAM_END;
}

// ------------------------------------------------------------
// Detect gzip magic numbers on a stream without consuming them
// ------------------------------------------------------------
static bool stream_has_gzip_magic(std::istream &in) {
    int c1 = in.get();
    if (c1 == EOF) {
        return false;
    }
    int c2 = in.get();
    if (c2 == EOF) {
        in.unget();
        return false;
    }
    bool isGz = (static_cast<unsigned char>(c1) == 0x1f && static_cast<unsigned char>(c2) == 0x8b);
    in.putback(static_cast<char>(c2));
    in.putback(static_cast<char>(c1));
    return isGz;
}

bool read_maybe_compressed(std::istream &in, std::string &out) {
    out.clear();
    if (stream_has_gzip_magic(in)) {
        return decompress_gzip_stream(in, out);
    }
    std::ostringstream oss;
    oss << in.rdbuf();
    out = oss.str();
    return true;
}

bool read_file_maybe_compressed(const std::string &path, std::string &out) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) {
        return false;
    }
    bool isGz = false;
    if (path.size() >= 3 && (path.compare(path.size() - 3, 3, ".gz") == 0)) {
        isGz = true;
    } else if (path.size() >= 4 && (path.compare(path.size() - 4, 4, ".bgz") == 0)) {
        isGz = true;
    } else if (path.size() >= 5 && (path.compare(path.size() - 5, 5, ".bgzf") == 0)) {
        isGz = true;
    }
    if (isGz || stream_has_gzip_magic(file)) {
        return decompress_gzip_stream(file, out);
    }
    std::ostringstream oss;
    oss << file.rdbuf();
    out = oss.str();
    return true;
}

// ------------------------------------------------------------
// StreamingGzipReader Implementation
// ------------------------------------------------------------

StreamingGzipReader::StreamingGzipReader(std::istream &in) : in_(in) {
    // Check for gzip magic bytes
    int c1 = in_.get();
    if (c1 == EOF) {
        eof_ = true;
        return;
    }
    int c2 = in_.get();
    if (c2 == EOF) {
        in_.unget();
        isCompressed_ = false;
        return;
    }

    isCompressed_ = (static_cast<unsigned char>(c1) == 0x1f &&
                     static_cast<unsigned char>(c2) == 0x8b);

    // Put bytes back
    in_.putback(static_cast<char>(c2));
    in_.putback(static_cast<char>(c1));

    if (isCompressed_) {
        if (!initZlib()) {
            error_ = true;
        }
    }
}

StreamingGzipReader::~StreamingGzipReader() {
    if (zstrm_) {
        z_stream *strm = static_cast<z_stream *>(zstrm_);
        inflateEnd(strm);
        delete strm;
        zstrm_ = nullptr;
    }
}

StreamingGzipReader::StreamingGzipReader(StreamingGzipReader &&other) noexcept
    : in_(other.in_),
      isCompressed_(other.isCompressed_),
      initialized_(other.initialized_),
      eof_(other.eof_),
      error_(other.error_),
      zstrm_(other.zstrm_),
      inBufPos_(other.inBufPos_),
      inBufLen_(other.inBufLen_),
      lineBuffer_(std::move(other.lineBuffer_)) {
    std::memcpy(inBuf_, other.inBuf_, CHUNK_SIZE);
    std::memcpy(outBuf_, other.outBuf_, CHUNK_SIZE);
    other.zstrm_ = nullptr;
    other.initialized_ = false;
}

StreamingGzipReader &StreamingGzipReader::operator=(StreamingGzipReader &&other) noexcept {
    if (this != &other) {
        // Clean up existing resources
        if (zstrm_) {
            z_stream *strm = static_cast<z_stream *>(zstrm_);
            inflateEnd(strm);
            delete strm;
        }

        isCompressed_ = other.isCompressed_;
        initialized_ = other.initialized_;
        eof_ = other.eof_;
        error_ = other.error_;
        zstrm_ = other.zstrm_;
        inBufPos_ = other.inBufPos_;
        inBufLen_ = other.inBufLen_;
        lineBuffer_ = std::move(other.lineBuffer_);
        std::memcpy(inBuf_, other.inBuf_, CHUNK_SIZE);
        std::memcpy(outBuf_, other.outBuf_, CHUNK_SIZE);

        other.zstrm_ = nullptr;
        other.initialized_ = false;
    }
    return *this;
}

bool StreamingGzipReader::initZlib() {
    z_stream *strm = new z_stream;
    std::memset(strm, 0, sizeof(z_stream));

    // 15 + 32 enables gzip decoding with automatic header detection
    if (inflateInit2(strm, 15 + 32) != Z_OK) {
        delete strm;
        return false;
    }

    zstrm_ = strm;
    initialized_ = true;
    return true;
}

int StreamingGzipReader::decompressChunk() {
    if (!zstrm_ || error_) {
        return -1;
    }

    z_stream *strm = static_cast<z_stream *>(zstrm_);

    // Read more input if needed
    if (strm->avail_in == 0 && !in_.eof()) {
        in_.read(inBuf_, CHUNK_SIZE);
        inBufLen_ = static_cast<size_t>(in_.gcount());
        if (inBufLen_ == 0 && in_.eof()) {
            eof_ = true;
            return 0;
        }
        strm->avail_in = static_cast<uInt>(inBufLen_);
        strm->next_in = reinterpret_cast<Bytef *>(inBuf_);
    }

    strm->avail_out = CHUNK_SIZE;
    strm->next_out = reinterpret_cast<Bytef *>(outBuf_);

    int ret = inflate(strm, Z_NO_FLUSH);

    if (ret == Z_STREAM_ERROR || ret == Z_NEED_DICT ||
        ret == Z_DATA_ERROR || ret == Z_MEM_ERROR) {
        error_ = true;
        return -1;
    }

    size_t have = CHUNK_SIZE - strm->avail_out;

    if (ret == Z_STREAM_END) {
        // Check for concatenated gzip streams (like BGZF)
        if (strm->avail_in > 0 || !in_.eof()) {
            // Reset for next stream
            inflateReset(strm);
        } else {
            eof_ = true;
        }
    }

    return static_cast<int>(have);
}

bool StreamingGzipReader::readUncompressed() {
    if (in_.eof()) {
        eof_ = true;
        return false;
    }

    in_.read(outBuf_, CHUNK_SIZE);
    size_t bytesRead = static_cast<size_t>(in_.gcount());

    if (bytesRead == 0) {
        eof_ = true;
        return false;
    }

    lineBuffer_.append(outBuf_, bytesRead);
    return true;
}

bool StreamingGzipReader::getline(std::string &line) {
    line.clear();

    if (error_) {
        return false;
    }

    while (true) {
        // Check if we have a complete line in the buffer
        size_t newlinePos = lineBuffer_.find('\n');
        if (newlinePos != std::string::npos) {
            // Extract the line (without newline)
            if (newlinePos > 0 && lineBuffer_[newlinePos - 1] == '\r') {
                // Handle CRLF
                line = lineBuffer_.substr(0, newlinePos - 1);
            } else {
                line = lineBuffer_.substr(0, newlinePos);
            }
            lineBuffer_.erase(0, newlinePos + 1);
            return true;
        }

        // Need more data
        if (eof_) {
            // Return remaining data as last line
            if (!lineBuffer_.empty()) {
                line = std::move(lineBuffer_);
                lineBuffer_.clear();
                return true;
            }
            return false;
        }

        if (isCompressed_) {
            int bytesDecompressed = decompressChunk();
            if (bytesDecompressed < 0) {
                return false;  // Error
            }
            if (bytesDecompressed > 0) {
                lineBuffer_.append(outBuf_, static_cast<size_t>(bytesDecompressed));
            }
        } else {
            if (!readUncompressed()) {
                // EOF reached, return remaining buffer if any
                if (!lineBuffer_.empty()) {
                    line = std::move(lineBuffer_);
                    lineBuffer_.clear();
                    return true;
                }
                return false;
            }
        }
    }
}

std::unique_ptr<StreamingGzipReader> make_streaming_reader(std::istream &in) {
    auto reader = std::make_unique<StreamingGzipReader>(in);
    if (reader->error()) {
        return nullptr;
    }
    return reader;
}

std::unique_ptr<StreamingGzipReader> make_streaming_reader(const std::string &path,
                                                            std::ifstream &fileStream) {
    fileStream.open(path, std::ios::binary);
    if (!fileStream.is_open()) {
        return nullptr;
    }
    return make_streaming_reader(fileStream);
}

} // namespace vcfx

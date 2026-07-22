#include "mls_project_clean.h"

#include <cerrno>
#include <clocale>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

[[noreturn]] void fail(const std::string& message) {
    const char* fixture = std::getenv("MLS_FIXTURE");
    if (fixture && *fixture)
        std::fprintf(stderr, "mls_clean_runner [%s]: %s\n", fixture, message.c_str());
    else
        std::fprintf(stderr, "mls_clean_runner: %s\n", message.c_str());
    std::exit(2);
}

const char* required_env(const char* name) {
    const char* value = std::getenv(name);
    if (!value || !*value) fail(std::string("missing environment variable ") + name);
    return value;
}

bool host_is_little_endian() {
    const uint16_t value = 1;
    return *reinterpret_cast<const uint8_t*>(&value) == 1;
}

uint32_t byte_swap32(uint32_t x) {
    return ((x & UINT32_C(0x000000ff)) << 24) |
           ((x & UINT32_C(0x0000ff00)) << 8) |
           ((x & UINT32_C(0x00ff0000)) >> 8) |
           ((x & UINT32_C(0xff000000)) >> 24);
}

size_t checked_value_count(size_t triples) {
    if (triples > std::numeric_limits<size_t>::max() / 3)
        fail("triple count overflows size_t");
    return triples * 3;
}

std::vector<float> read_raw_f32(const char* path, size_t triples) {
    const size_t values = checked_value_count(triples);
    if (values > std::numeric_limits<size_t>::max() / sizeof(float))
        fail("raw byte count overflows size_t");
    const size_t bytes = values * sizeof(float);
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) fail(std::string("cannot open raw input: ") + path);
    const std::streamoff length = input.tellg();
    if (length < 0 || static_cast<uint64_t>(length) != static_cast<uint64_t>(bytes))
        fail(std::string("raw input has wrong byte size: ") + path);
    input.seekg(0);
    std::vector<float> result(values);
    if (bytes != 0) input.read(reinterpret_cast<char*>(result.data()),
                               static_cast<std::streamsize>(bytes));
    if (!input) fail(std::string("failed to read raw input: ") + path);
    if (!host_is_little_endian()) {
        for (float& value : result) {
            uint32_t bits;
            std::memcpy(&bits, &value, sizeof(bits));
            bits = byte_swap32(bits);
            std::memcpy(&value, &bits, sizeof(bits));
        }
    }
    return result;
}

void write_raw_f32(const char* path, const std::vector<float>& values) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) fail(std::string("cannot open raw output: ") + path);
    if (host_is_little_endian()) {
        if (!values.empty())
            output.write(reinterpret_cast<const char*>(values.data()),
                         static_cast<std::streamsize>(values.size() * sizeof(float)));
    } else {
        for (float value : values) {
            uint32_t bits;
            std::memcpy(&bits, &value, sizeof(bits));
            bits = byte_swap32(bits);
            output.write(reinterpret_cast<const char*>(&bits), sizeof(bits));
        }
    }
    if (!output) fail(std::string("failed to write raw output: ") + path);
}

std::unordered_map<std::string, std::string> read_params(const char* path) {
    std::ifstream input(path);
    if (!input) fail(std::string("cannot open params: ") + path);
    std::unordered_map<std::string, std::string> params;
    std::string line;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const size_t equals = line.find('=');
        if (equals == std::string::npos || equals == 0)
            fail("malformed params line " + std::to_string(line_number));
        const std::string key = line.substr(0, equals);
        const std::string value = line.substr(equals + 1);
        if (!params.emplace(key, value).second)
            fail("duplicate params key: " + key);
    }
    if (!input.eof()) fail(std::string("failed to read params: ") + path);
    return params;
}

const std::string& param(const std::unordered_map<std::string, std::string>& params,
                         const char* key) {
    const auto it = params.find(key);
    if (it == params.end()) fail(std::string("missing params key: ") + key);
    return it->second;
}

size_t parse_size(const std::string& text, const char* key) {
    if (text.empty() || text[0] == '-') fail(std::string("invalid ") + key);
    errno = 0;
    char* end = nullptr;
    const unsigned long long value = std::strtoull(text.c_str(), &end, 10);
    if (errno || end == text.c_str() || *end != '\0' ||
        value > std::numeric_limits<size_t>::max())
        fail(std::string("invalid ") + key);
    return static_cast<size_t>(value);
}

float parse_float(const std::string& text, const char* key) {
    errno = 0;
    char* end = nullptr;
    const float value = std::strtof(text.c_str(), &end);
    if (errno || end == text.c_str() || *end != '\0' || !std::isfinite(value))
        fail(std::string("invalid ") + key);
    return value;
}

uint32_t parse_hex32(const std::string& text, const char* key) {
    errno = 0;
    char* end = nullptr;
    const unsigned long value = std::strtoul(text.c_str(), &end, 0);
    if (errno || end == text.c_str() || *end != '\0' || value > UINT32_MAX)
        fail(std::string("invalid ") + key);
    return static_cast<uint32_t>(value);
}

void parse_origin(const std::string& text, float origin[3]) {
    size_t begin = 0;
    for (int axis = 0; axis < 3; ++axis) {
        const size_t comma = (axis < 2) ? text.find(',', begin) : std::string::npos;
        if (axis < 2 && comma == std::string::npos) fail("invalid cell_origin");
        const size_t end = (axis < 2) ? comma : text.size();
        origin[axis] = parse_float(text.substr(begin, end - begin), "cell_origin");
        begin = end + 1;
    }
    if (begin < text.size() + 1) fail("invalid cell_origin");
}

} /* namespace */

int main() {
    if (!std::setlocale(LC_NUMERIC, "C")) fail("cannot select C numeric locale");
    const char* samples_path = required_env("MLS_SAMPLES_RAW");
    const char* queries_path = required_env("MLS_QUERIES_RAW");
    const char* params_path = required_env("MLS_PARAMS");
    const char* positions_path = required_env("MLS_POSITIONS_RAW");
    const char* normals_path = required_env("MLS_NORMALS_RAW");

    const auto params = read_params(params_path);
    if (param(params, "format") != "raw-little-endian-fp32-xyz-triples")
        fail("unsupported raw format");
    const size_t sample_count = parse_size(param(params, "sample_count"), "sample_count");
    const size_t query_count = parse_size(param(params, "query_count"), "query_count");
    const size_t iterations = parse_size(param(params, "iterations"), "iterations");
    if (iterations != 1) fail("iterations must equal 1 (one projection per call)");
    const float radius = parse_float(param(params, "radius"), "radius");
    if (!(radius > 0.0f)) fail("radius must be positive");
    const uint32_t declared_radius_bits = parse_hex32(param(params, "radius_bits"),
                                                      "radius_bits");
    uint32_t radius_bits;
    std::memcpy(&radius_bits, &radius, sizeof(radius_bits));
    if (radius_bits != declared_radius_bits) fail("radius and radius_bits disagree");
    float origin[3];
    parse_origin(param(params, "cell_origin"), origin);

    std::vector<float> samples = read_raw_f32(samples_path, sample_count);
    std::vector<float> queries = read_raw_f32(queries_path, query_count);
    std::vector<float> positions(checked_value_count(query_count));
    std::vector<float> normals(checked_value_count(query_count));

    /* Two cells per support radius gives bounded 5x5x5 candidate-cell lookup. */
    const float cell_width = radius * 0.5f;
    if (!(cell_width > 0.0f) || !std::isfinite(cell_width))
        fail("radius cannot produce a valid cell width");
    mls_clean_arena* arena = nullptr;
    mls_clean_status status = mls_clean_arena_create(samples.data(), sample_count,
                                                     cell_width, origin, &arena);
    if (status != MLS_CLEAN_STATUS_OK)
        fail(std::string("arena construction failed: ") + mls_clean_status_string(status));

    MLS_project_verts(arena, queries.data(), query_count, radius, origin,
                      positions.data(), normals.data());
    status = mls_clean_last_status();
    if (status != MLS_CLEAN_STATUS_OK) {
        mls_clean_arena_info info{};
        (void)mls_clean_arena_get_info(arena, &info);
        const std::string message = std::string("projection failed: ") +
            mls_clean_status_string(status) + " (problem_queries=" +
            std::to_string(info.last_problem_query_count) + ")";
        mls_clean_arena_destroy(arena);
        fail(message);
    }
    mls_clean_arena_destroy(arena);
    if (mls_clean_last_status() != MLS_CLEAN_STATUS_OK)
        fail("arena destruction failed");

    write_raw_f32(positions_path, positions);
    write_raw_f32(normals_path, normals);
    return 0;
}

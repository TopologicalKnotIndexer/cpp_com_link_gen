#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <chrono>
#include <condition_variable>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <optional>
#include <queue>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifndef DEBUG
#define DEBUG 0
#define CPP_COM_LINK_GEN_DEFINED_PD_DIAGRAM_DEBUG 1
#endif
#include "PdToDiagram2d.h"
#ifdef CPP_COM_LINK_GEN_DEFINED_PD_DIAGRAM_DEBUG
#undef DEBUG
#undef CPP_COM_LINK_GEN_DEFINED_PD_DIAGRAM_DEBUG
#endif
#ifdef SHOW_CERTAIN_DEBUG_MESSAGE
#undef SHOW_CERTAIN_DEBUG_MESSAGE
#endif
#ifdef SHOW_DEBUG_MESSAGE
#undef SHOW_DEBUG_MESSAGE
#endif
#ifdef THROW_EXCEPTION
#undef THROW_EXCEPTION
#endif
#ifdef DEFINE_EXCEPTION
#undef DEFINE_EXCEPTION
#endif
#ifdef PROCESS_EXCEPTION
#undef PROCESS_EXCEPTION
#endif

extern "C" {
char* cppkh_compute_pd_signed_variants_ex(const char* pd_code, const char* signs_text, int reorder_crossings);
const char* cppkh_last_error();
void cppkh_free(char* value);
}

namespace fs = std::filesystem;

namespace {

using Crossing = std::array<int, 4>;
using PDCode = std::vector<Crossing>;

constexpr const char* kComLinkGenVersion = "0.1.0";

struct Slot {
    size_t crossing = 0;
    int pos = 0;
};

struct LinkNameInfo {
    bool mirror = false;
    bool isLink = false;
    bool nonAlternating = false;
    int crossing = 0;
    int index = 0;
};

struct LinkRep {
    std::map<std::string, PDCode> definitions;
    std::vector<std::string> linkSet;
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> methods;
};

struct NormalizeResult {
    PDCode pd;
    std::unordered_map<int, int> oldToNew;
};

struct ConnectedSumResult {
    PDCode pd;
    std::unordered_map<int, int> aMap;
    std::unordered_map<int, int> bMap;
};

std::string trim(const std::string& s) {
    size_t first = 0;
    while (first < s.size() && std::isspace(static_cast<unsigned char>(s[first]))) ++first;
    size_t last = s.size();
    while (last > first && std::isspace(static_cast<unsigned char>(s[last - 1]))) --last;
    return s.substr(first, last - first);
}

bool startsWith(const std::string& value, const std::string& prefix) {
    return value.size() >= prefix.size() &&
           std::equal(prefix.begin(), prefix.end(), value.begin());
}

uint32_t leftRotate(uint32_t value, uint32_t amount) {
    return (value << amount) | (value >> (32U - amount));
}

std::string md5Hex(const std::string& input) {
    static const uint32_t shifts[64] = {
        7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22,
        5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20,
        4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
        6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21,
    };
    static const uint32_t constants[64] = {
        0xd76aa478U, 0xe8c7b756U, 0x242070dbU, 0xc1bdceeeU,
        0xf57c0fafU, 0x4787c62aU, 0xa8304613U, 0xfd469501U,
        0x698098d8U, 0x8b44f7afU, 0xffff5bb1U, 0x895cd7beU,
        0x6b901122U, 0xfd987193U, 0xa679438eU, 0x49b40821U,
        0xf61e2562U, 0xc040b340U, 0x265e5a51U, 0xe9b6c7aaU,
        0xd62f105dU, 0x02441453U, 0xd8a1e681U, 0xe7d3fbc8U,
        0x21e1cde6U, 0xc33707d6U, 0xf4d50d87U, 0x455a14edU,
        0xa9e3e905U, 0xfcefa3f8U, 0x676f02d9U, 0x8d2a4c8aU,
        0xfffa3942U, 0x8771f681U, 0x6d9d6122U, 0xfde5380cU,
        0xa4beea44U, 0x4bdecfa9U, 0xf6bb4b60U, 0xbebfbc70U,
        0x289b7ec6U, 0xeaa127faU, 0xd4ef3085U, 0x04881d05U,
        0xd9d4d039U, 0xe6db99e5U, 0x1fa27cf8U, 0xc4ac5665U,
        0xf4292244U, 0x432aff97U, 0xab9423a7U, 0xfc93a039U,
        0x655b59c3U, 0x8f0ccc92U, 0xffeff47dU, 0x85845dd1U,
        0x6fa87e4fU, 0xfe2ce6e0U, 0xa3014314U, 0x4e0811a1U,
        0xf7537e82U, 0xbd3af235U, 0x2ad7d2bbU, 0xeb86d391U,
    };

    std::vector<uint8_t> message(input.begin(), input.end());
    uint64_t bitLength = static_cast<uint64_t>(message.size()) * 8U;
    message.push_back(0x80U);
    while ((message.size() % 64U) != 56U) message.push_back(0U);
    for (int i = 0; i < 8; ++i) {
        message.push_back(static_cast<uint8_t>((bitLength >> (8 * i)) & 0xffU));
    }

    uint32_t a0 = 0x67452301U;
    uint32_t b0 = 0xefcdab89U;
    uint32_t c0 = 0x98badcfeU;
    uint32_t d0 = 0x10325476U;

    for (size_t offset = 0; offset < message.size(); offset += 64) {
        uint32_t words[16];
        for (int i = 0; i < 16; ++i) {
            size_t j = offset + static_cast<size_t>(i) * 4U;
            words[i] = static_cast<uint32_t>(message[j]) |
                       (static_cast<uint32_t>(message[j + 1]) << 8U) |
                       (static_cast<uint32_t>(message[j + 2]) << 16U) |
                       (static_cast<uint32_t>(message[j + 3]) << 24U);
        }

        uint32_t a = a0;
        uint32_t b = b0;
        uint32_t c = c0;
        uint32_t d = d0;

        for (uint32_t i = 0; i < 64; ++i) {
            uint32_t f = 0;
            uint32_t g = 0;
            if (i < 16) {
                f = (b & c) | ((~b) & d);
                g = i;
            } else if (i < 32) {
                f = (d & b) | ((~d) & c);
                g = (5U * i + 1U) % 16U;
            } else if (i < 48) {
                f = b ^ c ^ d;
                g = (3U * i + 5U) % 16U;
            } else {
                f = c ^ (b | (~d));
                g = (7U * i) % 16U;
            }

            uint32_t temp = d;
            d = c;
            c = b;
            b = b + leftRotate(a + f + constants[i] + words[g], shifts[i]);
            a = temp;
        }

        a0 += a;
        b0 += b;
        c0 += c;
        d0 += d;
    }

    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (uint32_t value : {a0, b0, c0, d0}) {
        for (int i = 0; i < 4; ++i) {
            out << std::setw(2) << ((value >> (8 * i)) & 0xffU);
        }
    }
    return out.str();
}

std::vector<int> parseIntegers(const std::string& text) {
    std::vector<int> values;
    size_t pos = 0;
    while (pos < text.size()) {
        while (pos < text.size() && text[pos] != '-' &&
               !std::isdigit(static_cast<unsigned char>(text[pos]))) {
            ++pos;
        }
        if (pos >= text.size()) break;
        int sign = 1;
        if (text[pos] == '-') {
            sign = -1;
            ++pos;
        }
        if (pos >= text.size() || !std::isdigit(static_cast<unsigned char>(text[pos]))) {
            throw std::runtime_error("invalid integer literal");
        }
        int value = 0;
        while (pos < text.size() && std::isdigit(static_cast<unsigned char>(text[pos]))) {
            value = value * 10 + (text[pos] - '0');
            ++pos;
        }
        values.push_back(sign * value);
    }
    return values;
}

PDCode parsePDCode(const std::string& text) {
    std::string payload = text;
    size_t marker = payload.find("PD_CODE:");
    if (marker != std::string::npos) payload = payload.substr(marker + 8);

    std::vector<int> values = parseIntegers(payload);
    if (values.empty()) return {};
    if (values.size() % 4 != 0) {
        throw std::runtime_error("PD code does not contain a multiple of four integers");
    }
    PDCode pd;
    pd.reserve(values.size() / 4);
    for (size_t i = 0; i < values.size(); i += 4) {
        pd.push_back({values[i], values[i + 1], values[i + 2], values[i + 3]});
    }
    return pd;
}

std::string formatPDCode(const PDCode& pd) {
    std::ostringstream out;
    out << "[";
    for (size_t i = 0; i < pd.size(); ++i) {
        if (i) out << ", ";
        out << "[";
        for (int j = 0; j < 4; ++j) {
            if (j) out << ", ";
            out << pd[i][j];
        }
        out << "]";
    }
    out << "]";
    return out.str();
}

std::string formatKnotTheoryPD(const PDCode& pd) {
    std::ostringstream out;
    out << "PD[";
    for (size_t i = 0; i < pd.size(); ++i) {
        if (i) out << ",";
        out << "X[";
        for (int j = 0; j < 4; ++j) {
            if (j) out << ",";
            out << pd[i][j];
        }
        out << "]";
    }
    out << "]";
    return out.str();
}

void validatePDCode(const PDCode& pd) {
    std::map<int, int> counts;
    for (const auto& crossing : pd) {
        for (int label : crossing) {
            if (label <= 0) throw std::runtime_error("PD labels must be positive integers");
            ++counts[label];
        }
    }
    for (const auto& item : counts) {
        if (item.second != 2) {
            throw std::runtime_error("each PD label must appear exactly twice");
        }
    }
}

void validateContiguousPDLabels(const PDCode& pd) {
    std::set<int> labels;
    for (const auto& crossing : pd) {
        for (int label : crossing) labels.insert(label);
    }
    int expected = static_cast<int>(labels.size());
    if (expected != static_cast<int>(pd.size()) * 2) {
        throw std::runtime_error("PD label count is inconsistent with crossing count");
    }
    for (int label = 1; label <= expected; ++label) {
        if (!labels.count(label)) {
            throw std::runtime_error("PD labels must be contiguous before cppkh computation");
        }
    }
}

int pairedPos(int pos) {
    if (pos < 0 || pos >= 4) throw std::runtime_error("invalid crossing slot");
    return (pos + 2) % 4;
}

std::map<int, std::vector<Slot>> labelSlots(const PDCode& pd) {
    std::map<int, std::vector<Slot>> slots;
    for (size_t i = 0; i < pd.size(); ++i) {
        for (int pos = 0; pos < 4; ++pos) {
            slots[pd[i][pos]].push_back({i, pos});
        }
    }
    return slots;
}

std::map<int, std::vector<int>> labelPairedNeighbors(const PDCode& pd) {
    std::map<int, std::vector<int>> neighbors;
    for (const auto& crossing : pd) {
        for (int pos = 0; pos < 4; ++pos) {
            neighbors[crossing[pos]].push_back(crossing[pairedPos(pos)]);
        }
    }
    return neighbors;
}

std::vector<std::vector<int>> componentsFromPDCode(const PDCode& pd) {
    validatePDCode(pd);
    std::map<int, std::set<int>> graph;
    for (const auto& crossing : pd) {
        graph[crossing[0]].insert(crossing[2]);
        graph[crossing[2]].insert(crossing[0]);
        graph[crossing[1]].insert(crossing[3]);
        graph[crossing[3]].insert(crossing[1]);
    }

    std::set<int> visited;
    std::vector<std::vector<int>> components;
    for (const auto& item : graph) {
        int start = item.first;
        if (visited.count(start)) continue;
        std::vector<int> component;
        std::vector<int> stack{start};
        visited.insert(start);
        while (!stack.empty()) {
            int now = stack.back();
            stack.pop_back();
            component.push_back(now);
            for (int nxt : graph[now]) {
                if (!visited.count(nxt)) {
                    visited.insert(nxt);
                    stack.push_back(nxt);
                }
            }
        }
        std::sort(component.begin(), component.end());
        components.push_back(component);
    }
    std::sort(components.begin(), components.end());
    return components;
}

std::vector<int> walkCycle(const std::map<int, std::set<int>>& graph, int start, int firstNext, size_t expected) {
    std::vector<int> cycle;
    cycle.reserve(expected);
    int prev = -1;
    int cur = start;
    int nxt = firstNext;
    cycle.push_back(start);
    for (;;) {
        if (nxt == start) break;
        cycle.push_back(nxt);
        if (cycle.size() > expected) throw std::runtime_error("component cycle did not close");
        const auto& options = graph.at(nxt);
        int candidate = -1;
        for (int item : options) {
            if (item != cur) {
                candidate = item;
                break;
            }
        }
        if (candidate < 0) candidate = start;
        prev = cur;
        cur = nxt;
        (void)prev;
        nxt = candidate;
    }
    if (cycle.size() != expected) throw std::runtime_error("component cycle length mismatch");
    return cycle;
}

std::vector<int> canonicalCycleForComponent(const PDCode& pd, const std::vector<int>& component) {
    if (component.empty()) return {};
    if (component.size() == 1) return component;
    if (component.size() == 2) return component;

    auto paired = labelPairedNeighbors(pd);
    std::set<int> componentSet(component.begin(), component.end());
    std::map<int, std::set<int>> graph;
    for (int label : component) {
        for (int nbr : paired[label]) {
            if (componentSet.count(nbr)) graph[label].insert(nbr);
        }
        if (graph[label].size() != 2) {
            throw std::runtime_error("component is not a simple PD cycle");
        }
    }

    int start = *std::min_element(component.begin(), component.end());
    std::vector<std::vector<int>> candidates;
    for (int firstNext : graph[start]) {
        candidates.push_back(walkCycle(graph, start, firstNext, component.size()));
    }
    if (candidates.size() != 2) throw std::runtime_error("expected two cycle orientations");
    return std::lexicographical_compare(
               candidates[1].begin(), candidates[1].end(),
               candidates[0].begin(), candidates[0].end())
               ? candidates[1]
               : candidates[0];
}

std::vector<std::vector<int>> canonicalCycles(const PDCode& pd) {
    std::vector<std::vector<int>> components = componentsFromPDCode(pd);
    std::vector<std::vector<int>> cycles;
    cycles.reserve(components.size());
    for (const auto& component : components) cycles.push_back(canonicalCycleForComponent(pd, component));
    std::sort(cycles.begin(), cycles.end());
    return cycles;
}

std::pair<std::unordered_map<int, int>, std::unordered_map<int, int>> nextPrevMaps(
    const std::vector<std::vector<int>>& cycles) {
    std::unordered_map<int, int> next;
    std::unordered_map<int, int> prev;
    for (const auto& cycle : cycles) {
        if (cycle.empty()) continue;
        for (size_t i = 0; i < cycle.size(); ++i) {
            int label = cycle[i];
            next[label] = cycle[(i + 1) % cycle.size()];
            prev[label] = cycle[(i + cycle.size() - 1) % cycle.size()];
        }
    }
    return {next, prev};
}

NormalizeResult normalizePDCode(const PDCode& input) {
    validatePDCode(input);
    if (input.empty()) return {input, {}};

    std::vector<std::vector<int>> oldCycles = canonicalCycles(input);
    std::unordered_map<int, int> oldToNew;
    std::vector<std::vector<int>> newCycles;
    int nextId = 1;
    for (const auto& cycle : oldCycles) {
        std::vector<int> newCycle;
        for (int oldLabel : cycle) {
            oldToNew[oldLabel] = nextId;
            newCycle.push_back(nextId);
            ++nextId;
        }
        newCycles.push_back(std::move(newCycle));
    }

    PDCode pd = input;
    for (auto& crossing : pd) {
        for (int& label : crossing) label = oldToNew.at(label);
    }

    auto maps = nextPrevMaps(newCycles);
    const auto& next = maps.first;
    for (auto& crossing : pd) {
        if (next.at(crossing[0]) == crossing[2]) {
            continue;
        } else if (next.at(crossing[2]) == crossing[0]) {
            crossing = {crossing[2], crossing[3], crossing[0], crossing[1]};
            continue;
        } else {
            throw std::runtime_error("crossing is inconsistent with component orientation");
        }
    }
    std::sort(pd.begin(), pd.end());
    validatePDCode(pd);
    return {pd, oldToNew};
}

Slot endpointForLabel(const PDCode& pd,
                      int label,
                      const std::unordered_map<int, int>& next,
                      const std::unordered_map<int, int>& prev,
                      bool wantAfter) {
    auto slots = labelSlots(pd);
    auto found = slots.find(label);
    if (found == slots.end() || found->second.size() != 2) {
        throw std::runtime_error("connected-sum label must appear exactly twice");
    }
    int nxt = next.at(label);
    int prv = prev.at(label);

    if (nxt == prv) {
        return wantAfter ? found->second[1] : found->second[0];
    }

    for (const Slot& slot : found->second) {
        int paired = pd[slot.crossing][pairedPos(slot.pos)];
        if (wantAfter && paired == nxt) return slot;
        if (!wantAfter && paired == prv) return slot;
    }
    throw std::runtime_error("could not locate oriented endpoint for label");
}

int maxLabel(const PDCode& pd) {
    int ans = 0;
    for (const auto& crossing : pd) {
        for (int label : crossing) ans = std::max(ans, label);
    }
    return ans;
}

PDCode offsetPDCode(PDCode pd, int offset) {
    for (auto& crossing : pd) {
        for (int& label : crossing) label += offset;
    }
    return pd;
}

ConnectedSumResult connectedSum(const PDCode& aInput, const PDCode& bInput, int aLabel, int bLabel) {
    validatePDCode(aInput);
    validatePDCode(bInput);
    if (aInput.empty()) {
        NormalizeResult normalized = normalizePDCode(bInput);
        ConnectedSumResult result;
        result.pd = normalized.pd;
        result.bMap = normalized.oldToNew;
        return result;
    }
    if (bInput.empty()) {
        NormalizeResult normalized = normalizePDCode(aInput);
        ConnectedSumResult result;
        result.pd = normalized.pd;
        result.aMap = normalized.oldToNew;
        return result;
    }

    int offset = maxLabel(aInput);
    PDCode a = aInput;
    PDCode b = offsetPDCode(bInput, offset);
    int bOffsetLabel = bLabel + offset;

    auto aMaps = nextPrevMaps(canonicalCycles(a));
    auto bMaps = nextPrevMaps(canonicalCycles(b));
    Slot aAfter = endpointForLabel(a, aLabel, aMaps.first, aMaps.second, true);
    Slot bAfter = endpointForLabel(b, bOffsetLabel, bMaps.first, bMaps.second, true);

    PDCode combined = a;
    combined.insert(combined.end(), b.begin(), b.end());
    combined[aAfter.crossing][aAfter.pos] = bOffsetLabel;
    combined[a.size() + bAfter.crossing][bAfter.pos] = aLabel;

    NormalizeResult normalized = normalizePDCode(combined);
    ConnectedSumResult result;
    result.pd = normalized.pd;
    for (const auto& item : normalized.oldToNew) {
        if (item.first <= offset) result.aMap[item.first] = item.second;
        else result.bMap[item.first - offset] = item.second;
    }
    return result;
}

std::string takeCppkhString(char* raw) {
    if (!raw) {
        const char* err = cppkh_last_error();
        throw std::runtime_error(err && *err ? err : "cppkh computation failed");
    }
    std::string value(raw);
    cppkh_free(raw);
    return value;
}

std::vector<std::string> splitLines(const std::string& text) {
    std::vector<std::string> lines;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(line);
    }
    return lines;
}

int parseInvariantFactor(const std::string& text) {
    std::string token = trim(text);
    if (token.empty()) throw std::runtime_error("empty invariant factor in Z[...]");
    size_t consumed = 0;
    int value = std::stoi(token, &consumed);
    if (consumed != token.size()) {
        throw std::runtime_error("bad invariant factor in Z[...]: " + token);
    }
    return value;
}

std::string sortedInvariantFactorsText(const std::string& text) {
    std::string body = trim(text);
    if (body.empty()) return "";

    std::vector<int> values;
    std::stringstream ss(body);
    std::string item;
    while (std::getline(ss, item, ',')) {
        values.push_back(parseInvariantFactor(item));
    }
    std::sort(values.begin(), values.end());

    std::ostringstream out;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) out << ",";
        out << values[i];
    }
    return out.str();
}

std::string normalizeKhovanovZInvariants(const std::string& text) {
    std::string out;
    size_t pos = 0;
    while (true) {
        size_t start = text.find("Z[", pos);
        if (start == std::string::npos) {
            out.append(text.substr(pos));
            break;
        }
        size_t end = text.find(']', start + 2);
        if (end == std::string::npos) {
            throw std::runtime_error("unterminated Z[...] in Khovanov text");
        }
        out.append(text.substr(pos, start - pos));
        out.append("Z[");
        out.append(sortedInvariantFactorsText(text.substr(start + 2, end - start - 2)));
        out.push_back(']');
        pos = end + 1;
    }
    return out;
}

int baseCrossingSign(const Crossing& crossing) {
    int b = crossing[1];
    int d = crossing[3];
    if (b - d == 1 || d - b > 1) return 1;
    if (d - b == 1 || b - d > 1) return -1;
    throw std::runtime_error("error finding crossing sign");
}

std::vector<int> baseCrossingSigns(const PDCode& pd) {
    std::vector<int> signs;
    signs.reserve(pd.size());
    for (const auto& crossing : pd) signs.push_back(baseCrossingSign(crossing));
    return signs;
}

void validateCppkhOrientedNumbering(const PDCode& pd) {
    validatePDCode(pd);
    if (pd.empty()) return;
    validateContiguousPDLabels(pd);

    std::vector<std::vector<int>> cycles = canonicalCycles(pd);
    auto maps = nextPrevMaps(cycles);
    const auto& next = maps.first;

    auto checkAdjacentPair = [&](int a, int b) {
        if (next.at(a) != b && next.at(b) != a) {
            throw std::runtime_error("PD strand endpoints are not adjacent in component orientation");
        }
    };

    for (const auto& crossing : pd) {
        checkAdjacentPair(crossing[0], crossing[2]);
        checkAdjacentPair(crossing[1], crossing[3]);
        (void)baseCrossingSign(crossing);
    }
}

PDCode preparePDForCppkh(const PDCode& pd) {
    PDCode normalized = normalizePDCode(pd).pd;
    validateCppkhOrientedNumbering(normalized);
    return normalized;
}

std::mutex& cppkhMutex() {
    static std::mutex mutex;
    return mutex;
}

std::unordered_map<int, size_t> labelComponentIndex(const std::vector<std::vector<int>>& components) {
    std::unordered_map<int, size_t> index;
    for (size_t i = 0; i < components.size(); ++i) {
        for (int label : components[i]) index[label] = i;
    }
    return index;
}

std::string formatCrossingSigns(const std::vector<int>& signs) {
    if (signs.empty()) return "[]";
    std::ostringstream out;
    for (size_t i = 0; i < signs.size(); ++i) {
        if (i) out << ' ';
        out << signs[i];
    }
    return out.str();
}

std::string orientationSignsDocument(
    const PDCode& pd,
    const std::vector<std::vector<int>>& components,
    uint64_t orientationCount) {
    std::vector<int> baseSigns = baseCrossingSigns(pd);
    std::unordered_map<int, size_t> componentOf = labelComponentIndex(components);
    std::ostringstream document;

    for (uint64_t mask = 0; mask < orientationCount; ++mask) {
        std::vector<int> signs = baseSigns;
        for (size_t i = 0; i < pd.size(); ++i) {
            const auto& crossing = pd[i];
            size_t strandA = componentOf.at(crossing[0]);
            size_t strandAPaired = componentOf.at(crossing[2]);
            size_t strandB = componentOf.at(crossing[1]);
            size_t strandBPaired = componentOf.at(crossing[3]);
            if (strandA != strandAPaired || strandB != strandBPaired) {
                throw std::runtime_error("crossing strands are inconsistent with components");
            }
            bool reverseA = ((mask >> strandA) & 1ULL) != 0;
            bool reverseB = ((mask >> strandB) & 1ULL) != 0;
            if (strandA != strandB && reverseA != reverseB) signs[i] = -signs[i];
        }
        document << formatCrossingSigns(signs) << "\n";
    }
    return document.str();
}

std::string computeKhovanovSinglePD(const PDCode& pd) {
    PDCode oriented = preparePDForCppkh(pd);
    std::string pdText = formatKnotTheoryPD(oriented);
    std::string signsText = formatCrossingSigns(baseCrossingSigns(oriented)) + "\n";

    std::string rawOutput;
    {
        std::lock_guard<std::mutex> lock(cppkhMutex());
        rawOutput = takeCppkhString(cppkh_compute_pd_signed_variants_ex(pdText.c_str(), signsText.c_str(), 1));
    }

    std::vector<std::string> lines;
    for (std::string line : splitLines(rawOutput)) {
        line = trim(line);
        if (!line.empty()) lines.push_back(normalizeKhovanovZInvariants(line));
    }
    if (lines.size() != 1) {
        std::ostringstream err;
        err << "single PD Khovanov computation returned " << lines.size()
            << " non-empty line(s), expected exactly one";
        throw std::runtime_error(err.str());
    }
    return lines.front();
}

std::vector<std::string> computeKhovanovAllOrientations(const PDCode& pd) {
    PDCode oriented = preparePDForCppkh(pd);
    std::vector<std::vector<int>> components = componentsFromPDCode(oriented);
    if (components.size() >= 63) throw std::runtime_error("too many components to enumerate orientations");
    uint64_t orientationCount = components.empty() ? 1ULL : (1ULL << components.size());
    uint64_t maxDistinct = components.empty() ? 1ULL : (1ULL << (components.size() - 1));

    std::string pdText = formatKnotTheoryPD(oriented);
    std::string signsText = orientationSignsDocument(oriented, components, orientationCount);

    std::string rawOutput;
    {
        std::lock_guard<std::mutex> lock(cppkhMutex());
        rawOutput = takeCppkhString(cppkh_compute_pd_signed_variants_ex(pdText.c_str(), signsText.c_str(), 1));
    }

    std::set<std::string> unique;
    for (std::string line : splitLines(rawOutput)) {
        line = trim(line);
        if (!line.empty()) unique.insert(normalizeKhovanovZInvariants(line));
    }

    if (unique.size() > maxDistinct) {
        std::ostringstream err;
        err << "orientation enumeration produced " << unique.size()
            << " distinct Khovanov homology values for a link with "
            << components.size() << " component(s); expected at most "
            << maxDistinct << ". This usually indicates an invalid PD code, "
            << "connected-sum bug, or component-orientation bug.";
        throw std::runtime_error(err.str());
    }
    return std::vector<std::string>(unique.begin(), unique.end());
}

LinkNameInfo parseLinkName(const std::string& name) {
    static const std::regex re(R"(^(m?)([LK])(\d+)([an])(\d+)$)");
    std::smatch match;
    if (!std::regex_match(name, match, re)) {
        throw std::runtime_error("invalid link name: " + name);
    }
    LinkNameInfo info;
    info.mirror = match[1].str() == "m";
    info.isLink = match[2].str() == "L";
    info.crossing = std::stoi(match[3].str());
    info.nonAlternating = match[4].str() == "n";
    info.index = std::stoi(match[5].str());
    return info;
}

bool linkNameLess(const std::string& a, const std::string& b) {
    LinkNameInfo x = parseLinkName(a);
    LinkNameInfo y = parseLinkName(b);
    return std::tie(x.crossing, x.isLink, x.nonAlternating, x.index, x.mirror) <
           std::tie(y.crossing, y.isLink, y.nonAlternating, y.index, y.mirror);
}

fs::path sourceRoot() {
#ifdef CPP_COM_LINK_GEN_SOURCE_DIR
    return fs::path(CPP_COM_LINK_GEN_SOURCE_DIR);
#else
    return fs::current_path();
#endif
}

fs::path defaultDataRoot() {
    fs::path cwdData = fs::current_path() / "data";
    if (fs::exists(cwdData / "prime_link_knot_10")) return cwdData;
    fs::path sourceData = sourceRoot() / "data";
    return sourceData;
}

std::set<std::string> loadAmphicheiral(const fs::path& dataRoot) {
    std::set<std::string> ans;
    fs::path dir = dataRoot / "prime_link_knot_10" / "amphicheiral";
    if (!fs::is_directory(dir)) throw std::runtime_error("missing amphicheiral data: " + dir.string());
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        std::ifstream in(entry.path());
        if (!in) throw std::runtime_error("cannot read " + entry.path().string());
        std::string line;
        while (std::getline(in, line)) {
            line = trim(line);
            if (!line.empty()) ans.insert(line);
        }
    }
    return ans;
}

std::map<std::string, PDCode> loadPrimePDTable(const fs::path& dataRoot) {
    std::set<std::string> amphicheiral = loadAmphicheiral(dataRoot);
    fs::path dir = dataRoot / "prime_link_knot_10" / "pd_code";
    if (!fs::is_directory(dir)) throw std::runtime_error("missing PD table data: " + dir.string());

    std::map<std::string, PDCode> table;
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        std::ifstream in(entry.path());
        if (!in) throw std::runtime_error("cannot read " + entry.path().string());
        std::string line;
        while (std::getline(in, line)) {
            line = trim(line);
            if (line.empty()) continue;
            size_t colon = line.find(':');
            if (colon == std::string::npos) throw std::runtime_error("bad data line: " + line);
            std::string name = trim(line.substr(0, colon));
            if (startsWith(name, "m") && amphicheiral.count(name.substr(1))) continue;
            table[name] = parsePDCode(line.substr(colon + 1));
            validatePDCode(table[name]);
        }
    }
    return table;
}

std::vector<std::string> sortedPrimeNames(const std::map<std::string, PDCode>& table) {
    std::vector<std::string> names;
    names.reserve(table.size());
    for (const auto& item : table) names.push_back(item.first);
    std::sort(names.begin(), names.end(), linkNameLess);
    return names;
}

std::vector<std::vector<std::string>> allPrimeCombinations(const std::vector<std::string>& names, int totalCrossing) {
    std::vector<std::vector<std::string>> solutions;
    std::vector<std::string> current;
    std::function<void(size_t, int)> dfs = [&](size_t lastPos, int curCross) {
        if (curCross > totalCrossing) return;
        if (!current.empty()) solutions.push_back(current);
        if (curCross >= totalCrossing) return;
        for (size_t i = lastPos; i < names.size(); ++i) {
            int crossing = parseLinkName(names[i]).crossing;
            if (curCross + crossing > totalCrossing) continue;
            current.push_back(names[i]);
            dfs(i, curCross + crossing);
            current.pop_back();
        }
    };
    dfs(0, 0);
    return solutions;
}

std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> possibleEdges(const std::vector<int>& counts) {
    std::vector<std::pair<int, int>> nodes;
    for (size_t i = 0; i < counts.size(); ++i) {
        for (int j = 1; j <= counts[i]; ++j) nodes.push_back({static_cast<int>(i + 1), j});
    }

    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> edges;
    for (size_t i = 0; i < nodes.size(); ++i) {
        for (size_t j = i + 1; j < nodes.size(); ++j) {
            if (nodes[i].first != nodes[j].first) edges.push_back({nodes[i], nodes[j]});
        }
    }
    return edges;
}

bool edgeSetConnects(const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>& edges, int groups) {
    std::vector<std::vector<int>> graph(groups + 1);
    for (const auto& edge : edges) {
        int a = edge.first.first;
        int b = edge.second.first;
        if (a == b) continue;
        graph[a].push_back(b);
        graph[b].push_back(a);
    }
    std::vector<int> stack{1};
    std::vector<char> visited(groups + 1, 0);
    visited[1] = 1;
    while (!stack.empty()) {
        int now = stack.back();
        stack.pop_back();
        for (int nxt : graph[now]) {
            if (!visited[nxt]) {
                visited[nxt] = 1;
                stack.push_back(nxt);
            }
        }
    }
    for (int i = 1; i <= groups; ++i) {
        if (!visited[i]) return false;
    }
    return true;
}

std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> edgeCombinations(
    const std::vector<int>& counts,
    int chooseCount) {
    auto edges = possibleEdges(counts);
    std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> out;
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> current;
    std::function<void(size_t)> dfs = [&](size_t begin) {
        if (static_cast<int>(current.size()) == chooseCount) {
            if (edgeSetConnects(current, static_cast<int>(counts.size()))) out.push_back(current);
            return;
        }
        for (size_t i = begin; i < edges.size(); ++i) {
            current.push_back(edges[i]);
            dfs(i + 1);
            current.pop_back();
        }
    };
    dfs(0);
    return out;
}

std::string makeLinkRepresentation(const std::vector<std::string>& solution,
                                   const std::map<std::string, PDCode>& table,
                                   const std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>& edges,
                                   int totalCrossing,
                                   int maxPrimeCnt) {
    std::ostringstream out;
    out << "// automatically generated by com-link-gen-10 (version:" << kComLinkGenVersion << ")\n";
    out << "// - command: com_link_gen(" << totalCrossing << ", " << maxPrimeCnt << ")\n";

    std::vector<std::string> unique = solution;
    std::sort(unique.begin(), unique.end(), linkNameLess);
    unique.erase(std::unique(unique.begin(), unique.end()), unique.end());
    for (const std::string& name : unique) {
        out << name << ": " << formatPDCode(table.at(name)) << "\n";
    }
    out << "[";
    for (size_t i = 0; i < solution.size(); ++i) {
        if (i) out << ", ";
        out << solution[i];
    }
    out << "]\n";
    for (const auto& edge : edges) {
        out << "L[" << edge.first.first << ", " << edge.first.second << "]#"
            << "L[" << edge.second.first << ", " << edge.second.second << "]\n";
    }
    return out.str();
}

std::vector<std::string> generateAllRepresentations(const fs::path& dataRoot, int totalCrossing, int maxPrimeCnt) {
    if (totalCrossing >= 11) throw std::runtime_error("total_crs >= 11 is not supported");
    if (totalCrossing <= 1) throw std::runtime_error("total_crs must be greater than 1");
    if (maxPrimeCnt <= 0) throw std::runtime_error("max_prime_cnt must be positive");

    auto table = loadPrimePDTable(dataRoot);
    auto names = sortedPrimeNames(table);
    auto solutions = allPrimeCombinations(names, totalCrossing);

    std::vector<std::string> output;
    for (const auto& solution : solutions) {
        if (static_cast<int>(solution.size()) > maxPrimeCnt) continue;
        std::vector<int> componentCounts;
        for (const std::string& name : solution) {
            componentCounts.push_back(static_cast<int>(componentsFromPDCode(table.at(name)).size()));
        }
        auto edgeSets = edgeCombinations(componentCounts, static_cast<int>(solution.size()) - 1);
        for (const auto& edges : edgeSets) {
            output.push_back(makeLinkRepresentation(solution, table, edges, totalCrossing, maxPrimeCnt));
        }
    }
    return output;
}

std::vector<std::string> parseNameList(const std::string& line) {
    std::string s = trim(line);
    if (s.size() < 2 || s.front() != '[' || s.back() != ']') {
        throw std::runtime_error("invalid link set line: " + line);
    }
    s = s.substr(1, s.size() - 2);
    std::vector<std::string> names;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        item = trim(item);
        if (!item.empty()) names.push_back(item);
    }
    return names;
}

LinkRep parseLinkRep(const std::string& text) {
    LinkRep rep;
    std::istringstream in(text);
    std::string line;
    int linkSetLines = 0;
    static const std::regex methodRe(R"(L\[(\d+)\s*,\s*(\d+)\]\s*#\s*L\[(\d+)\s*,\s*(\d+)\])");

    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || startsWith(line, "//")) continue;
        if (line.find(':') != std::string::npos) {
            size_t colon = line.find(':');
            std::string name = trim(line.substr(0, colon));
            rep.definitions[name] = parsePDCode(line.substr(colon + 1));
            validatePDCode(rep.definitions[name]);
        } else if (line.find('#') != std::string::npos) {
            std::smatch match;
            if (!std::regex_search(line, match, methodRe)) {
                throw std::runtime_error("invalid link method line: " + line);
            }
            rep.methods.push_back({
                {std::stoi(match[1].str()), std::stoi(match[2].str())},
                {std::stoi(match[3].str()), std::stoi(match[4].str())},
            });
        } else if (startsWith(line, "[") && !line.empty() && line.back() == ']') {
            ++linkSetLines;
            rep.linkSet = parseNameList(line);
        } else {
            throw std::runtime_error("unrecognized link representation line: " + line);
        }
    }
    if (linkSetLines != 1) throw std::runtime_error("link representation must contain exactly one link set");
    return rep;
}

PDCode linkRepToPDCode(const std::string& text) {
    LinkRep rep = parseLinkRep(text);
    std::vector<std::optional<PDCode>> linkDict;
    std::vector<std::vector<std::pair<int, int>>> repState;

    for (const std::string& name : rep.linkSet) {
        auto found = rep.definitions.find(name);
        if (found == rep.definitions.end()) throw std::runtime_error("undefined link name: " + name);
        linkDict.push_back(found->second);
        std::vector<std::pair<int, int>> states;
        for (const auto& component : componentsFromPDCode(found->second)) {
            states.push_back({static_cast<int>(linkDict.size()) - 1, component.front()});
        }
        repState.push_back(std::move(states));
    }

    if (linkDict.size() == 1) {
        return normalizePDCode(*linkDict.front()).pd;
    }

    for (const auto& method : rep.methods) {
        int aId = method.first.first - 1;
        int aComp = method.first.second - 1;
        int bId = method.second.first - 1;
        int bComp = method.second.second - 1;
        if (aId < 0 || bId < 0 ||
            aId >= static_cast<int>(repState.size()) ||
            bId >= static_cast<int>(repState.size()) ||
            aComp < 0 || bComp < 0 ||
            aComp >= static_cast<int>(repState[aId].size()) ||
            bComp >= static_cast<int>(repState[bId].size())) {
            throw std::runtime_error("link method references an out-of-range component");
        }

        auto [nowA, labelA] = repState[aId][aComp];
        auto [nowB, labelB] = repState[bId][bComp];
        if (nowA == nowB) throw std::runtime_error("link method connects an already merged cluster");
        if (!linkDict[nowA].has_value() || !linkDict[nowB].has_value()) {
            throw std::runtime_error("internal cluster state is inconsistent");
        }

        PDCode pdA = *linkDict[nowA];
        PDCode pdB = *linkDict[nowB];
        linkDict[nowA].reset();
        linkDict[nowB].reset();
        ConnectedSumResult summed = connectedSum(pdA, pdB, labelA, labelB);
        linkDict[nowA] = summed.pd;

        for (auto& perOriginal : repState) {
            for (auto& state : perOriginal) {
                if (state.first == nowA) {
                    state.second = summed.aMap.at(state.second);
                } else if (state.first == nowB) {
                    state.first = nowA;
                    state.second = summed.bMap.at(state.second);
                }
            }
        }
    }

    std::vector<PDCode> remaining;
    for (const auto& item : linkDict) {
        if (item.has_value()) remaining.push_back(*item);
    }
    if (remaining.size() != 1) {
        throw std::runtime_error("link representation did not merge to one connected cluster");
    }
    return normalizePDCode(remaining.front()).pd;
}

fs::path generatedDir(const fs::path& dataRoot, int totalCrossing, int maxPrimeCnt) {
    std::ostringstream name;
    name << "com_link_gen_10-v" << kComLinkGenVersion
         << "-com_link_gen-" << totalCrossing << "-" << maxPrimeCnt;
    return dataRoot / name.str();
}

template <typename Fn>
void parallelFor(size_t count, int jobs, Fn&& fn) {
    if (jobs <= 0) throw std::runtime_error("jobs/process_count must be positive");
    if (count == 0) return;
    int workerCount = std::min<int>(jobs, static_cast<int>(count));
    std::atomic<size_t> next{0};
    std::mutex errorMutex;
    std::exception_ptr firstError = nullptr;
    std::vector<std::thread> workers;
    workers.reserve(workerCount);
    for (int w = 0; w < workerCount; ++w) {
        workers.emplace_back([&]() {
            for (;;) {
                size_t i = next.fetch_add(1);
                if (i >= count) break;
                try {
                    fn(i);
                } catch (...) {
                    std::lock_guard<std::mutex> lock(errorMutex);
                    if (!firstError) firstError = std::current_exception();
                    break;
                }
            }
        });
    }
    for (auto& worker : workers) worker.join();
    if (firstError) std::rethrow_exception(firstError);
}

std::string zeroPaddedIndex(size_t index) {
    std::ostringstream out;
    out.width(7);
    out.fill('0');
    out << index;
    return out.str();
}

fs::path generateAllFiles(const fs::path& dataRoot, int totalCrossing, int maxPrimeCnt, int jobs) {
    fs::path outDir = generatedDir(dataRoot, totalCrossing, maxPrimeCnt);
    fs::create_directories(outDir);
    auto reps = generateAllRepresentations(dataRoot, totalCrossing, maxPrimeCnt);
    std::cout << "Generating " << reps.size() << " files into " << outDir << "\n";

    parallelFor(reps.size(), jobs, [&](size_t i) {
        PDCode pd = preparePDForCppkh(linkRepToPDCode(reps[i]));
        fs::path file = outDir / (zeroPaddedIndex(i + 1) + ".txt");
        std::ofstream out(file, std::ios::binary);
        if (!out) throw std::runtime_error("cannot write " + file.string());
        out << "// PD_CODE: " << formatPDCode(pd) << "\n" << reps[i];
    });
    return outDir;
}

std::vector<fs::path> indexedGeneratedFiles(const fs::path& dir) {
    if (!fs::is_directory(dir)) throw std::runtime_error("not a directory: " + dir.string());
    std::vector<std::pair<int, fs::path>> indexed;
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (!entry.is_regular_file() || entry.path().extension() != ".txt") continue;
        std::string stem = entry.path().stem().string();
        if (stem.empty() || !std::all_of(stem.begin(), stem.end(), [](char c) {
                return std::isdigit(static_cast<unsigned char>(c));
            })) {
            continue;
        }
        indexed.push_back({std::stoi(stem), entry.path()});
    }
    std::sort(indexed.begin(), indexed.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });
    std::vector<fs::path> files;
    for (const auto& item : indexed) files.push_back(item.second);
    return files;
}

std::string readFile(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot read " + path.string());
    std::ostringstream out;
    out << in.rdbuf();
    return out.str();
}

void writeFile(const fs::path& path, const std::string& content) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("cannot write " + path.string());
    out << content;
}

bool hasTxtExtension(const fs::path& path) {
    std::string ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return ext == ".txt";
}

std::vector<fs::path> sortedTxtFiles(const fs::path& dir) {
    if (!fs::is_directory(dir)) throw std::runtime_error("not a directory: " + dir.string());
    std::vector<fs::path> files;
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.is_regular_file() && hasTxtExtension(entry.path())) files.push_back(entry.path());
    }
    std::sort(files.begin(), files.end());
    return files;
}

std::vector<fs::path> sortedClusterLinkFiles(const fs::path& clusterDir) {
    if (!fs::is_directory(clusterDir)) throw std::runtime_error("not a directory: " + clusterDir.string());
    std::vector<fs::path> files;
    for (const auto& entry : fs::recursive_directory_iterator(clusterDir)) {
        if (!entry.is_regular_file() || !hasTxtExtension(entry.path())) continue;
        if (entry.path().filename() == "khovanov.txt") continue;
        files.push_back(entry.path());
    }
    std::sort(files.begin(), files.end());
    return files;
}

std::string removeAsciiSpaces(std::string value) {
    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
    return value;
}

std::string khovanovKeyFromContent(const std::string& content, const fs::path& source) {
    std::vector<std::string> lines;
    std::istringstream in(content);
    std::string line;
    while (std::getline(in, line)) {
        if (line.find("KHOVANOV") == std::string::npos) continue;
        size_t colon = line.find(':');
        if (colon == std::string::npos) throw std::runtime_error("bad KHOVANOV line in " + source.string());
        lines.push_back(normalizeKhovanovZInvariants(trim(line.substr(colon + 1))));
    }
    if (lines.empty()) throw std::runtime_error("KHOVANOV header not found in " + source.string());
    std::sort(lines.begin(), lines.end());

    std::ostringstream out;
    for (size_t i = 0; i < lines.size(); ++i) {
        if (i) out << "#";
        out << lines[i];
    }
    return removeAsciiSpaces(out.str());
}

std::string pdCodeKeyFromContent(const std::string& content, const fs::path& source) {
    std::istringstream in(content);
    std::string line;
    while (std::getline(in, line)) {
        size_t marker = line.find("PD_CODE:");
        if (marker == std::string::npos) continue;
        return formatPDCode(parsePDCode(line.substr(marker + 8)));
    }
    throw std::runtime_error("PD_CODE header not found in " + source.string());
}

std::string khovanovTextFromKey(std::string key) {
    std::replace(key.begin(), key.end(), '#', '\n');
    if (!key.empty() && key.back() != '\n') key.push_back('\n');
    return key;
}

fs::path defaultClusterDirForDataDir(const fs::path& dataDir) {
    fs::path parent = dataDir.has_parent_path() ? dataDir.parent_path() : fs::current_path();
    return parent / "cluster";
}

fs::path defaultReClusterDirForClusterDir(const fs::path& clusterDir) {
    fs::path parent = clusterDir.has_parent_path() ? clusterDir.parent_path() : fs::current_path();
    return parent / "re_cluster";
}

struct ClusterSummary {
    size_t inputFiles = 0;
    size_t uniquePd = 0;
    size_t ignoredDuplicatePd = 0;
    size_t khovanovClasses = 0;
    std::map<int, int> classSizeHistogram;
};

void printHistogram(const std::map<int, int>& histogram) {
    std::cout << "{";
    bool first = true;
    for (const auto& item : histogram) {
        if (!first) std::cout << ", ";
        first = false;
        std::cout << item.first << ": " << item.second;
    }
    std::cout << "}\n";
}

ClusterSummary classifyByKhovanov(const fs::path& dataDir, const fs::path& clusterDir, bool dryRun) {
    std::vector<fs::path> files = sortedTxtFiles(dataDir);
    if (!dryRun) {
        if (fs::exists(clusterDir)) fs::remove_all(clusterDir);
        fs::create_directories(clusterDir);
    }

    std::map<std::string, std::vector<fs::path>> khovanovToFiles;
    std::unordered_map<std::string, std::string> pdToKhovanov;
    ClusterSummary summary;
    summary.inputFiles = files.size();

    for (const fs::path& file : files) {
        std::string content = readFile(file);
        std::string khovanov = khovanovKeyFromContent(content, file);
        std::string pdCode = pdCodeKeyFromContent(content, file);

        auto [it, inserted] = pdToKhovanov.emplace(pdCode, khovanov);
        if (!inserted) {
            if (it->second != khovanov) {
                throw std::runtime_error("same PD_CODE has different Khovanov values: " + file.string());
            }
            ++summary.ignoredDuplicatePd;
            continue;
        }

        khovanovToFiles[khovanov].push_back(file);
        if (!dryRun) {
            fs::path folder = clusterDir / md5Hex(khovanov);
            fs::create_directories(folder);
            writeFile(folder / "khovanov.txt", khovanovTextFromKey(khovanov));
            fs::copy_file(file, folder / file.filename(), fs::copy_options::overwrite_existing);
        }
    }

    summary.uniquePd = pdToKhovanov.size();
    summary.khovanovClasses = khovanovToFiles.size();
    for (const auto& item : khovanovToFiles) {
        ++summary.classSizeHistogram[static_cast<int>(item.second.size())];
    }

    std::cout << "cluster_dir: " << clusterDir << "\n";
    printHistogram(summary.classSizeHistogram);
    std::cout << "ignore_pd_code: " << summary.ignoredDuplicatePd << "\n";
    std::cout << "unique_pd_code: " << summary.uniquePd << "\n";
    std::cout << "khovanov_classes: " << summary.khovanovClasses << "\n";
    return summary;
}

std::string xmlEscape(const std::string& text) {
    std::string out;
    for (char ch : text) {
        switch (ch) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            case '\'': out += "&apos;"; break;
            default: out.push_back(ch); break;
        }
    }
    return out;
}

struct SvgPoint {
    double x = 0.0;
    double y = 0.0;
};

std::string svgNumber(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(1) << value;
    return out.str();
}

constexpr int kDiagramTop = 1 << 0;
constexpr int kDiagramRight = 1 << 1;
constexpr int kDiagramBottom = 1 << 2;
constexpr int kDiagramLeft = 1 << 3;

struct DiagramMatrixBounds {
    int minRow = 0;
    int minCol = 0;
    int maxRow = -1;
    int maxCol = -1;

    bool empty() const {
        return maxRow < minRow || maxCol < minCol;
    }

    int rows() const {
        return empty() ? 0 : maxRow - minRow + 1;
    }

    int cols() const {
        return empty() ? 0 : maxCol - minCol + 1;
    }
};

struct DiagramLayoutCandidate {
    IntMatrix matrix;
    double score = 0.0;
    unsigned int seed = 0;
    int borderSocket = -1;

    DiagramLayoutCandidate(IntMatrix matrixIn,
                           double scoreIn,
                           unsigned int seedIn,
                           int borderSocketIn)
        : matrix(std::move(matrixIn)),
          score(scoreIn),
          seed(seedIn),
          borderSocket(borderSocketIn) {}
};

PDCode renumberPdLabelsForDiagram(const PDCode& pd) {
    validatePDCode(pd);
    std::map<int, int> labelMap;
    for (const Crossing& crossing : pd) {
        for (int label : crossing) labelMap.emplace(label, 0);
    }

    int nextLabel = 1;
    for (auto& item : labelMap) item.second = nextLabel++;

    PDCode normalized = pd;
    for (Crossing& crossing : normalized) {
        for (int& label : crossing) label = labelMap.at(label);
    }
    return normalized;
}

std::vector<int> diagramBorderSocketCandidates(const PDCode& pd) {
    std::vector<int> candidates{-1};
    if (pd.empty()) return candidates;

    std::map<int, std::set<int>> graph;
    for (const Crossing& crossing : pd) {
        graph[crossing[0]].insert(crossing[2]);
        graph[crossing[2]].insert(crossing[0]);
        graph[crossing[1]].insert(crossing[3]);
        graph[crossing[3]].insert(crossing[1]);
    }

    std::set<int> visited;
    for (const auto& item : graph) {
        const int start = item.first;
        if (visited.count(start)) continue;

        int representative = start;
        std::vector<int> stack{start};
        visited.insert(start);
        while (!stack.empty()) {
            const int now = stack.back();
            stack.pop_back();
            representative = std::min(representative, now);
            for (int next : graph[now]) {
                if (!visited.count(next)) {
                    visited.insert(next);
                    stack.push_back(next);
                }
            }
        }
        candidates.push_back(representative);
    }

    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    if (candidates.size() > 5) candidates.resize(5);
    if (std::find(candidates.begin(), candidates.end(), -1) == candidates.end()) {
        candidates.insert(candidates.begin(), -1);
    }
    return candidates;
}

DiagramMatrixBounds diagramMatrixBounds(const IntMatrix& matrix) {
    DiagramMatrixBounds bounds;
    for (int row = 0; row < matrix.getRowCnt(); ++row) {
        for (int col = 0; col < matrix.getColCnt(); ++col) {
            if (matrix.getPos(row, col) == 0) continue;
            if (bounds.empty()) {
                bounds.minRow = bounds.maxRow = row;
                bounds.minCol = bounds.maxCol = col;
            } else {
                bounds.minRow = std::min(bounds.minRow, row);
                bounds.maxRow = std::max(bounds.maxRow, row);
                bounds.minCol = std::min(bounds.minCol, col);
                bounds.maxCol = std::max(bounds.maxCol, col);
            }
        }
    }
    return bounds;
}

int diagramBitCount4(int mask) {
    int count = 0;
    for (int bit : {kDiagramTop, kDiagramRight, kDiagramBottom, kDiagramLeft}) {
        if (mask & bit) ++count;
    }
    return count;
}

std::vector<int> diagramMaskDirections(int mask) {
    std::vector<int> dirs;
    for (int bit : {kDiagramTop, kDiagramRight, kDiagramBottom, kDiagramLeft}) {
        if (mask & bit) dirs.push_back(bit);
    }
    return dirs;
}

int diagramOppositeDirection(int direction) {
    switch (direction) {
        case kDiagramTop: return kDiagramBottom;
        case kDiagramRight: return kDiagramLeft;
        case kDiagramBottom: return kDiagramTop;
        case kDiagramLeft: return kDiagramRight;
        default: return 0;
    }
}

int diagramDirectionIndex(int direction) {
    switch (direction) {
        case kDiagramTop: return 0;
        case kDiagramRight: return 1;
        case kDiagramBottom: return 2;
        case kDiagramLeft: return 3;
        default: return -1;
    }
}

std::pair<int, int> diagramStepCell(int row, int col, int direction) {
    switch (direction) {
        case kDiagramTop: return {row - 1, col};
        case kDiagramRight: return {row, col + 1};
        case kDiagramBottom: return {row + 1, col};
        case kDiagramLeft: return {row, col - 1};
        default: return {row, col};
    }
}

int diagramLineMask(const IntMatrix& matrix, int row, int col) {
    const int val = matrix.getPos(row, col);
    if (val <= 0) return 0;

    int mask = 0;
    const int top = matrix.getPos(row - 1, col);
    const int right = matrix.getPos(row, col + 1);
    const int bottom = matrix.getPos(row + 1, col);
    const int left = matrix.getPos(row, col - 1);
    if (top == val || top < 0) mask |= kDiagramTop;
    if (right == val || right < 0) mask |= kDiagramRight;
    if (bottom == val || bottom < 0) mask |= kDiagramBottom;
    if (left == val || left < 0) mask |= kDiagramLeft;
    return mask;
}

bool isDiagramStraightMask(int mask) {
    return mask == (kDiagramTop | kDiagramBottom) ||
           mask == (kDiagramLeft | kDiagramRight);
}

bool isDiagramCornerMask(int mask) {
    return mask == (kDiagramTop | kDiagramRight) ||
           mask == (kDiagramRight | kDiagramBottom) ||
           mask == (kDiagramBottom | kDiagramLeft) ||
           mask == (kDiagramLeft | kDiagramTop);
}

double scoreDiagramMatrix(const IntMatrix& matrix) {
    const DiagramMatrixBounds bounds = diagramMatrixBounds(matrix);
    if (bounds.empty()) return 0.0;

    const int rows = bounds.rows();
    const int cols = bounds.cols();
    const int area = rows * cols;
    int nonZero = 0;
    int turns = 0;
    int unsupported = 0;
    int endpointPenalty = 0;

    for (int row = bounds.minRow; row <= bounds.maxRow; ++row) {
        for (int col = bounds.minCol; col <= bounds.maxCol; ++col) {
            const int val = matrix.getPos(row, col);
            if (val == 0) continue;
            ++nonZero;
            if (val < 0) continue;

            const int mask = diagramLineMask(matrix, row, col);
            const int degree = diagramBitCount4(mask);
            if (isDiagramCornerMask(mask)) ++turns;
            if (!isDiagramCornerMask(mask) && !isDiagramStraightMask(mask)) ++unsupported;
            if (degree != 2) endpointPenalty += std::abs(degree - 2);
        }
    }

    const int blanks = area - nonZero;
    const int imbalance = std::abs(rows - cols);
    return static_cast<double>(area) * 1000.0 +
           static_cast<double>(blanks) * 25.0 +
           static_cast<double>(turns) * 12.0 +
           static_cast<double>(imbalance) * 5.0 +
           static_cast<double>(unsupported) * 100000.0 +
           static_cast<double>(endpointPenalty) * 50000.0;
}

IntMatrix buildOptimizedDiagramMatrix(const PDCode& normalizedPd) {
    constexpr unsigned int kSeedStart = 42;
    constexpr int kSeedAttemptsPerBorder = 64;
    constexpr int kMinimumAttemptsAfterSuccess = 16;
    constexpr auto kOptimizeTimeBudget = std::chrono::milliseconds(2000);

    const std::string pdInput = formatPDCode(normalizedPd);
    const std::vector<int> borderSockets = diagramBorderSocketCandidates(normalizedPd);
    PdToDiagram2d converter;
    std::optional<DiagramLayoutCandidate> best;
    std::string lastError;
    int attempts = 0;
    const auto started = std::chrono::steady_clock::now();

    for (int borderSocket : borderSockets) {
        for (int seedOffset = 0; seedOffset < kSeedAttemptsPerBorder; ++seedOffset) {
            const unsigned int seed = kSeedStart + static_cast<unsigned int>(seedOffset);
            ++attempts;
            try {
                std::stringstream input(pdInput);
                auto layout = converter.tryConvertOnce(seed, borderSocket, input);
                IntMatrix matrix = std::get<1>(layout);
                const double score = scoreDiagramMatrix(matrix);
                if (!best || score < best->score) {
                    best.emplace(std::move(matrix), score, seed, borderSocket);
                }
            } catch (const std::exception& error) {
                lastError = error.what();
            } catch (...) {
                lastError = "unknown layout error";
            }

            if (best && attempts >= kMinimumAttemptsAfterSuccess &&
                std::chrono::steady_clock::now() - started >= kOptimizeTimeBudget) {
                return best->matrix;
            }
        }
    }

    if (!best) {
        std::string message = "pd-code-to-diagram could not lay out this PD code";
        if (!lastError.empty()) message += ": " + lastError;
        throw std::runtime_error(message);
    }
    return best->matrix;
}

void appendSvgLine(std::ostringstream& svg,
                   double x1,
                   double y1,
                   double x2,
                   double y2) {
    svg << "<line class=\"strand\" x1=\"" << svgNumber(x1)
        << "\" y1=\"" << svgNumber(y1)
        << "\" x2=\"" << svgNumber(x2)
        << "\" y2=\"" << svgNumber(y2) << "\"/>\n";
}

SvgPoint diagramPortPoint(double x, double y, double tile, int direction) {
    const double midX = x + tile / 2.0;
    const double midY = y + tile / 2.0;
    switch (direction) {
        case kDiagramTop: return {midX, y};
        case kDiagramRight: return {x + tile, midY};
        case kDiagramBottom: return {midX, y + tile};
        case kDiagramLeft: return {x, midY};
        default: return {midX, midY};
    }
}

int diagramArcSweepFlag(int entryDirection, int exitDirection) {
    const int entryIndex = diagramDirectionIndex(entryDirection);
    const int exitIndex = diagramDirectionIndex(exitDirection);
    if (entryIndex < 0 || exitIndex < 0) return 0;
    return exitIndex == (entryIndex + 1) % 4 ? 0 : 1;
}

void appendSvgCornerArc(std::ostringstream& svg,
                        double x,
                        double y,
                        double tile,
                        int entryDirection,
                        int exitDirection) {
    const SvgPoint start = diagramPortPoint(x, y, tile, entryDirection);
    const SvgPoint end = diagramPortPoint(x, y, tile, exitDirection);
    const double radius = tile / 2.0;
    const int sweep = diagramArcSweepFlag(entryDirection, exitDirection);
    svg << "<path class=\"strand\" d=\"M " << svgNumber(start.x) << " " << svgNumber(start.y)
        << " A " << svgNumber(radius) << " " << svgNumber(radius)
        << " 0 0 " << sweep
        << " " << svgNumber(end.x) << " " << svgNumber(end.y) << "\"/>\n";
}

void appendSvgRegularTile(std::ostringstream& svg,
                          const IntMatrix& matrix,
                          int row,
                          int col,
                          double x,
                          double y,
                          double tile) {
    const double midX = x + tile / 2.0;
    const double midY = y + tile / 2.0;
    const int mask = diagramLineMask(matrix, row, col);

    switch (mask) {
        case kDiagramTop | kDiagramRight:
            appendSvgCornerArc(svg, x, y, tile, kDiagramTop, kDiagramRight);
            return;
        case kDiagramRight | kDiagramBottom:
            appendSvgCornerArc(svg, x, y, tile, kDiagramRight, kDiagramBottom);
            return;
        case kDiagramBottom | kDiagramLeft:
            appendSvgCornerArc(svg, x, y, tile, kDiagramBottom, kDiagramLeft);
            return;
        case kDiagramLeft | kDiagramTop:
            appendSvgCornerArc(svg, x, y, tile, kDiagramLeft, kDiagramTop);
            return;
        case kDiagramTop | kDiagramBottom:
            appendSvgLine(svg, midX, y, midX, y + tile);
            return;
        case kDiagramLeft | kDiagramRight:
            appendSvgLine(svg, x, midY, x + tile, midY);
            return;
        default:
            break;
    }

    if (mask & kDiagramTop) appendSvgLine(svg, midX, midY, midX, y);
    if (mask & kDiagramRight) appendSvgLine(svg, midX, midY, x + tile, midY);
    if (mask & kDiagramBottom) appendSvgLine(svg, midX, midY, midX, y + tile);
    if (mask & kDiagramLeft) appendSvgLine(svg, midX, midY, x, midY);
}

bool isTraceableDiagramCell(const IntMatrix& matrix, int row, int col) {
    if (matrix.getPos(row, col) <= 0) return false;
    const int mask = diagramLineMask(matrix, row, col);
    return diagramBitCount4(mask) == 2 &&
           (isDiagramStraightMask(mask) || isDiagramCornerMask(mask));
}

bool sameTraceableDiagramArcNeighbor(const IntMatrix& matrix,
                                     int row,
                                     int col,
                                     int direction) {
    const int value = matrix.getPos(row, col);
    if (value <= 0) return false;
    const auto next = diagramStepCell(row, col, direction);
    if (matrix.getPos(next.first, next.second) != value) return false;
    return isTraceableDiagramCell(matrix, next.first, next.second);
}

int otherTraceDirection(const IntMatrix& matrix, int row, int col, int entryDirection) {
    const std::vector<int> dirs = diagramMaskDirections(diagramLineMask(matrix, row, col));
    if (dirs.size() != 2) return 0;
    if (dirs[0] == entryDirection) return dirs[1];
    if (dirs[1] == entryDirection) return dirs[0];
    return 0;
}

struct DiagramTraceCursor {
    int row = 0;
    int col = 0;
    int entryDirection = 0;
};

DiagramTraceCursor rewindDiagramTraceStart(const IntMatrix& matrix,
                                           int row,
                                           int col,
                                           int entryDirection) {
    DiagramTraceCursor cursor{row, col, entryDirection};
    const int guardLimit = std::max(4, matrix.getRowCnt() * matrix.getColCnt() + 4);
    std::set<std::tuple<int, int, int>> seen;

    for (int guard = 0; guard < guardLimit; ++guard) {
        const auto state = std::make_tuple(cursor.row, cursor.col, cursor.entryDirection);
        if (!seen.insert(state).second) break;
        if (!sameTraceableDiagramArcNeighbor(matrix, cursor.row, cursor.col, cursor.entryDirection)) break;

        const auto prev = diagramStepCell(cursor.row, cursor.col, cursor.entryDirection);
        const int connectedSide = diagramOppositeDirection(cursor.entryDirection);
        const int previousEntry = otherTraceDirection(matrix, prev.first, prev.second, connectedSide);
        if (previousEntry == 0) break;
        cursor = DiagramTraceCursor{prev.first, prev.second, previousEntry};
    }

    return cursor;
}

void appendSvgTraceSegment(std::ostringstream& path,
                           double x,
                           double y,
                           double tile,
                           int entryDirection,
                           int exitDirection) {
    const SvgPoint end = diagramPortPoint(x, y, tile, exitDirection);
    if ((entryDirection == kDiagramTop && exitDirection == kDiagramBottom) ||
        (entryDirection == kDiagramBottom && exitDirection == kDiagramTop) ||
        (entryDirection == kDiagramLeft && exitDirection == kDiagramRight) ||
        (entryDirection == kDiagramRight && exitDirection == kDiagramLeft)) {
        path << " L " << svgNumber(end.x) << " " << svgNumber(end.y);
        return;
    }

    const double radius = tile / 2.0;
    const int sweep = diagramArcSweepFlag(entryDirection, exitDirection);
    path << " A " << svgNumber(radius) << " " << svgNumber(radius)
         << " 0 0 " << sweep
         << " " << svgNumber(end.x) << " " << svgNumber(end.y);
}

void appendSvgRegularTracedPaths(std::ostringstream& svg,
                                 const IntMatrix& matrix,
                                 const DiagramMatrixBounds& bounds,
                                 double padding,
                                 double tile,
                                 std::set<std::pair<int, int>>& tracedCells) {
    for (int row = bounds.minRow; row <= bounds.maxRow; ++row) {
        for (int col = bounds.minCol; col <= bounds.maxCol; ++col) {
            if (!isTraceableDiagramCell(matrix, row, col)) continue;
            if (tracedCells.count({row, col})) continue;

            const std::vector<int> dirs = diagramMaskDirections(diagramLineMask(matrix, row, col));
            if (dirs.size() != 2) continue;
            DiagramTraceCursor cursor = rewindDiagramTraceStart(matrix, row, col, dirs[0]);
            if (tracedCells.count({cursor.row, cursor.col})) {
                cursor = rewindDiagramTraceStart(matrix, row, col, dirs[1]);
            }
            if (tracedCells.count({cursor.row, cursor.col})) continue;

            const double startX = padding + static_cast<double>(cursor.col - bounds.minCol) * tile;
            const double startY = padding + static_cast<double>(cursor.row - bounds.minRow) * tile;
            const SvgPoint start = diagramPortPoint(startX, startY, tile, cursor.entryDirection);

            std::ostringstream path;
            path << "M " << svgNumber(start.x) << " " << svgNumber(start.y);
            bool wroteSegment = false;

            const int guardLimit = std::max(4, matrix.getRowCnt() * matrix.getColCnt() + 4);
            for (int guard = 0; guard < guardLimit; ++guard) {
                if (!isTraceableDiagramCell(matrix, cursor.row, cursor.col)) break;
                if (tracedCells.count({cursor.row, cursor.col})) break;
                tracedCells.insert({cursor.row, cursor.col});

                const int exitDirection = otherTraceDirection(matrix, cursor.row, cursor.col, cursor.entryDirection);
                if (exitDirection == 0) break;

                const double x = padding + static_cast<double>(cursor.col - bounds.minCol) * tile;
                const double y = padding + static_cast<double>(cursor.row - bounds.minRow) * tile;
                appendSvgTraceSegment(path, x, y, tile, cursor.entryDirection, exitDirection);
                wroteSegment = true;

                if (!sameTraceableDiagramArcNeighbor(matrix, cursor.row, cursor.col, exitDirection)) break;
                const auto next = diagramStepCell(cursor.row, cursor.col, exitDirection);
                if (tracedCells.count({next.first, next.second})) break;
                cursor = DiagramTraceCursor{
                    next.first,
                    next.second,
                    diagramOppositeDirection(exitDirection),
                };
            }

            if (wroteSegment) {
                svg << "<path class=\"strand\" d=\"" << path.str() << "\"/>\n";
            }
        }
    }
}

void appendSvgCrossingTile(std::ostringstream& svg,
                           int crossingValue,
                           double x,
                           double y,
                           double tile) {
    const double midX = x + tile / 2.0;
    const double midY = y + tile / 2.0;
    const double gapRadius = 6.0;

    if (crossingValue == -1) {
        appendSvgLine(svg, midX, y, midX, y + tile);
        svg << "<circle class=\"gap\" cx=\"" << svgNumber(midX)
            << "\" cy=\"" << svgNumber(midY)
            << "\" r=\"" << svgNumber(gapRadius) << "\"/>\n";
        appendSvgLine(svg, x, midY, x + tile, midY);
    } else if (crossingValue == -2) {
        appendSvgLine(svg, x, midY, x + tile, midY);
        svg << "<circle class=\"gap\" cx=\"" << svgNumber(midX)
            << "\" cy=\"" << svgNumber(midY)
            << "\" r=\"" << svgNumber(gapRadius) << "\"/>\n";
        appendSvgLine(svg, midX, y, midX, y + tile);
    }
}

void appendSvgArcLabel(std::ostringstream& svg,
                       int value,
                       double x,
                       double y,
                       const char* anchor) {
    if (value <= 0) return;
    svg << "<text class=\"arc-label\" x=\"" << svgNumber(x)
        << "\" y=\"" << svgNumber(y)
        << "\" text-anchor=\"" << anchor << "\">"
        << value << "</text>\n";
}

void appendSvgCrossingLabels(std::ostringstream& svg,
                             const IntMatrix& matrix,
                             int row,
                             int col,
                             double x,
                             double y,
                             double tile) {
    constexpr double margin = 3.0;
    constexpr double fontBaseline = 8.0;

    appendSvgArcLabel(svg, matrix.getPos(row - 1, col), x + tile - margin, y - margin, "end");
    appendSvgArcLabel(svg, matrix.getPos(row, col + 1), x + tile + margin, y + margin + fontBaseline, "start");
    appendSvgArcLabel(svg, matrix.getPos(row + 1, col), x + margin, y + tile + margin + fontBaseline, "start");
    appendSvgArcLabel(svg, matrix.getPos(row, col - 1), x - margin, y + tile - margin, "end");
}

std::string renderUnknotSvg(const std::string& title) {
    constexpr int width = 180;
    constexpr int height = 150;
    std::ostringstream svg;
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 " << width << " " << height
        << "\" role=\"img\" aria-label=\"" << xmlEscape(title) << "\">\n";
    svg << "<title>" << xmlEscape(title) << "</title>\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    svg << "<circle cx=\"90\" cy=\"75\" r=\"48\" fill=\"none\" stroke=\"#111827\""
        << " stroke-width=\"4\"/>\n";
    svg << "</svg>\n";
    return svg.str();
}

std::string renderDiagramErrorSvg(const std::string& title, const std::string& error) {
    constexpr int width = 520;
    constexpr int height = 140;
    std::ostringstream svg;
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 " << width << " " << height
        << "\" role=\"img\" aria-label=\"" << xmlEscape(title) << "\">\n";
    svg << "<title>" << xmlEscape(title) << "</title>\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    svg << "<text x=\"20\" y=\"38\" font-family=\"Arial,DejaVu Sans,sans-serif\""
        << " font-size=\"16\" font-weight=\"700\" fill=\"#991b1b\">PD diagram layout failed</text>\n";
    svg << "<text x=\"20\" y=\"70\" font-family=\"Arial,DejaVu Sans,sans-serif\""
        << " font-size=\"12\" fill=\"#374151\">" << xmlEscape(error).substr(0, 220) << "</text>\n";
    svg << "</svg>\n";
    return svg.str();
}

std::string renderPdMatrixSvg(const IntMatrix& matrix, const std::string& title) {
    const DiagramMatrixBounds bounds = diagramMatrixBounds(matrix);
    if (bounds.empty()) return renderUnknotSvg(title);

    constexpr double tile = 30.0;
    constexpr double padding = 14.0;
    const int rows = bounds.rows();
    const int cols = bounds.cols();
    const int width = static_cast<int>(cols * tile + padding * 2.0);
    const int height = static_cast<int>(rows * tile + padding * 2.0);

    std::ostringstream svg;
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 " << width << " " << height
        << "\" role=\"img\" aria-label=\"" << xmlEscape(title) << "\">\n";
    svg << "<title>" << xmlEscape(title) << "</title>\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    svg << "<style>"
        << ".strand{fill:none;stroke:#111827;stroke-width:4;stroke-linecap:butt;"
        << "stroke-linejoin:round;shape-rendering:geometricPrecision}"
        << ".gap{fill:white;stroke:white;stroke-width:0}"
        << ".arc-label{font-family:Arial,DejaVu Sans,sans-serif;font-size:9px;"
        << "font-weight:700;fill:#dc2626;stroke:white;stroke-width:3px;"
        << "paint-order:stroke fill;stroke-linejoin:round}"
        << "</style>\n";

    std::set<std::pair<int, int>> tracedCells;
    appendSvgRegularTracedPaths(svg, matrix, bounds, padding, tile, tracedCells);

    for (int row = bounds.minRow; row <= bounds.maxRow; ++row) {
        for (int col = bounds.minCol; col <= bounds.maxCol; ++col) {
            const int value = matrix.getPos(row, col);
            if (value <= 0) continue;
            if (tracedCells.count({row, col})) continue;
            const double x = padding + static_cast<double>(col - bounds.minCol) * tile;
            const double y = padding + static_cast<double>(row - bounds.minRow) * tile;
            appendSvgRegularTile(svg, matrix, row, col, x, y, tile);
        }
    }

    for (int row = bounds.minRow; row <= bounds.maxRow; ++row) {
        for (int col = bounds.minCol; col <= bounds.maxCol; ++col) {
            const int value = matrix.getPos(row, col);
            if (value != -1 && value != -2) continue;
            const double x = padding + static_cast<double>(col - bounds.minCol) * tile;
            const double y = padding + static_cast<double>(row - bounds.minRow) * tile;
            appendSvgCrossingTile(svg, value, x, y, tile);
        }
    }

    for (int row = bounds.minRow; row <= bounds.maxRow; ++row) {
        for (int col = bounds.minCol; col <= bounds.maxCol; ++col) {
            const int value = matrix.getPos(row, col);
            if (value != -1 && value != -2) continue;
            const double x = padding + static_cast<double>(col - bounds.minCol) * tile;
            const double y = padding + static_cast<double>(row - bounds.minRow) * tile;
            appendSvgCrossingLabels(svg, matrix, row, col, x, y, tile);
        }
    }

    svg << "</svg>\n";
    return svg.str();
}

std::string renderPDCodeSvg(const PDCode& pd, const std::string& title) {
    try {
        if (pd.empty()) return renderUnknotSvg(title);
        PDCode diagramPd = renumberPdLabelsForDiagram(pd);
        const IntMatrix matrix = buildOptimizedDiagramMatrix(diagramPd);
        return renderPdMatrixSvg(matrix, title);
    } catch (const std::exception& error) {
        return renderDiagramErrorSvg(title, error.what());
    }
}

void generateDiagramSvgForFile(const fs::path& txtPath, bool force) {
    fs::path outPath = txtPath;
    outPath.replace_extension(".svg");
    if (fs::exists(outPath) && !force) return;

    std::string content = readFile(txtPath);
    PDCode pd = parsePDCode(pdCodeKeyFromContent(content, txtPath));
    writeFile(outPath, renderPDCodeSvg(pd, txtPath.filename().string()));
}

size_t generateClusterDiagrams(const fs::path& clusterDir, int jobs, bool force) {
    std::vector<fs::path> files = sortedClusterLinkFiles(clusterDir);
    parallelFor(files.size(), jobs, [&](size_t i) {
        generateDiagramSvgForFile(files[i], force);
    });
    std::cout << "diagram_svg: " << files.size() << "\n";
    return files.size();
}

int countClusterLinkFiles(const fs::path& dir) {
    if (!fs::is_directory(dir)) throw std::runtime_error("not a directory: " + dir.string());
    int count = 0;
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (!entry.is_regular_file() || !hasTxtExtension(entry.path())) continue;
        if (entry.path().filename() == "khovanov.txt") continue;
        ++count;
    }
    return count;
}

std::map<int, int> reClusterByClassSize(const fs::path& clusterDir, const fs::path& outDir) {
    if (!fs::is_directory(clusterDir)) throw std::runtime_error("not a directory: " + clusterDir.string());
    if (fs::exists(outDir)) fs::remove_all(outDir);
    fs::create_directories(outDir);

    std::map<int, int> histogram;
    std::vector<fs::path> folders;
    for (const auto& entry : fs::directory_iterator(clusterDir)) {
        if (entry.is_directory()) folders.push_back(entry.path());
    }
    std::sort(folders.begin(), folders.end());

    for (const fs::path& folder : folders) {
        int linkCount = countClusterLinkFiles(folder);
        ++histogram[linkCount];
        std::ostringstream bucketName;
        bucketName << std::setw(3) << std::setfill('0') << linkCount;
        fs::path bucket = outDir / bucketName.str();
        fs::create_directories(bucket);
        fs::copy(folder, bucket / folder.filename(),
                 fs::copy_options::recursive | fs::copy_options::overwrite_existing);
    }

    std::cout << "re_cluster_dir: " << outDir << "\n";
    printHistogram(histogram);
    return histogram;
}

void runPostprocess(const fs::path& dataDir,
                    const fs::path& clusterDir,
                    const fs::path& reClusterDir,
                    int jobs,
                    bool diagrams,
                    bool forceDiagrams) {
    classifyByKhovanov(dataDir, clusterDir, false);
    if (diagrams) generateClusterDiagrams(clusterDir, jobs, forceDiagrams);
    reClusterByClassSize(clusterDir, reClusterDir);
}

void processOneFile(const fs::path& file) {
    if (!fs::is_regular_file(file)) throw std::runtime_error("not a file: " + file.string());
    std::string oldContent = readFile(file);
    if (oldContent.find("KHOVANOV:") != std::string::npos) return;

    std::istringstream in(oldContent);
    std::string line;
    std::optional<PDCode> pd;
    while (std::getline(in, line)) {
        size_t marker = line.find("PD_CODE:");
        if (marker != std::string::npos) {
            pd = parsePDCode(line.substr(marker + 8));
            break;
        }
    }
    if (!pd.has_value()) throw std::runtime_error("PD_CODE header not found in " + file.string());
    std::vector<std::string> khovanov = computeKhovanovAllOrientations(*pd);
    std::ostringstream prefix;
    for (const std::string& item : khovanov) {
        prefix << "// KHOVANOV: " << item << "\n";
    }
    writeFile(file, prefix.str() + oldContent);
}

void processKhovanovParallel(const fs::path& dir, int jobs) {
    auto files = indexedGeneratedFiles(dir);
    std::cout << "Processing " << files.size() << " files in " << dir << "\n";
    parallelFor(files.size(), jobs, [&](size_t i) {
        processOneFile(files[i]);
    });
}

void processKhovanovModulo(const fs::path& dir, int mod, int res) {
    if (mod <= 0) throw std::runtime_error("mod must be positive");
    if (res < 0 || res >= mod) throw std::runtime_error("res must satisfy 0 <= res < mod");
    auto files = indexedGeneratedFiles(dir);
    for (const fs::path& file : files) {
        int index = std::stoi(file.stem().string());
        if (index % mod == res) processOneFile(file);
    }
}

std::optional<std::string> takeOption(std::vector<std::string>& args, const std::string& name) {
    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == name) {
            if (i + 1 >= args.size()) throw std::runtime_error("missing value for " + name);
            std::string value = args[i + 1];
            args.erase(args.begin() + static_cast<std::ptrdiff_t>(i),
                       args.begin() + static_cast<std::ptrdiff_t>(i + 2));
            return value;
        }
    }
    return std::nullopt;
}

int optionInt(std::vector<std::string>& args, const std::string& name, int fallback) {
    auto value = takeOption(args, name);
    return value ? std::stoi(*value) : fallback;
}

fs::path optionPath(std::vector<std::string>& args, const std::string& name, const fs::path& fallback) {
    auto value = takeOption(args, name);
    return value ? fs::path(*value) : fallback;
}

bool takeFlag(std::vector<std::string>& args, const std::string& name) {
    auto it = std::find(args.begin(), args.end(), name);
    if (it == args.end()) return false;
    args.erase(it);
    return true;
}

void usage() {
    std::cout
        << "Usage:\n"
        << "  cpp_com_link_gen generate --total-crs 10 --max-prime-cnt 3 --jobs 8 [--data-root data]\n"
        << "  cpp_com_link_gen khovanov --dir DIR --jobs 8\n"
        << "  cpp_com_link_gen postprocess --dir DIR --jobs 8 [--cluster-dir DIR] [--re-cluster-dir DIR] [--no-diagrams]\n"
        << "  cpp_com_link_gen classify --dir DIR [--cluster-dir DIR] [--dry-run]\n"
        << "  cpp_com_link_gen diagrams --cluster-dir DIR --jobs 8 [--force]\n"
        << "  cpp_com_link_gen re-cluster --cluster-dir DIR [--out-dir DIR]\n"
        << "  cpp_com_link_gen all --total-crs 10 --max-prime-cnt 3 --jobs 8 [--data-root data]\n"
        << "  cpp_com_link_gen default --jobs 8 [--total-crs 10 --max-prime-cnt 3]\n"
        << "  cpp_com_link_gen process-one FILE\n"
        << "  cpp_com_link_gen legacy --dir DIR --mod M --res R\n"
        << "  cpp_com_link_gen pd --file LINK_REP.txt\n"
        << "  cpp_com_link_gen svg (--pd \"[[1,5,2,4],...]\" | --file generated.txt) --out diagram.svg\n"
        << "  cpp_com_link_gen kh --pd \"[[1,5,2,4],...]\"\n"
        << "  cpp_com_link_gen kh-all-orientations --pd \"[[1,5,2,4],...]\"\n";
}

int runCommand(std::vector<std::string> args) {
    if (args.empty()) {
        int jobs = 1;
        std::cout << "process_count>>>";
        std::cin >> jobs;
        processKhovanovParallel(generatedDir(defaultDataRoot(), 10, 3), jobs);
        return 0;
    }
    std::string command = args.front();
    args.erase(args.begin());
    if (command == "--help" || command == "-h" || command == "help") {
        usage();
        return 0;
    }

    fs::path dataRoot = optionPath(args, "--data-root", defaultDataRoot());

    if (command == "generate") {
        int total = optionInt(args, "--total-crs", 10);
        int maxPrime = optionInt(args, "--max-prime-cnt", 3);
        int jobs = optionInt(args, "--jobs", 1);
        fs::path out = generateAllFiles(dataRoot, total, maxPrime, jobs);
        std::cout << out << "\n";
    } else if (command == "khovanov") {
        fs::path dir = optionPath(args, "--dir", generatedDir(dataRoot, 10, 3));
        int jobs = optionInt(args, "--jobs", 1);
        processKhovanovParallel(dir, jobs);
    } else if (command == "postprocess") {
        fs::path dir = optionPath(args, "--dir", generatedDir(dataRoot, 10, 3));
        int jobs = optionInt(args, "--jobs", 1);
        fs::path clusterDir = optionPath(args, "--cluster-dir", defaultClusterDirForDataDir(dir));
        fs::path reClusterDir = optionPath(args, "--re-cluster-dir", defaultReClusterDirForClusterDir(clusterDir));
        bool noDiagrams = takeFlag(args, "--no-diagrams");
        bool forceDiagrams = takeFlag(args, "--force-diagrams") || takeFlag(args, "--force");
        runPostprocess(dir, clusterDir, reClusterDir, jobs, !noDiagrams, forceDiagrams);
    } else if (command == "classify") {
        fs::path dir = optionPath(args, "--dir", generatedDir(dataRoot, 10, 3));
        fs::path clusterDir = optionPath(args, "--cluster-dir", defaultClusterDirForDataDir(dir));
        bool dryRun = takeFlag(args, "--dry-run");
        classifyByKhovanov(dir, clusterDir, dryRun);
    } else if (command == "diagrams") {
        fs::path clusterDir = optionPath(
            args, "--cluster-dir", defaultClusterDirForDataDir(generatedDir(dataRoot, 10, 3)));
        int jobs = optionInt(args, "--jobs", 1);
        bool force = takeFlag(args, "--force");
        generateClusterDiagrams(clusterDir, jobs, force);
    } else if (command == "re-cluster" || command == "recluster") {
        fs::path clusterDir = optionPath(
            args, "--cluster-dir", defaultClusterDirForDataDir(generatedDir(dataRoot, 10, 3)));
        fs::path outDir = optionPath(args, "--out-dir", defaultReClusterDirForClusterDir(clusterDir));
        reClusterByClassSize(clusterDir, outDir);
    } else if (command == "default") {
        int total = optionInt(args, "--total-crs", 10);
        int maxPrime = optionInt(args, "--max-prime-cnt", 3);
        int jobs = optionInt(args, "--jobs", 1);
        processKhovanovParallel(generatedDir(dataRoot, total, maxPrime), jobs);
    } else if (command == "all") {
        int total = optionInt(args, "--total-crs", 10);
        int maxPrime = optionInt(args, "--max-prime-cnt", 3);
        int jobs = optionInt(args, "--jobs", 1);
        fs::path out = generateAllFiles(dataRoot, total, maxPrime, jobs);
        processKhovanovParallel(out, jobs);
    } else if (command == "process-one") {
        if (args.empty()) throw std::runtime_error("process-one needs a file path");
        fs::path file = args.front();
        args.erase(args.begin());
        processOneFile(file);
    } else if (command == "legacy") {
        fs::path dir = optionPath(args, "--dir", generatedDir(dataRoot, 10, 3));
        int mod = optionInt(args, "--mod", 1);
        int res = optionInt(args, "--res", 0);
        processKhovanovModulo(dir, mod, res);
    } else if (command == "pd") {
        fs::path file = optionPath(args, "--file", {});
        if (file.empty()) throw std::runtime_error("pd needs --file");
        std::cout << formatPDCode(linkRepToPDCode(readFile(file))) << "\n";
    } else if (command == "svg") {
        auto pdText = takeOption(args, "--pd");
        fs::path file = optionPath(args, "--file", {});
        fs::path out = optionPath(args, "--out", {});
        if (!pdText && file.empty()) throw std::runtime_error("svg needs --pd or --file");
        if (pdText && !file.empty()) throw std::runtime_error("svg accepts only one of --pd and --file");
        if (out.empty()) throw std::runtime_error("svg needs --out");

        PDCode pd = pdText ? parsePDCode(*pdText)
                           : parsePDCode(pdCodeKeyFromContent(readFile(file), file));
        const std::string title = file.empty() ? out.filename().string() : file.filename().string();
        writeFile(out, renderPDCodeSvg(pd, title));
    } else if (command == "kh") {
        auto pdText = takeOption(args, "--pd");
        if (!pdText) throw std::runtime_error("kh needs --pd");
        std::cout << computeKhovanovSinglePD(parsePDCode(*pdText)) << "\n";
    } else if (command == "kh-all-orientations") {
        auto pdText = takeOption(args, "--pd");
        if (!pdText) throw std::runtime_error("kh-all-orientations needs --pd");
        for (const std::string& item : computeKhovanovAllOrientations(parsePDCode(*pdText))) {
            std::cout << item << "\n";
        }
    } else {
        throw std::runtime_error("unknown command: " + command);
    }

    if (!args.empty()) {
        throw std::runtime_error("unused argument: " + args.front());
    }
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    try {
        std::vector<std::string> args;
        for (int i = 1; i < argc; ++i) args.emplace_back(argv[i]);
        return runCommand(std::move(args));
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}

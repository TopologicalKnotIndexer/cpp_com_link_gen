#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <condition_variable>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
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
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

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
        }
        if (next.at(crossing[2]) == crossing[0]) {
            crossing = {crossing[2], crossing[3], crossing[0], crossing[1]};
            continue;
        }
        throw std::runtime_error("crossing is inconsistent with component orientation");
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

std::vector<std::string> computeKhovanovAllOrientations(const PDCode& pd) {
    validatePDCode(pd);
    std::vector<std::vector<int>> components = componentsFromPDCode(pd);
    if (components.size() >= 63) throw std::runtime_error("too many components to enumerate orientations");
    uint64_t orientationCount = components.empty() ? 1ULL : (1ULL << components.size());
    uint64_t maxDistinct = components.empty() ? 1ULL : (1ULL << (components.size() - 1));

    std::string pdText = formatKnotTheoryPD(pd);
    std::string signsText = orientationSignsDocument(pd, components, orientationCount);

    static std::mutex cppkhMutex;
    std::string rawOutput;
    {
        std::lock_guard<std::mutex> lock(cppkhMutex);
        rawOutput = takeCppkhString(cppkh_compute_pd_signed_variants_ex(pdText.c_str(), signsText.c_str(), 1));
    }

    std::set<std::string> unique;
    for (std::string line : splitLines(rawOutput)) {
        line = trim(line);
        if (!line.empty()) unique.insert(line);
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
        PDCode pd = linkRepToPDCode(reps[i]);
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

void usage() {
    std::cout
        << "Usage:\n"
        << "  cpp_com_link_gen generate --total-crs 10 --max-prime-cnt 3 --jobs 8 [--data-root data]\n"
        << "  cpp_com_link_gen khovanov --dir DIR --jobs 8\n"
        << "  cpp_com_link_gen all --total-crs 10 --max-prime-cnt 3 --jobs 8 [--data-root data]\n"
        << "  cpp_com_link_gen default --jobs 8 [--total-crs 10 --max-prime-cnt 3]\n"
        << "  cpp_com_link_gen process-one FILE\n"
        << "  cpp_com_link_gen legacy --dir DIR --mod M --res R\n"
        << "  cpp_com_link_gen pd --file LINK_REP.txt\n"
        << "  cpp_com_link_gen kh --pd \"[[1,5,2,4],...]\"\n";
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
    } else if (command == "kh") {
        auto pdText = takeOption(args, "--pd");
        if (!pdText) throw std::runtime_error("kh needs --pd");
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

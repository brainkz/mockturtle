#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <tuple>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <stack>
#include <array>
#include <fmt/format.h>
#include <thread>

typedef unsigned short US;
typedef unsigned int   UI;
typedef unsigned long  UL;
typedef const unsigned short CUS;
typedef const unsigned int   CUI;
typedef const unsigned long  CUL;
typedef unsigned long long ULL;

constexpr UI NUM_VARS = 4u;
typedef unsigned int TT;
// typedef kitty::static_truth_table<NUM_IN> TT;

constexpr US fDFF   = 0;
constexpr US fNOT   = 1;
constexpr US fMERGE = 2;
constexpr US fOR    = 3;
constexpr US fAND   = 4;
constexpr US fXOR   = 5;
constexpr US fOR3   = 6;
constexpr US fAND3  = 7;
constexpr US fMAJ3  = 8;
constexpr US fCB    = 9;
constexpr US fSPL   = 10;
constexpr US fPI    = 11;
constexpr US fNOFUNC= 99;

const UI kNumThreads = 100;


// constexpr bool accel_cost = true;
#define accel_cost false

constexpr std::array<UI,12> COSTS = {7, 9, 8, 8, 8, 7, 11, 11, 11, 8, 7, 0};
std::unordered_map<US, std::string> F2STR { 
    {fDFF   , "DFF"},
    {fNOT   , "NOT"},
    {fMERGE , "MRG"},
    {fOR    , "OR "},
    {fAND   , "AND"},
    {fXOR   , "XOR"},
    {fOR3   , "OR3"},
    {fAND3  , "AND3"},
    {fMAJ3  , "MAJ3"},
    {fCB    , "CB "},
    {fSPL   , "SPL"},
    {fPI    , "PI "},
    {fNOFUNC, "N/A"}
    }; 

constexpr UL NUM_TT = (1 << (1 << NUM_VARS));
constexpr TT ONES = NUM_TT - 1;
constexpr UL INF = 0xFFFFFF;

// Hash combiner
template <typename T>
static void hash_combine(std::size_t& seed, const T& val) {
    seed ^= std::hash<T>{}(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

ULL calculate_hash(TT func, US last_func, UI cost, UI depth, bool xorable, std::vector<ULL> parent_hashes = {})
{
    std::size_t seed = 0;
    hash_combine(seed, func);
    hash_combine(seed, last_func);
    hash_combine(seed, cost);
    hash_combine(seed, depth);
    hash_combine(seed, xorable);
    for (const auto parent_hash : parent_hashes) {
        hash_combine(seed, parent_hash);
    }
    return seed;
};

class Node {
public:
    TT func = 0;
    US last_func = fNOFUNC;
    UI cost = INF;
    UI depth = INF;
    bool xorable = false;
    std::vector<ULL> parent_hashes;
    UI lvl = INF;
    ULL hash = 0;

    #if accel_cost
    std::unordered_map<ULL, UI> predecessor_count;
    std::unordered_set<ULL> non_splittable_pred;
    #endif
    Node() = default;
    Node(const Node& other) : func(other.func), last_func(other.last_func), cost(other.cost), depth(other.depth), xorable(other.xorable), parent_hashes(other.parent_hashes), hash(other.hash), lvl(other.lvl) {}

    Node& operator=(const Node& other) {
        if (this != &other) {
            func = other.func;
            last_func = other.last_func;
            cost = other.cost;
            depth = other.depth;
            xorable = other.xorable;
            parent_hashes = other.parent_hashes;
            hash = other.hash;
            lvl = other.lvl;
        }
        return *this;
    }

    Node(TT _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes(_parent_hashes), lvl(depth / 3)
    {
        // Calculate hash based on the hashes of parent_hashes and specified fields
        hash = calculate_hash(func, last_func, cost, depth, xorable, parent_hashes);
        // hash = calculate_hash();
    }

    Node(TT _func, US _last_func, UI _cost, UI _depth, bool _xorable)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes{}, lvl(depth / 3)
    {
        // Calculate hash based on the hashes of parent_hashes and specified fields
        // It is assumed the node has no parents
        hash = calculate_hash(func, last_func, cost, depth, xorable);
    }

    Node(TT _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes, ULL _hash)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes(_parent_hashes), lvl(depth / 3), hash(_hash)
    {
        // Hash is precalculated based on the hashes of parent_hashes and specified fields
    }

    ~Node() { }

    // Equality operator
    bool operator==(const Node& other) const {
        return hash == other.hash;
    }

    ULL get_hash() {
        return hash;
    }

    std::string to_str() const
    {
        // std::string s;
        // for (auto & ph : parent_hashes)
        // {
        //     Node & ref = GNM[ph];
        //     s += fmt::format(" {} ", ref.func);
        // }
        return fmt::format("{:016b}{} | {} | c{}, d{}, l{} |", func, (xorable?'x':' '), F2STR[(US)last_func], cost, depth, lvl);
        // return fmt::format("{:04x}{} | {} | {}, {} | ", func, (xorable?'x':' '), F2STR[(US)last_func], cost, depth);
        // return "Func: " + std::to_string(func) + "|Last: " + std::to_string(last_func) + "|Cost: " + std::to_string(cost) + "|Depth: " + std::to_string(depth) + "|X: " + std::to_string(xorable);
    }

private:

};

namespace std {
    template<> struct hash<Node> {
        std::size_t operator()(const Node& node) const {
            return std::hash<ULL>{}(node.hash);
        }
    };
}


#pragma region write_output

void write_csv_gnm(const std::unordered_map<ULL, Node>& gnm, const std::string& filename) {
    // Open output file
    std::ofstream outfile(filename);

    // Write header row to CSV file
    outfile << "Hash,Func,Last Func,Cost,Depth,Xorable,Parent Hashes,Lvl" << std::endl;

    // Write data to CSV file
    for (const auto& [hash, node] : gnm) {
        outfile << hash << ",";
        outfile << node.func << ",";
        outfile << node.last_func << ",";
        outfile << node.cost << ",";
        outfile << node.depth << ",";
        outfile << node.xorable << ",";
        for (const auto& parent_hash : node.parent_hashes) {
            outfile << parent_hash << "|";
        }
        outfile << ",";
        outfile << node.lvl << std::endl;
    }

    // Close output file
    outfile.close();
}

void write_csv_arr(const std::array<ULL, NUM_TT>& arr_hashes, const std::string& filename) {
    // Open output file
    std::ofstream outfile(filename);

    // Write header row to CSV file
    outfile << "Hash" << std::endl;

    // Write data to CSV file
    for (const auto& hash : arr_hashes) {
        outfile << hash << std::endl;
    }

    // Close output file
    outfile.close();
}

std::unordered_map<ULL, Node> read_csv_gnm(const std::string& filename) {
    // Open input file
    std::ifstream infile(filename);
    
    std::unordered_map<ULL, Node> gnm;

    // Parse CSV file and populate GNM variable
    std::string line;
    std::getline(infile, line);  // skip header row
    while (std::getline(infile, line)) {
        std::stringstream ss;
        ss.str(line);
        std::string field;
        ULL hash;
        Node node;
        std::getline(ss, field, ',');
        std::stringstream(field) >> hash;
        std::getline(ss, field, ',');
        std::stringstream(field) >> node.func;
        std::getline(ss, field, ',');
        std::stringstream(field) >> node.last_func;
        std::getline(ss, field, ',');
        std::stringstream(field) >> node.cost;
        std::getline(ss, field, ',');
        std::stringstream(field) >> node.depth;
        std::getline(ss, field, ',');
        std::stringstream(field) >> node.xorable;
        std::getline(ss, field, ',');
        std::stringstream parent_hashes_ss(field);
        while (std::getline(parent_hashes_ss, field, '|')) {
            ULL parent_hash;
            std::stringstream(field) >> parent_hash;
            node.parent_hashes.push_back(parent_hash);
        }
        std::getline(ss, field, ',');
        std::stringstream(field) >> node.lvl;
        gnm[hash] = node;
    }

    // Close input file
    infile.close();
    return gnm;
}

std::array<ULL, NUM_TT> read_csv_arr(const std::string& filename) {
    // Open input file
    std::ifstream infile(filename);
    std::array<ULL, NUM_TT> arr_hashes;

    // Parse CSV file and populate array
    std::string line;
    std::getline(infile, line);  // skip header row
    UI index = 0;
    while (std::getline(infile, line) && index < NUM_TT) {
        std::stringstream ss(line);
        std::string field;
        std::getline(ss, field, ',');
        std::stringstream(field) >> arr_hashes[index++];
    }

    // Close input file
    infile.close();
    return arr_hashes;
}

#pragma endregion

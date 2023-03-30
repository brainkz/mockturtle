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
#include <kitty/kitty.hpp>


typedef unsigned short US;
typedef unsigned int   UI;
typedef unsigned long  UL;
typedef const unsigned short CUS;
typedef const unsigned int   CUI;
typedef const unsigned long  CUL;
typedef unsigned long long ULL;

bool DEBUG = false;

constexpr UI NUM_VARS = 4u;
typedef kitty::static_truth_table<NUM_VARS> staticTT;

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

constexpr UI kNumThreads = 100;

const std::string GENLIB_PHASE = "UNKNOWN";
constexpr float GENLIB_INPUT_LOAD = 1;
constexpr float GENLIB_MAX_LOAD = 999;
constexpr float GENLIB_RISE_BLOCK_DELAY   = 0.025;
constexpr float GENLIB_RISE_FANOUT_DELAY  = 0.025;
constexpr float GENLIB_FALL_BLOCK_DELAY   = 0.025;
constexpr float GENLIB_FALL_FANOUT_DELAY  = 0.025;

std::vector<UI> PI {0x5555, 0x3333, 0x0F0F, 0x00FF};

std::unordered_map<UI, std::string> eqn_dict;
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
std::unordered_map<US, std::string> PI2LETTER { 
    {0x00FF, "d"},
    {0x0F0F, "c"},
    {0x3333, "b"},
    {0x5555, "a"},
    {0x0000, "0"},
    {0xFFFF, "1"},
    }; 

constexpr UL NUM_TT = (1 << (1 << NUM_VARS));
constexpr UI ONES = NUM_TT - 1;
constexpr UL INF = 0xFFFFFF;

// Hash combiner
template <typename T>
static void hash_combine(std::size_t& seed, const T& val) {
    seed ^= std::hash<T>{}(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

ULL calculate_hash(UI func, US last_func, UI cost, UI depth, bool xorable, std::vector<ULL> parent_hashes = {})
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

// inline void apply_swaps(staticTT & tt, const std::vector<std::pair<uint8_t,uint8_t>> & swaps)
// {
//     for (auto [i,j] : swaps)
//     {
//         kitty::swap_inplace( tt, i, j );
//     }
// }
inline UI apply_swaps(UI func, const std::vector<std::pair<uint8_t,uint8_t>> & swaps)
{
    std::vector<ULL> words = {func};
    staticTT tt;
    kitty::create_from_words(tt, words.begin(), words.end());
    for (auto [i,j] : swaps)
    {
        kitty::swap_inplace( tt, i, j );
    }
    return tt._bits;
}

class Node {
public:
    UI func = 0;
    US last_func = fNOFUNC;
    UI cost = INF;
    UI depth = INF;
    bool xorable = false;
    std::vector<ULL> parent_hashes;
    UI lvl = INF;
    ULL hash = 0;

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

    Node(UI _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes(_parent_hashes), lvl(depth / 3)
    {
        // Calculate hash based on the hashes of parent_hashes and specified fields
        hash = calculate_hash(func, last_func, cost, depth, xorable, parent_hashes);
        // hash = calculate_hash();
    }

    Node(UI _func, US _last_func, UI _cost, UI _depth, bool _xorable)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes{}, lvl(depth / 3)
    {
        // Calculate hash based on the hashes of parent_hashes and specified fields
        // It is assumed the node has no parents
        hash = calculate_hash(func, last_func, cost, depth, xorable);
    }

    Node(UI _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes, ULL _hash)
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
        return fmt::format("{:>20}: {:016b}{} | {} | c{}, d{}, l{} |", hash, func, (xorable?'x':' '), F2STR[(US)last_func], cost, depth, lvl);
        // return fmt::format("{:04x}{} | {} | {}, {} | ", func, (xorable?'x':' '), F2STR[(US)last_func], cost, depth);
        // return "Func: " + std::to_string(func) + "|Last: " + std::to_string(last_func) + "|Cost: " + std::to_string(cost) + "|Depth: " + std::to_string(depth) + "|X: " + std::to_string(xorable);
    }

    // std::string genlib_eqn(std::unordered_map<ULL, Node> & nodemap, std::vector<std::string> & symbol) const
    std::string genlib_eqn(std::unordered_map<ULL, Node> & nodemap, std::vector<std::string> & symbol, std::vector<UI> & indices) const
    {
        // fmt::print("\t\tAccessing genlib_eqn: {}\n", to_str());
        if (last_func == fPI)
        {
            // std::vector<std::string> _symbol_tmp = { "d", "c", "b", "a" };
            // std::vector<std::string> _symbol = { "d", "c", "b", "a" };
            // for (auto idx : indices)
            // {
            //     _symbol.push_back();
            // }
            switch (func)
            {
                case 0x5555: return symbol[indices[3]];
                case 0x3333: return symbol[indices[2]];
                case 0x0F0F: return symbol[indices[1]];
                case 0x00FF: return symbol[indices[0]];
                // case 0x5555: return _symbol[indices[3]];
                // case 0x3333: return _symbol[indices[2]];
                // case 0x0F0F: return _symbol[indices[1]];
                // case 0x00FF: return _symbol[indices[0]];
                // case 0x5555: return "a";
                // case 0x3333: return "b";
                // case 0x0F0F: return "c";
                // case 0x00FF: return "d";
                case 0x0000: return "0";
                case 0xFFFF: return "1";
            }
        }
        else if (last_func == fDFF)
        {
            assert(parent_hashes.size() == 1);
            return fmt::format("{}", nodemap[parent_hashes.back()].genlib_eqn(nodemap, symbol, indices));
        }
        else if (last_func == fNOT)
        {
            assert(parent_hashes.size() == 1);
            return fmt::format("!{}", nodemap[parent_hashes.back()].genlib_eqn(nodemap, symbol, indices));
        }
        else if (last_func == fCB)
        {
            assert(parent_hashes.size() == 2);
            return fmt::format("({0}+{1})", nodemap[parent_hashes.front()].genlib_eqn(nodemap, symbol, indices), nodemap[parent_hashes.back()].genlib_eqn(nodemap, symbol, indices));
        }
        else if (last_func == fOR || last_func == fMERGE)
        {
            assert(parent_hashes.size() == 2);
            return fmt::format("({0}|{1})", nodemap[parent_hashes.front()].genlib_eqn(nodemap, symbol, indices), nodemap[parent_hashes.back()].genlib_eqn(nodemap, symbol, indices));
        }
        else if (last_func == fAND)
        {
            assert(parent_hashes.size() == 2);
            return fmt::format("({0}&{1})", nodemap[parent_hashes.front()].genlib_eqn(nodemap, symbol, indices), nodemap[parent_hashes.back()].genlib_eqn(nodemap, symbol, indices));
        }
        else if (last_func == fXOR)
        {
            assert(parent_hashes.size() == 2);
            return fmt::format("({0}^{1})", nodemap[parent_hashes.front()].genlib_eqn(nodemap, symbol, indices), nodemap[parent_hashes.back()].genlib_eqn(nodemap, symbol, indices));
        }
        else
        {
            if (DEBUG) {fmt::print("Unsupported function {}", to_str());}
            return "";
        }
    }

    // std::string to_genlib(std::unordered_map<ULL, Node> & nodemap, const std::vector<UI> & levels, const std::vector<UI> PI_funcs, const std::vector<std::pair<uint8_t,uint8_t>> & swaps = {}) const
    // {
    //     std::vector<UI> pis;
    //     std::string str = fmt::format("GATE 0x{:04x}_{} {} O={};\n", func, fmt::join(levels, ""), cost,  genlib_eqn(nodemap, pis));

    //     for (auto & pi : pis)
    //     {
    //         UL idx = std::find(PI_funcs.begin(), PI_funcs.end(), pi) - PI_funcs.begin();
    //         auto true_lvl = (depth + 1) / 3;
    //         std::string line = fmt::format("\tPIN {} {} {} {} {:d} {:0.3f} {:d} {:0.3f}\n", 
    //         PI2LETTER[pi], GENLIB_PHASE, GENLIB_INPUT_LOAD, GENLIB_MAX_LOAD, 
    //         true_lvl - levels[idx], GENLIB_RISE_FANOUT_DELAY, true_lvl - levels[idx], GENLIB_FALL_FANOUT_DELAY);
    //         str.append(line);
    //         // if (DEBUG) {fmt::print(line);}
    //     }
    //     // if (DEBUG) {fmt::print(str);}
    //     return str;
    // }

    std::string to_stack(std::unordered_map<ULL, Node> & nodemap, std::vector<std::string> & symbol)
    {
        if (last_func == fPI)
        {
            assert(func == 0x5555 || func == 0x3333 || func == 0x0F0F || func == 0x00FF || func == 0 || func == 0xFFFF);
            switch (func)
            {
                case 0x5555: return symbol[0];
                case 0x3333: return symbol[1];
                case 0x0F0F: return symbol[2];
                case 0x00FF: return symbol[3];
                case 0x0000: return "0";
                case 0xFFFF: return "1";
            }
        }
        else if (last_func == fDFF)
        {
            assert(parent_hashes.size() == 1);
            return fmt::format("DFF({})", nodemap[parent_hashes.back()].to_stack(nodemap, symbol));
        }
        else if (last_func == fNOT)
        {
            assert(parent_hashes.size() == 1);
            return fmt::format("NOT({})", nodemap[parent_hashes.back()].to_stack(nodemap, symbol));
        }
        else if (last_func == fCB || last_func == fMERGE)
        {
            assert(parent_hashes.size() == 2);
            return fmt::format("CB({0}, {1})", nodemap[parent_hashes.front()].to_stack(nodemap, symbol), nodemap[parent_hashes.back()].to_stack(nodemap, symbol));
        }
        else if ( last_func == fOR )
        {
            assert(parent_hashes.size() == 2);
            return fmt::format("OR({0}, {1})", nodemap[parent_hashes.front()].to_stack(nodemap, symbol), nodemap[parent_hashes.back()].to_stack(nodemap, symbol));
        }
        else if (last_func == fAND)
        {
            assert(parent_hashes.size() == 2);
            return fmt::format("AND({0}, {1})", nodemap[parent_hashes.front()].to_stack(nodemap, symbol), nodemap[parent_hashes.back()].to_stack(nodemap, symbol));
        }
        else if (last_func == fXOR)
        {
            assert(parent_hashes.size() == 2);
            return fmt::format("XOR({0}, {1})", nodemap[parent_hashes.front()].to_stack(nodemap, symbol), nodemap[parent_hashes.back()].to_stack(nodemap, symbol));
        }
        else
        {
            // if (DEBUG) {
            fmt::print("Unsupported function {}", to_str());
                // }
            return "";
        }
    }

    std::tuple<bool, UI> redundancy_check(std::unordered_map<ULL, Node> & GNM)
    {
        // std::vector<UI>   pi_funcs     {0x00FF, 0x0F0F, 0x3333, 0x5555};
        std::vector<UI>   pi_funcs     {0x5555, 0x3333, 0x0F0F, 0x00FF};
        std::vector<bool> has_dff      {false, false, false, false};
        std::vector<bool> has_other    {false, false, false, false};
        std::vector<bool> is_reached   {false, false, false, false};
        std::vector<ULL> stack = parent_hashes;
        if (DEBUG) {fmt::print("\tAnalyzing stack for node {}\n", to_str());}
        while (!stack.empty())
        {
            ULL n_hash = stack.back();
            stack.pop_back();
            Node& n = GNM.at(n_hash); // [n_hash];
            if (DEBUG) {fmt::print("\t\tAnalyzing node: {}\n", n.to_str());}
            auto it = std::find(pi_funcs.begin(), pi_funcs.end(), n.func);
            if (it != pi_funcs.end()) // if the function is a PI
            {
                auto idx = it - pi_funcs.begin();
                if (DEBUG) {fmt::print("\t\tFound PI at idx {} for {:04x}\n", idx, n.func);}
                // if (DEBUG) {fmt::print("\t\tLast_func == {}\n", n.last_func);}
                // if (DEBUG) {fmt::print("\t\fDFF == {}\n", fDFF);}
                if (n.last_func == fDFF) 
                {
                    has_dff[idx] = true;
                }
                else if (n.last_func == fPI)
                {
                    is_reached[idx] = true;
                }
                else 
                {
                    has_other[idx] = true;
                }
                if (DEBUG) {fmt::print("\t\t\tNew has_dff:\t{}\n", fmt::join(has_dff, "\t"));}
                if (DEBUG) {fmt::print("\t\t\tNew has_other:\t{}\n", fmt::join(has_other, "\t"));}
                if (DEBUG) {fmt::print("\t\t\tNew is_reached:\t{}\n", fmt::join(is_reached, "\t"));}
            }
            /*
                if (n.parent_hashes.size() > 1)
                {
                    bool all_dff = true;
                    for (auto phash : n.parent_hashes)
                    {
                        Node & p = GNM.at(phash);
                        if (p.last_func != fDFF) 
                        {
                            all_dff = false;
                            break;
                        }
                    }
                    if (all_dff)
                    {
                        return std::make_tuple(false, 0);
                    }
                }
            */
            stack.insert(stack.end(), n.parent_hashes.begin(), n.parent_hashes.end());
        }   

        bool is_valid = true;
        UI support_size = NUM_VARS;

        if (DEBUG) {fmt::print("\t\tAnalyzing vectors\n");}
        for (auto i = 0u; i < NUM_VARS; ++i)
        {
            if (DEBUG) {fmt::print("\t\tPI : {}\t func:{:04x} \t has_dff: {} | has_other: {}| is_reached: {}\n", i, pi_funcs[i], has_dff[i],  has_other[i], is_reached[i]);   }
            if (has_other[i])
            {
                if (DEBUG) {fmt::print("\t\t\t OK PI\n");}
                continue;
            }
            else if ( has_dff[i] ) 
            {
                assert(is_reached[i]); 
                if (DEBUG) {fmt::print("\t\t\t Violating PI\n");}
                return std::make_tuple(false, 0);
            }
            else if ( ~has_dff[i] ) 
            {             
                if (is_reached[i]) 
                {
                    if (DEBUG) {fmt::print("\t\t\t OK PI\n");}
                    continue;
                }
                else
                {
                    if (DEBUG) {fmt::print("\t\t\t Redundant PI\n");}
                    support_size--; 
                    continue;
                }
            }
        }
        return std::make_tuple(true, support_size);
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


std::vector<UI> gen_pi_func(UI nvars)
{
    if (nvars == 1)
    {
        return {1};
    }
    UI power = (1 << (1 << (nvars - 1)));
    UI factor = power + 1;
    std::vector<UI> out = {power - 1};
    for (auto n : gen_pi_func(nvars - 1))
    {
        out.push_back(factor * n);
    }
    return out;
}

#pragma region write_output

void write_csv_gnm(const std::unordered_map<ULL, Node>& gnm, const std::string& filename) {
    // Open output file
    std::ofstream outfile(filename);

    // Write header row to CSV file
    outfile << "Hash,Func,Last Func,Cost,Depth,Xorable,Parent Hashes,Lvl" << std::endl;

    // Write data to CSV file
    for (const auto& [hash, n] : gnm) {
        std::string str = fmt::format("{0},{1},{2},{3},{4},{5:d},{6},{7}", 
                                        hash, n.func, n.last_func, n.cost, n.depth, n.xorable, 
                                        fmt::join(n.parent_hashes, "|"), n.lvl);
        outfile << str << std::endl;
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
        outfile << fmt::format("{0}", hash) << std::endl;
    }
    // Close output file
    outfile.close();
}

std::unordered_map<ULL, Node> read_csv_gnm(const std::string& filename) 
{
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
        node.hash = hash;
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


bool is_good(ULL hash, Node & node, std::unordered_map<ULL, bool> & status, std::unordered_map<ULL, Node>& all_hashes)
{
    bool good = true;
    for (auto & p_hash : node.parent_hashes)
    {
        Node & p_node = all_hashes[p_hash];
        if (status.find(p_hash) == status.end()) // if a parent status is unknown
        {
            status[p_hash] = is_good(p_hash, p_node, status, all_hashes);
        }
        good &= status[p_hash];
    }
    return good;
}

std::unordered_map<ULL, bool> subset_of_pi(std::vector<ULL>& pi, std::unordered_map<ULL, Node>& all_hashes)
{
    std::unordered_map<ULL, bool> status;
    for (ULL hash : pi)
    {
        status[hash] = true;
    }

    for (auto & [hash, node] : all_hashes)
    {
        status[hash] = is_good(hash, node, status, all_hashes);
    }
    return status;
}

#pragma endregion

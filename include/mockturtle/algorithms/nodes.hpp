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
#include <limits>
#include <kitty/kitty.hpp>

bool DEBUG = false;

constexpr uint8_t NUM_VARS = 4u;
typedef kitty::static_truth_table<NUM_VARS> TT;

template <typename number>
inline TT num2tt( number num )
{
    TT tt;
    kitty::create_from_words(tt, &num, &num+1);
    return tt;
}

constexpr uint8_t fDFF   = 0;
constexpr uint8_t fNOT   = 1;
constexpr uint8_t fMERGE = 2;
constexpr uint8_t fOR    = 3;
constexpr uint8_t fAND   = 4;
constexpr uint8_t fXOR   = 5;
constexpr uint8_t fOR3   = 6;
constexpr uint8_t fAND3  = 7;
constexpr uint8_t fMAJ3  = 8;
constexpr uint8_t fCB    = 9;
constexpr uint8_t fSPL   = 10;
constexpr uint8_t fPI    = 11;
constexpr uint8_t fNOFUNC= 99;

constexpr uint8_t kNumThreads = 100;

const std::string GENLIB_PHASE = "UNKNOWN";
constexpr float GENLIB_INPUT_LOAD = 1;
constexpr float GENLIB_MAX_LOAD = 999;
constexpr float GENLIB_RISE_BLOCK_DELAY   = 0.025;
constexpr float GENLIB_RISE_FANOUT_DELAY  = 0.025;
constexpr float GENLIB_FALL_BLOCK_DELAY   = 0.025;
constexpr float GENLIB_FALL_FANOUT_DELAY  = 0.025;


const TT const0 =  num2tt(0);
const TT const1 = ~const0;
const std::vector<TT> PI {  num2tt(0x5555), 
                            num2tt(0x3333), 
                            num2tt(0x0F0F), 
                            num2tt(0x00FF)};

// constexpr bool accel_cost = true;
#define accel_cost false

constexpr std::array<uint32_t,12> COSTS = {7, 9, 8, 8, 8, 7, 11, 11, 11, 8, 7, 0};
std::unordered_map<uint8_t, std::string> F2STR { 
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
std::unordered_map<uint32_t, std::string> PI2LETTER { 
    {0x00FF, "d"},
    {0x0F0F, "c"},
    {0x3333, "b"},
    {0x5555, "a"},
    {0x0000, "0"},
    {0xFFFF, "1"},
    }; 

constexpr uint64_t NUM_TT = (1 << (1 << NUM_VARS));
constexpr uint64_t INF64 = std::numeric_limits<uint64_t>::max();
constexpr uint32_t INF32 = std::numeric_limits<uint32_t>::max();
constexpr uint16_t INF16 = std::numeric_limits<uint16_t>::max();
constexpr uint8_t  INF8  = std::numeric_limits<uint8_t>::max();

// Hash combiner
template <typename T>
static void hash_combine(std::size_t& seed, const T& val) {
    seed ^= std::hash<T>{}(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

inline bool is_pi(TT func)
{
    return kitty::is_const0(func) || kitty::is_const0(~func) || func._bits == 0x5555 || func._bits == 0x3333 || func._bits == 0x0F0F || func._bits == 0x00FF;
}

inline bool is_const(TT func)
{
    return kitty::is_const0(func) || kitty::is_const0(~func);
}


uint64_t calculate_hash(TT func, uint8_t last_func, uint32_t cost, uint8_t depth, bool xorable, std::vector<uint64_t> parent_hashes = {})
{
    std::size_t seed = 0;
    hash_combine(seed, func._bits);
    hash_combine(seed, last_func);
    hash_combine(seed, cost);
    hash_combine(seed, depth);
    hash_combine(seed, xorable);
    for (const auto parent_hash : parent_hashes) {
        hash_combine(seed, parent_hash);
    }
    return seed;
};

class Node 
{
public:
    TT func = num2tt(0);
    uint8_t last_func = fNOFUNC;
    uint32_t cost = INF32;
    uint8_t depth = INF8;
    bool xorable = false;
    std::vector<uint64_t> parent_hashes;
    uint8_t lvl = INF8;
    uint64_t hash = 0;
    std::vector<uint64_t> pi_hashes;

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

    Node(TT _func, uint8_t _last_func, uint32_t _cost, uint8_t _depth, bool _xorable, std::vector<uint64_t> _parent_hashes)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes(_parent_hashes), lvl(depth / 3)
    {
        // Calculate hash based on the hashes of parent_hashes and specified fields
        hash = calculate_hash(func, last_func, cost, depth, xorable, parent_hashes);
        // hash = calculate_hash();
    }

    Node(TT _func, uint8_t _last_func, uint32_t _cost, uint8_t _depth, bool _xorable)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes{}, lvl(depth / 3)
    {
        // Calculate hash based on the hashes of parent_hashes and specified fields
        // It is assumed the node has no parents
        hash = calculate_hash(func, last_func, cost, depth, xorable);
    }

    Node(TT _func, uint8_t _last_func, uint32_t _cost, uint8_t _depth, bool _xorable, std::vector<uint64_t> _parent_hashes, uint64_t _hash)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes(_parent_hashes), lvl(depth / 3), hash(_hash)
    {}

    ~Node() { }

    bool operator==(const Node& other) const {
        return hash == other.hash;
    }

    uint64_t get_hash() 
    {
        return hash;
    }

    std::string to_str() const
    {
        return fmt::format("{0:>20}: 0x{1:04x}=0b{1:016b}{2} | {3} | c{4}, d{5}, l{6} |", hash, func._bits, (xorable?'x':' '), F2STR[last_func], cost, depth, lvl);
    }

    std::string genlib_eqn(std::unordered_map<uint64_t, Node> & nodemap, std::vector<std::string> & symbol, std::vector<uint8_t> & indices) const
    {
        // fmt::print("\t\tAccessing genlib_eqn: {}\n", to_str());
        if (last_func == fPI)
        {
            switch (func._bits)
            {
                case 0x5555: return symbol[3];
                case 0x3333: return symbol[2];
                case 0x0F0F: return symbol[1];
                case 0x00FF: return symbol[0];
                /*
                    case 0x5555: return symbol[indices[3]];
                    case 0x3333: return symbol[indices[2]];
                    case 0x0F0F: return symbol[indices[1]];
                    case 0x00FF: return symbol[indices[0]];
                    case 0x5555: return _symbol[indices[3]];
                    case 0x3333: return _symbol[indices[2]];
                    case 0x0F0F: return _symbol[indices[1]];
                    case 0x00FF: return _symbol[indices[0]];
                    case 0x5555: return "a";
                    case 0x3333: return "b";
                    case 0x0F0F: return "c";
                    case 0x00FF: return "d"; */
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

    std::string to_stack(std::unordered_map<uint64_t, Node> & nodemap, std::vector<std::string> & symbol)
    {
        if (last_func == fPI)
        {
            assert(is_pi(func));
            switch (func._bits)
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

    std::tuple<bool, uint32_t> redundancy_check(std::unordered_map<uint64_t, Node> & GNM)
    {
        std::vector<bool> has_dff      {false, false, false, false};
        std::vector<bool> has_other    {false, false, false, false};
        std::vector<bool> is_reached   {false, false, false, false};
        std::vector<uint64_t> stack = parent_hashes;
        if (DEBUG) {fmt::print("\tAnalyzing stack for node {}\n", to_str());}
        while (!stack.empty())
        {
            uint64_t n_hash = stack.back();
            stack.pop_back();
            Node& n = GNM.at(n_hash); // [n_hash];
            if (DEBUG) {fmt::print("\t\tAnalyzing node: {}\n", n.to_str());}
            auto it = std::find(PI.begin(), PI.end(), n.func);
            if (it != PI.end()) // if the function is a PI
            {
                auto idx = it - PI.begin();
                if (DEBUG) {fmt::print("\t\tFound PI at idx {} for {:04x}\n", idx, n.func._bits);}
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
            stack.insert(stack.end(), n.parent_hashes.begin(), n.parent_hashes.end());
        }   

        bool is_valid = true;
        uint32_t support_size = NUM_VARS;

        if (DEBUG) {fmt::print("\t\tAnalyzing vectors\n");}
        for (auto i = 0u; i < NUM_VARS; ++i)
        {
            if (DEBUG) {fmt::print("\t\tPI : {}\t func:{:04x} \t has_dff: {} | has_other: {}| is_reached: {}\n", i, PI[i]._bits, has_dff[i],  has_other[i], is_reached[i]);   }
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
            return std::hash<uint64_t>{}(node.hash);
        }
    };
}

std::vector<uint64_t> gen_pi(uint32_t nvars)
{
    if (nvars == 1)
    {
        return { 1 };
    }
    uint64_t power = (1 << (1 << (nvars - 1)));
    uint64_t factor = power + 1;
    std::vector<uint64_t> out = {power - 1};
    for (uint64_t n : gen_pi(nvars - 1))
    {
        out.push_back(factor * n);
    }
    return out;
}

#pragma region write_output

void write_csv_gnm(const std::unordered_map<uint64_t, Node>& gnm, const std::string& filename) {
    // Open output file
    std::ofstream outfile(filename);

    // Write header row to CSV file
    outfile << "Hash,Func,Last Func,Cost,Depth,Xorable,Parent Hashes,Lvl" << std::endl;

    // Write data to CSV file
    for (const auto& [hash, n] : gnm) {
        std::string str = fmt::format("{0},{1},{2},{3},{4},{5:d},{6},{7}", 
                                        hash, n.func._bits, n.last_func, n.cost, n.depth, n.xorable, 
                                        fmt::join(n.parent_hashes, "|"), n.lvl);
        outfile << str << std::endl;
    }

    // Close output file
    outfile.close();
}

void write_csv_arr(const std::array<uint64_t, NUM_TT>& arr_hashes, const std::string& filename) {
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

std::unordered_map<uint64_t, Node> read_csv_gnm(const std::string& filename) 
{
    // Open input file
    std::ifstream infile(filename);
    
    std::unordered_map<uint64_t, Node> gnm;

    // Parse CSV file and populate GNM variable
    std::string line;
    std::getline(infile, line);  // skip header row
    while (std::getline(infile, line)) {
        std::stringstream ss;
        ss.str(line);
        std::string field;
        uint64_t hash;
        uint32_t func;
        Node node;
        std::getline(ss, field, ',');
        std::stringstream(field) >> hash;
        node.hash = hash;
        std::getline(ss, field, ',');
        std::stringstream(field) >> func;
        node.func = num2tt(func);
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
        while (std::getline(parent_hashes_ss, field, '|')) 
        {
            uint64_t parent_hash;
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

std::array<uint64_t, NUM_TT> read_csv_arr(const std::string& filename) {
    // Open input file
    std::ifstream infile(filename);
    std::array<uint64_t, NUM_TT> arr_hashes;

    // Parse CSV file and populate array
    std::string line;
    std::getline(infile, line);  // skip header row
    uint32_t index = 0;
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


bool is_good(uint64_t hash, Node & node, std::unordered_map<uint64_t, bool> & status, std::unordered_map<uint64_t, Node>& all_hashes)
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

std::unordered_map<uint64_t, bool> subset_of_pi(std::vector<uint64_t>& pi, std::unordered_map<uint64_t, Node>& all_hashes)
{
    std::unordered_map<uint64_t, bool> status;
    for (uint64_t hash : pi)
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

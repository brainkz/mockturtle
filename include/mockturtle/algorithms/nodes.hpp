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
typedef kitty::static_truth_table<4u> TT4;
typedef kitty::static_truth_table<3u> TT3;
typedef kitty::static_truth_table<2u> TT2;
typedef kitty::static_truth_table<1u> TT1;

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


// constexpr bool accel_cost = true;
#define accel_cost false

#define VECTOR_CONTAINS(vec, value) (std::find(vec.begin(), vec.end(), value) != vec.end())


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

std::array<TT4, NUM_VARS> PI_WORDS = {{0x5555, 0x3333, 0x0F0F, 0x00FF}};

std::unordered_map<uint16_t, std::string> PI2LETTER { 
    {0x00FF, "d"},
    {0x0F0F, "c"},
    {0x3333, "b"},
    {0x5555, "a"},
    {0x0000, "0"},
    {0xFFFF, "1"},
    }; 

std::unordered_map<uint16_t, uint16_t> dummy_map {
    {0x00FF, 0x00FF},
    {0x0F0F, 0x0F0F},
    {0x3333, 0x3333},
    {0x5555, 0x5555}
    };

std::unordered_map<uint16_t, std::string> FN2LETTER { 
    {0x00FF, "d"},
    {0x0F0F, "c"},
    {0x3333, "b"},
    {0x5555, "a"},
    {0x0F, "c"},
    {0x33, "b"},
    {0x55, "a"},
    {0x3, "b"},
    {0x5, "a"},
    {0x1, "a"},
    {0x0000, "CONST0"},
    {0xFFFF, "CONST1"},
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

// TODO: find P representative and minimize functions in support
// TODO: permute the delays accordingly

// TODO: in next function, test whether one delay pattern dominates the other. 
// TODO: Perhaps, a delay object will help


#if false
    // TODO : Reduce support 
    // kitty::static_truth_table<4u> 



    std::vector<uint8_t> min_base_delays( UI& func , std::array<US,NUM_VARS> delays)
    {
        std::vector<uint8_t> support;
        kitty::static_truth_table<NUM_VARS> tt;
        tt._bits = func;

        auto k = 0u;
        for ( auto i = 0u; i < tt.num_vars(); ++i )
        {
            if ( !has_var( tt, i ) )
            {
            continue;
            }
            if ( k < i )
            {
            swap_inplace( tt, k, i );
            }
            support.push_back( i );
            ++k;
        }

        return support;
    }
#endif 

template <typename TT>
bool _tt_gt (const std::pair<TT, uint8_t>& a, const std::pair<TT, uint8_t>& b) { return ~(a.first < b.first); }

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
        return fmt::format("{:016b}{} | {} | c{}, d{}, l{} |", func, (xorable?'x':' '), F2STR[(US)last_func], cost, depth, lvl);
        // return fmt::format("{:04x}{} | {} | {}, {} | ", func, (xorable?'x':' '), F2STR[(US)last_func], cost, depth);
        // return "Func: " + std::to_string(func) + "|Last: " + std::to_string(last_func) + "|Cost: " + std::to_string(cost) + "|Depth: " + std::to_string(depth) + "|X: " + std::to_string(xorable);
    }

    std::string genlib_eqn(std::unordered_map<ULL, Node> & nodemap, std::vector<UI> & pis) const
    {
        // if (DEBUG) {fmt::print("\t\tAccessing genlib_eqn: {}\n", to_str());}
        if (last_func == fPI)
        {
            if (std::find(pis.begin(), pis.end(), func) == pis.end())
            {
                pis.push_back(func);
            }
            return PI2LETTER[func];
        }
        else if (last_func == fDFF)
        {
            assert(parent_hashes.size() == 1);
            return fmt::format("{}", nodemap[parent_hashes.back()].genlib_eqn(nodemap, pis));
        }
        else if (last_func == fNOT)
        {
            assert(parent_hashes.size() == 1);
            return fmt::format("!{}", nodemap[parent_hashes.back()].genlib_eqn(nodemap, pis));
        }
        else if (last_func == fCB)
        {
            assert(parent_hashes.size() == 2);
            return fmt::format("({0}+{1})", nodemap[parent_hashes.front()].genlib_eqn(nodemap, pis), nodemap[parent_hashes.back()].genlib_eqn(nodemap, pis));
        }
        else if (last_func == fOR || last_func == fMERGE)
        {
            assert(parent_hashes.size() == 2);
            return fmt::format("({0}|{1})", nodemap[parent_hashes.front()].genlib_eqn(nodemap, pis), nodemap[parent_hashes.back()].genlib_eqn(nodemap, pis));
        }
        else if (last_func == fAND)
        {
            assert(parent_hashes.size() == 2);
            return fmt::format("({0}&{1})", nodemap[parent_hashes.front()].genlib_eqn(nodemap, pis), nodemap[parent_hashes.back()].genlib_eqn(nodemap, pis));
        }
        else if (last_func == fXOR)
        {
            assert(parent_hashes.size() == 2);
            // return fmt::format("(!({0})*({1})+({0})*!({1}))", nodemap[parent_hashes.front()].genlib_eqn(nodemap, pis), nodemap[parent_hashes.back()].genlib_eqn(nodemap, pis));
            return fmt::format("({0}^{1})", nodemap[parent_hashes.front()].genlib_eqn(nodemap, pis), nodemap[parent_hashes.back()].genlib_eqn(nodemap, pis));
        }
        else
        {
            if (DEBUG) {fmt::print("Unsupported function {}", to_str());}
            return "";
        }
    }

    std::string to_genlib(std::unordered_map<ULL, Node> & nodemap, const std::vector<UI> & levels, const std::vector<UI> PI_funcs) const
    {
        std::vector<UI> pis;
        std::string str = fmt::format("GATE 0x{:04x}_{} {} O={};\n", func, fmt::join(levels, ""), cost,  genlib_eqn(nodemap, pis));
        for (auto & pi : pis)
        {
            UL idx = std::find(PI_funcs.begin(), PI_funcs.end(), pi) - PI_funcs.begin();
            auto true_lvl = (depth + 1) / 3;
            std::string line = fmt::format("\tPIN {} {} {} {} {:d} {:0.3f} {:d} {:0.3f}\n", 
            PI2LETTER[pi], GENLIB_PHASE, GENLIB_INPUT_LOAD, GENLIB_MAX_LOAD, 
            true_lvl - levels[idx], GENLIB_RISE_FANOUT_DELAY, true_lvl - levels[idx], GENLIB_FALL_FANOUT_DELAY);
            str.append(line);
            // if (DEBUG) {fmt::print(line);}
        }
        // if (DEBUG) {fmt::print(str);}
        return str;
    }

    // std::string to_genlib(std::unordered_map<ULL, Node> & nodemap, const std::vector<UI> & levels, std::unordered_map<UI, std::string> pi2symbol) const
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

    std::string to_stack(std::unordered_map<ULL, Node> & nodemap, const std::unordered_map<uint16_t, uint16_t> & pi_map = dummy_map) const
    {
        if (last_func == fPI)
        {
            // auto fn = pi_map.at(func);
            // return FN2LETTER[fn];
            return PI2LETTER[func];
        }
        else if (last_func == fDFF)
        {
            assert(parent_hashes.size() == 1);
            return fmt::format("DFF({})", nodemap[parent_hashes.back()].to_stack(nodemap, pi_map));
        }
        else if (last_func == fNOT)
        {
            assert(parent_hashes.size() == 1);
            return fmt::format("NOT({})", nodemap[parent_hashes.back()].to_stack(nodemap, pi_map));
        }
        else if (last_func == fCB || last_func == fMERGE)
        {
            assert(parent_hashes.size() == 2);
            return fmt::format("CB({0}, {1})", nodemap[parent_hashes.front()].to_stack(nodemap, pi_map), nodemap[parent_hashes.back()].to_stack(nodemap, pi_map));
        }
        else if ( last_func == fOR )
        {
            assert(parent_hashes.size() == 2);
            return fmt::format("OR({0}, {1})", nodemap[parent_hashes.front()].to_stack(nodemap, pi_map), nodemap[parent_hashes.back()].to_stack(nodemap, pi_map));
        }
        else if (last_func == fAND)
        {
            assert(parent_hashes.size() == 2);
            return fmt::format("AND({0}, {1})", nodemap[parent_hashes.front()].to_stack(nodemap, pi_map), nodemap[parent_hashes.back()].to_stack(nodemap, pi_map));
        }
        else if (last_func == fXOR)
        {
            assert(parent_hashes.size() == 2);
            return fmt::format("XOR({0}, {1})", nodemap[parent_hashes.front()].to_stack(nodemap, pi_map), nodemap[parent_hashes.back()].to_stack(nodemap, pi_map));
        }
        else
        {
            if (DEBUG) {fmt::print("Unsupported function {}", to_str());}
            return "";
        }
    }

    auto process_nodes(std::unordered_map<ULL, Node> & GNM)
    {
        std::vector<std::pair<ULL,US>> stack { std::make_pair(hash, 0) };
        std::vector<ULL> seen;
        std::array<US,NUM_VARS> delays {{ 0xFF, 0xFF, 0xFF, 0xFF }};
        bool valid = true;
        while (!stack.empty())
        {
            auto [nhash, nlvl] = stack.back();
            stack.pop_back();
            
            if (VECTOR_CONTAINS(seen, nhash))
            {
                continue;
            }
            else
            {
                seen.push_back( nhash );
            }

            Node& n = GNM[nhash];
            if (n.last_func == fNOFUNC) 
            {
                valid = false;
                continue;
            }
            else if (n.last_func == fPI)
            {
                switch (n.func)
                {
                    case 0x5555: delays[0] = nlvl; break;
                    case 0x3333: delays[1] = nlvl; break;
                    case 0x0F0F: delays[2] = nlvl; break;
                    case 0x00FF: delays[3] = nlvl; break;
                }
                continue;
            }

            US new_lvl = nlvl + ((n.last_func == fDFF) | (n.last_func == fNOT) | (n.last_func == fXOR));

            for (ULL phash : n.parent_hashes)
            {
                stack.push_back( std::make_pair(phash, new_lvl) );
            }
        }

    }
    /*
        std::tuple<bool, UI> process_nodes(std::unordered_map<ULL, Node> & GNM)
        {
            std::array<UI,  4> pi_funcs     {{0x5555, 0x3333, 0x0F0F, 0x00FF}};
            std::array<bool,4> has_dff      {{false, false, false, false}};
            std::array<bool,4> has_other    {{false, false, false, false}};
            std::array<bool,4> is_reached   {{false, false, false, false}};
            
            std::vector<ULL> stack = { hash };
            if (DEBUG) {fmt::print("\tAnalyzing stack for node {}\n", to_str());}
            while (!stack.empty())
            {
                ULL n_hash = stack.back();
                stack.pop_back();
                Node& n = GNM[n_hash];
                if (n.last_func == fNOFUNC) continue;
                if (DEBUG) {fmt::print("\t\tAnalyzing node: {}\n", n.to_str());}
                auto it = std::find(pi_funcs.begin(), pi_funcs.end(), n.func);
                if (it != pi_funcs.end()) // if the function is a PI
                {
                    auto idx = it - pi_funcs.begin();
                    if (DEBUG) {fmt::print("\t\tFound PI at idx {} for {:04x}\n", idx, n.func);}
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
    */
    std::tuple<bool, UI> redundancy_check(std::unordered_map<ULL, Node> & GNM)
    {
        // std::vector<UI>   pi_funcs     {0x00FF, 0x0F0F, 0x3333, 0x5555};
        std::array<UI,  4> pi_funcs     {0x5555, 0x3333, 0x0F0F, 0x00FF};
        std::array<bool,4> has_dff      {false, false, false, false};
        std::array<bool,4> has_other    {false, false, false, false};
        std::array<bool,4> is_reached   {false, false, false, false};
        
        std::vector<ULL> stack = parent_hashes;
        if (DEBUG) {fmt::print("\tAnalyzing stack for node {}\n", to_str());}
        while (!stack.empty())
        {
            ULL n_hash = stack.back();
            stack.pop_back();
            Node& n = GNM[n_hash];
            if (n.last_func == fNOFUNC) continue;
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

// 1. Reduce support if redundant variables are found. 
// 2. Determine delays.
// std::pair<kitty::static_truth_table<3u>,std::vector<std::pair<kitty::static_truth_table<3u>, uint8_t>>> reduce_tt3( TT4 & tt, TT4 & p_mask, std::vector<std::pair<TT4, uint8_t>> pi_delays )
template <unsigned N>
std::pair<kitty::static_truth_table<N>, std::vector<std::pair<kitty::static_truth_table<N>, uint8_t>>> reduce_tt( TT4 &tt, TT4 &p_mask, const std::vector<std::pair<TT4, uint8_t>> pi_delays) 
{
    kitty::static_truth_table<N> tt_reduced;
    tt_reduced._bits = 0;
    std::vector<std::pair<kitty::static_truth_table<N>, uint8_t>> new_pi_delays;
    for (auto [pi, delay] : pi_delays)
    {
        kitty::static_truth_table<N> zero_tt;
        zero_tt._bits = 0;
        new_pi_delays.push_back(std::make_pair(zero_tt, delay));
    }

    for (uint8_t i = 0u, k = 0u; i < (1 << NUM_VARS); ++i )
    {
        if ( kitty::get_bit(p_mask, i) )
        {
            if (kitty::get_bit(tt, i))
            {
                kitty::set_bit( tt_reduced, k );
            }
            for (auto j = 0u; j < new_pi_delays.size(); ++j)
            {
                if (kitty::get_bit(pi_delays[j].first, i))
                {
                    kitty::set_bit( new_pi_delays[j].first, k );
                }
            }
            k++;
        }
    }
    return std::make_pair(tt_reduced, new_pi_delays);
}

template <typename TT>
std::tuple<TT, std::vector<std::pair<TT, uint8_t>>> canonize_p( const TT& tt, std::vector<std::pair<TT, uint8_t>> pi_delays )
{
    const auto num_vars = tt.num_vars();

    /* Special case for n = 0 (const 0/1) or n = 1 (PI) */
    if ( num_vars == 0 || num_vars == 1 )
    {
        return std::make_tuple( tt, pi_delays );
    }

    assert( num_vars >= 2 && num_vars <= 7 );

    auto t1 = tt;
    auto tmin = t1;
    std::vector<std::pair<TT, uint8_t>> new_pi_delays = pi_delays;

    const auto& swaps = kitty::detail::swaps[num_vars - 2u];

    int best_swap = -1;
    for ( std::size_t i = 0; i < swaps.size(); ++i )
    {
        const auto pos = swaps[i];
        kitty::swap_adjacent_inplace( t1, pos );
        for (auto j = 0u; j < pi_delays.size(); ++j)
        {
            kitty::swap_adjacent_inplace( pi_delays[j].first, pos );
        }

        if ( t1 < tmin )
        {
            best_swap = static_cast<int>( i );
            tmin = t1;
            new_pi_delays = pi_delays;
        }
    }
    return std::make_tuple( tmin, new_pi_delays );
}

/* Returns indices of support, and p_mask 
    // std::tuple<std::vector<uint8_t>, std::vector<TT4>, TT4, uint8_t> reduce_support(uint64_t hash, std::unordered_map<ULL, Node> & hashmap) */
std::tuple<std::vector<uint8_t>, TT4> get_support(uint64_t hash, std::unordered_map<ULL, Node> & hashmap)
{
    Node & n = hashmap[hash];
    TT4 tt;
    tt._bits = n.func;
    std::vector<uint8_t> support_idx;
    // std::vector<TT4> redundant;
    TT4 p_mask; // mask describing which variables to take
    p_mask._bits = 0xFFFF;
    // uint8_t p_mask_weight = (1 << NUM_VARS);
    for (uint8_t i = 0u; i < NUM_VARS; ++i)
    {
        if (kitty::has_var(tt, i))
        {
            support_idx.push_back(i);
        }
        else
        {
            // redundant.push_back(PI_WORDS[i]);
            p_mask &= PI_WORDS[i];
            // p_mask_weight >>= 1;
        }
    }
    return std::make_tuple( support_idx, p_mask ) ;
    // return std::make_tuple( support_idx, redundant, p_mask, p_mask_weight ) ;
}


/* Returns the functions of PIs and their correponding delays*/
std::pair<std::vector<std::pair<TT4, uint8_t>>, bool> get_delays(uint64_t hash, std::unordered_map<ULL, Node> & hashmap)
{
    std::vector<std::pair<uint64_t, uint8_t>> stack { std::make_pair(hash , 0) };
    std::vector<uint64_t> seen;
    std::vector<std::pair<TT4, uint8_t>> pi_delays;
    while (!stack.empty())
    {
        auto [h, delay] = stack.back();
        Node & n = hashmap[h];
        // fmt::print("\tStack: analyzing {}: {}\n", h, n.to_str());
        stack.pop_back();
        if (std::find( seen.begin(), seen.end(), h) != seen.end() )
        {
            continue;
        }
        seen.push_back(h);

        // Node & n = hashmap[h];
        if (n.last_func == fNOFUNC)
        {
            return std::make_pair( pi_delays , false );
        }
        else if (n.last_func == fPI)
        {
            TT4 tt;
            tt._bits = n.func;
            pi_delays.push_back( std::make_pair(tt, delay) );
            continue;
        }

        uint8_t new_delay = delay + (n.last_func == fDFF || n.last_func == fNOT || n.last_func == fXOR);
        for (uint64_t phash : n.parent_hashes)        
        {
            stack.push_back( std::make_pair( phash, new_delay ) );
        }
    }
    std::sort(pi_delays.begin(), pi_delays.end(), [](const std::pair<TT4, uint8_t>& a, const std::pair<TT4, uint8_t>& b) { return ~(a.first < b.first); });
    return std::make_pair( pi_delays , true );
}

// template <unsigned N>
// std::string to_genlib(std::unordered_map<ULL, Node> & nodemap, const std::vector<std::pair<TT4, uint8_t>> pi_delays, const std::vector<std::pair<kitty::static_truth_table<N>, uint8_t>> new_pi_delays) const
// {
//     std::vector<UI> pis;
//     std::string str = fmt::format("GATE 0x{:04x}_{} {} O={};\n\t#{}\n", func, fmt::join(levels, ""), cost,  genlib_eqn(nodemap, pis), to_stack(nodemap));

//     for (auto & pi : pis)
//     {
//         UL idx = std::find(PI_funcs.begin(), PI_funcs.end(), pi) - PI_funcs.begin();
//         std::string line = fmt::format("\tPIN {} {} {} {} {:d} {:0.3f} {:d} {:0.3f}\n", 
//         PI2LETTER[pi], GENLIB_PHASE, GENLIB_INPUT_LOAD, GENLIB_MAX_LOAD, 
//         true_lvl - levels[idx], GENLIB_RISE_FANOUT_DELAY, true_lvl - levels[idx], GENLIB_FALL_FANOUT_DELAY);
//         str.append(line);
//         // if (DEBUG) {fmt::print(line);}
//     }
//     // if (DEBUG) {fmt::print(str);}
//     return str;
// }

std::tuple<bool, uint8_t, uint16_t, std::vector<uint8_t>, std::unordered_map<uint16_t, uint16_t>> process_node(const uint64_t hash, std::unordered_map<ULL, Node> & hashmap) // , const std::vector<UI> & levels
{
    // std::vector<std::pair<TT4, uint8_t>> 
    auto [pi_delays, status] = get_delays(hash, hashmap);
    if (!status)
    {
        return std::make_tuple( false, 0, 0, std::vector<uint8_t>{}, std::unordered_map<uint16_t, uint16_t>{} );
    }
    
    auto [support_idx, p_mask] = get_support(hash, hashmap);
    Node & n = hashmap[hash];

    TT4 tt;
    tt._bits = n.func;

    std::unordered_map<uint16_t, uint16_t> pi_map;

    if (support_idx.size() == 4)
    {
        auto [tmin, new_pi_delays] = canonize_p( tt, pi_delays );
        // auto [tmin, _0, perm] = kitty::exact_p_canonization( tt );
        for (auto i = 0u; i < 4; ++i)
        {
            pi_map[ pi_delays[i].first._bits ] = new_pi_delays[i].first._bits;
        }

        auto sorted_delays = new_pi_delays; // delays sorted in descending order
        std::sort(sorted_delays.begin(), sorted_delays.end(), [](const std::pair<TT4, uint8_t>& a, const std::pair<TT4, uint8_t>& b) { return ~(a.first < b.first); });

        std::vector<uint8_t> delay_pattern;
        for (auto [pi, delay] : sorted_delays )
        {
            delay_pattern.push_back(delay);
        }
        return std::make_tuple( true, support_idx.size(), tmin._bits, delay_pattern, pi_map );

        #if false
            fmt::print("{}\n", n.to_stack(hashmap));
            fmt::print("Original TT: {0:04x}={0:016b}\n", n.func);
            fmt::print("Original PI: {:04x}, {:04x}, {:04x}, {:04x}\n", pi_delays[0].first._bits, pi_delays[1].first._bits, pi_delays[2].first._bits, pi_delays[3].first._bits);
            fmt::print("Original Delays: {}, {}, {}, {}\n", pi_delays[0].second, pi_delays[1].second, pi_delays[2].second, pi_delays[3].second);
            fmt::print("Canonized TT: {}={}\n", kitty::to_hex(tmin), kitty::to_binary(tmin));
            fmt::print("Canonized PI: {:04x}, {:04x}, {:04x}, {:04x}\n\n", new_pi_delays[0].first._bits, new_pi_delays[1].first._bits, new_pi_delays[2].first._bits, new_pi_delays[3].first._bits);
            fmt::print("Canonized Delays: {}, {}, {}, {}\n\n", new_pi_delays[0].second, new_pi_delays[1].second, new_pi_delays[2].second, new_pi_delays[3].second);
        #endif
    }
    if (support_idx.size() == 3)
    {
        auto [tt_reduced, reduced_pi_delays] = reduce_tt<3>( tt, p_mask, pi_delays );
        auto [tmin, new_pi_delays] = canonize_p( tt_reduced, reduced_pi_delays );
        for (auto i = 0u; i < 3; ++i)
        {
            pi_map[ pi_delays[i].first._bits ] = new_pi_delays[i].first._bits;
        }

        auto sorted_delays = new_pi_delays; // delays sorted in descending order
        std::sort(sorted_delays.begin(), sorted_delays.end(), [](const std::pair<TT3, uint8_t>& a, const std::pair<TT3, uint8_t>& b) { return ~(a.first < b.first); });

        std::vector<uint8_t> delay_pattern;
        for (auto [pi, delay] : sorted_delays )
        {
            delay_pattern.push_back(delay);
        }
        return std::make_tuple( true, support_idx.size(), tmin._bits, delay_pattern, pi_map );

        // std::tuple<TT, int, std::vector<std::pair<TT4, uint8_t>>>
        // auto [tmin, new_pi_delays] = canonize_p( tt, pi_delays );
        #if false
            fmt::print("{}\n", n.to_stack(hashmap));
            fmt::print("Original TT: {0:04x}={0:016b}\n", n.func);
            fmt::print("Original PI: {:04x}, {:04x}, {:04x}\n", pi_delays[0].first._bits, pi_delays[1].first._bits, pi_delays[2].first._bits);
            fmt::print("Original Delays: {}, {}, {}\n", pi_delays[0].second, pi_delays[1].second, pi_delays[2].second);
            
            fmt::print("Reduced TT: {}={}\n", kitty::to_hex(tt_reduced), kitty::to_binary(tt_reduced));
            fmt::print("Reduced PI: {:02x}, {:02x}, {:02x}\n", reduced_pi_delays[0].first._bits, reduced_pi_delays[1].first._bits, reduced_pi_delays[2].first._bits);
            fmt::print("Reduced Delays: {}, {}, {}\n\n", reduced_pi_delays[0].second, reduced_pi_delays[1].second, reduced_pi_delays[2].second);
            
            fmt::print("Canonized TT: {}={}\n", kitty::to_hex(tmin), kitty::to_binary(tmin));
            fmt::print("Canonized PI: {:02x}, {:02x}, {:02x}\n", new_pi_delays[0].first._bits, new_pi_delays[1].first._bits, new_pi_delays[2].first._bits);
            fmt::print("Canonized Delays: {}, {}, {}\n\n", new_pi_delays[0].second, new_pi_delays[1].second, new_pi_delays[2].second);
        #endif
    }
    if (support_idx.size() == 2)
    {
        auto [tt_reduced, reduced_pi_delays] = reduce_tt<2>( tt, p_mask, pi_delays );
        auto [tmin, new_pi_delays] = canonize_p( tt_reduced, reduced_pi_delays );
        for (auto i = 0u; i < 2; ++i)
        {
            pi_map[ pi_delays[i].first._bits ] = new_pi_delays[i].first._bits;
        }
        
        auto sorted_delays = new_pi_delays; // delays sorted in descending order
        std::sort(sorted_delays.begin(), sorted_delays.end(), [](const std::pair<TT2, uint8_t>& a, const std::pair<TT2, uint8_t>& b) { return ~(a.first < b.first); });
        std::vector<uint8_t> delay_pattern;
        for (auto [pi, delay] : sorted_delays )
        {
            delay_pattern.push_back(delay);
        }
        return std::make_tuple( true, support_idx.size(), tmin._bits, delay_pattern, pi_map );

        // std::tuple<TT, int, std::vector<std::pair<TT4, uint8_t>>>
        // auto [tmin, new_pi_delays] = canonize_p( tt, pi_delays );
        #if false
            fmt::print("{}\n", n.to_stack(hashmap));
            fmt::print("Original TT: {0:04x}={0:016b}\n", n.func);
            fmt::print("Original PI: {:04x}, {:04x}\n", pi_delays[0].first._bits, pi_delays[1].first._bits);
            fmt::print("Original Delays: {}, {}\n", pi_delays[0].second, pi_delays[1].second);
            
            fmt::print("Reduced TT: {}={}\n", kitty::to_hex(tt_reduced), kitty::to_binary(tt_reduced));
            fmt::print("Reduced PI: {:01x}, {:01x}\n", reduced_pi_delays[0].first._bits, reduced_pi_delays[1].first._bits);
            fmt::print("Reduced Delays: {}, {}\n\n", reduced_pi_delays[0].second, reduced_pi_delays[1].second);
            
            fmt::print("Canonized TT: {}={}\n", kitty::to_hex(tmin), kitty::to_binary(tmin));
            fmt::print("Canonized PI: {:01x}, {:01x}\n", new_pi_delays[0].first._bits, new_pi_delays[1].first._bits);
            fmt::print("Canonized Delays: {}, {}\n\n", new_pi_delays[0].second, new_pi_delays[1].second);
        #endif
    }
    if (support_idx.size() == 1)
    {
        auto [tmin, new_pi_delays] = reduce_tt<1>( tt, p_mask, pi_delays );
        pi_map[ pi_delays[0].first._bits ] = new_pi_delays[0].first._bits;

        auto sorted_delays = new_pi_delays; // cost_hash sorted in descending order
        std::sort(sorted_delays.begin(), sorted_delays.end(), [](const std::pair<TT1, uint8_t>& a, const std::pair<TT1, uint8_t>& b) { return ~(a.first < b.first); });
        std::vector<uint8_t> delay_pattern;
        for (auto [pi, delay] : sorted_delays )
        {
            delay_pattern.push_back(delay);
        }
        return std::make_tuple( true, support_idx.size(), tmin._bits, delay_pattern, pi_map );

        // std::tuple<TT, int, std::vector<std::pair<TT4, uint8_t>>>
        // auto [tmin, new_pi_delays] = canonize_p( tt, pi_delays );
        #if false
            fmt::print("LEN_PI_DELAYS: {}\n", pi_delays.size());
            fmt::print("{}\n", n.to_stack(hashmap));
            fmt::print("Original TT: {0:04x}={0:016b}\n", n.func);
            fmt::print("Original PI: {:016b}\n", pi_delays[0].first._bits);
            fmt::print("Original Delays: {}\n", pi_delays[0].second);
            
            fmt::print("Canonized TT: {}={}\n", kitty::to_hex(tmin), kitty::to_binary(tmin));
            fmt::print("Canonized PI: {:02b}\n", new_pi_delays[0].first._bits);
            fmt::print("Canonized Delays: {}\n\n", new_pi_delays[0].second);
        #endif
    }
    return std::make_tuple( false, 0, 0, std::vector<uint8_t>{}, std::unordered_map<uint16_t, uint16_t>{} );
}

    // if (support_idx.size() == 4u)
    // {
    //     auto[tt_p, best_swap] = canonize_p(tt);
    // }
    // else if (support_idx.size() == 3u)
    // {
    //     uint16_t reduced_tt = reduce_tt ( tt, p_mask, hamming_weight );

    // }

    // std::vector<uint64_t> stack { hash };



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

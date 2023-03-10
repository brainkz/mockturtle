#include <iostream>
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
#include <algorithm>
#include <execution>

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

const UI kNumThreads = 10;


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
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes(_parent_hashes), lvl((depth + 1) / 3)
    {
        // Calculate hash based on the hashes of parent_hashes and specified fields
        hash = calculate_hash(func, last_func, cost, depth, xorable, parent_hashes);
        // hash = calculate_hash();
    }

    Node(TT _func, US _last_func, UI _cost, UI _depth, bool _xorable)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes{}, lvl((depth + 1) / 3)
    {
        // Calculate hash based on the hashes of parent_hashes and specified fields
        // It is assumed the node has no parents
        hash = calculate_hash(func, last_func, cost, depth, xorable);
    }

    Node(TT _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes, ULL _hash)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes(_parent_hashes), lvl((depth + 1) / 3), hash(_hash)
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
        return fmt::format("{:016b}{} | {} | {}, {} | ", func, (xorable?'x':' '), F2STR[(US)last_func], cost, depth);
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

std::unordered_map<ULL, Node> GNM;
std::array<Node, NUM_TT> GEA;
std::array<Node, NUM_TT> GEX;

#if accel_cost 
    std::tuple<std::unordered_map<ULL, UI>, std::unordered_set<ULL>, UI> pred_count(ULL hash)
    {
        std::stack<ULL> stack;
        stack.push(hash);

        std::unordered_map<ULL, UI> ct_spl;
        std::unordered_set<ULL> non_splittable_nodes;
        UI gate_cost = 0;

        // fmt::print("\t\t\tInitial cost {}\n", gate_cost);
        while (!stack.empty())
        {
            ULL n_hash = stack.top();
            stack.pop();
            Node& n = GNM[n_hash];
            // fmt::print("\t\t\tProcessing node {}\n", n.to_str());
            ct_spl[n_hash]++;
            if (ct_spl[n_hash] == 1)
            {
                // fmt::print("\t\t\tFirst time, adding the cost of {} ({}). Total cost is {}\n", F2STR[n.last_func], COSTS[n.last_func], gate_cost);
                gate_cost += COSTS[n.last_func];
                for (const auto& p_hash : n.parent_hashes)
                {
                    Node& p = GNM[p_hash];
                    stack.push(p_hash);
                    // fmt::print("\t\t\tAdding parent node {} to stack\n", p.to_str());
                }
            }
            else
            {
                gate_cost += COSTS[fSPL];
            }
            if (n.last_func == fAND || n.last_func == fOR)
            {
                for (const auto& p_hash : n.parent_hashes)
                {
                    non_splittable_nodes.emplace(p_hash);
                }
            }
        }

        return std::make_tuple(ct_spl, non_splittable_nodes, gate_cost);
    }

    std::tuple<std::unordered_map<ULL, UI>, std::unordered_set<ULL>, UI> combine(ULL h1, ULL h2)
    {
        Node& n1 = GNM[h1];
        Node& n2 = GNM[h2];

        std::stack<ULL> stack;
        stack.push(h2);

        std::unordered_map<ULL, UI> ct_spl;
        std::unordered_set<ULL> non_splittable_nodes;
        UI gate_cost = 0;

        // fmt::print("\t\t\tInitial cost {}\n", gate_cost);
        while (!stack.empty())
        {
            ULL n_hash = stack.top();
            stack.pop();
            Node& n = GNM[n_hash];
            // fmt::print("\t\t\tProcessing node {}\n", n.to_str());
            ct_spl[n_hash]++;
            if (ct_spl[n_hash] == 1)
            {
                // fmt::print("\t\t\tFirst time, adding the cost of {} ({}). Total cost is {}\n", F2STR[n.last_func], COSTS[n.last_func], gate_cost);
                gate_cost += COSTS[n.last_func];
                for (const auto& p_hash : n.parent_hashes)
                {
                    Node& p = GNM[p_hash];
                    stack.push(p_hash);
                    // fmt::print("\t\t\tAdding parent node {} to stack\n", p.to_str());
                }
            }
            else
            {
                gate_cost += COSTS[fSPL];
            }
            if (n.last_func == fAND || n.last_func == fOR)
            {
                for (const auto& p_hash : n.parent_hashes)
                {
                    non_splittable_nodes.emplace(p_hash);
                }
            }
        }

        return std::make_tuple(ct_spl, non_splittable_nodes, gate_cost);
    }



#endif

std::tuple<ULL, bool> force_create_node(TT _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes, ULL _hash)
{
    auto it = GNM.find(_hash);
    bool node_is_new = (it == GNM.end());
    if (node_is_new) //the node is new
    {
        GNM[_hash] = Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
        std::make_tuple(_hash, true);
    }
    return std::make_tuple(_hash, false); 
}

std::tuple<ULL, bool> force_create_node(TT _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes = {})
{
    ULL _hash = calculate_hash(_func, _last_func, _cost, _depth, _xorable, _parent_hashes);
    return force_create_node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
}

std::tuple<ULL, bool> create_node(TT _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes, ULL _hash, std::unordered_map<ULL, Node>& hash_map = GNM)
{
    auto it = hash_map.find(_hash);
    bool node_is_new = (it == hash_map.end());
    if (node_is_new) //the node is new
    {
        fmt::print("\t\t\t\tChecking for creation {}\n", Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash).to_str());
        UI lvl = (_depth + 1) / 3;
        Node& best_any = GEA[_func];
        fmt::print("\t\t\t\tConsidering {}: {}, {}\n", _func,          lvl,         _cost);
        fmt::print("\t\t\t\tBest A      {}: {}, {}\n", _func, best_any.lvl, best_any.cost);

        if (_xorable)
        {
            fmt::print("\t\t\t\tFunction is xorable\n");
            Node& best_xor = GEX[_func];
            if (std::tie(lvl, _cost ) <= std::tie(best_any.lvl, best_any.cost))
            {
                hash_map[_hash] = Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
                GEA[_func] = hash_map[_hash];
                GEX[_func] = hash_map[_hash];
                fmt::print("\t\t\t\tNode is better than any existing function\n");
                return std::make_tuple(_hash, true);
            }
            else if (std::tie(lvl, _cost ) <= std::tie(best_xor.lvl, best_xor.cost))
            {
                hash_map[_hash] = Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
                GEX[_func] = hash_map[_hash];
                fmt::print("\t\t\t\tNode is better than any existing xorable function\n");
                return std::make_tuple(_hash, true);
            }
        }
        else 
        {
            fmt::print("\t\t\t\tFunction is NOT xorable\n");
            if (std::tie(lvl, _cost ) <= std::tie(best_any.lvl, best_any.cost))
            {
                hash_map[_hash] = Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
                GEA[_func] = hash_map[_hash];
                fmt::print("\t\t\t\tNode is better than any existing function\n");
                return std::make_tuple(_hash, true);
            }
        }
    }
    else{
        fmt::print("\t\t\t\tNode already exists\n");
    }
    fmt::print("\t\t\t\tNode is not created\n");
    return std::make_tuple(_hash,  false); 
    
}

std::tuple<ULL, bool> create_node_no_upd(TT _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes, ULL _hash, std::unordered_map<ULL, Node>& hash_map = GNM)
{
    auto it = hash_map.find(_hash);
    bool node_is_new = (it == hash_map.end());
    if (node_is_new) //the node is new
    {
        fmt::print("\t\t\t\tChecking for creation {}\n", Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash).to_str());
        UI lvl = (_depth + 1) / 3;
        Node& best_any = GEA[_func];
        fmt::print("\t\t\t\tConsidering {}: {}, {}\n", _func,          lvl,         _cost);
        fmt::print("\t\t\t\tBest A      {}: {}, {}\n", _func, best_any.lvl, best_any.cost);

        if (_xorable)
        {
            fmt::print("\t\t\t\tFunction is xorable\n");
            Node& best_xor = GEX[_func];
            if (std::tie(lvl, _cost ) <= std::tie(best_any.lvl, best_any.cost))
            {
                hash_map[_hash] = Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
                fmt::print("\t\t\t\tNode is better than any existing function\n");
                return std::make_tuple(_hash, true);
            }
            else if (std::tie(lvl, _cost ) <= std::tie(best_xor.lvl, best_xor.cost))
            {
                hash_map[_hash] = Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
                fmt::print("\t\t\t\tNode is better than any existing xorable function\n");
                return std::make_tuple(_hash, true);
            }
        }
        else 
        {
            fmt::print("\t\t\t\tFunction is NOT xorable\n");
            if (std::tie(lvl, _cost ) <= std::tie(best_any.lvl, best_any.cost))
            {
                hash_map[_hash] = Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
                fmt::print("\t\t\t\tNode is better than any existing function\n");
                return std::make_tuple(_hash, true);
            }
        }
    }
    else{
        fmt::print("\t\t\t\tNode already exists\n");
    }
    fmt::print("\t\t\t\tNode is not created\n");
    return std::make_tuple(_hash,  false); 
}

std::tuple<ULL, bool> create_node(TT _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes = {}, std::unordered_map<ULL, Node>& hash_map = GNM)
{
    ULL _hash = calculate_hash(_func, _last_func, _cost, _depth, _xorable, _parent_hashes);
    return create_node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash, hash_map);
}

std::tuple<ULL, bool> create_node_no_upd(TT _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes = {}, std::unordered_map<ULL, Node>& hash_map = GNM)
{
    ULL _hash = calculate_hash(_func, _last_func, _cost, _depth, _xorable, _parent_hashes);
    return create_node_no_upd(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash, hash_map);
}

void register_nodes(std::unordered_map<ULL, Node>& hash_map)
{
    // for (auto & [hash, node] : hash_map)
    for (auto it_tmp = hash_map.begin(); it_tmp != hash_map.end(); it_tmp++)
    {
        // Move an item from map1 to map2
        // auto it = map1.find(2);
        // if (it != map1.end()) {
        //     map2.insert(std::move(*it));
        //     map1.erase(it);
        // }
        auto [hash, node] = *it_tmp;
        auto it_gnm = GNM.find(hash);
        bool node_is_new = (it_gnm == GNM.end());
        if (node_is_new) //the node is new
        {
            Node& best_any = GEA[node.func];

            if (node.xorable)
            {
                Node& best_xor = GEX[node.func];
                if (std::tie(node.lvl, node.cost) <= std::tie(best_any.lvl, best_any.cost))
                {
                    GNM.insert(std::move(*it_tmp)); hash_map.erase(it_tmp); //move the node
                    GEA[node.func] = GNM[hash];
                    GEX[node.func] = GNM[hash];
                    // fmt::print("\t\t\t\tNode is better than any existing function\n");
                }
                else if (std::tie(node.lvl, node.cost) <= std::tie(best_xor.lvl, best_xor.cost))
                {
                    GNM.insert(std::move(*it_tmp)); hash_map.erase(it_tmp); //move the node
                    GEX[node.func] = GNM[hash];
                    // fmt::print("\t\t\t\tNode is better than any existing xorable function\n");
                }
            }
            else 
            {
                // fmt::print("\t\t\t\tFunction is NOT xorable\n");
                if (std::tie(node.lvl, node.cost) <= std::tie(best_any.lvl, best_any.cost))
                {
                    GNM.insert(std::move(*it_tmp)); hash_map.erase(it_tmp); //move the node
                    GEA[node.func] = GNM[hash];
                    // fmt::print("\t\t\t\tNode is better than any existing function\n");
                }
            }
        }
        // else{
        //     fmt::print("\t\t\t\tNode already exists\n");
        // }
        // fmt::print("\t\t\t\tNode is not created\n");

    }
}

std::vector<TT> gen_pi_func(UI nvars)
{
    if (nvars == 1)
    {
        return {1};
    }
    TT power = (1 << (1 << (nvars - 1)));
    TT factor = power + 1;
    std::vector<TT> out = {power - 1};
    for (auto n : gen_pi_func(nvars - 1))
    {
        out.push_back(factor * n);
    }
    return out;
}

std::vector<ULL> select_depth(US min_depth, US max_depth)
{
    std::vector<ULL> out;
    for (auto & [key, n] : GNM)
    {
        if (n.depth >= min_depth && n.depth <= max_depth && n.func != 0 && n.func != ONES){
            out.push_back(key);
        }
    }
    return out;
}

UI node_cost(const Node& n1, const Node& n2, UI gate_cost)
{
    std::stack<ULL> stack;
    stack.push(n1.hash);
    stack.push(n2.hash);

    std::unordered_map<ULL, UI> ct_spl;
    std::unordered_set<ULL> non_splittable_nodes;

    // fmt::print("\t\t\tInitial cost {}\n", gate_cost);
    while (!stack.empty())
    {
        ULL n_hash = stack.top();
        stack.pop();
        Node& n = GNM[n_hash];
        // fmt::print("\t\t\tProcessing node {}\n", n.to_str());
        ct_spl[n_hash]++;
        if (ct_spl[n_hash] == 1)
        {
            gate_cost += COSTS[n.last_func];
            // fmt::print("\t\t\tFirst time, adding the cost of {} ({}). Total cost is {}\n", F2STR[n.last_func], COSTS[n.last_func], gate_cost);
            for (const auto& p_hash : n.parent_hashes)
            {
                Node& p = GNM[p_hash];
                stack.push(p_hash);
                // fmt::print("\t\t\tAdding parent node {} to stack\n", p.to_str());
            }
        }
        else
        {
            gate_cost += COSTS[fSPL];
            // fmt::print("\t\t\tNot first time, adding the cost of SPL ({}). Total cost is {}\n", COSTS[n.last_func], gate_cost);
        }
        if (n.last_func == fAND || n.last_func == fOR)
        {
            for (const auto& p_hash : n.parent_hashes)
            {
                non_splittable_nodes.emplace(p_hash);
            }
        }
    }

    for (const ULL & n_hash : non_splittable_nodes)
    {
        auto count = ct_spl[n_hash];
        if (ct_spl[n_hash] > 1)
        {
            Node& n = GNM[n_hash];
            auto last_func = n.last_func;
            gate_cost += COSTS[last_func] * (count - 1); // duplicate a gate
            if (last_func == fXOR)
            {
                gate_cost += (count - 1) * COSTS[fSPL]; // need to do twice more splittings for duplicating an XOR gate
            }
        }
    }
    return gate_cost;
}

std::tuple<ULL, bool> check_cb(Node& ni, Node& nj, std::vector<ULL>& fresh_nodes)
{
    fmt::print("\t\tA: {}\n", ni.to_str());
    fmt::print("\t\tB: {}\n", nj.to_str());
    TT func = ni.func | nj.func;
    if (func == ni.func || func == nj.func || func == 0 || func == ONES)
    {
        Node& ref = GEX[func];
        return std::make_tuple(ref.hash, false);
    }
    bool xorable = (func == (ni.func ^ nj.func));
    US last_func = fCB;
    US depth = ni.lvl * 3 + 1;
    UI cost = node_cost(ni, nj, COSTS[last_func]);
    auto [nhash, added] = create_node(func, fCB, cost, depth, xorable, {ni.hash, nj.hash});
    if (added)
    {
        Node& n =  GNM[nhash];
        fresh_nodes.push_back(nhash);
        fmt::print("\t\tC: {}\n", n .to_str());
    }
    return std::make_tuple(nhash, added);
}

std::tuple<ULL, bool> check_xor(Node& ni, Node& nj)
{
    fmt::print("\t\tA: {}\n", ni.to_str());
    fmt::print("\t\tB: {}\n", nj.to_str());
    TT func = ni.func ^ nj.func;
    // if (func == 0 || func == ONES)
    // {
    //     Node& ref = GEX[func];
    //     return std::make_tuple(ref.hash, false);
    // }
    US depth = ni.lvl * 3 + 2;
    UI cost = node_cost(ni, nj, COSTS[fXOR]);
    auto [nhash, added] = create_node(func, fXOR, cost, depth, true, {ni.hash, nj.hash});
    if (added)
    {
        Node& n =  GNM[nhash];
        fmt::print("\t\tC: {}\n", n .to_str());
    }
    return std::make_tuple(nhash, added);
}

std::tuple<ULL, bool, ULL, bool> check_and_or(Node& ni, Node& nj)
{
    fmt::print("\t\tA: {}\n", ni.to_str());
    fmt::print("\t\tB: {}\n", nj.to_str());
    US depth = ni.lvl * 3 + 3;
    TT func_and = ni.func & nj.func;
    TT func_or  = ni.func | nj.func;
    fmt::print("\t\t\t&: {}\n", func_and);
    fmt::print("\t\t\t&: {}\n", func_or );
    // if (func == 0 || func == ONES)
    // {
    //     Node& ref = GEX[func];
    //     return std::make_tuple(ref.hash, false);
    // }
    UI cost_and = node_cost(ni, nj, COSTS[fAND]);
    UI cost_or  = cost_and - COSTS[fAND] + COSTS[fOR];
    auto [hash_and, added_and] = create_node(func_and, fAND, cost_and, depth, true, {ni.hash, nj.hash});
    auto [hash_or , added_or ] = create_node(func_or , fOR , cost_or , depth, true, {ni.hash, nj.hash});
    if (added_and)
    {
        Node& n =  GNM[hash_and];
        fmt::print("\t\t&: {}\n", n .to_str());
    }
    if (added_or )
    {
        Node& n =  GNM[hash_or ];
        fmt::print("\t\t|: {}\n", n .to_str());
    }
    return std::make_tuple(hash_and, added_and, hash_or, added_or);
}

void remove_dominated()
{
    std::unordered_set<ULL> protected_nodes;
    std::stack<ULL> stack;
    for (auto & [hash, node] : GNM)
    {
        if (node.last_func == fDFF) continue;
        Node& best = node.xorable? GEX[node.func] : GEA[node.func];
        if (node.lvl < best.lvl || (node.lvl == best.lvl && node.cost <= best.cost) )
        {
            fmt::print("\t\t\t\tProtecting node {} (best is {}) [{} {}]\n", node.to_str(), best.to_str(), node.hash, best.hash);
            stack.push(hash);
        }
    } 

    while (!stack.empty())
    {
        ULL hash = stack.top();
        stack.pop();
        
        Node& node = GNM[hash];
        if (protected_nodes.find(hash) == protected_nodes.end())
        {
            for (ULL p_hash : node.parent_hashes)
            {
                stack.push(p_hash);
            }
        }
        protected_nodes.insert(hash);
    }
    fmt::print("\t\t\t\tProtecting {} nodes out of {}\n", protected_nodes.size(), GNM.size());

    std::unordered_set<ULL> to_be_removed;

    for (auto & [hash, node] : GNM)
    {
        if (node.last_func == fDFF) continue;
        if (protected_nodes.find(hash) != protected_nodes.end()) continue;
        Node& best = node.xorable? GEX[node.func] : GEA[node.func];
        if (node.lvl > best.lvl || (node.lvl == best.lvl && node.cost > best.cost) )
        {
            to_be_removed.insert(hash);
        }
    }

    for (auto hash : to_be_removed)
    {
        fmt::print("\t\t\t\tDestroying dominated node {}\n", GNM[hash].to_str());
        GNM.erase(hash);
    }   
}

void cb_generation(US lvl)
{
    std::vector<ULL> new_nodes = select_depth(lvl * 3, lvl * 3 + 1);
    std::vector<ULL> old_nodes;
    bool go_on;
    do {
        go_on = false;
        std::vector<ULL> fresh_nodes;

        // first, combine new nodes with old nodes
        fmt::print("\tProcessing old {} nodes with new {} nodes\n", old_nodes.size(), new_nodes.size());
        #pragma omp parallel for num_threads(10)
        for (UL k = 0u; k < new_nodes.size() * old_nodes.size(); k++)
        {
            UL i = k / old_nodes.size();
            UL j = k % old_nodes.size();
            fmt::print("\tCB ({}, {}, {}):\n", i, j, k);
            Node& ni =  GNM[new_nodes[i]];
            Node& nj =  GNM[old_nodes[j]];
            auto [nhash, added] = check_cb(ni, nj, fresh_nodes);
            go_on |= added;
        }
        // next, combine new nodes 
        fmt::print("\tProcessing new {} nodes\n", new_nodes.size());
        #pragma omp parallel for num_threads(10)
        for (UL k = 0u, n = new_nodes.size(); k < n * (n - 1) / 2; k++)
        {
            UL i = n - 2 - UL(sqrt(4*n*(n-1) - 8*k - 7)/2.0 - 0.5);
            UL j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
            fmt::print("\tCB ({}, {}, {}):\n", i, j, k);
            Node& ni =  GNM[new_nodes[i]];
            Node& nj =  GNM[new_nodes[j]];
            auto [nhash, added] = check_cb(ni, nj, fresh_nodes);
            go_on |= added;
        }
        old_nodes.insert(old_nodes.end(), new_nodes.begin(), new_nodes.end());
        new_nodes.clear();
        new_nodes.insert(new_nodes.end(), fresh_nodes.begin(), fresh_nodes.end());
        fresh_nodes.clear();
        remove_dominated();
    } while(go_on);
}

#if false

void cb_generation(US lvl)
{
    std::vector<ULL> new_nodes = select_depth(lvl * 3, lvl * 3 + 1);
    std::vector<ULL> old_nodes;
    bool go_on;
    do {
        go_on = false;
        std::vector<ULL> fresh_nodes;

        // first, combine new nodes with old nodes
        fmt::print("\tProcessing old {} nodes with new {} nodes\n", old_nodes.size(), new_nodes.size());
        std::for_each(
            std::execution::par,
            new_nodes.begin(),
            new_nodes.end(),
            [](auto&& item)
            {
                //do stuff with item
            }
        );


        #pragma omp parallel for num_threads(10)
        for (UL k = 0u; k < new_nodes.size() * old_nodes.size(); k++)
        {
            UL i = k / old_nodes.size();
            UL j = k % old_nodes.size();
            fmt::print("\tCB ({}, {}, {}):\n", i, j, k);
            Node& ni =  GNM[new_nodes[i]];
            Node& nj =  GNM[old_nodes[j]];
            auto [nhash, added] = check_cb(ni, nj, fresh_nodes);
            go_on |= added;
        }
        // next, combine new nodes 
        fmt::print("\tProcessing new {} nodes\n", new_nodes.size());
        #pragma omp parallel for num_threads(10)
        for (UL k = 0u, n = new_nodes.size(); k < n * (n - 1) / 2; k++)
        {
            UL i = n - 2 - UL(sqrt(4*n*(n-1) - 8*k - 7)/2.0 - 0.5);
            UL j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
            fmt::print("\tCB ({}, {}, {}):\n", i, j, k);
            Node& ni =  GNM[new_nodes[i]];
            Node& nj =  GNM[new_nodes[j]];
            auto [nhash, added] = check_cb(ni, nj, fresh_nodes);
            go_on |= added;
        }
        old_nodes.insert(old_nodes.end(), new_nodes.begin(), new_nodes.end());
        new_nodes.clear();
        new_nodes.insert(new_nodes.end(), fresh_nodes.begin(), fresh_nodes.end());
        fresh_nodes.clear();
        remove_dominated();
    } while(go_on);
}

std::unordered_map<ULL, Node> cb_generation_chunk(std::vector<ULL>& nodes)
{
    std::unordered_map<ULL, Node> tmp_map;

    for (UL k = 0u, n = nodes.size(); k < n * (n - 1) / 2; k++)
    {
        UL i = n - 2 - UL(sqrt(4*n*(n-1) - 8*k - 7)/2.0 - 0.5);
        UL j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
        fmt::print("\tCB ({}, {}, {}):\n", i, j, k);
        Node& ni =  GNM[nodes[i]];
        Node& nj =  GNM[nodes[j]];

        TT func = ni.func | nj.func;
        if (func == ni.func || func == nj.func || func == 0 || func == ONES)
        {
            continue;
        }
        bool xorable = (func == (ni.func ^ nj.func));
        US last_func = fCB;
        US depth = ni.lvl * 3 + 1;
        UI cost = node_cost(ni, nj, COSTS[last_func]);
        create_node_no_upd(func, fCB, cost, depth, xorable, {ni.hash, nj.hash}, tmp_map);
    }
    return tmp_map;
}

void cb_generation_parallel(US lvl)
{
    // Create a vector of threads
    bool go_on;
    do {
        go_on = false;
        std::vector<std::thread> threads;
        std::vector<ULL> nodes = select_depth(lvl * 3, lvl * 3 + 1);
        // Determine the size of each chunk
        const int chunk_size = nodes.size() / kNumThreads;  
        std::vector<std::unordered_map<ULL, Node>> results(kNumThreads);

        // Start the threads to process each chunk of data
        auto begin = nodes.begin();
        auto end = begin + chunk_size;
        for (int i = 0; i < kNumThreads; ++i) {
            if (i == kNumThreads - 1) {
            end = nodes.end();
            }
            std::vector<ULL> tmp_nodes;
            tmp_nodes.insert(tmp_nodes.end(), begin, end);
            std::thread t = [&tmp_nodes, results, i](){
                results[i] = cb_generation_chunk(tmp_nodes);
            };

            threads.emplace_back([&tmp_nodes, results, i](){
                results[i] = cb_generation_chunk(tmp_nodes);
            });
            begin = end;
            end = begin + chunk_size;
        }



        // Join the threads
        for (auto& t : threads) {
            t.join();
        }

    } while(go_on);
}
#endif 

void as_generation(US lvl)
{
    std::vector<ULL> nodes = select_depth(lvl * 3, lvl * 3 + 1);

    US tgt_depth = lvl * 3 + 2;
    fmt::print("\tDFF/NOT {} nodes\n", nodes.size());
    for (ULL hash : nodes)
    {
        Node& ni =  GNM[hash];
        force_create_node(ni.func, fDFF, ni.cost + COSTS[fDFF], tgt_depth, true, {hash});
        create_node(ni.func ^ ONES, fNOT, ni.cost + COSTS[fNOT], tgt_depth, true, {hash});
    }
    // xor-combine the nodes
    fmt::print("\tXOR {} nodes\n", nodes.size());
    #pragma omp parallel for num_threads(10)
    for (UL k = 0u, n = nodes.size(); k < n * (n - 1) / 2; k++)
    {
        UL i = n - 2 - UL(sqrt(4*n*(n-1) - 8*k - 7)/2.0 - 0.5);
        UL j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
        fmt::print("\tCB ({}, {}, {}):\n", i, j, k);
        Node& ni =  GNM[nodes[i]];
        Node& nj =  GNM[nodes[j]];
        if (!ni.xorable || !nj.xorable) continue;
        auto [nhash, added] = check_xor(ni, nj);
    }
    remove_dominated();
}
void sa_generation(US lvl)
{
    std::vector<ULL> nodes = select_depth(lvl * 3 + 2, lvl * 3 + 2);

    US tgt_depth = lvl * 3 + 3;
    // xor-combine the nodes
    fmt::print("\tAND/OR {} nodes\n", nodes.size());
    #pragma omp parallel for num_threads(10)
    for (UL k = 0u, n = nodes.size(); k < n * (n - 1) / 2; k++)
    {
        UL i = n - 2 - UL(sqrt(4*n*(n-1) - 8*k - 7)/2.0 - 0.5);
        UL j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
        fmt::print("\tAND/OR ({}, {}, {}):\n", i, j, k);
        Node& ni =  GNM[nodes[i]];
        Node& nj =  GNM[nodes[j]];
        if (!ni.xorable || !nj.xorable) continue;
        auto [hash_and, added_and, hash_or, added_or] = check_and_or(ni, nj);
    }
    remove_dominated();
}

int main() {

    for (TT func : gen_pi_func(NUM_VARS))
    {
        create_node(func, fPI, 0, 0, true);
    }
    create_node(   0, fPI, 0, 0, true);
    create_node(ONES, fPI, 0, 0, true);
    

    std::for_each(
        std::execution::par,
        foo.begin(),
        foo.end(),
        [](auto&& item)
        {
            //do stuff with item
        }
    );
    
    US depth = 1;
    for (US lvl = 0; lvl < 3; lvl ++)
    {
        cb_generation(lvl);
        as_generation(lvl);
        sa_generation(lvl);
    }
    return 0;
}
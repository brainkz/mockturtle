#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <tuple>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <array>
#include <fmt/format.h>
#include <thread>
#include "mockturtle/algorithms/nodes.hpp"
#include <algorithm>
#include <locale>

#include <iomanip>
#include <ctime>

#define continue_if(condition) if (condition) continue
#define break_if(condition) if (condition) break
#define shortcut_cost(init_cost) if (init_cost > limit) return std::make_tuple(INF32, false)

#define COMPUTE_IJ(k, N, i, j) \
    uint64_t i = N - 2 - static_cast<uint64_t>(std::sqrt(4*N*(N - 1) - 8*k - 7) / 2.0 - 0.5); \
    uint64_t j = k + i + 1 - N * (N - 1) / 2 + (N - i) * ((N - i) - 1) / 2;

inline void unique_merge(std::vector<uint64_t>& a, std::vector<uint64_t>& b, std::vector<uint64_t>& c)
{
    std::vector<uint64_t> c(a.size() + b.size());
    auto it = std::set_union(a.begin(), a.end(), b.begin(), b.end(), c.begin());
    c.resize(it - c.begin());
}

inline void unique_insert(std::vector<uint64_t>& v, const uint64_t& element) {
    if (std::find(v.begin(), v.end(), element) == v.end()) {
        v.push_back(element);
    }
}

std::string get_current_time_formatted() 
{
    std::time_t t = std::time(nullptr);
    std::tm* local_time = std::localtime(&t);

    std::stringstream ss;
    ss << std::put_time(local_time, "%Y%m%d_%H%M%S");
    return ss.str();
}

constexpr uint64_t MAX_SIZE = 20'000'000'000;
constexpr uint64_t NUM_THREADS = 100;

inline void join_threads(std::vector<std::thread>& threads)
{
    for (auto& t : threads) 
    {
        t.join();
    }
}
inline void check_threads(std::vector<std::thread>& threads)
{
    for (auto& t : threads) 
    {
        assert(!t.joinable());
    }
}

// std::unordered_map<std::array<uint8_t, 4>, std::unordered_map<uint64_t, Node>> GNM_global;
std::unordered_map<uint64_t, Node> GNM_global;
std::unordered_map<uint64_t, Node> GNM;
std::array<uint64_t, NUM_TT> GEA;
std::array<uint64_t, NUM_TT> GEX;
std::unordered_map<std::array<uint8_t, NUM_VARS>, std::array<uint64_t, NUM_TT>> GEA_global;
std::unordered_map<std::array<uint8_t, NUM_VARS>, std::array<uint64_t, NUM_TT>> GEX_global;
std::unordered_map<uint64_t, std::vector<uint64_t>> NODE_PIS; //maps each node to hashes of its PI-s
std::unordered_set<uint64_t> PI_HASHES;

std::array<std::array<uint8_t, 4>, 35> sets_of_levels {{
    {{0,0,0,0}}, 
    {{0,0,0,1}}, 
    {{0,0,0,2}}, 
    {{0,0,0,3}}, 
    {{0,0,0,4}}, 
    {{0,0,1,1}}, 
    {{0,0,1,2}}, 
    {{0,0,1,3}}, 
    {{0,0,1,4}}, 
    {{0,0,2,2}}, 
    {{0,0,2,3}}, 
    {{0,0,2,4}}, 
    {{0,0,3,3}}, 
    {{0,0,3,4}}, 
    {{0,0,4,4}}, 
    {{0,1,1,1}}, 
    {{0,1,1,2}}, 
    {{0,1,1,3}}, 
    {{0,1,1,4}}, 
    {{0,1,2,2}}, 
    {{0,1,2,3}}, 
    {{0,1,2,4}}, 
    {{0,1,3,3}}, 
    {{0,1,3,4}}, 
    {{0,1,4,4}}, 
    {{0,2,2,2}}, 
    {{0,2,2,3}}, 
    {{0,2,2,4}}, 
    {{0,2,3,3}}, 
    {{0,2,3,4}}, 
    {{0,2,4,4}}, 
    {{0,3,3,3}}, 
    {{0,3,3,4}}, 
    {{0,3,4,4}}, 
    {{0,4,4,4}}
}};

inline uint64_t ceil(uint64_t x, uint64_t y)
{
    return 1 + ((x - 1) / y);
}


inline Node& get_gea(TT func) 
{
    uint64_t hash = GEA[func._bits];
    return GNM[hash];
}
inline Node& get_gex(TT func)
{
    uint64_t hash = GEX[func._bits];
    return GNM[hash];
}

inline void check_node(const uint64_t hash)
{
    const Node & n = GNM[hash]; //.at(hash); // [hash]; //
    fmt::print("Checking node {} | num_parents = {} | {:L}\n", n.to_str(), n.parent_hashes.size(), hash);
    if (n.last_func == fNOFUNC)
    {
        assert(n.parent_hashes.empty());
        return;
    }
    else if (is_const(n.func))
    {
        assert(n.last_func == fPI);
    }
    assert(n.cost != INF32);
    assert(n.depth != INF8);
}

inline void check_nodes(const std::vector<uint64_t>& nodes)
{
    if (nodes.empty()) 
    {
        return;
    }
    for (const uint64_t hash : nodes)
    {
        check_node(hash);
    }
}

inline void check_GNM()
{
    uint64_t violation_ct = 0u;
    for (auto [hash, n] : GNM)
    {
        fmt::print("Checking {} \n", n.to_str());
        if (n.last_func == fNOFUNC)
        {
            violation_ct++;
        }
        else
        {
            assert(kitty::is_const0(n.func) || n.last_func == fPI);
            assert(n.cost != INF32);
            assert(n.depth != INF8);
        }
    }
}

inline void print_GNM()
{
    for (auto [hash, n] : GNM)
    {
        fmt::print("{}\n", n.to_str());
    }
}

std::tuple<uint64_t, bool> force_create_node(TT _func, uint8_t _last_func, uint32_t _cost, uint8_t _depth, bool _xorable, std::vector<uint64_t> _parent_hashes, uint64_t _hash)
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

std::tuple<uint64_t, bool> force_create_node(TT _func, uint8_t _last_func, uint32_t _cost, uint8_t _depth, bool _xorable, std::vector<uint64_t> _parent_hashes = {})
{
    uint64_t _hash = calculate_hash(_func, _last_func, _cost, _depth, _xorable, _parent_hashes);
    return force_create_node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
}

std::tuple<uint64_t, bool> create_node(TT _func, uint8_t _last_func, uint32_t _cost, uint8_t _depth, bool _xorable, std::vector<uint64_t> _parent_hashes, uint64_t _hash)
//, std::unordered_map<uint64_t, Node>& hash_map = GNM)
{
    auto it = GNM.find(_hash);
    bool node_is_new = (it == GNM.end());
    if (node_is_new) //the node is new
    {
        // fmt::print("\t\t\t\tChecking for creation {}\n", Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash).to_str());
        uint8_t lvl = _depth / 3;
        const Node& best_any = get_gea(_func);
        // fmt::print("\t\t\t\tConsidering {}: {}, {}\n", _func,          lvl,         _cost);
        // fmt::print("\t\t\t\tBest A      {}: {}, {}\n", _func, best_any.lvl, best_any.cost);

        if (_xorable)
        {
            // const Node& best_xor = GEX[_func];
            const Node& best_xor = get_gex(_func);    
            // fmt::print("\t\t\t\tBest X      {}: {}, {}\n", _func, best_xor.lvl, best_xor.cost);
            if (std::tie(lvl, _cost ) <= std::tie(best_any.lvl, best_any.cost))
            {
                GNM[_hash] = Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
                GEA[_func._bits] = _hash; //hash_map[_hash];
                GEX[_func._bits] = _hash; //hash_map[_hash];
                // fmt::print("\t\t\t\tNode is better than any existing function\n");
                return std::make_tuple(_hash, true);
            }
            else if (std::tie(lvl, _cost ) <= std::tie(best_xor.lvl, best_xor.cost))
            {
                GNM[_hash] = Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
                GEX[_func._bits] = _hash; //hash_map[_hash];
                // fmt::print("\t\t\t\tNode is better than any existing xorable function\n");
                return std::make_tuple(_hash, true);
            }
        }
        else 
        {
            // fmt::print("\t\t\t\tFunction is NOT xorable\n");
            if (std::tie(lvl, _cost ) <= std::tie(best_any.lvl, best_any.cost))
            {
                GNM[_hash] = Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
                GEA[_func._bits] = _hash; //hash_map[_hash];
                // fmt::print("\t\t\t\tNode is better than any existing function\n");
                return std::make_tuple(_hash, true);
            }
        }
    }
    else{
        // fmt::print("\t\t\t\tNode already exists\n");
    }
    // fmt::print("\t\t\t\tNode is not created\n");
    return std::make_tuple(_hash,  false); 
}

std::tuple<uint64_t, bool> create_node(TT _func, uint8_t _last_func, uint32_t _cost, uint8_t _depth, bool _xorable, std::vector<uint64_t> _parent_hashes = {})
{
    uint64_t _hash = calculate_hash(_func, _last_func, _cost, _depth, _xorable, _parent_hashes);
    return create_node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
}

std::tuple<uint64_t, bool> create_node_no_upd(TT _func, uint8_t _last_func, uint32_t _cost, uint8_t _depth, bool _xorable, std::vector<uint64_t> _parent_hashes, uint64_t _hash, std::unordered_map<uint64_t, Node>& hash_map = GNM)
{
    auto it = hash_map.find(_hash);
    bool node_is_new = (it == hash_map.end());
    if (node_is_new) //the node is new
    {
        // fmt::print("\t\t\t\tChecking for creation {}\n", Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash).to_str());
        uint8_t lvl = _depth / 3;
        const Node& best_any = get_gea(_func); // GEA[_func];
        // fmt::print("\t\t\t\tConsidering {}: {}, {}\n", _func,          lvl,         _cost);
        // fmt::print("\t\t\t\tBest A      {}: {}, {}\n", _func, best_any.lvl, best_any.cost);

        if (_xorable)
        {
            // fmt::print("\t\t\t\tFunction is xorable\n");
            const Node& best_xor = get_gex(_func); // GEX[_func];
            if (std::tie(lvl, _cost ) <= std::tie(best_any.lvl, best_any.cost))
            {
                hash_map[_hash] = Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
                // fmt::print("\t\t\t\tNode is better than any existing function\n");
                return std::make_tuple(_hash, true);
            }
            else if (std::tie(lvl, _cost ) <= std::tie(best_xor.lvl, best_xor.cost))
            {
                hash_map[_hash] = Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
                // fmt::print("\t\t\t\tNode is better than any existing xorable function\n");
                return std::make_tuple(_hash, true);
            }
        }
        else 
        {
            // fmt::print("\t\t\t\tFunction is NOT xorable\n");
            if (std::tie(lvl, _cost ) <= std::tie(best_any.lvl, best_any.cost))
            {
                hash_map[_hash] = Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
                // fmt::print("\t\t\t\tNode is better than any existing function\n");
                return std::make_tuple(_hash, true);
            }
        }
    }
    else{
        // fmt::print("\t\t\t\tNode already exists\n");
    }
    // fmt::print("\t\t\t\tNode is not created\n");
    return std::make_tuple(_hash,  false); 
}

std::tuple<uint64_t, bool> create_node_no_upd(TT _func, uint8_t _last_func, uint32_t _cost, uint8_t _depth, bool _xorable, std::vector<uint64_t> _parent_hashes = {}, std::unordered_map<uint64_t, Node>& hash_map = GNM)
{
    uint64_t _hash = calculate_hash(_func, _last_func, _cost, _depth, _xorable, _parent_hashes);
    return create_node_no_upd(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash, hash_map);
}

void register_nodes(std::unordered_map<uint64_t, Node>& hash_map)
{
    // for (auto & [hash, node] : hash_map)
    for (auto it_tmp = hash_map.begin(); it_tmp != hash_map.end(); it_tmp++)
    {
        auto [hash, node] = *it_tmp;
        auto it_gnm = GNM.find(hash);
        bool node_is_new = (it_gnm == GNM.end());
        if (node_is_new) //the node is new
        {
            const Node& best_any = get_gea(node.func); //GEA[node.func];

            if (node.xorable)
            {
                const Node& best_xor = get_gex(node.func);// GEX[node.func];
                if (std::tie(node.lvl, node.cost) <= std::tie(best_any.lvl, best_any.cost))
                {
                    GNM.insert(std::move(*it_tmp)); hash_map.erase(it_tmp); //move the node
                    GEA[node.func._bits] = hash; //GNM[hash];
                    GEX[node.func._bits] = hash; //GNM[hash];
                    // fmt::print("\t\t\t\tNode is better than any existing function\n");
                }
                else if (std::tie(node.lvl, node.cost) <= std::tie(best_xor.lvl, best_xor.cost))
                {
                    GNM.insert(std::move(*it_tmp)); hash_map.erase(it_tmp); //move the node
                    GEX[node.func._bits] = hash; //GNM[hash];
                    // fmt::print("\t\t\t\tNode is better than any existing xorable function\n");
                }
            }
            else 
            {
                // fmt::print("\t\t\t\tFunction is NOT xorable\n");
                if (std::tie(node.lvl, node.cost) <= std::tie(best_any.lvl, best_any.cost))
                {
                    GNM.insert(std::move(*it_tmp)); hash_map.erase(it_tmp); //move the node
                    GEA[node.func._bits] = hash; //GNM[hash];
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

std::vector<uint64_t> select_depth(uint8_t min_depth, uint8_t max_depth)
{
    std::vector<uint64_t> out;
    for (auto & [key, n] : GNM)
    {
        if (n.depth >= min_depth && n.depth <= max_depth && !kitty::is_const0(n.func) && !kitty::is_const0(~n.func))
        {
            out.push_back(key);
        }
    }
    return out;
}

std::vector<std::pair<uint64_t, uint64_t>> divide_range(uint64_t N, uint64_t k) 
{
    std::vector<std::pair<uint64_t, uint64_t>> ranges(k);
    uint64_t base_length = N / k;
    uint64_t extra = N % k;

    uint64_t start = 0;
    for (uint64_t i = 0; i < k; ++i) {
        uint64_t end = start + base_length - 1;
        if (extra > 0) {
            end += 1;
            extra -= 1;
        }
        ranges[i] = std::make_pair(start, end + 1);
        start = end + 1;
    }

    return ranges;
}

uint32_t node_cost(const Node& n1, const Node& n2, uint32_t gate_cost)
{
    std::vector<uint64_t> stack;
    stack.push_back(n1.hash);
    stack.push_back(n2.hash);

    std::unordered_map<uint64_t, uint8_t> ct_spl;
    std::unordered_set<uint64_t> non_splittable_nodes;

    // fmt::print("\t\t\tInitial cost {}\n", gate_cost);
    while (!stack.empty())
    {
        const uint64_t n_hash = stack.back();
        stack.pop_back();
        const Node& n = GNM[n_hash]; //.at(n_hash); // [n_hash];
        // fmt::print("\t\t\tProcessing node {}\n", n.to_str());
        ct_spl[n_hash]++;
        if (ct_spl[n_hash] == 1)
        {
            gate_cost += COSTS[n.last_func];
            // fmt::print("\t\t\tFirst time, adding the cost of {} ({}). Total cost is {}\n", F2STR[n.last_func], COSTS[n.last_func], gate_cost);
            for (const uint64_t & p_hash : n.parent_hashes)
            {
                stack.push_back(p_hash); // Node& p = GNM[p_hash]; fmt::print("\t\t\tAdding parent node {} to stack\n", p.to_str());
            }
            if (n.last_func == fAND || n.last_func == fOR)
            {
                for (const uint64_t & p_hash : n.parent_hashes)
                {
                    non_splittable_nodes.emplace(p_hash);
                }
            }
        }
        else
        {
            gate_cost += COSTS[fSPL]; // fmt::print("\t\t\tNot first time, adding the cost of SPL ({}). Total cost is {}\n", COSTS[n.last_func], gate_cost);
        }
    }

    for (const uint64_t & n_hash : non_splittable_nodes)
    {
        uint8_t count = ct_spl[n_hash];
        if (ct_spl[n_hash] > 1)
        {
            const Node & n = GNM[n_hash]; //.at(n_hash); //[n_hash];
            uint8_t last_func = n.last_func;
            gate_cost += COSTS[last_func] * (count - 1); // duplicate a gate
            if (last_func == fXOR)
            {
                gate_cost += (count - 1) * COSTS[fSPL]; // need to do twice more splittings for duplicating an XOR gate
            }
        }
    }
    return gate_cost;
}

std::tuple<uint32_t, std::vector<std::pair<uint64_t, uint8_t>>, std::vector<uint64_t>> node_cost_single(const uint64_t h1, uint32_t gate_cost, const uint8_t depth)
{
    std::vector<uint64_t> stack { h1 };
    stack.reserve(2 * depth);

    std::vector<std::pair<uint64_t, uint8_t>> ct_spl;
    ct_spl.reserve(2 * depth);
    std::vector<uint64_t> non_splittable_nodes;

    while (!stack.empty())
    {
        const uint64_t n_hash = stack.back();
        stack.pop_back();
        const Node& n = GNM[n_hash]; //.at(n_hash); //[n_hash];

        auto it = std::find_if(ct_spl.begin(), ct_spl.end(), [n_hash](const std::pair<uint64_t, uint8_t>& p) { return p.first == n_hash; });
        if (it == ct_spl.end())
        {
            ct_spl.emplace_back(n_hash, 1);
            gate_cost += COSTS[n.last_func];

            for (const uint64_t p_hash : n.parent_hashes)
            {
                stack.push_back(p_hash);
            }
            if (n.last_func == fAND || n.last_func == fOR)
            {
                for (const uint64_t p_hash : n.parent_hashes)
                {
                    if (std::find(non_splittable_nodes.begin(), non_splittable_nodes.end(), p_hash) == non_splittable_nodes.end())
                    {
                        non_splittable_nodes.push_back(p_hash);
                    }
                }
            }
        }
        else
        {
            it->second++;
            gate_cost += COSTS[fSPL];
        }
    }

    return std::make_tuple(gate_cost, ct_spl, non_splittable_nodes);
}

// assumes precomputed gate_cost, ct_spl, and non_splittable_nodes
std::tuple<uint32_t, bool> node_cost(uint64_t h2, uint32_t init_cost, std::vector<std::pair<uint64_t, uint8_t>> ct_spl, std::vector<uint64_t> non_splittable_nodes, const uint32_t limit, const uint8_t depth, bool verbose=false)
{
    std::vector<uint64_t> stack { h2 };
    stack.reserve( 2 * depth ); 
    if (verbose) fmt::print("\t\t\tProcessing node {}, limit = {}\n", GNM[h2].to_str(), limit);
    shortcut_cost(init_cost);
    while (!stack.empty())
    {
        const uint64_t n_hash = stack.back();
        stack.pop_back();
        // if (verbose) fmt::print("\t\t\tProcessing node {}\n", n.to_str());

        auto it = std::find_if(ct_spl.begin(), ct_spl.end(), [n_hash](const std::pair<uint64_t, uint8_t>& p) { return p.first == n_hash; });
        if (it == ct_spl.end())
        {
            ct_spl.emplace_back(n_hash, 1);
            it = std::prev(ct_spl.end());
        }
        else
        {
            it->second++;
        }

        if (verbose) fmt::print("\t\t\tIncremented the ct_spl for {}\n", n_hash);
        if (it->second == 1)
        {
            // assert(GNM.find(n_hash) != GNM.end());
            const Node & n = GNM[n_hash]; // .at(n_hash); //[n_hash];
            if (verbose) fmt::print("\t\t\tAccessing the cost of {}\n", n.last_func);
            init_cost += COSTS[n.last_func];
            shortcut_cost(init_cost);     
            if (verbose) fmt::print("\t\t\tFirst time, adding the cost of {} ({}). Total cost is {}\n", F2STR[n.last_func], COSTS[n.last_func], init_cost);
            for (const uint64_t p_hash : n.parent_hashes)
            {
                // assert(GNM.find(p_hash) != GNM.end());
                const Node& p = GNM[p_hash]; // .at(p_hash); //[p_hash];
                stack.push_back(p_hash);
                if (verbose) fmt::print("\t\t\tAdding parent node {} to stack\n", p.to_str());
            }
            if (n.last_func == fAND || n.last_func == fOR)
            {
                for (const uint64_t p_hash : n.parent_hashes)
                {
                    if (std::find(non_splittable_nodes.begin(), non_splittable_nodes.end(), p_hash) == non_splittable_nodes.end())
                    {
                        non_splittable_nodes.push_back(p_hash);
                    }
                }
            }
        }
        else
        {
            init_cost += COSTS[fSPL];
            shortcut_cost(init_cost);  
            if (verbose) fmt::print("\t\t\tNot first time, adding the cost of SPL ({}). Total cost is {}\n", COSTS[fSPL], init_cost);
        }
    }

    for (const uint64_t & n_hash : non_splittable_nodes)
    {
        auto it = std::find_if(ct_spl.begin(), ct_spl.end(), [n_hash](const std::pair<uint64_t, uint8_t>& p) { return p.first == n_hash; });
        if (it != ct_spl.end() && it->second > 1)
        {
            const Node& n = GNM[n_hash]; // .at(n_hash); //[n_hash];
            auto last_func = n.last_func;
            init_cost += COSTS[last_func] * (it->second - 1); // duplicate a gate
            if (last_func == fXOR)
            {
                init_cost += (it->second - 1) * COSTS[fSPL]; // need to do twice more splittings for duplicating an XOR gate
            }
            shortcut_cost(init_cost);  
        }
    }
    if (verbose) fmt::print("final cost is {} and status is {}\n", init_cost, true);
    return std::make_tuple(init_cost, true);
}

void thread_old_new(const std::vector<uint64_t>& hashes_i, const std::vector<uint64_t>& hashes_j, std::vector<std::vector<TT>>& funcs, std::vector<std::vector<bool>>& xorables, std::vector<std::vector<uint32_t>>& costs, const uint64_t start_row, const uint64_t end_row, const uint32_t gate_cost, const uint64_t offset, const uint8_t tgt_lvl) 
{
    // const int M = hashes_A.size();
    const int N = hashes_j.size();
    fmt::print("\t\tProcessing old - new pairs between #{:L} and #{:L} (out of total {:L})\n", start_row, end_row, hashes_i.size());
    for (uint64_t i = start_row; i < end_row; ++i) {
        const uint64_t hi = hashes_i[i];
        const Node & ni = GNM[hi]; // .at(hi); // [hi];
        auto [init_cost, ct_spl, non_splittable_nodes] = node_cost_single(hi, gate_cost, tgt_lvl*3);
        // fmt::print("Preliminary calculations: init_cost = {}, #predecessors = {}, #non-spllittables = {}\n", init_cost, ct_spl.size(), non_splittable_nodes.size());
        for (int j = 0; j < N; ++j) 
        {
            const uint64_t hj = hashes_j[j];
            const Node & nj = GNM[hj]; // .at(hj); // [hj];

            // fmt::print("\tCost calculation: ni = {}, nj = {}\n", ni.to_str(), nj.to_str());

            const TT func = ni.func | nj.func;
            continue_if (kitty::is_const0(func) || kitty::is_const0(~func) || func == ni.func || func == nj.func);
            const bool xorable = ((ni.func ^ nj.func) == func) & ni.xorable & nj.xorable;
            const Node & best_n = xorable ? get_gex(func) : get_gea(func);
            if (best_n.lvl >= tgt_lvl)
            {
                auto [cost, status] = node_cost(hj, init_cost, ct_spl, non_splittable_nodes, best_n.cost, tgt_lvl*3); // costs[i][j] = node_cost(ni, nj, gate_cost);
                // fmt::print("\tFinal calculation: cost = {}, status = {}\n", cost, status);
                funcs[i-offset][j] = func;
                xorables[i-offset][j] = xorable;
                costs[i-offset][j] = cost;
            }
        }
    }
}

void thread_old_new_wrapper(const std::vector<uint64_t>& old_nodes, const std::vector<uint64_t>& new_nodes, const uint8_t depth, std::vector<uint64_t>& fresh_nodes, uint64_t num_threads = NUM_THREADS)
{
    const uint64_t M = old_nodes.size();
    if (M == 0)
    {
        return;
    }
    const uint64_t N = new_nodes.size();
    const uint8_t tgt_lvl = depth / 3;
    uint64_t chunk_size = M / num_threads;
    if (chunk_size == 0)
    {
        num_threads = 1;
        chunk_size = M;
    }

    if (M * N < MAX_SIZE)
    {

        std::vector<std::vector<TT>> funcs(M, std::vector<TT>(N, const0));
        std::vector<std::vector<bool>> xorables(M, std::vector<bool>(N, false));
        std::vector<std::vector<uint32_t>> costs(M, std::vector<uint32_t>(N, INF32));

        // fmt::print("\t\tCalculating func, xorable, and cost for {} old nodes and {} new nodes. Threads: {}; Chunk size: {} \n", M, N, num_threads, chunk_size);
        std::vector<std::thread> threads;
        // threads.reserve(num_threads);
        for (uint64_t i = 0; i < num_threads; i++) {
            const uint64_t start_row = i * chunk_size;
            const uint64_t end_row = (i == num_threads - 1) ? M : (i + 1) * chunk_size;
            threads.push_back(
                std::thread(thread_old_new, std::ref(old_nodes), std::ref(new_nodes), std::ref(funcs), std::ref(xorables), std::ref(costs), start_row, end_row, COSTS[fCB], 0, tgt_lvl)
            );
        }
        join_threads(threads);
        check_threads(threads);
        for (uint64_t i = 0u; i < M; i++)
        {
            // fmt::print("Combining {} out of {} ({:f}\%) \n", (i+1), old_nodes.size(), 100 * (float)(i+1) / (float)old_nodes.size());
            uint64_t hi = old_nodes[i];
            for (uint64_t j = 0u; j < N; j++)
            {
                uint64_t hj = new_nodes[j];
                uint32_t cost = costs[i][j];
                continue_if (cost >= INF32);
                auto [nhash, added] = create_node(funcs[i][j], fCB, cost, depth, xorables[i][j], {hi, hj});
                if (added)
                {
                    fresh_nodes.push_back(nhash);
                }
            }
        }
    }
    else //need to split into chunks
    {
        // case 1 - row size is not too large - should in principle always be the case for 4-input UI-s
        // In this case, divide the combinations in one dimension - split by rows
        fmt::print("\t\tToo many combinations between {:L} old nodes and {:L} new nodes to fit into memory (total {:L}). Splitting... \n", M, N, M*N);
        if (N < MAX_SIZE / NUM_THREADS)
        {
            uint64_t nrows_per_thread = MAX_SIZE / N / NUM_THREADS;
            // fmt::print("\t\tProcessing {} rows per thread with {} threads... \n", nrows_per_thread, NUM_THREADS);
            uint64_t Mlocal = nrows_per_thread * NUM_THREADS;

            for (uint64_t start = 0u; start < M; start += Mlocal)
            {
                uint64_t end = std::min(start + Mlocal, M);
                std::vector<std::vector<TT>> funcs(Mlocal, std::vector<TT>(N, const0));
                std::vector<std::vector<bool>> xorables(Mlocal, std::vector<bool>(N, false));
                std::vector<std::vector<uint32_t>> costs(Mlocal, std::vector<uint32_t>(N, INF32));

                std::vector<std::thread> threads;
                // threads.reserve(num_threads);
                for (uint64_t i = 0u, start_row = start; start_row < end; i++, start_row += nrows_per_thread)
                {
                    uint64_t end_row = std::min(start_row + nrows_per_thread, end);
                    threads.push_back(
                        std::thread(thread_old_new, std::ref(old_nodes), std::ref(new_nodes), std::ref(funcs), std::ref(xorables), std::ref(costs), start_row, end_row, COSTS[fCB], start, tgt_lvl)
                    );
                }
                join_threads(threads);
                check_threads(threads);
                fmt::print("Combining results\n");
                for (auto row = start; row < end; row++)
                {
                    // fmt::print("Combining {} out of {} ({:f}\%) \n", (i+1), old_nodes.size(), 100 * (float)(i+1) / (float)old_nodes.size());
                    uint64_t i = row - start;
                    uint64_t hi = old_nodes[row];
                    // const Node& ni =  GNM[hi];
                    for (uint64_t j = 0u; j < N; j++)
                    {
                        uint64_t hj = new_nodes[j];
                        uint32_t cost = costs[i][j];
                        TT func = funcs[i][j];
                        continue_if (cost >= INF32) ;
                        // const Node& nj =  GNM[hj];
                        auto [nhash, added] = create_node(func, fCB, cost, depth, xorables[i][j], {hi, hj});
                        if (added)
                        {
                            fresh_nodes.push_back(nhash);
                        }
                    }
                }
            }
        }
    }
}

void threaded_new_new(const std::vector<uint64_t>& hashes, std::vector<TT>& funcs, std::vector<bool>& xorables, std::vector<uint32_t>& costs, const uint64_t start_k, const uint64_t end_k, const uint32_t gate_cost, const uint64_t offset, const uint8_t tgt_lvl) 
{
    // fmt::print("\t\tFunction is called\n");
    const uint64_t N = hashes.size();
    fmt::print("\t\tCB: Iterating between {:L} and {:L}\n", start_k, end_k);
    uint64_t old_i;
    uint32_t init_cost;
    std::vector<std::pair<uint64_t, uint8_t>> ct_spl;
    std::vector<uint64_t> non_splittable_nodes;
    // fmt::print("\t\tProcessing new - new combinations between #{} and #{}\n", start_k, end_k);
    for (uint64_t k = start_k; k < end_k; k++)
    {   
        COMPUTE_IJ(k, N, i, j);
        const uint64_t hi = hashes[i]; // fmt::print("\t\t hash_i={}\n", hi);
        const Node & ni = GNM[hi]; // .at(hi); //[hi]; // fmt::print("\t\t Retrieved ni={}\n", ni.to_str());
        if (i != old_i || k == start_k) // update [init_cost, ct_spl, non_splittable_nodes] if the row has changed
        {
            // fmt::print("\tUpdating\n");
            auto [init_cost_tmp, ct_spl_tmp, non_splittable_nodes_tmp] = node_cost_single(hi, gate_cost, tgt_lvl*3);
            init_cost = init_cost_tmp;
            ct_spl = ct_spl_tmp;
            non_splittable_nodes = non_splittable_nodes_tmp;
        }

        // fmt::print("\t\tk={}, i={}, j={}\n", k, i, j);
        const uint64_t hj = hashes[j]; // fmt::print("\t\t hash_j={}\n", hj);
        const Node & nj = GNM[hj]; // .at(hj); // [hj]; // fmt::print("\t\t Retrieved nj={}\n", nj.to_str());

        // fmt::print("\t\t\tCalculating cost for k {}:\n\t\t\t\t#{}: {}\n\t\t\t\t#{}: {}\n", k, i, ni.to_str(), j, nj.to_str());

        const TT func = ni.func | nj.func;
        continue_if (kitty::is_const0(func) || kitty::is_const0(~func) || func == ni.func || func == nj.func);
        // bool xorable = (ni.func ^ nj.func) == func; // fmt::print("\t\t Xorable status = {}\n", xorables[k]);
        const bool xorable = ((ni.func ^ nj.func) == func) & ni.xorable & nj.xorable;
        
        const Node & best_n = xorable ? get_gex(func) : get_gea(func);
        if (best_n.lvl >= tgt_lvl)
        {
            auto [cost, status] = node_cost(hj, init_cost, ct_spl, non_splittable_nodes, best_n.cost, tgt_lvl*3); // costs[i][j] = node_cost(ni, nj, gate_cost);
            // fmt::print("\t\t\t\t Result: cost {}, success {}\n", cost, status);
            if (status)
            {
                costs[k - offset] = cost;
                xorables[k - offset] = xorable;
                funcs[k - offset] = func; // fmt::print("\t\t Wrote func = {0:016b} = {0:04x}\n",  funcs[k]);
            }
        }
        old_i = i;
    }
}


void threaded_new_new_wrapper(const std::vector<uint64_t>& new_nodes, const uint8_t depth, std::vector<uint64_t>& fresh_nodes, uint64_t num_threads = NUM_THREADS)
{
    const uint64_t N = new_nodes.size();
    if (N <= 1)
    {
        return;
    }
    const uint64_t Ncombs = N * (N - 1) / 2;
    const uint8_t tgt_lvl = depth / 3;

    if (Ncombs < MAX_SIZE)
    {
        uint64_t chunk_size = Ncombs / (num_threads - 1);
        if (chunk_size == 0)
        {
            num_threads = 1;
            chunk_size = Ncombs;
        }
        
        std::vector<TT> funcs(Ncombs, num2tt(0));
        std::vector<bool> xorables(Ncombs, false);
        std::vector<uint32_t> costs(Ncombs, INF32);

        std::vector<std::thread> threads;
        // threads.reserve(num_threads);
        for (uint64_t start_k = 0u, i = 0u; start_k < Ncombs; start_k += chunk_size, i++)
        {
            uint64_t end_k = std::min(start_k + chunk_size, Ncombs);
            // fmt::print("\t\t\tcreating thread from {} to {} (out of {}) \n", start_k, end_k, Ncombs);
            threads.push_back(
                std::thread(threaded_new_new, std::ref(new_nodes), std::ref(funcs), std::ref(xorables), std::ref(costs), start_k, end_k, COSTS[fCB], 0, tgt_lvl)
            );
        }
        join_threads(threads);
        check_threads(threads);

        fmt::print("\t\tCombining new - new pairs\n");

        for (uint64_t k = 0u; k < Ncombs; k++)
        {
            COMPUTE_IJ(k, N, i, j);
            uint32_t cost = costs[k];
            // fmt::print("\tCombining: {} (#{}) and {} (#{}), cost {}, func {:016b}, xorable {}\n", GNM[new_nodes[i]].to_str(), i, GNM[new_nodes[j]].to_str(), j, cost, funcs[k], xorables[k]);
            // fmt::print("\tBest any l-{}, c-{}, best xorable l-{}, c-{}\n", get_gea(funcs[k]).lvl, get_gea(funcs[k]).cost, get_gex(funcs[k]).lvl, get_gex(funcs[k]).cost);
            continue_if(cost >= INF32);
            auto [nhash, added] = create_node(funcs[k], fCB, costs[k], depth, xorables[k], {new_nodes[i], new_nodes[j]});
            if (added)
            {
                // fmt::print("\t\t added node\n");
                fresh_nodes.push_back(nhash);
            }
        }
    }
    else //too many combinations to fit into memory
    {
        fmt::print("\t\tToo many pairs among {:L} new nodes to fit into memory (total {:L}). Splitting... \n", N, Ncombs);

        uint64_t NUM_MACRO_CHUNKS = ceil(Ncombs, MAX_SIZE);
        fmt::print("\t\tInto {} chunks \n", NUM_MACRO_CHUNKS);
        for (auto [START, END]  : divide_range(Ncombs, NUM_MACRO_CHUNKS))
        {
            uint64_t chunk_size = END - START;
            std::vector<TT> funcs(chunk_size, num2tt(0));
            std::vector<bool> xorables(chunk_size, false);
            std::vector<uint32_t> costs(chunk_size, INF32);
            std::vector<std::thread> threads;
            threads.reserve(num_threads);

            for (auto [start, end] : divide_range(chunk_size, NUM_THREADS))
            {
                uint64_t start_k = start + START;
                uint64_t end_k   = end   + START;
                threads.push_back(
                    std::thread(threaded_new_new, std::ref(new_nodes), std::ref(funcs), std::ref(xorables), std::ref(costs), start_k, end_k, COSTS[fCB], START, tgt_lvl)
                );
            }
            join_threads(threads);
            check_threads(threads);

            fmt::print("\t\tCombining NEW-NEW CB results\n");
            for (uint64_t k = START; k < END; k++)
            {
                COMPUTE_IJ(k, N, i, j);
                uint32_t cost = costs[k - START];
                TT func = funcs[k-START];
                continue_if (cost >= INF32);
                auto [nhash, added] = create_node(func, fCB, cost, depth, xorables[k-START], {new_nodes[i], new_nodes[j]});
                if (added)
                {
                    fresh_nodes.push_back(nhash);
                }
            }
        }
    }
}

void threaded_xor(const std::vector<uint64_t>& hashes, std::vector<TT>& funcs, std::vector<uint32_t>& costs, const uint64_t start_k, const uint64_t end_k, const uint32_t gate_cost, const uint64_t offset, const uint8_t tgt_lvl) 
{
    // fmt::print("\t\tFunction is called\n");
    const uint64_t N = hashes.size();
    fmt::print("\t\tXOR: Iterating between {:L} and {:L}\n", start_k, end_k);
    uint64_t old_i;
    uint32_t init_cost;
    std::vector<std::pair<uint64_t, uint8_t>> ct_spl;
    std::vector<uint64_t> non_splittable_nodes;
    for (uint64_t k = start_k; k < end_k; k++)
    {   
        COMPUTE_IJ(k, N, i, j);
        const uint64_t hi = hashes[i];     // fmt::print("\t\t hash_i={}\n", hi);
        const Node & ni = GNM[hi]; // .at(hi); // [hi];    // fmt::print("\t\t Retrieved ni={}\n", ni.to_str());

        if (i != old_i || k == start_k) // update [init_cost, ct_spl, non_splittable_nodes] if the row has changed
        {
            // fmt::print("\tUpdating\n");
            auto [init_cost_tmp, ct_spl_tmp, non_splittable_nodes_tmp] = node_cost_single(hi, gate_cost, tgt_lvl*3);
            init_cost = init_cost_tmp;
            ct_spl = ct_spl_tmp;
            non_splittable_nodes = non_splittable_nodes_tmp;
        }

        // fmt::print("\t\tk={}, i={}, j={}\n", k, i, j);
        const uint64_t hj = hashes[j];     // fmt::print("\t\t hash_j={}\n", hj);
        const Node & nj = GNM[hj]; // .at(hj); // [hj];    // fmt::print("\t\t Retrieved nj={}\n", nj.to_str());

        TT func = ni.func ^ nj.func;
        const Node & best_n = get_gex(func);
        if (best_n.lvl >= tgt_lvl) 
        {
            auto [cost, status] = node_cost(hj, init_cost, ct_spl, non_splittable_nodes, best_n.cost, tgt_lvl*3); // costs[i][j] = node_cost(ni, nj, gate_cost);
            if (status)
            {
                funcs[k - offset] = func;  // fmt::print("\t\t Wrote func = {0:016b} = {0:04x}\n",  funcs[k]);
                costs[k - offset] = cost;   // fmt::print("\t\t Cost = {}\n", costs[k]);
            }
        }
        old_i = i;
    }
}

void threaded_xor_wrapper(const std::vector<uint64_t>& nodes, const uint8_t depth, uint64_t num_threads = NUM_THREADS)
{
    // fmt::print("\t\tEntered wrapper\n");
    std::vector<uint64_t> valid_nodes;
    for (auto hash : nodes)
    {   
        if (GNM[hash].xorable) // [hash]
        {
            valid_nodes.push_back(hash);
        }
    }
    
    const uint64_t N = valid_nodes.size();
    if (N <= 1) return;

    const uint64_t Ncombs = N * (N - 1) / 2;
    const uint8_t tgt_lvl = depth / 3;
    // TODO : bring splitting to XOR and AND/OR
    if (Ncombs < MAX_SIZE)
    {
        uint64_t chunk_size = Ncombs / num_threads;
        if (chunk_size == 0)
        {
            num_threads = 1;
            chunk_size = Ncombs;
        }
        
        std::vector<TT> funcs(Ncombs);
        std::vector<uint32_t> costs(Ncombs);

        // fmt::print("\t\tCalculating func and cost for {} nodes. Threads: {}; Chunk size: {} \n", N, num_threads, chunk_size);
        std::vector<std::thread> threads;
        // threads.reserve(num_threads);
        // Launch threads to square elements of the vector
        for (uint64_t i = 0; i < num_threads; ++i) {
            const int start_k = i * chunk_size;
            const int end_k = (i == num_threads - 1) ? Ncombs : (i + 1) * chunk_size;
            // fmt::print("\t\tChunk #{}: Iterating k from {} to {} \n", i, start_k, end_k);
            threads.push_back(
                std::thread(threaded_xor, std::ref(valid_nodes), std::ref(funcs), std::ref(costs), start_k, end_k, COSTS[fXOR], 0, tgt_lvl)
            );
        }
        join_threads(threads);
        check_threads(threads);

        for (uint64_t k = 0u; k < Ncombs; k++)
        {
            COMPUTE_IJ(k, N, i, j);
            continue_if (costs[k] >= INF32);
            create_node(funcs[k], fXOR, costs[k], depth, true, {valid_nodes[i], valid_nodes[j]}); //auto [nhash, added] = 
        }
    }
    else //too many combinations to fit into memory
    {
        fmt::print("\t\tToo many pairs among {:L} xorable nodes to fit into memory (total {:L}). Splitting... \n", N, Ncombs);

        uint64_t NUM_MACRO_CHUNKS = ceil(Ncombs, MAX_SIZE);
        for (auto [START, END]  : divide_range(Ncombs, NUM_MACRO_CHUNKS))
        {
            uint64_t chunk_size = END - START;
            std::vector<TT> funcs(chunk_size, num2tt(0));
            std::vector<uint32_t> costs(chunk_size, INF32);
            std::vector<std::thread> threads;
            threads.reserve(num_threads);
            for (auto [start, end] : divide_range(chunk_size, NUM_THREADS))
            {
                uint64_t start_k = start + START;
                uint64_t end_k   = end   + START;
                threads.push_back(
                    std::thread(threaded_xor, std::ref(valid_nodes), std::ref(funcs), std::ref(costs), start_k, end_k, COSTS[fXOR], START, tgt_lvl)
                );
            }
            join_threads(threads);
            check_threads(threads);
            fmt::print("\t\tCombining XOR results\n");
            for (uint64_t k = START; k < END; k++)
            {
                COMPUTE_IJ(k, N, i, j);
                uint32_t cost = costs[k - START];
                TT func = funcs[k-START];
                continue_if (cost >= INF32);
                create_node(func, fXOR, cost, depth, true, {valid_nodes[i], valid_nodes[j]});
            }
        }
    }
}

void threaded_and_or(const std::vector<uint64_t>& hashes, std::vector<TT>& funcs_and, std::vector<uint32_t>& costs_and, std::vector<TT>& funcs_or, std::vector<uint32_t>& costs_or, const uint64_t start_k, const uint64_t end_k, const uint32_t gate_cost_and, const uint64_t offset, const uint8_t tgt_lvl) 
{
    // fmt::print("\t\tFunction is called\n");
    const uint64_t N = hashes.size();
    fmt::print("\t\tAND/OR: Iterating between {:L} and {:L}\n", start_k, end_k);
    uint64_t old_i;
    uint32_t init_cost;
    std::vector<std::pair<uint64_t, uint8_t>> ct_spl;
    std::vector<uint64_t> non_splittable_nodes;

    for (uint64_t k = start_k; k < end_k; k++)
    {   
        COMPUTE_IJ(k, N, i, j);
        const uint64_t hi = hashes[i]; // fmt::print("\t\t hash_i={}\n", hi);
        const Node & ni = GNM[hi]; // fmt::print("\t\t Retrieved ni={}\n", ni.to_str());
        if (i != old_i || k == start_k) // update [init_cost, ct_spl, non_splittable_nodes] if the row has changed
        {
            //             fmt::print("\tUpdating\n");
            auto [init_cost_tmp, ct_spl_tmp, non_splittable_nodes_tmp] = node_cost_single(hi, gate_cost_and, tgt_lvl*3);
            init_cost = init_cost_tmp;
            ct_spl = ct_spl_tmp;
            non_splittable_nodes = non_splittable_nodes_tmp;
        }

        const uint64_t hj = hashes[j]; // fmt::print("\t\t hash_j={}\n", hj);
        const Node & nj = GNM[hj]; // fmt::print("\t\t Retrieved nj={}\n", nj.to_str());

        const TT func_and = ni.func & nj.func;
        const uint64_t best_and_hash = GEX[func_and._bits];
        uint8_t best_and_lvl;
        uint32_t best_and_cost;
        if (!(kitty::is_const0(func_and) || !kitty::is_const0(~func_and) || func_and == ni.func || func_and == nj.func))
        {
            const Node & best_and = GNM[best_and_hash]; // [best_and_hash];
            best_and_lvl = best_and.lvl;
            best_and_cost = best_and.cost;
        }
        else
        {
            best_and_lvl = 0;
            best_and_cost = 0;
        }
        const TT func_or  = ni.func | nj.func;
        const uint64_t best_or_hash  = GEX[func_or._bits];
        uint8_t best_or_lvl;
        uint32_t best_or_cost;
        if (!(kitty::is_const0(func_or) || !kitty::is_const0(~func_or) || func_or == ni.func || func_or == nj.func))
        {
            const Node & best_or = GNM[best_or_hash];
            best_or_lvl = best_or.lvl;
            best_or_cost = best_or.cost;
        }
        else
        {
            best_or_lvl = 0;
            best_or_cost = 0;
        }

        if (best_and_lvl >= tgt_lvl && best_or_lvl >= tgt_lvl) 
        {
            uint8_t limit = std::max(best_and_cost, best_or_cost);
            auto [cost_and, status] = node_cost(hj, init_cost, ct_spl, non_splittable_nodes, limit, tgt_lvl*3);
            if (status)
            {
                uint32_t cost_or = cost_and - COSTS[fAND] + COSTS[fOR]; 
                funcs_and[k - offset] = func_and;  
                costs_and[k - offset] = cost_and;   
                funcs_or[k - offset] = func_or;  
                costs_or[k - offset] = cost_or; 
            }  
        }
        else if (best_and_lvl >= tgt_lvl)
        {
            auto [cost_and, status] = node_cost(hj, init_cost, ct_spl, non_splittable_nodes, best_and_cost, tgt_lvl*3);
            if (status)
            {
                funcs_and[k - offset] = func_and; 
                costs_and[k - offset] = cost_and; 
            }
        }
        else if (best_or_lvl >= tgt_lvl) 
        {
            auto [cost_or, status] = node_cost(hj, init_cost, ct_spl, non_splittable_nodes, best_or_cost, tgt_lvl*3); // costs[i][j] = node_cost(ni, nj, gate_cost);
            if (status)
            {
                funcs_or[k - offset] = func_or;  
                costs_or[k - offset] = cost_or;  
            }
        }
        old_i = i;
    }
}

void threaded_and_or_wrapper(const std::vector<uint64_t>& nodes, const uint8_t depth, uint8_t num_threads = NUM_THREADS)
{
    const uint64_t N = nodes.size();
    if (N <= 1) return;

    const uint64_t Ncombs = N * (N - 1) / 2;
    const uint8_t tgt_lvl = depth / 3;
    if (Ncombs < MAX_SIZE)
    {
        uint64_t chunk_size = Ncombs / num_threads;
        if (chunk_size == 0)
        {
            num_threads = 1;
            chunk_size = Ncombs;
        }
        
        std::vector<TT> funcs_and(Ncombs);
        std::vector<uint32_t> costs_and(Ncombs);
        std::vector<TT> funcs_or(Ncombs);
        std::vector<uint32_t> costs_or(Ncombs);

        std::vector<std::thread> threads;
        // threads.reserve(num_threads);
        for (uint64_t i = 0; i < num_threads; ++i) {
            const int start_k = i * chunk_size;
            const int end_k = (i == num_threads - 1) ? Ncombs : (i + 1) * chunk_size;
            // fmt::print("\t\tChunk #{}: Iterating k from {} to {} \n", i, start_k, end_k);
            threads.push_back(
                std::thread(threaded_and_or, std::ref(nodes), std::ref(funcs_and), std::ref(costs_and), std::ref(funcs_or), std::ref(costs_or), start_k, end_k, COSTS[fAND], 0, tgt_lvl)
            );
        }
        join_threads(threads);
        check_threads(threads);

        for (uint64_t k = 0u; k < Ncombs; k++)
        {
            COMPUTE_IJ(k, N, i, j);
            if (costs_and[k] <= INF32)
            {
                create_node(funcs_and[k], fAND, costs_and[k], depth, true, {nodes[i], nodes[j]});// auto [nhash_and, added_and]
            }
            if (costs_or[k] <= INF32)
            {
                create_node(funcs_or[k],   fOR, costs_or[k] , depth, true, {nodes[i], nodes[j]});// auto [nhash_or , added_or ]
            }
        }
    }
    else
    {
        fmt::print("\t\tToo many pairs among {:L} AS nodes to fit into memory (total {:L}). Splitting... \n", N, Ncombs);
        uint64_t NUM_MACRO_CHUNKS = ceil(Ncombs, MAX_SIZE);
        for (auto [START, END]  : divide_range(Ncombs, NUM_MACRO_CHUNKS))
        {
            uint64_t chunk_size = END - START;
            std::vector<TT> funcs_and(chunk_size, num2tt(0));
            std::vector<uint32_t> costs_and(chunk_size, INF32);
            std::vector<TT> funcs_or(chunk_size, num2tt(0));
            std::vector<uint32_t> costs_or(chunk_size, INF32);
            std::vector<std::thread> threads;
            threads.reserve(num_threads);

            for (auto [start, end] : divide_range(chunk_size, NUM_THREADS))
            {
                uint64_t start_k = start + START;
                uint64_t end_k   = end   + START;
                threads.push_back(
                    std::thread(threaded_and_or, std::ref(nodes), std::ref(funcs_and), std::ref(costs_and), std::ref(funcs_or), std::ref(costs_or), start_k, end_k, COSTS[fAND], START, tgt_lvl)
                );
            }
            join_threads(threads);
            check_threads(threads);

            fmt::print("\t\tCombining AND/OR results\n");
            for (uint64_t k = START; k < END; k++)
            {
                COMPUTE_IJ(k, N, i, j);
                if (costs_and[k-START] <= INF32)
                {
                    create_node(funcs_and[k-START], fAND, costs_and[k-START], depth, true, {nodes[i], nodes[j]});
                }
                if (costs_or[k-START] <= INF32)
                {
                    create_node(funcs_or[k-START],  fOR , costs_or[k-START] , depth, true, {nodes[i], nodes[j]});
                }
            }
        }
    }
}

std::unordered_set<uint64_t> remove_dominated()
{
    std::unordered_set<uint64_t> protected_nodes;
    std::vector<uint64_t> stack;
    for (const auto & [hash, node] : GNM)
    {
        continue_if (node.last_func == fDFF);
        const Node& best = node.xorable? get_gex(node.func) : get_gea(node.func); //GEX[node.func] : GEA[node.func];
        if (node.lvl < best.lvl || (node.lvl == best.lvl && node.cost <= best.cost) )
        {
            // fmt::print("\t\t\t\tProtecting node {} (best is {}) [{} {}]\n", node.to_str(), best.to_str(), node.hash, best.hash);
            stack.push_back(hash);
        }
    } 

    while (!stack.empty())
    {
        const uint64_t hash = stack.back();
        stack.pop_back();
        
        const Node& node = GNM[hash]; // [hash];
        if (protected_nodes.find(hash) == protected_nodes.end())
        {
            for (const uint64_t p_hash : node.parent_hashes)
            {
                stack.push_back(p_hash);
            }
        }
        protected_nodes.insert(hash);
    }
    // fmt::print("\t\t\t\tProtecting {} nodes out of {}\n", protected_nodes.size(), GNM.size());

    std::unordered_set<uint64_t> to_be_removed;

    for (const auto & [hash, node] : GNM)
    {
        continue_if (node.last_func == fDFF);
        continue_if (protected_nodes.find(hash) != protected_nodes.end());
        // const Node& best = node.xorable? GEX[node.func] : GEA[node.func];
        const Node& best = node.xorable? get_gex(node.func) : get_gea(node.func); //GEX[node.func] : GEA[node.func];
        if (node.lvl > best.lvl || (node.lvl == best.lvl && node.cost > best.cost) )
        {
            to_be_removed.insert(hash);
        }
    }

    for (const auto hash : to_be_removed)
    {
        // fmt::print("\t\t\t\tDestroying dominated node {}\n", GNM[hash].to_str());
        GNM.erase(hash);
    }   
    return to_be_removed;
}

void cb_generation(uint8_t lvl)
{
    const uint8_t min_lvl = (lvl == 0) ? 0 : lvl * 3 - 1;
    std::vector<uint64_t> new_nodes = select_depth(min_lvl, lvl * 3 + 1);
    fmt::print("\tChecking new nodes after selection between depths {:L} and {:L}...\n", lvl * 3, lvl * 3 + 1);
    // check_nodes(new_nodes);

    uint8_t depth = lvl * 3 + 1;
    std::vector<uint64_t> old_nodes;
    bool go_on;
    std::vector<uint64_t> fresh_nodes;
    uint8_t iteration = 0;
    do {
        // first, combine new nodes with old nodes
        fmt::print("\tProcessing {:L} old nodes with {:L} new nodes\n", old_nodes.size(), new_nodes.size());
        thread_old_new_wrapper(old_nodes, new_nodes, depth, fresh_nodes);
        fmt::print("\tProcessing pair among {:L} new nodes\n", old_nodes.size(), new_nodes.size());
        threaded_new_new_wrapper(new_nodes, depth, fresh_nodes);
        std::unordered_set<uint64_t> removed_nodes = remove_dominated();
        for (const uint64_t new_hash : new_nodes)
        {
            if (removed_nodes.find(new_hash) == removed_nodes.end())
            {
                old_nodes.push_back(new_hash);
            }
        }
        new_nodes.clear();
        for (const uint64_t fresh_hash : fresh_nodes)
        {
            if (removed_nodes.find(fresh_hash) == removed_nodes.end())
            {
                new_nodes.push_back(fresh_hash);
            }
        }
        go_on = fresh_nodes.size() > 0;
        fresh_nodes.clear(); 
        iteration++;
    } while(go_on);
}

void as_generation(uint8_t lvl)
{
    uint8_t min_lvl = (lvl == 0) ? 0 : lvl * 3 - 1;
    std::vector<uint64_t> nodes = select_depth(min_lvl, lvl * 3 + 1);

    uint8_t tgt_depth = lvl * 3 + 2;
    fmt::print("\t{} : DFF/NOT {:L} nodes\n", lvl, nodes.size());
    for (const uint64_t hash : nodes)
    {
        const Node& ni =  GNM[hash];
        force_create_node(ni.func, fDFF, ni.cost + COSTS[fDFF], tgt_depth, true, {hash});
        create_node(~ni.func, fNOT, ni.cost + COSTS[fNOT], tgt_depth, true, {hash});
    }
    // xor-combine the nodes
    fmt::print("\t{}: XOR {:L} nodes\n", lvl, nodes.size());
    threaded_xor_wrapper(nodes, tgt_depth);
    remove_dominated();
}

void sa_generation(uint8_t lvl)
{
    std::vector<uint64_t> nodes = select_depth(lvl * 3 + 2, lvl * 3 + 2);
    uint8_t tgt_depth = lvl * 3 + 3;

    // and/or-combine the nodes
    fmt::print("\t{}: AND/OR {:L} nodes\n", lvl, nodes.size());
    threaded_and_or_wrapper(nodes, tgt_depth);
    remove_dominated();
}

inline bool is_done(std::array<uint64_t, NUM_TT> & GEX)  
{
    // TODO: modify to keep track of best levels considering other sets of levels
    #pragma vector
    for (auto i = 0u; i < NUM_TT; i++)
    {
        if (GEX[i] == 0)
        {
            return false;
        }
    }
    return true;
}

inline uint32_t count_done()  
{
    #pragma vector
    uint32_t ctr = 0;
    for (auto i = 0u; i < NUM_TT; i++)
    {
        ctr += (GNM[GEX[i]].cost < INF32);
    }
    return ctr;
}

bool is_subset(const std::vector<uint64_t>& v1, const std::vector<uint64_t>& v2) {
    // Check if all elements of v1 are in v2
    return std::all_of(v1.begin(), v1.end(), [&v2](const int& val) {
        return std::find(v2.begin(), v2.end(), val) != v2.end();
    });
}

std::vector<uint64_t> get_pi_recursive(uint64_t hash, std::unordered_map<uint64_t, Node>& nodemap, std::unordered_set<uint64_t>& all_hashes)
{
    auto it = NODE_PIS.find(hash);
    if (it == NODE_PIS.end())
    {
        Node & n = nodemap[hash]; 
        if (n.last_func == fPI) // the node is itself PI, insert hash. No continuation is needed
        {
            std::vector<uint64_t> pis = { hash };
            NODE_PIS[hash] = pis;
            PI_HASHES.emplace(hash);
            all_hashes.erase(hash);
            return pis;
        }
        else if (n.last_func == fNOFUNC)
        {
            // Don't throw, since this may interrupt large scale optimization
            fmt::print("\tWARNING: encountered invalid node.");
            all_hashes.erase(hash);
            return {};
        }
        else // in all other cases, just call the function once again
        {
            std::vector<uint64_t> pis;
            for (uint64_t phash : n.parent_hashes)
            {
                pis.emplace_back(get_pi_recursive(phash, nodemap, all_hashes));
            }
            NODE_PIS[hash] = pis;
            all_hashes.erase(hash);
            return pis;
        }
    }
    all_hashes.erase(hash);
    return it->second;
}

std::unordered_set<uint64_t> get_hashes(const std::unordered_map<uint64_t, Node>& hashmap) 
{
    std::unordered_set<uint64_t> hashes;
    for (const auto& [hash, value] : hashmap) 
    {
        Node & n = GNM_global[hash];
        if (n.last_func != fNOFUNC && NODE_PIS.find(hash) == NODE_PIS.end())
        {
            hashes.emplace(hash);
        }
    }
    return hashes;
}

void update_NODE_PIS()
{
    std::unordered_set<uint64_t> all_hashes = get_hashes(GNM_global);
    while (!all_hashes.empty())
    {
        get_pi_recursive(*all_hashes.begin(), GNM_global, all_hashes);
    }
}

void get_nodemap(std::array<uint8_t, NUM_VARS>& levels, std::unordered_map<uint64_t, Node>& hashmap, std::array<uint64_t, NUM_TT>& local_GEA, std::array<uint64_t, NUM_TT>& local_GEX) 
// TODO: also UPDATE GEX and GEA
{
    update_NODE_PIS();
    std::vector<uint64_t> valid_pi_hashes;
    for (const uint64_t pi_hash : PI_HASHES)
    {
        const Node & pi = GNM_global[pi_hash];
        uint8_t pi_idx = INF8;
        switch (pi.func._bits)
        {
            case 0x5555: pi_idx = 0;
            case 0x3333: pi_idx = 1;
            case 0x0F0F: pi_idx = 2;
            case 0x00FF: pi_idx = 3;
        }
        if (levels[pi_idx] == pi.lvl)
        {
            valid_pi_hashes.push_back(pi.hash);
        }
    }
    
    for (auto & [hash, pis] : NODE_PIS)
    {
        if ( std::includes(valid_pi_hashes.begin(), valid_pi_hashes.end(), pis.begin(), pis.end()) )
        {
            hashmap[hash] = GNM_global[hash];
        }
        else 
        {
            hashmap.erase(hash);
        }
    }
}

int main() 
{
    std::locale::global(std::locale("en_US.UTF-8"));

    //outer dimension - function
    //middle dimension - set of delays
    //inner dimension - delay of each input
    // std::vector<std::vector<std::vector<uint8_t>>> best_delays;

    for (std::array<uint8_t, 4> levels : sets_of_levels) 
    {
        std::fill(GEA.begin(), GEA.end(), 0);
        std::fill(GEX.begin(), GEX.end(), 0);
        // IMPORTANT: update GNM_global, GEA_global, GEX_global after a set of levels finished processing
        // IMPORTANT: create local GNM, GEX, and GEA at each iteration
        // IMPORTANT: calculate delay for each node in GNM and mark the node hopeless if the node has strictily worse delay than the existing node in GNM_global
        // IMPORTANT: consider hopeless nodes completed to terminate the iteration early
        // IMPORTANT: make sure only four PIs (+2 constants) reside in GNM at each iteration

        // IMPORTANT: also create functions to shift the existing nodes one level (depth += 3) forward. 
        // IMPORTANT: add these nodes to GNM_global and don't forget to update GEA/GEX
        uint8_t i = 0u;
        for (uint64_t func : gen_pi(NUM_VARS))
        {
            create_node(num2tt(func), fPI, 0, levels[i++]*3 + 1, true);
        }

        std::string level_prefix = fmt::format("exploration_test/x3_{}", fmt::join(levels, ""));

        fmt::print("Analyzing levels: {}\n", level_prefix);
        fmt::print("After PI:\n");
        print_GNM();
        create_node(const0, fPI, 0, 0, true);
        create_node(const1, fPI, 0, 0, true);
        // fmt::print("After const:\n");
        // print_GNM();
        
        // for (uint8_t lvl = 0; lvl < 50; lvl ++)
        for (uint8_t lvl = 0; lvl < 50; lvl ++)
        {
            uint64_t start_n_TT = count_done();
            fmt::print("Processing lvl {}: CB\n", lvl);
            // fmt::print("Checking integrity of GNM before CB\n");
            // check_GNM();
            cb_generation(lvl);
            // fmt::print("Checking integrity of GNM after CB\n");
            // check_GNM();
            fmt::print("\n\n\t\nCompleted {}\n\n\n", count_done());
            // write_csv_gnm(GNM, fmt::format("{}_gnm_cb_{}_{}.csv", level_prefix, lvl, get_current_time_formatted()));
            // write_csv_arr(GEA, fmt::format("{}_gea_cb_{}_{}.csv", level_prefix, lvl, get_current_time_formatted()));
            // write_csv_arr(GEX, fmt::format("{}_gex_cb_{}_{}.csv", level_prefix, lvl, get_current_time_formatted()));
            break_if (is_done(GEX) && lvl >= levels.back()) ;

            // break_if (lvl == 1); 

            fmt::print("Processing lvl {}: AS\n", lvl);
            // fmt::print("Checking integrity of GNM before AS\n");
            // check_GNM();
            as_generation(lvl);
            // fmt::print("Checking integrity of GNM after AS\n");
            // check_GNM();
            fmt::print("\n\n\t\nCompleted {}\n\n\n", count_done());
            // write_csv_gnm(GNM, fmt::format("{}_gnm_as_{}_{}.csv", level_prefix, lvl, get_current_time_formatted()));
            // write_csv_arr(GEA, fmt::format("{}_gea_as_{}_{}.csv", level_prefix, lvl, get_current_time_formatted()));
            // write_csv_arr(GEX, fmt::format("{}_gex_as_{}_{}.csv", level_prefix, lvl, get_current_time_formatted()));
            break_if (is_done(GEX) && lvl >= levels.back()) ;

            fmt::print("Processing lvl {}: SA\n", lvl);
            // fmt::print("Checking integrity of GNM before SA\n");
            // check_GNM();
            sa_generation(lvl);
            // fmt::print("Checking integrity of GNM after SA\n");
            // check_GNM();
            fmt::print("\n\n\t\nCompleted {}\n\n\n", count_done());
            // write_csv_gnm(GNM, fmt::format("{}_gnm_sa_{}_{}.csv", level_prefix, lvl, get_current_time_formatted()));
            // write_csv_arr(GEA, fmt::format("{}_gea_sa_{}_{}.csv", level_prefix, lvl, get_current_time_formatted()));
            // write_csv_arr(GEX, fmt::format("{}_gex_sa_{}_{}.csv", level_prefix, lvl, get_current_time_formatted()));
            break_if (is_done(GEX) && lvl >= levels.back()) ;

            uint64_t end_n_TT = count_done();
            break_if (start_n_TT == end_n_TT && lvl >= levels.back());
        }

        // Write data to CSV files
        write_csv_gnm(GNM, fmt::format("{}_gnm.csv", level_prefix));
        write_csv_arr(GEA, fmt::format("{}_gea.csv", level_prefix));
        write_csv_arr(GEX, fmt::format("{}_gex.csv", level_prefix));

        std::unordered_map<uint64_t, Node> GNM_new;
        std::array<uint64_t, NUM_TT> GEA_new;
        std::array<uint64_t, NUM_TT> GEX_new;

        // Read data from CSV files
        GNM_new = read_csv_gnm(fmt::format("{}_gnm.csv", level_prefix));
        GEA_new = read_csv_arr(fmt::format("{}_gea.csv", level_prefix));
        GEX_new = read_csv_arr(fmt::format("{}_gex.csv", level_prefix));

        std::cout << GNM.size() << " " << GNM_new.size() << " " << (GNM == GNM_new) << std::endl;
        std::cout << GEA.size() << " " << GEA_new.size() << " " << (GEA == GEA_new) << std::endl;
        std::cout << GEX.size() << " " << GEX_new.size() << " " << (GEX == GEX_new) << std::endl;

        for (const auto & [hash, n1] : GNM)
        {
            const Node& n2 = GNM_new[hash];
            bool scalar_fields_equal = (n1.func == n2.func && n1.last_func == n2.last_func && n1.cost == n2.cost && n1.depth == n2.depth && n1.xorable == n2.xorable);
            std::vector<uint64_t> p1, p2;
            p1.insert(p1.end(), n1.parent_hashes.begin(), n1.parent_hashes.end());
            p2.insert(p2.end(), n2.parent_hashes.begin(), n2.parent_hashes.end());

            std::sort(p1.begin(), p1.end());
            std::sort(p2.begin(), p2.end());
            bool parents_equal = (p1 == p2);

            if (!(scalar_fields_equal && parents_equal))
            {
                fmt::print("Hash: {}\n\n", hash);
                fmt::print("\tfunc\t: {} - {}\n", n1.func._bits, n2.func._bits);
                fmt::print("\tlast_func: {} - {}\n", n1.func._bits, n2.func._bits);
                fmt::print("\tcost\t: {} - {}\n", n1.cost, n2.cost);
                fmt::print("\tdepth\t: {} - {}\n", n1.depth, n2.depth);
                fmt::print("\txorable\t: {} - {}\n", n1.xorable, n2.xorable);
                for (const uint64_t phash : n1.parent_hashes)
                {
                    fmt::print("\t\tParents 1\t: {}\n", phash);
                }
                for (const uint64_t phash : n2.parent_hashes)
                {
                    fmt::print("\t\tParents 2\t: {}\n", phash);
                }
            }
        }

        GNM_global.insert(GNM.begin(), GNM.end());
        GEA_global[levels] = GEA;
        GEX_global[levels] = GEX;
    }
    return 0;
}

/*
nohup ./experiments/explore_sfq_vec > 20230329_vec/x3.log 2>&1 &
*/
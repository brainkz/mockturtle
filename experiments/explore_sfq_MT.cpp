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
#include "mockturtle/algorithms/nodes.hpp"
#include <algorithm>
// #include <execution>

constexpr US THREADING_MODE = 0; //0 threading, 1 batch processing, 2 serial
constexpr UL MAX_SIZE = 200'000'000;
constexpr UI NUM_THREADS = 100;

inline void join_threads(std::vector<std::thread>& threads, UL maxthreads = INF)
{
    auto i = 0u;
    for (auto& t : threads) {
        if (i == maxthreads)
        {
            break;
        }
        if (t.joinable()) {
            t.join();
        }
    }
}

std::unordered_map<ULL, Node> GNM;
// std::unordered_map<ULL, std::unordered_map<ULL,UI>> GN_CT;
// std::unordered_map<ULL, std::unordered_set<ULL>> GN_CT;
// std::array<Node, NUM_TT> GEA;
// std::array<Node, NUM_TT> GEX;
std::array<ULL, NUM_TT> GEA;
std::array<ULL, NUM_TT> GEX;

template<typename dtype_x, typename dtype_y>
inline dtype_y ceil(dtype_x x, dtype_y y)
{
    return 1 + ((x - 1) / y);
}


inline Node& get_gea(UI func) 
{
    ULL hash = GEA[func];
    return GNM[hash];
}
inline Node& get_gex(UI func)
{
    ULL hash = GEX[func];
    return GNM[hash];
}

inline void check_node(ULL hash)
{
    Node & n = GNM[hash];
    fmt::print("Checking node {} | num_parents = {} | {:L}\n", n.to_str(), n.parent_hashes.size(), hash);
    if (n.last_func == fNOFUNC)
    {
        assert(n.parent_hashes.empty());
        return;
    }
    else if (n.func == 0 || n.func == ONES)
    {
        assert(n.last_func == fPI);
    }
    assert(n.cost != INF);
    assert(n.depth != INF);
}

inline void check_nodes(const std::vector<ULL>& nodes)
{
    if (nodes.empty()) 
    {
        return;
    }
    for (ULL hash : nodes)
    {
        check_node(hash);
    }
}

inline void check_GNM()
{
    UI violation_ct = 0u;
    for (auto [hash, n] : GNM)
    {
        fmt::print("Checking {} \n", n.to_str());
        if (n.last_func == fNOFUNC)
        {
            violation_ct++;
        }
        else
        {
            assert(n.func != 0 || n.last_func == fPI);
            assert(n.cost != INF);
            assert(n.depth != INF);
        }
        // if (violation_ct > 1)
        // {
        //     assert(false);
        // }
    }
}

inline void print_GNM()
{
    for (auto [hash, n] : GNM)
    {
        fmt::print("{}\n", n.to_str());
    }
}

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

std::tuple<ULL, bool> create_node(TT _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes, ULL _hash)
//, std::unordered_map<ULL, Node>& hash_map = GNM)
{
    auto it = GNM.find(_hash);
    bool node_is_new = (it == GNM.end());
    if (node_is_new) //the node is new
    {
        // fmt::print("\t\t\t\tChecking for creation {}\n", Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash).to_str());
        UI lvl = _depth / 3;
        Node& best_any = get_gea(_func);
        // fmt::print("\t\t\t\tConsidering {}: {}, {}\n", _func,          lvl,         _cost);
        // fmt::print("\t\t\t\tBest A      {}: {}, {}\n", _func, best_any.lvl, best_any.cost);

        if (_xorable)
        {
            // Node& best_xor = GEX[_func];
            Node& best_xor = get_gex(_func);    
            // fmt::print("\t\t\t\tBest X      {}: {}, {}\n", _func, best_xor.lvl, best_xor.cost);
            if (std::tie(lvl, _cost ) <= std::tie(best_any.lvl, best_any.cost))
            {
                GNM[_hash] = Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
                GEA[_func] = _hash; //hash_map[_hash];
                GEX[_func] = _hash; //hash_map[_hash];
                // fmt::print("\t\t\t\tNode is better than any existing function\n");
                return std::make_tuple(_hash, true);
            }
            else if (std::tie(lvl, _cost ) <= std::tie(best_xor.lvl, best_xor.cost))
            {
                GNM[_hash] = Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
                GEX[_func] = _hash; //hash_map[_hash];
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
                GEA[_func] = _hash; //hash_map[_hash];
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


std::tuple<ULL, bool> create_node(TT _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes = {})
{
    ULL _hash = calculate_hash(_func, _last_func, _cost, _depth, _xorable, _parent_hashes);
    return create_node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash);
}

std::tuple<ULL, bool> create_node_no_upd(TT _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes, ULL _hash, std::unordered_map<ULL, Node>& hash_map = GNM)
{
    auto it = hash_map.find(_hash);
    bool node_is_new = (it == hash_map.end());
    if (node_is_new) //the node is new
    {
        // fmt::print("\t\t\t\tChecking for creation {}\n", Node(_func, _last_func, _cost, _depth, _xorable, _parent_hashes, _hash).to_str());
        UI lvl = _depth / 3;
        Node& best_any = get_gea(_func); // GEA[_func];
        // fmt::print("\t\t\t\tConsidering {}: {}, {}\n", _func,          lvl,         _cost);
        // fmt::print("\t\t\t\tBest A      {}: {}, {}\n", _func, best_any.lvl, best_any.cost);

        if (_xorable)
        {
            // fmt::print("\t\t\t\tFunction is xorable\n");
            Node& best_xor = get_gex(_func); // GEX[_func];
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
            Node& best_any = get_gea(node.func); //GEA[node.func];

            if (node.xorable)
            {
                Node& best_xor = get_gex(node.func);// GEX[node.func];
                if (std::tie(node.lvl, node.cost) <= std::tie(best_any.lvl, best_any.cost))
                {
                    GNM.insert(std::move(*it_tmp)); hash_map.erase(it_tmp); //move the node
                    GEA[node.func] = hash; //GNM[hash];
                    GEX[node.func] = hash; //GNM[hash];
                    // fmt::print("\t\t\t\tNode is better than any existing function\n");
                }
                else if (std::tie(node.lvl, node.cost) <= std::tie(best_xor.lvl, best_xor.cost))
                {
                    GNM.insert(std::move(*it_tmp)); hash_map.erase(it_tmp); //move the node
                    GEX[node.func] = hash; //GNM[hash];
                    // fmt::print("\t\t\t\tNode is better than any existing xorable function\n");
                }
            }
            else 
            {
                // fmt::print("\t\t\t\tFunction is NOT xorable\n");
                if (std::tie(node.lvl, node.cost) <= std::tie(best_any.lvl, best_any.cost))
                {
                    GNM.insert(std::move(*it_tmp)); hash_map.erase(it_tmp); //move the node
                    GEA[node.func] = hash; //GNM[hash];
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
                // Node& p = GNM[p_hash];
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

std::tuple<UI, std::unordered_map<ULL, UI>, std::unordered_set<ULL>> node_cost_single(const ULL h1, UI gate_cost)
{
    std::stack<ULL> stack;
    stack.push(h1);

    std::unordered_map<ULL, UI> ct_spl;
    std::unordered_set<ULL> non_splittable_nodes;

    // fmt::print("\t\tPrecalculating {}\n", GNM[h1].to_str());
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
                // Node& p = GNM[p_hash];
                stack.push(p_hash);
                // fmt::print("\t\t\tAdding parent node {} to stack\n", p.to_str());
            }
        }
        else
        {
            gate_cost += COSTS[fSPL];
            // fmt::print("\t\t\tNot first time, adding the cost of SPL ({}). Total cost is {}\n", COSTS[fSPL], gate_cost);
        }
        if (n.last_func == fAND || n.last_func == fOR)
        {
            for (const auto& p_hash : n.parent_hashes)
            {
                non_splittable_nodes.emplace(p_hash);
            }
        }
    }
    // fmt::print("\t\tFinished precalculating for {} : cost = {}\n", GNM[h1].to_str(), gate_cost);
    return std::make_tuple(gate_cost, ct_spl, non_splittable_nodes);
}

// assumes precomputed gate_cost, ct_spl, and non_splittable_nodes
std::tuple<UI, bool> node_cost(ULL h2, UI init_cost, std::unordered_map<ULL, UI> ct_spl, std::unordered_set<ULL> non_splittable_nodes, UI limit)
{
    std::stack<ULL> stack;
    stack.push(h2);
    // fmt::print("\t\t\tProcessing node {}, limit = {}\n", GNM[h2].to_str(), limit);
    if (init_cost > limit) return std::make_tuple(INF, false);
    while (!stack.empty())
    {
        ULL n_hash = stack.top();
        stack.pop();
        Node& n = GNM[n_hash];
        // fmt::print("\t\t\tProcessing node {}\n", n.to_str());
        ct_spl[n_hash]++;
        if (ct_spl[n_hash] == 1)
        {
            init_cost += COSTS[n.last_func];
            if (init_cost > limit) return std::make_tuple(INF, false);      
            // fmt::print("\t\t\tFirst time, adding the cost of {} ({}). Total cost is {}\n", F2STR[n.last_func], COSTS[n.last_func], init_cost);
            for (const auto& p_hash : n.parent_hashes)
            {
                // Node& p = GNM[p_hash];
                stack.push(p_hash);
                // fmt::print("\t\t\tAdding parent node {} to stack\n", p.to_str());
            }
        }
        else
        {
            init_cost += COSTS[fSPL];
            if (init_cost > limit) return std::make_tuple(INF, false);
            // fmt::print("\t\t\tNot first time, adding the cost of SPL ({}). Total cost is {}\n", COSTS[fSPL], init_cost);
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
            init_cost += COSTS[last_func] * (count - 1); // duplicate a gate
            if (init_cost > limit) return std::make_tuple(INF, false);
            if (last_func == fXOR)
            {
                init_cost += (count - 1) * COSTS[fSPL]; // need to do twice more splittings for duplicating an XOR gate
                if (init_cost > limit) return std::make_tuple(INF, false);
            }
        }
    }
    // fmt::print("final cost is {} and status is {}\n", init_cost, true);
    return std::make_tuple(init_cost, true);
}

void thread_old_new(const std::vector<ULL>& hashes_i, const std::vector<ULL>& hashes_j, std::vector<std::vector<UI>>& funcs, std::vector<std::vector<bool>>& xorables, std::vector<std::vector<UI>>& costs, const ULL start_row, const ULL end_row, const UI gate_cost, const ULL offset, const UI tgt_lvl) 
{
    // const int M = hashes_A.size();
    const int N = hashes_j.size();
    fmt::print("\t\tProcessing old - new pairs between #{:L} and #{:L} (out of total {:L})\n", start_row, end_row, hashes_i.size());
    for (ULL i = start_row; i < end_row; ++i) {
        ULL hi = hashes_i[i];
        Node & ni = GNM[hi];
        auto [init_cost, ct_spl, non_splittable_nodes] = node_cost_single(hi, gate_cost);
        // fmt::print("Preliminary calculations: init_cost = {}, #predecessors = {}, #non-spllittables = {}\n", init_cost, ct_spl.size(), non_splittable_nodes.size());
        for (int j = 0; j < N; ++j) 
        {
            ULL hj = hashes_j[j];
            Node & nj = GNM[hj];

            // fmt::print("\tCost calculation: ni = {}, nj = {}\n", ni.to_str(), nj.to_str());

            TT func = ni.func | nj.func;
            bool xorable = (ni.func ^ nj.func) == func;
            Node & best_n = xorable ? get_gex(func) : get_gea(func);
            if (best_n.lvl >= tgt_lvl)
            {
                auto [cost, status] = node_cost(hj, init_cost, ct_spl, non_splittable_nodes, best_n.cost); // costs[i][j] = node_cost(ni, nj, gate_cost);
                // fmt::print("\tFinal calculation: cost = {}, status = {}\n", cost, status);
                funcs[i-offset][j] = func;
                xorables[i-offset][j] = xorable;
                costs[i-offset][j] = cost;
            }
        }
    }
}

void thread_old_new_wrapper(const std::vector<ULL>& old_nodes, const std::vector<ULL>& new_nodes, const UI depth, std::vector<ULL>& fresh_nodes, UI num_threads = NUM_THREADS)
{
    const ULL M = old_nodes.size();
    if (M == 0)
    {
        return;
    }
    const ULL N = new_nodes.size();
    const UI tgt_lvl = depth / 3;
    ULL chunk_size = M / num_threads;
    if (chunk_size == 0)
    {
        num_threads = 1;
        chunk_size = M;
    }

    if (M * N < MAX_SIZE)
    {

        std::vector<std::vector<UI>> funcs(M, std::vector<UI>(N, 0));
        std::vector<std::vector<bool>> xorables(M, std::vector<bool>(N, false));
        std::vector<std::vector<UI>> costs(M, std::vector<UI>(N, INF));

        // fmt::print("\t\tCalculating func, xorable, and cost for {} old nodes and {} new nodes. Threads: {}; Chunk size: {} \n", M, N, num_threads, chunk_size);
        std::vector<std::thread> threads;
        threads.reserve(num_threads);
        for (UI i = 0; i < num_threads; i++) {
            const ULL start_row = i * chunk_size;
            const ULL end_row = (i == num_threads - 1) ? M : (i + 1) * chunk_size;
            threads.push_back(
                std::thread(thread_old_new, std::ref(old_nodes), std::ref(new_nodes), std::ref(funcs), std::ref(xorables), std::ref(costs), start_row, end_row, COSTS[fCB], 0, tgt_lvl)
            );
        }
        join_threads(threads);
        for (auto i = 0u; i < M; i++)
        {
            // fmt::print("Combining {} out of {} ({:f}\%) \n", (i+1), old_nodes.size(), 100 * (float)(i+1) / (float)old_nodes.size());
            ULL hi = old_nodes[i];
            for (auto j = 0u; j < N; j++)
            {
                ULL hj = new_nodes[j];
                UI cost = costs[i][j];
                if (cost >= INF) continue;
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
        // case 1 - row size is not too large - should in principle always be the case for 4-input TT-s
        // In this case, divide the combinations in one dimension - split by rows
        fmt::print("\t\tToo many combinations between {:L} old nodes and {:L} new nodes to fit into memory (total {:L}). Splitting... \n", M, N, M*N);
        if (N < MAX_SIZE / NUM_THREADS)
        {
            ULL nrows_per_thread = MAX_SIZE / N / NUM_THREADS;
            // fmt::print("\t\tProcessing {} rows per thread with {} threads... \n", nrows_per_thread, NUM_THREADS);
            ULL Mlocal = nrows_per_thread * NUM_THREADS;

            // std::vector<UI> end_indices;
            for (auto start = 0u; start < M; start += Mlocal)
            {
                auto end = std::min(start + Mlocal, M);
                std::vector<std::vector<UI>> funcs(Mlocal, std::vector<UI>(N, 0));
                std::vector<std::vector<bool>> xorables(Mlocal, std::vector<bool>(N, false));
                std::vector<std::vector<UI>> costs(Mlocal, std::vector<UI>(N, INF));

                std::vector<std::thread> threads;
                threads.reserve(num_threads);
                for (auto i = 0u, start_row = start; start_row < end; i++, start_row += nrows_per_thread)
                {
                    auto end_row = std::min(start_row + nrows_per_thread, end);
                    threads.push_back(
                        std::thread(thread_old_new, std::ref(old_nodes), std::ref(new_nodes), std::ref(funcs), std::ref(xorables), std::ref(costs), start_row, end_row, COSTS[fCB], start, tgt_lvl)
                    );
                }
                join_threads(threads);
                for (auto row = start; row < end; row++)
                {
                    // fmt::print("Combining {} out of {} ({:f}\%) \n", (i+1), old_nodes.size(), 100 * (float)(i+1) / (float)old_nodes.size());
                    auto i = row - start;
                    ULL hi = old_nodes[row];
                    // Node& ni =  GNM[hi];
                    for (auto j = 0u; j < N; j++)
                    {
                        ULL hj = new_nodes[j];
                        UI cost = costs[i][j];
                        if (cost >= INF) continue;
                        // Node& nj =  GNM[hj];
                        auto [nhash, added] = create_node(funcs[i][j], fCB, costs[i][j], depth, xorables[i][j], {hi, hj});
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

void threaded_new_new(const std::vector<ULL>& hashes, std::vector<UI>& funcs, std::vector<bool>& xorables, std::vector<UI>& costs, const ULL start_k, const ULL end_k, const UI gate_cost, const ULL offset, const UI tgt_lvl) 
{
    // fmt::print("\t\tFunction is called\n");
    const ULL N = hashes.size();
    fmt::print("\t\tCB: Iterating between {:L} and {:L}\n", start_k, end_k);
    ULL old_i;
    UI init_cost;
    std::unordered_map<ULL, UI> ct_spl;
    std::unordered_set<ULL> non_splittable_nodes;
    // fmt::print("\t\tProcessing new - new combinations between #{} and #{}\n", start_k, end_k);
    for (ULL k = start_k; k < end_k; k++)
    {   
        ULL i = N - 2 - ULL(sqrt(4*N*(N - 1) - 8*k - 7)/2.0 - 0.5);
        ULL hi = hashes[i]; // fmt::print("\t\t hash_i={}\n", hi);
        Node & ni = GNM[hi]; // fmt::print("\t\t Retrieved ni={}\n", ni.to_str());
        if (i != old_i || k == start_k) // update [init_cost, ct_spl, non_splittable_nodes] if the row has changed
        {
            // fmt::print("\tUpdating\n");
            auto [init_cost_tmp, ct_spl_tmp, non_splittable_nodes_tmp] = node_cost_single(hi, gate_cost);
            init_cost = init_cost_tmp;
            ct_spl = ct_spl_tmp;
            non_splittable_nodes = non_splittable_nodes_tmp;
        }

        ULL j = k + i + 1 - N * (N - 1) / 2 + (N - i)*((N - i) - 1) / 2; // fmt::print("\t\tk={}, i={}, j={}\n", k, i, j);
        ULL hj = hashes[j]; // fmt::print("\t\t hash_j={}\n", hj);
        Node & nj = GNM[hj]; // fmt::print("\t\t Retrieved nj={}\n", nj.to_str());

        // fmt::print("\t\t\tCalculating cost for k {}:\n\t\t\t\t#{}: {}\n\t\t\t\t#{}: {}\n", k, i, ni.to_str(), j, nj.to_str());

        TT func = ni.func | nj.func;
        bool xorable = (ni.func ^ nj.func) == func; // fmt::print("\t\t Xorable status = {}\n", xorables[k]);
        Node & best_n = xorable ? get_gex(func) : get_gea(func);
        if (best_n.lvl >= tgt_lvl)
        {
            auto [cost, status] = node_cost(hj, init_cost, ct_spl, non_splittable_nodes, best_n.cost); // costs[i][j] = node_cost(ni, nj, gate_cost);
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


void threaded_new_new_wrapper(const std::vector<ULL>& new_nodes, const UI depth, std::vector<ULL>& fresh_nodes, UI num_threads = NUM_THREADS)
{
    const ULL N = new_nodes.size();
    // check_nodes(new_nodes);

    if (N <= 1)
    {
        return;
    }
    const ULL Ncombs = N * (N - 1) / 2;
    const UI tgt_lvl = depth / 3;

    if (Ncombs < MAX_SIZE)
    {
        ULL chunk_size = Ncombs / (num_threads - 1);
        if (chunk_size == 0)
        {
            num_threads = 1;
            chunk_size = Ncombs;
        }
        
        std::vector<UI> funcs(Ncombs, 0);
        std::vector<bool> xorables(Ncombs, false);
        std::vector<UI> costs(Ncombs, INF);

        // fmt::print("\t\tCalculating func, xorable, and cost for {} new nodes. Pairs : {}; Threads: {}; Chunk size: {} \n", N, Ncombs, num_threads, chunk_size);

        // auto ibefore = 0u;
        // for (auto hash : new_nodes)
        // {
        //     fmt::print("BEFORE THREADS {}: {} {}\n", ++ibefore, GNM[hash].to_str(), hash);
        // }
        
        std::vector<std::thread> threads;
        threads.reserve(num_threads);
        for (ULL start_k = 0u, i = 0u; start_k < Ncombs; start_k += chunk_size, i++)
        {
            ULL end_k = std::min(start_k + chunk_size, Ncombs);
            // fmt::print("\t\t\tcreating thread from {} to {} (out of {}) \n", start_k, end_k, Ncombs);
            threads.push_back(
                std::thread(threaded_new_new, std::ref(new_nodes), std::ref(funcs), std::ref(xorables), std::ref(costs), start_k, end_k, COSTS[fCB], 0, tgt_lvl)
            );
        }
        join_threads(threads);

        // auto iafter = 0u;
        // for (auto hash : new_nodes)
        // {
        //     fmt::print("AFTER THREADS {}: {} {}\n", ++iafter, GNM[hash].to_str(), hash);
        // }

        fmt::print("\t\tCombining new - new pairs\n");


        for (ULL k = 0u; k < Ncombs; k++)
        {
            ULL i = N - 2 - ULL(sqrt(4*N*(N - 1) - 8*k - 7)/2.0 - 0.5);       // ULL hi = new_nodes[i]; // Node & ni = GNM[hi];
            ULL j = k + i + 1 - N * (N - 1) / 2 + (N - i)*((N - i) - 1) / 2; // ULL hj = new_nodes[j]; // Node & nj = GNM[hj];
            UI cost = costs[k];
            // fmt::print("\tCombining: {} (#{}) and {} (#{}), cost {}, func {:016b}, xorable {}\n", GNM[new_nodes[i]].to_str(), i, GNM[new_nodes[j]].to_str(), j, cost, funcs[k], xorables[k]);
            // fmt::print("\tBest any l-{}, c-{}, best xorable l-{}, c-{}\n", get_gea(funcs[k]).lvl, get_gea(funcs[k]).cost, get_gex(funcs[k]).lvl, get_gex(funcs[k]).cost);
            if (cost >= INF) continue;
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
        ULL Nchunks = ceil(Ncombs, MAX_SIZE);
        ULL k_per_thread = Ncombs / Nchunks / NUM_THREADS; 
        k_per_thread += (k_per_thread == 0);
        ULL Ncombs_chunk = k_per_thread * NUM_THREADS;
        std::vector<UI> funcs(Ncombs_chunk, 0);
        std::vector<bool> xorables(Ncombs_chunk, false);
        std::vector<UI> costs(Ncombs_chunk, INF);

        std::vector<std::thread> threads;
        threads.reserve(num_threads);
        ULL offset = 0u;
        for (ULL start_k = 0u; start_k < Ncombs; start_k += k_per_thread)
        {
            ULL end_k = std::min(start_k + k_per_thread, Ncombs);
            threads.push_back(
                std::thread(threaded_new_new, std::ref(new_nodes), std::ref(funcs), std::ref(xorables), std::ref(costs), start_k, end_k, COSTS[fCB], offset, tgt_lvl)
            );
            if (threads.size() == num_threads || end_k == Ncombs)
            {
                join_threads(threads);
                fmt::print("\t\tCombining new - new combinations\n");
                for (ULL k = offset; k < end_k; k++)
                {
                    ULL i = N - 2 - ULL(sqrt(4*N*(N - 1) - 8*k - 7)/2.0 - 0.5);       // ULL hi = new_nodes[i]; // Node & ni = GNM[hi];
                    ULL j = k + i + 1 - N * (N - 1) / 2 + (N - i)*((N - i) - 1) / 2; // ULL hj = new_nodes[j]; // Node & nj = GNM[hj];
                    UI cost = costs[k - offset];
                    if (cost >= INF) continue;
                    auto [nhash, added] = create_node(funcs[k - offset], fCB, costs[k - offset], depth, xorables[k - offset], {new_nodes[i], new_nodes[j]});
                    if (added)
                    {
                        // check_node(nhash);
                        fresh_nodes.push_back(nhash);
                    }
                }

                // Reset the variables for the next set of threads
                threads.clear();
                offset = end_k;
                std::fill(funcs.begin(), funcs.end(), 0u);
                std::fill(xorables.begin(), xorables.end(), false);
                std::fill(costs.begin(), costs.end(), INF);
            }
        }
    }
}

void threaded_xor(const std::vector<ULL>& hashes, std::vector<UI>& funcs, std::vector<UI>& costs, const ULL start_k, const ULL end_k, const UI gate_cost, const ULL offset, const UI tgt_lvl) 
{
    // fmt::print("\t\tFunction is called\n");
    const ULL N = hashes.size();
    fmt::print("\t\tXOR: Iterating between {:L} and {:L}\n", start_k, end_k);
    ULL old_i;
    UI init_cost;
    std::unordered_map<ULL, UI> ct_spl;
    std::unordered_set<ULL> non_splittable_nodes;
    for (ULL k = start_k; k < end_k; k++)
    {   
        ULL i = N - 2 - UL(sqrt(4*N*(N - 1) - 8*k - 7)/2.0 - 0.5);
        ULL hi = hashes[i];     // fmt::print("\t\t hash_i={}\n", hi);
        Node & ni = GNM[hi];    // fmt::print("\t\t Retrieved ni={}\n", ni.to_str());

        if (i != old_i || k == start_k) // update [init_cost, ct_spl, non_splittable_nodes] if the row has changed
        {
            // fmt::print("\tUpdating\n");
            auto [init_cost_tmp, ct_spl_tmp, non_splittable_nodes_tmp] = node_cost_single(hi, gate_cost);
            init_cost = init_cost_tmp;
            ct_spl = ct_spl_tmp;
            non_splittable_nodes = non_splittable_nodes_tmp;
        }

        ULL j = k + i + 1 - N * (N - 1) / 2 + (N - i)*((N - i) - 1) / 2;    // fmt::print("\t\tk={}, i={}, j={}\n", k, i, j);
        ULL hj = hashes[j];     // fmt::print("\t\t hash_j={}\n", hj);
        Node & nj = GNM[hj];    // fmt::print("\t\t Retrieved nj={}\n", nj.to_str());

        TT func = ni.func ^ nj.func;
        Node & best_n = get_gex(func);
        if (best_n.lvl >= tgt_lvl) 
        {
            auto [cost, status] = node_cost(hj, init_cost, ct_spl, non_splittable_nodes, best_n.cost); // costs[i][j] = node_cost(ni, nj, gate_cost);
            if (status)
            {
                funcs[k - offset] = func;  // fmt::print("\t\t Wrote func = {0:016b} = {0:04x}\n",  funcs[k]);
                costs[k - offset] = cost;   // fmt::print("\t\t Cost = {}\n", costs[k]);
            }
        }
        old_i = i;
    }
}

void threaded_xor_wrapper(const std::vector<ULL>& nodes, const UI depth, UI num_threads = NUM_THREADS)
{
    // fmt::print("\t\tEntered wrapper\n");
    std::vector<ULL> valid_nodes;
    for (auto hash : nodes)
    {   
        if (GNM[hash].xorable)
        {
            valid_nodes.push_back(hash);
        }
    }
    
    const ULL N = valid_nodes.size();
    if (N <= 1) return;

    const ULL Ncombs = N * (N - 1) / 2;
    const UI tgt_lvl = depth / 3;
    // TODO : bring splitting to XOR and AND/OR
    if (Ncombs < MAX_SIZE)
    {
        ULL chunk_size = Ncombs / num_threads;
        if (chunk_size == 0)
        {
            num_threads = 1;
            chunk_size = Ncombs;
        }
        
        std::vector<UI> funcs(Ncombs);
        std::vector<UI> costs(Ncombs);

        // fmt::print("\t\tCalculating func and cost for {} nodes. Threads: {}; Chunk size: {} \n", N, num_threads, chunk_size);
        std::vector<std::thread> threads;
        threads.reserve(num_threads);
        // Launch threads to square elements of the vector
        for (UI i = 0; i < num_threads; ++i) {
            const int start_k = i * chunk_size;
            const int end_k = (i == num_threads - 1) ? Ncombs : (i + 1) * chunk_size;
            // fmt::print("\t\tChunk #{}: Iterating k from {} to {} \n", i, start_k, end_k);
            threads.push_back(
                std::thread(threaded_xor, std::ref(valid_nodes), std::ref(funcs), std::ref(costs), start_k, end_k, COSTS[fXOR], 0, tgt_lvl)
            );
        }
        join_threads(threads);

        for (ULL k = 0u; k < Ncombs; k++)
        {
            ULL i = N - 2 - ULL(sqrt(4*N*(N - 1) - 8*k - 7)/2.0 - 0.5);           // ULL hi = new_nodes[i];  // Node & ni = GNM[hi];
            ULL j = k + i + 1 - N * (N - 1) / 2 + (N - i)*((N - i) - 1) / 2;     // ULL hj = new_nodes[j];  // Node & nj = GNM[hj];
            if (costs[k] >= INF) continue;
            create_node(funcs[k], fXOR, costs[k], depth, true, {valid_nodes[i], valid_nodes[j]}); //auto [nhash, added] = 
        }
    }
    else //too many combinations to fit into memory
    {
        fmt::print("\t\tToo many pairs among {:L} xorable nodes to fit into memory (total {:L}). Splitting... \n", N, Ncombs);
        ULL Nchunks = ceil(Ncombs, MAX_SIZE);
        ULL k_per_thread = Ncombs / Nchunks / NUM_THREADS; 
        k_per_thread += (k_per_thread == 0);
        ULL Ncombs_chunk = k_per_thread * NUM_THREADS;
        std::vector<UI> funcs(Ncombs_chunk, 0);
        std::vector<UI> costs(Ncombs_chunk, INF);

        std::vector<std::thread> threads;
        threads.reserve(num_threads);
        ULL offset = 0u;
        for (ULL start_k = 0u; start_k < Ncombs; start_k += k_per_thread)
        {
            ULL end_k = std::min(start_k + k_per_thread, Ncombs);
            threads.push_back(
                std::thread(threaded_xor, std::ref(valid_nodes), std::ref(funcs), std::ref(costs), start_k, end_k, COSTS[fXOR], offset, tgt_lvl)
            );
            if (threads.size() == num_threads || end_k == Ncombs)
            {
                join_threads(threads);
                fmt::print("\t\tCombining XOR results\n");
                for (ULL k = offset; k < end_k; k++)
                {
                    ULL i = N - 2 - ULL(sqrt(4*N*(N - 1) - 8*k - 7)/2.0 - 0.5);       // ULL hi = new_nodes[i]; // Node & ni = GNM[hi];
                    ULL j = k + i + 1 - N * (N - 1) / 2 + (N - i)*((N - i) - 1) / 2; // ULL hj = new_nodes[j]; // Node & nj = GNM[hj];
                    UI cost = costs[k - offset];
                    if (cost >= INF) continue;
                    create_node(funcs[k - offset], fXOR, costs[k - offset], depth, true, {valid_nodes[i], valid_nodes[j]});
                }

                // Reset the variables for the next set of threads
                // std::vector<std::thread> threads(num_threads);
                threads.clear();
                offset = end_k;
                std::fill(funcs.begin(), funcs.end(), 0u);
                std::fill(costs.begin(), costs.end(), INF);
            }
        }
    }
}

void threaded_and_or(const std::vector<ULL>& hashes, std::vector<UI>& funcs_and, std::vector<UI>& costs_and, std::vector<UI>& funcs_or, std::vector<UI>& costs_or, const ULL start_k, const ULL end_k, const UI gate_cost_and, const ULL offset, const UI tgt_lvl) 
{
    // fmt::print("\t\tFunction is called\n");
    const ULL N = hashes.size();
    fmt::print("\t\tAND/OR: Iterating between {:L} and {:L}\n", start_k, end_k);
    ULL old_i;
    UI init_cost;
    std::unordered_map<ULL, UI> ct_spl;
    std::unordered_set<ULL> non_splittable_nodes;
    for (ULL k = start_k; k < end_k; k++)
    {   
        ULL i = N - 2 - ULL(sqrt(4*N*(N - 1) - 8*k - 7)/2.0 - 0.5);
        ULL hi = hashes[i]; // fmt::print("\t\t hash_i={}\n", hi);
        Node & ni = GNM[hi]; // fmt::print("\t\t Retrieved ni={}\n", ni.to_str());
        if (i != old_i || k == start_k) // update [init_cost, ct_spl, non_splittable_nodes] if the row has changed
        {
            // fmt::print("\tUpdating\n");
            auto [init_cost_tmp, ct_spl_tmp, non_splittable_nodes_tmp] = node_cost_single(hi, gate_cost_and);
            init_cost = init_cost_tmp;
            ct_spl = ct_spl_tmp;
            non_splittable_nodes = non_splittable_nodes_tmp;
        }

        ULL j = k + i + 1 - N * (N - 1) / 2 + (N - i)*((N - i) - 1) / 2; // fmt::print("\t\tk={}, i={}, j={}\n", k, i, j);
        ULL hj = hashes[j]; // fmt::print("\t\t hash_j={}\n", hj);
        Node & nj = GNM[hj]; // fmt::print("\t\t Retrieved nj={}\n", nj.to_str());

        TT func_and = ni.func & nj.func;
        TT func_or  = ni.func | nj.func;

        Node & best_and = get_gex(func_and);
        Node & best_or  = get_gex(func_or );
        if (best_and.lvl >= tgt_lvl && best_or.lvl >= tgt_lvl) 
        {
            UI limit = std::max(best_and.cost, best_or.cost);
            auto [cost_and, status] = node_cost(hj, init_cost, ct_spl, non_splittable_nodes, limit); // costs[i][j] = node_cost(ni, nj, gate_cost);
            if (status)
            {
                UI cost_or = cost_and - COSTS[fAND] + COSTS[fOR];
                funcs_and[k - offset] = func_and;  // fmt::print("\t\t Wrote func = {0:016b} = {0:04x}\n",  funcs[k]);
                costs_and[k - offset] = cost_and;   // fmt::print("\t\t Cost = {}\n", costs[k]);
                funcs_or[k - offset] = func_or;  // fmt::print("\t\t Wrote func = {0:016b} = {0:04x}\n",  funcs[k]);
                costs_or[k - offset] = cost_or;   // fmt::print("\t\t Cost = {}\n", costs[k]);
            }
        }
        else if (best_and.lvl >= tgt_lvl)
        {
            auto [cost_and, status] = node_cost(hj, init_cost, ct_spl, non_splittable_nodes, best_and.cost); // costs[i][j] = node_cost(ni, nj, gate_cost);
            if (status)
            {
                funcs_and[k - offset] = func_and;  // fmt::print("\t\t Wrote func = {0:016b} = {0:04x}\n",  funcs[k]);
                costs_and[k - offset] = cost_and;   // fmt::print("\t\t Cost = {}\n", costs[k]);
            }
        }
        else if (best_or.lvl >= tgt_lvl) 
        {
            auto [cost_or, status] = node_cost(hj, init_cost, ct_spl, non_splittable_nodes, best_or.cost); // costs[i][j] = node_cost(ni, nj, gate_cost);
            if (status)
            {
                funcs_or[k - offset] = func_or;  // fmt::print("\t\t Wrote func = {0:016b} = {0:04x}\n",  funcs[k]);
                costs_or[k - offset] = cost_or;   // fmt::print("\t\t Cost = {}\n", costs[k]);
            }
        }
        old_i = i;
    }
}

void threaded_and_or_wrapper(const std::vector<ULL>& nodes, const UI depth, UI num_threads = NUM_THREADS)
{
    const ULL N = nodes.size();
    if (N <= 1) return;

    const ULL Ncombs = N * (N - 1) / 2;
    const UI tgt_lvl = depth / 3;
    // TODO : bring splitting to XOR and AND/OR
    if (Ncombs < MAX_SIZE)
    {
        ULL chunk_size = Ncombs / num_threads;
        if (chunk_size == 0)
        {
            num_threads = 1;
            chunk_size = Ncombs;
        }
        
        std::vector<UI> funcs_and(Ncombs);
        std::vector<UI> costs_and(Ncombs);
        std::vector<UI> funcs_or(Ncombs);
        std::vector<UI> costs_or(Ncombs);

        // fmt::print("\t\tCalculating func, xorable, and cost for {} new nodes. Threads: {}; Chunk size: {} \n", N, num_threads, chunk_size);
        std::vector<std::thread> threads;
        threads.reserve(num_threads);
        // Launch threads to square elements of the vector
        for (UI i = 0; i < num_threads; ++i) {
            const int start_k = i * chunk_size;
            const int end_k = (i == num_threads - 1) ? Ncombs : (i + 1) * chunk_size;
            // fmt::print("\t\tChunk #{}: Iterating k from {} to {} \n", i, start_k, end_k);
            threads.push_back(
                std::thread(threaded_and_or, std::ref(nodes), std::ref(funcs_and), std::ref(costs_and), std::ref(funcs_or), std::ref(costs_or), start_k, end_k, COSTS[fAND], 0, tgt_lvl)
            );
        }
        join_threads(threads);

        for (ULL k = 0u; k < Ncombs; k++)
        {
            ULL i = N - 2 - ULL(sqrt(4*N*(N - 1) - 8*k - 7)/2.0 - 0.5);           // ULL hi = new_nodes[i];  // Node & ni = GNM[hi];
            ULL j = k + i + 1 - N * (N - 1) / 2 + (N - i)*((N - i) - 1) / 2;     // ULL hj = new_nodes[j];  // Node & nj = GNM[hj];
            if (costs_and[k] <= INF)
            {
                auto [nhash_and, added_and] = create_node(funcs_and[k], fAND, costs_and[k], depth, true, {nodes[i], nodes[j]});//
            }
            if (costs_or[k] <= INF)
            {
                auto [nhash_or , added_or ] = create_node(funcs_or[k],   fOR, costs_or[k] , depth, true, {nodes[i], nodes[j]});//
            }
        }
    }
    else
    {
        fmt::print("\t\tToo many pairs among {:L} AS nodes to fit into memory (total {:L}). Splitting... \n", N, Ncombs);
        ULL Nchunks = ceil(Ncombs, MAX_SIZE);
        ULL k_per_thread = Ncombs / Nchunks / NUM_THREADS; 
        k_per_thread += (k_per_thread == 0);
        ULL Ncombs_chunk = k_per_thread * NUM_THREADS;
        std::vector<UI> funcs_and(Ncombs_chunk, 0);
        std::vector<UI> costs_and(Ncombs_chunk, INF);
        std::vector<UI> funcs_or(Ncombs_chunk, INF);
        std::vector<UI> costs_or(Ncombs_chunk, INF);

        std::vector<std::thread> threads(num_threads);
        threads.reserve(num_threads);
        ULL offset = 0u;
        for (ULL start_k = 0u; start_k < Ncombs; start_k += k_per_thread)
        {
            ULL end_k = std::min(start_k + k_per_thread, Ncombs);
            threads.push_back(
                std::thread(threaded_and_or, std::ref(nodes), std::ref(funcs_and), std::ref(costs_and), std::ref(funcs_or), std::ref(costs_or), start_k, end_k, COSTS[fAND], offset, tgt_lvl)
            );
            if (threads.size() == num_threads || end_k == Ncombs)
            {
                join_threads(threads);
                fmt::print("\t\tCombining AND/OR results\n");
                for (ULL k = offset; k < end_k; k++)
                {
                    ULL i = N - 2 - ULL(sqrt(4*N*(N - 1) - 8*k - 7)/2.0 - 0.5);       // ULL hi = new_nodes[i]; // Node & ni = GNM[hi];
                    ULL j = k + i + 1 - N * (N - 1) / 2 + (N - i)*((N - i) - 1) / 2; // ULL hj = new_nodes[j]; // Node & nj = GNM[hj];
                    // create_node(funcs[k - offset], fXOR, costs[k - offset], depth, true, {valid_nodes[i], valid_nodes[j]});
                    if (costs_and[k] <= INF)
                    {
                        create_node(funcs_and[k-offset], fAND, costs_and[k-offset], depth, true, {nodes[i], nodes[j]});//auto [nhash_and, added_and]
                    }
                    if (costs_or[k] <= INF)
                    {
                        create_node(funcs_or[k-offset],  fOR , costs_or[k-offset] , depth, true, {nodes[i], nodes[j]});//auto [nhash_or , added_or ]
                    }
                }

                // Reset the variables for the next set of threads
                // std::vector<std::thread> threads(num_threads);
                threads.clear();
                offset = end_k;
                std::fill(funcs_and.begin(), funcs_and.end(), 0u);
                std::fill(costs_and.begin(), costs_and.end(), INF);
                std::fill(funcs_or.begin(), funcs_or.end(), 0u);
                std::fill(costs_or.begin(), costs_or.end(), INF);
            }
        }
    }
}

std::tuple<ULL, bool> check_cb(Node& ni, Node& nj, std::vector<ULL>& fresh_nodes)
{
    // fmt::print("\t\tA: {}\n", ni.to_str());
    // fmt::print("\t\tB: {}\n", nj.to_str());
    TT func = ni.func | nj.func;
    if (func == ni.func || func == nj.func || func == 0 || func == ONES)
    {
        // Node& ref = GEX[func];
        return std::make_tuple(GEX[func], false);
    }
    bool xorable = (func == (ni.func ^ nj.func));
    US depth = ni.lvl * 3 + 1;
    UI cost = node_cost(ni, nj, COSTS[fCB]);
    auto [nhash, added] = create_node(func, fCB, cost, depth, xorable, {ni.hash, nj.hash});
    if (added)
    {
        // Node& n =  GNM[nhash];
        fresh_nodes.push_back(nhash);
        // fmt::print("\t\tC: {}\n", n .to_str());
    }
    return std::make_tuple(nhash, added);
}

std::tuple<ULL, bool> check_xor(Node& ni, Node& nj)
{
    // fmt::print("\t\tA: {}\n", ni.to_str());
    // fmt::print("\t\tB: {}\n", nj.to_str());
    TT func = ni.func ^ nj.func;
    // if (func == 0 || func == ONES)
    // {
    //     Node& ref = GEX[func];
    //     return std::make_tuple(ref.hash, false);
    // }
    US depth = ni.lvl * 3 + 2;
    UI cost = node_cost(ni, nj, COSTS[fXOR]);
    auto [nhash, added] = create_node(func, fXOR, cost, depth, true, {ni.hash, nj.hash}); // 
    // if (added)
    // {
    //     Node& n =  GNM[nhash];
    //     fmt::print("\t\tC: {}\n", n .to_str());
    // }
    return std::make_tuple(nhash, added);
}

std::tuple<ULL, bool, ULL, bool> check_and_or(Node& ni, Node& nj)
{
    // fmt::print("\t\tA: {}\n", ni.to_str());
    // fmt::print("\t\tB: {}\n", nj.to_str());
    US depth = ni.lvl * 3 + 3;
    TT func_and = ni.func & nj.func;
    TT func_or  = ni.func | nj.func;
    // fmt::print("\t\t\t&: {}\n", func_and);
    // fmt::print("\t\t\t&: {}\n", func_or );
    // if (func == 0 || func == ONES)
    // {
    //     Node& ref = GEX[func];
    //     return std::make_tuple(ref.hash, false);
    // }
    UI cost_and = node_cost(ni, nj, COSTS[fAND]);
    UI cost_or  = cost_and - COSTS[fAND] + COSTS[fOR];
    auto [hash_and, added_and] = create_node(func_and, fAND, cost_and, depth, true, {ni.hash, nj.hash});
    auto [hash_or , added_or ] = create_node(func_or , fOR , cost_or , depth, true, {ni.hash, nj.hash});
    // if (added_and)
    // {
    //     Node& n =  GNM[hash_and];
    //     fmt::print("\t\t&: {}\n", n .to_str());
    // }
    // if (added_or )
    // {
    //     Node& n =  GNM[hash_or ];
    //     fmt::print("\t\t|: {}\n", n .to_str());
    // }
    return std::make_tuple(hash_and, added_and, hash_or, added_or);
}

std::unordered_set<ULL> remove_dominated()
{
    std::unordered_set<ULL> protected_nodes;
    std::stack<ULL> stack;
    for (auto & [hash, node] : GNM)
    {
        if (node.last_func == fDFF) continue;
        Node& best = node.xorable? get_gex(node.func) : get_gea(node.func); //GEX[node.func] : GEA[node.func];
        if (node.lvl < best.lvl || (node.lvl == best.lvl && node.cost <= best.cost) )
        {
            // fmt::print("\t\t\t\tProtecting node {} (best is {}) [{} {}]\n", node.to_str(), best.to_str(), node.hash, best.hash);
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
    // fmt::print("\t\t\t\tProtecting {} nodes out of {}\n", protected_nodes.size(), GNM.size());

    std::unordered_set<ULL> to_be_removed;

    for (auto & [hash, node] : GNM)
    {
        if (node.last_func == fDFF) continue;
        if (protected_nodes.find(hash) != protected_nodes.end()) continue;
        // Node& best = node.xorable? GEX[node.func] : GEA[node.func];
        Node& best = node.xorable? get_gex(node.func) : get_gea(node.func); //GEX[node.func] : GEA[node.func];
        if (node.lvl > best.lvl || (node.lvl == best.lvl && node.cost > best.cost) )
        {
            to_be_removed.insert(hash);
        }
    }

    for (auto hash : to_be_removed)
    {
        // fmt::print("\t\t\t\tDestroying dominated node {}\n", GNM[hash].to_str());
        GNM.erase(hash);
    }   
    return to_be_removed;
}

void cb_generation(US lvl)
{
    std::vector<ULL> new_nodes = select_depth(lvl * 3, lvl * 3 + 1);
    fmt::print("\tChecking new nodes after selection between levels {:L} and {:L}...\n", lvl * 3, lvl * 3 + 1);
    // check_nodes(new_nodes);

    US depth = lvl * 3 + 1;
    std::vector<ULL> old_nodes;
    bool go_on;
    std::vector<ULL> fresh_nodes;
    UI iteration = 0;
    do {
    #if (THREADING_MODE == 0)
        // first, combine new nodes with old nodes
        fmt::print("\tProcessing {:L} old nodes with {:L} new nodes\n", old_nodes.size(), new_nodes.size());
        thread_old_new_wrapper(old_nodes, new_nodes, depth, fresh_nodes);
        
        // fmt::print("\t{}: Checking old nodes after old-new...\n", iteration);
        // check_nodes(old_nodes);
        // fmt::print("\t{}: Checking fresh nodes after old-new...\n", iteration);
        // check_nodes(fresh_nodes);
        // fmt::print("\t{}: Checking new nodes after old-new...\n", iteration);
        // check_nodes(new_nodes);
        // fmt::print("\t{}: Checking GNM after old-new...\n", iteration);
        // check_GNM();

        write_csv_gnm(GNM, fmt::format("output_gnm_cb_{}_{}.csv", lvl, iteration));
        write_csv_arr(GEA, fmt::format("output_gea_cb_{}_{}.csv", lvl, iteration));
        write_csv_arr(GEX, fmt::format("output_gex_cb_{}_{}.csv", lvl, iteration));

        // next, combine new nodes 
        // fmt::print("\tProcessing pairs among {} new nodes\n", new_nodes.size());
        // check_GNM();
        // fmt::print("{}: GNM size before new-new : {}\n", iteration, GNM.size());   
        threaded_new_new_wrapper(new_nodes, depth, fresh_nodes);
        // fmt::print("{}: GNM size after new-new : {}\n", iteration, GNM.size());   
        // check_GNM();

        write_csv_gnm(GNM, fmt::format("output_gnm_cb_{}_{}.csv", lvl, iteration));
        write_csv_arr(GEA, fmt::format("output_gea_cb_{}_{}.csv", lvl, iteration));
        write_csv_arr(GEX, fmt::format("output_gex_cb_{}_{}.csv", lvl, iteration));

    #elif (THREADING_MODE == 1)
            // first, combine new nodes with old nodes
            fmt::print("\tProcessing old {} nodes with new {} nodes\n", old_nodes.size(), new_nodes.size());
            #pragma omp parallel for num_threads(10)
            for (ULL k = 0u; k < new_nodes.size() * old_nodes.size(); k++)
            {
                ULL i = k / old_nodes.size();
                ULL j = k % old_nodes.size();
                fmt::print("\tCB ({}, {}, {}):\n", i, j, k);
                Node& ni =  GNM[new_nodes[i]];
                Node& nj =  GNM[old_nodes[j]];
                auto [nhash, added] = check_cb(ni, nj, fresh_nodes);
                go_on |= added;
            }
            // next, combine new nodes 
            fmt::print("\tProcessing new {} nodes\n", new_nodes.size());
            #pragma omp parallel for num_threads(10)
            for (ULL k = 0u, n = new_nodes.size(); k < n * (n - 1) / 2; k++)
            {
                ULL i = n - 2 - ULL(sqrt(4*n*(n-1) - 8*k - 7)/2.0 - 0.5);
                ULL j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
                fmt::print("\tCB ({}, {}, {}):\n", i, j, k);
                Node& ni =  GNM[new_nodes[i]];
                Node& nj =  GNM[new_nodes[j]];
                auto [nhash, added] = check_cb(ni, nj, fresh_nodes);
                go_on |= added;
            }
    #elif (THREADING_MODE == 2)
        // first, combine new nodes with old nodes
        fmt::print("\tProcessing {} old nodes with {} new nodes\n", old_nodes.size(), new_nodes.size());
        auto i = 0u;
        for (auto hi : old_nodes)
        {
            fmt::print("\nCompleted {} out of {} ({:f}\%) \r", i++, old_nodes.size(), 100 * (float)i / (float)old_nodes.size());
            Node& ni =  GNM[hi];
            auto [gate_cost, ct_spl, non_splittable_nodes] = node_cost_single(hi, COSTS[fCB]);
            for (auto hj : new_nodes)
            {
                Node& nj =  GNM[hj];
                TT func = ni.func | nj.func;
                if (func == ni.func || func == nj.func || func == 0 || func == ONES) continue;
                // if (func == ni.func || func == nj.func || func == 0 || func == ONES) continue;
                bool xorable = (func == (ni.func ^ nj.func));
                //if inputs are already more expensive, skip
                UI limit = (xorable ? GNM[GEX[func]].cost : GNM[GEA[func]].cost);
                // fmt::print("calculating cb cost for {} and {} up to limit {}\n", ni.to_str(), nj.to_str(), INF);
                auto [cost, status] = node_cost(hj, gate_cost, ct_spl, non_splittable_nodes, limit);
                // fmt::print("final cost is {} and status is {}\n", cost, status);
                if (!status) continue;
                auto [nhash, added] = create_node(func, fCB, cost, depth, xorable, {ni.hash, nj.hash});
                if (added)
                {
                    Node& n =  GNM[nhash];
                    fresh_nodes.push_back(nhash);
                    // fmt::print("\t\tC: {}\n", n .to_str());
                    go_on = true;
                }
            }
        }

        // next, combine new nodes 
        for (auto i = 0u; i < new_nodes.size(); i++)
        {
            fmt::print("\nCompleted {} out of {} ({:f}\%) \r", i, new_nodes.size(), 100 * (float)i / (float)new_nodes.size());
            ULL hi = new_nodes[i];
            Node& ni =  GNM[hi];
            auto [gate_cost, ct_spl, non_splittable_nodes] = node_cost_single(hi, COSTS[fCB]);
            for (auto j = i+1; j < new_nodes.size(); j++)
            {
                ULL hj = new_nodes[j];
                Node& nj =  GNM[hj];
                TT func = ni.func | nj.func;
                if (func == ni.func || func == nj.func || func == 0 || func == ONES) continue;
                bool xorable = (func == (ni.func ^ nj.func));

                UI limit = (xorable ? GNM[GEX[func]].cost : GNM[GEA[func]].cost);
                // fmt::print("calculating cb cost for {} and {} up to limit {}\n", ni.to_str(), nj.to_str(), INF);
                auto [cost, status] = node_cost(hj, gate_cost, ct_spl, non_splittable_nodes, limit);
                // fmt::print("1 final cost is {} and status is {}\n", cost, status);
                if (!status) continue;
                // fmt::print("2 final cost is {} and status is {}\n", cost, status);
                // UI cost = node_cost(hj, gate_cost, ct_spl, non_splittable_nodes);
                
                auto [nhash, added] = create_node(func, fCB, cost, depth, xorable, {ni.hash, nj.hash});
                if (added)
                {
                    Node& n =  GNM[nhash];
                    fresh_nodes.push_back(nhash);
                    // fmt::print("\t\tC: {}\n", n .to_str());
                }
            }
        }
    #endif

        // fmt::print("\t{}: Checking old nodes after new-new...\n", iteration);
        // check_nodes(old_nodes);
        // fmt::print("\t{}: Checking fresh nodes after new-new...\n", iteration);
        // check_nodes(fresh_nodes);
        // fmt::print("\t{}: Checking new nodes after new-new...\n", iteration);
        // check_nodes(new_nodes);
        
        std::unordered_set<ULL> removed_nodes = remove_dominated();
        for (auto new_hash : new_nodes)
        {
            if (removed_nodes.find(new_hash) == removed_nodes.end())
            {
                old_nodes.push_back(new_hash);
            }
        }
        new_nodes.clear();
        for (auto fresh_hash : fresh_nodes)
        {
            if (removed_nodes.find(fresh_hash) == removed_nodes.end())
            {
                new_nodes.push_back(fresh_hash);
            }
        }
        // new_nodes.clear();
        // old_nodes.insert(old_nodes.end(), new_nodes.begin(), new_nodes.end());
        // new_nodes.clear();
        // new_nodes.insert(new_nodes.end(), fresh_nodes.begin(), fresh_nodes.end());
        go_on = fresh_nodes.size() > 0;
        fresh_nodes.clear(); 
        iteration++;
        // fmt::print("GOON VAR IS {}", go_on);
    } while(go_on);
}

void as_generation(US lvl)
{
    std::vector<ULL> nodes = select_depth(lvl * 3, lvl * 3 + 1);

    US tgt_depth = lvl * 3 + 2;
    fmt::print("\t{} : DFF/NOT {:L} nodes\n", lvl, nodes.size());
    for (ULL hash : nodes)
    {
        Node& ni =  GNM[hash];
        force_create_node(ni.func, fDFF, ni.cost + COSTS[fDFF], tgt_depth, true, {hash});
        create_node(ni.func ^ ONES, fNOT, ni.cost + COSTS[fNOT], tgt_depth, true, {hash});
    }
    // xor-combine the nodes
    
    fmt::print("\t{}: XOR {:L} nodes\n", lvl, nodes.size());
    #if (THREADING_MODE == 0)   
        threaded_xor_wrapper(nodes, tgt_depth);
    #elif (THREADING_MODE == 2)
        #pragma omp parallel for num_threads(10)
        for (ULL k = 0u, n = nodes.size(); k < n * (n - 1) / 2; k++)
        {
            ULL i = n - 2 - ULL(sqrt(4*n*(n-1) - 8*k - 7)/2.0 - 0.5);
            ULL j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
            fmt::print("\tXOR ({}, {}, {}):\n", i, j, k);
            Node& ni =  GNM[nodes[i]];
            Node& nj =  GNM[nodes[j]];
            if (!ni.xorable || !nj.xorable) continue;
            auto [nhash, added] = check_xor(ni, nj);
        }
    #elif (THREADING_MODE == 1)
        for (auto i = 0u; i < nodes.size(); i++)
        {
            fmt::print("\nCompleted {} out of {} ({:f}\%) \r", i, nodes.size(), 100 * (float)i / (float)nodes.size());
            ULL hi = nodes[i];
            Node& ni =  GNM[hi];
            if (!ni.xorable) continue;
            auto [gate_cost, ct_spl, non_splittable_nodes] = node_cost_single(hi, COSTS[fXOR]);
            for (auto j = i+1; j < nodes.size(); j++)
            {
                ULL hj = nodes[j];
                Node& nj =  GNM[hj];
                if (!nj.xorable) continue;

                // fmt::print("\t\tA: {}\n", ni.to_str());
                // fmt::print("\t\tB: {}\n", nj.to_str());
                TT func = ni.func ^ nj.func;
                // UI limit = std::min(GNM[GEX[func]].cost);
                auto [cost,status] = node_cost(hj, gate_cost, ct_spl, non_splittable_nodes, GNM[GEX[func]].cost);
                if (!status) continue;
                auto [nhash, added] = create_node(func, fXOR, cost, tgt_depth, true, {ni.hash, nj.hash});
                if (added)
                {
                    Node& n =  GNM[nhash];
                    // fmt::print("\t\tC: {}\n", n .to_str());
                }
            }
        }

    #endif
    remove_dominated();
}
void sa_generation(US lvl)
{
    std::vector<ULL> nodes = select_depth(lvl * 3 + 2, lvl * 3 + 2);
    US tgt_depth = lvl * 3 + 3;

    // US tgt_depth = lvl * 3 + 3;
    // xor-combine the nodes
    fmt::print("\t{}: AND/OR {:L} nodes\n", lvl, nodes.size());
    #if (THREADING_MODE == 0)
        threaded_and_or_wrapper(nodes, tgt_depth);
    #elif (THREADING_MODE == 2)
        #pragma omp parallel for num_threads(10)
        for (ULL k = 0u, n = nodes.size(); k < n * (n - 1) / 2; k++)
        {
            ULL i = n - 2 - ULL(sqrt(4*n*(n-1) - 8*k - 7)/2.0 - 0.5);
            ULL j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
            fmt::print("\tAND/OR ({}, {}, {}):\n", i, j, k);
            Node& ni =  GNM[nodes[i]];
            Node& nj =  GNM[nodes[j]];
            if (!ni.xorable || !nj.xorable) continue;
            auto [hash_and, added_and, hash_or, added_or] = check_and_or(ni, nj);
        }
    #elif (THREADING_MODE == 1)
        // ULL total_combs = nodes.size() * (nodes.size() - 1) / 2;
        for (auto i = 0u; i < nodes.size(); i++)
        {
            fmt::print("\nCompleted {} out of {} ({:f}\%) \r", i, nodes.size(), 100 * (float)i / (float)nodes.size());
            ULL hi = nodes[i];
            Node& ni =  GNM[hi];
            auto [gate_cost_and, ct_spl, non_splittable_nodes] = node_cost_single(hi, COSTS[fAND]);
            
            for (auto j = i+1; j < nodes.size(); j++)
            {
                ULL hj = nodes[j];
                Node& nj =  GNM[hj];

                TT func_and = ni.func & nj.func;
                TT func_or  = ni.func | nj.func;
                UI limit = std::max(GNM[GEX[func_and]].cost, GNM[GEX[func_or]].cost);

                auto [cost_and, status] = node_cost(hj, gate_cost_and, ct_spl, non_splittable_nodes, INF);
                UI cost_or = cost_and - COSTS[fAND] + COSTS[fOR];
                if (!status) continue;

                auto [hash_and, added_and] = create_node(func_and, fAND, cost_and, tgt_depth, true, {ni.hash, nj.hash});
                auto [hash_or , added_or ] = create_node(func_or , fOR , cost_or , tgt_depth, true, {ni.hash, nj.hash});
                // if (added_and)
                // {
                //     Node& n =  GNM[hash_and];
                //     fmt::print("\t\t&: {}\n", n .to_str());
                // }
                // if (added_or )
                // {
                //     Node& n =  GNM[hash_or ];
                //     fmt::print("\t\t|: {}\n", n .to_str());
                // }
            }
        }
    #endif 
    remove_dominated();
}

inline bool is_done()  
{
    #pragma vector
    for (auto i = 0u; i < NUM_TT; i++)
    {
        // if (GEX[i].cost == INF || GEX[i].depth == INF)
        if (get_gex(i).cost == INF || get_gex(i).depth == INF)
        {
            return false;
        }
    }
    return true;
}

inline UI count_done()  
{
    #pragma vector
    UI ctr = 0;
    for (auto i = 0u; i < NUM_TT; i++)
    {
        ctr += (get_gex(i).cost < INF);
        // if (GEX[i].cost == INF || GEX[i].depth == INF)
        // if (get_gex(i).cost == INF || get_gex(i).depth == INF)
        // {
        //     return false;
        // }
    }
    return ctr;
}

bool is_subset(const std::vector<ULL>& v1, const std::vector<ULL>& v2) {
    // Check if all elements of v1 are in v2
    return std::all_of(v1.begin(), v1.end(), [&v2](const int& val) {
        return std::find(v2.begin(), v2.end(), val) != v2.end();
    });
}

int main() {
    // fmt::print("Start:\n");
    // print_GNM();
    for (TT func : gen_pi_func(NUM_VARS))
    {
        create_node(func, fPI, 0, 0, true);
    }
    // fmt::print("After PI:\n");
    // print_GNM();
    create_node(   0, fPI, 0, 0, true);
    create_node(ONES, fPI, 0, 0, true);
    // fmt::print("After const:\n");
    // print_GNM();
    
    for (US lvl = 0; lvl < 5; lvl ++)
    {
        fmt::print("Processing lvl {}: CB\n", lvl);
        // fmt::print("Checking integrity of GNM before CB\n");
        // check_GNM();
        cb_generation(lvl);
        // fmt::print("Checking integrity of GNM after CB\n");
        // check_GNM();
        fmt::print("\n\n\t\nCompleted {}\n\n\n", count_done());
    write_csv_gnm(GNM, fmt::format("output_gnm_cb_{}.csv", lvl));
    write_csv_arr(GEA, fmt::format("output_gea_cb_{}.csv", lvl));
    write_csv_arr(GEX, fmt::format("output_gex_cb_{}.csv", lvl));
    //     if (is_done()) break;
        fmt::print("Processing lvl {}: AS\n", lvl);
        // fmt::print("Checking integrity of GNM before AS\n");
        // check_GNM();
        as_generation(lvl);
        // fmt::print("Checking integrity of GNM after AS\n");
        // check_GNM();
        fmt::print("\n\n\t\nCompleted {}\n\n\n", count_done());
    write_csv_gnm(GNM, fmt::format("output_gnm_as_{}.csv", lvl));
    write_csv_arr(GEA, fmt::format("output_gea_as_{}.csv", lvl));
    write_csv_arr(GEX, fmt::format("output_gex_as_{}.csv", lvl));
        if (is_done()) break;
        fmt::print("Processing lvl {}: SA\n", lvl);
        // fmt::print("Checking integrity of GNM before SA\n");
        // check_GNM();
        sa_generation(lvl);
        // fmt::print("Checking integrity of GNM after SA\n");
        // check_GNM();
        fmt::print("\n\n\t\nCompleted {}\n\n\n", count_done());
    write_csv_gnm(GNM, fmt::format("output_gnm_sa_{}.csv", lvl));
    write_csv_arr(GEA, fmt::format("output_gea_sa_{}.csv", lvl));
    write_csv_arr(GEX, fmt::format("output_gex_sa_{}.csv", lvl));
        if (is_done()) break;
    }

    // Write data to CSV files
    // write_csv_gnm(GNM, "output_gnm.csv");
    // write_csv_arr(GEA, "output_gea.csv");
    // write_csv_arr(GEX, "output_gex.csv");

    std::unordered_map<ULL, Node> GNM_new;
    std::array<ULL, NUM_TT> GEA_new;
    std::array<ULL, NUM_TT> GEX_new;

    // Read data from CSV files
    GNM_new = read_csv_gnm("output_gnm.csv");
    GEA_new = read_csv_arr("output_gea.csv");
    GEX_new = read_csv_arr("output_gex.csv");

    std::cout << GNM.size() << " " << GNM_new.size() << " " << (GNM == GNM_new) << std::endl;
    std::cout << GEA.size() << " " << GEA_new.size() << " " << (GEA == GEA_new) << std::endl;
    std::cout << GEX.size() << " " << GEX_new.size() << " " << (GEX == GEX_new) << std::endl;

    for (auto & [hash, n1] : GNM)
    {
        Node& n2 = GNM_new[hash];
        bool scalar_fields_equal = (n1.func == n2.func && n1.last_func == n2.last_func && n1.cost == n2.cost && n1.depth == n2.depth && n1.xorable == n2.xorable);
        std::vector<ULL> p1, p2;
        p1.insert(p1.end(), n1.parent_hashes.begin(), n1.parent_hashes.end());
        p2.insert(p2.end(), n2.parent_hashes.begin(), n2.parent_hashes.end());

        std::sort(p1.begin(), p1.end());
        std::sort(p2.begin(), p2.end());
        bool parents_equal = (p1 == p2);

        if (scalar_fields_equal && parents_equal)
        {
            continue;
        }
        else
        {
            fmt::print("Hash: {}\n\n", hash);
            fmt::print("\tfunc\t: {} - {}\n", n1.func, n2.func);
            fmt::print("\tlast_func: {} - {}\n", n1.func, n2.func);
            fmt::print("\tcost\t: {} - {}\n", n1.cost, n2.cost);
            fmt::print("\tdepth\t: {} - {}\n", n1.depth, n2.depth);
            fmt::print("\txorable\t: {} - {}\n", n1.xorable, n2.xorable);
            for (auto phash : n1.parent_hashes)
            {
                fmt::print("\t\tParents 1\t: {}\n", phash);
            }
            for (auto phash : n2.parent_hashes)
            {
                fmt::print("\t\tParents 2\t: {}\n", phash);
            }
        }
    }

    return 0;
}

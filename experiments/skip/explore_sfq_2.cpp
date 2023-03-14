#include <iostream>
#include <vector>
#include <tuple>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <stack>
#include <array>

typedef unsigned short US;
typedef unsigned int   UI;
typedef unsigned long  UL;
typedef const unsigned short CUS;
typedef const unsigned int   CUI;
typedef const unsigned long  CUL;
typedef unsigned long long ULL;
typedef UI TT;

constexpr UI NUM_VARS = 4u;
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

constexpr UI COSTS[] = {7, 9, 8, 8, 8, 7, 11, 11, 11, 8, 7};

constexpr UL NUM_TT = (1 << (1 << NUM_VARS));
constexpr TT ONES = NUM_TT - 1;
constexpr UL INF = 0xFFFFFF;

class Node {
public:
    UI func = 0;
    US last_func = fNOFUNC;
    UI cost = INF;
    UI depth = INF;
    bool xorable = false;
    std::vector<Node> parents;
    UI lvl = INF;
    ULL hash = 0;

    Node() = default;
    Node(const Node& other) : func(other.func), last_func(other.last_func), cost(other.cost), depth(other.depth), xorable(other.xorable), parents(other.parents), hash(other.hash), lvl(other.lvl) {}
    Node& operator=(const Node& other) {
    if (this != &other) {
        func = other.func;
        last_func = other.last_func;
        cost = other.cost;
        depth = other.depth;
        xorable = other.xorable;
        parents = other.parents;
        hash = other.hash;
        lvl = other.lvl;
        }
        return *this;
    }


    Node(UI _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<Node> _parents)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parents(_parents), lvl((depth + 1) / 3)
    {
        // Calculate hash based on the hashes of parents and specified fields
        std::size_t seed = 0;
        hash_combine(seed, func);
        hash_combine(seed, last_func);
        hash_combine(seed, cost);
        hash_combine(seed, depth);
        hash_combine(seed, xorable);
        for (const auto& parent : parents) {
            hash_combine(seed, parent.hash);
        }
        hash = seed;
    }

    Node(UI _func, US _last_func, UI _cost, UI _depth, bool _xorable)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parents{}, lvl((depth + 1) / 3)
    {
        // Calculate hash based on the hashes of parents and specified fields
        std::size_t seed = 0;
        hash_combine(seed, func);
        hash_combine(seed, last_func);
        hash_combine(seed, cost);
        hash_combine(seed, depth);
        hash_combine(seed, xorable);
        hash = seed;
    }

    Node(UI _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<Node> _parents, ULL _hash)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parents(_parents), lvl((depth + 1) / 3), hash(_hash)
    {
        // Hash is precalculated based on the hashes of parents and specified fields
    }

    // Equality operator
    bool operator==(const Node& other) const {
        return hash == other.hash;
    }

    ULL get_hash() {
        return hash;
    }

    std::string to_str() const
    {
        return "Func: " + std::to_string(func) + "|Last: " + std::to_string(last_func) + "|Cost: " + std::to_string(cost) + "|Depth: " + std::to_string(depth) + "|X: " + std::to_string(xorable);
    }

private:
    // Hash combiner
    template <typename T>
    static void hash_combine(std::size_t& seed, const T& val) {
        seed ^= std::hash<T>{}(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
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

class NodeWrapper {
public:
    NodeWrapper(Node node)
        : node_(node)
    {
        // Check if a node with the same hash already exists in the map
        auto it = node_map_.find(node.hash);
        if (it != node_map_.end()) {
            // Node already exists in map, use existing node
            node_ = it->second;
        } else {
            // Node doesn't exist in map, create new node and add to map
            node_map_[node.hash] = node_;
            Node & best_any = earliest_any[node_.func];
            Node & best_xor = earliest_xor[node_.func];
            // If the new node is better than any existing one
            if (std::tie(node_.lvl, node_.cost) < std::tie(best_any.lvl, best_any.cost))
            {
                earliest_any[node_.func] = node_;
                earliest_xor[node_.func] = node_;
            }
            // If the new node is better than any existing xorable one
            else if (std::tie(node_.lvl, node_.cost) < std::tie(best_xor.lvl, best_xor.cost))
            {
                earliest_xor[node_.func] = node_;
            }
        }
    }

    NodeWrapper(UI func, US last_func, UI cost, UI depth, bool xorable, std::vector<Node>& parents)
    {
        // Calculate hash based on the specified fields
        std::size_t seed = 0;
        hash_combine(seed, func);
        hash_combine(seed, last_func);
        hash_combine(seed, cost);
        hash_combine(seed, depth);
        hash_combine(seed, xorable);
        for (const auto& parent : parents) {
            hash_combine(seed, parent.hash);
        }
        ULL hash = seed;

        // Check if a node with the same hash already exists in the map
        auto it = node_map_.find(hash);
        if (it != node_map_.end()) {
            // Node already exists in map, use existing node
            node_ = it->second;
        } else {
            // Node doesn't exist in map, create new node and add to map
            node_ = Node(func, last_func, cost, depth, xorable, parents);
            node_map_[hash] = node_;
            Node & best_any = earliest_any[node_.func];
            Node & best_xor = earliest_xor[node_.func];
            // If the new node is better than any existing one
            if (std::tie(node_.lvl, node_.cost) < std::tie(best_any.lvl, best_any.cost))
            {
                earliest_any[node_.func] = node_;
                earliest_xor[node_.func] = node_;
            }
            // If the new node is better than any existing xorable one
            else if (std::tie(node_.lvl, node_.cost) < std::tie(best_xor.lvl, best_xor.cost))
            {
                earliest_xor[node_.func] = node_;
            }
        }
    }

    NodeWrapper(UI func, US last_func, UI cost, UI depth, bool xorable, std::vector<NodeWrapper>& parents)
    {
        // Calculate hash based on the specified fields
        std::size_t seed = 0;
        hash_combine(seed, func);
        hash_combine(seed, last_func);
        hash_combine(seed, cost);
        hash_combine(seed, depth);
        hash_combine(seed, xorable);
        for (const auto& parent : parents) {
            hash_combine(seed, parent.hash());
        }
        ULL hash = seed;

        // Check if a node with the same hash already exists in the map
        auto it = node_map_.find(hash);
        if (it != node_map_.end()) {
            // Node already exists in map, use existing node
            node_ = it->second;
        } else {
            // Node doesn't exist in map, create new node and add to map
            std::vector<Node> node_parents;
            for (const auto& parent : parents) {
                node_parents.push_back(parent.node());
            }
            node_ = Node(func, last_func, cost, depth, xorable, node_parents);
            node_map_[hash] = node_;
            Node & best_any = earliest_any[node_.func];
            Node & best_xor = earliest_xor[node_.func];
            // If the new node is better than any existing one
            if (std::tie(node_.lvl, node_.cost) < std::tie(best_any.lvl, best_any.cost))
            {
                earliest_any[node_.func] = node_;
                earliest_xor[node_.func] = node_;
            }
            // If the new node is better than any existing xorable one
            else if (std::tie(node_.lvl, node_.cost) < std::tie(best_xor.lvl, best_xor.cost))
            {
                earliest_xor[node_.func] = node_;
            }
        }
    }

    NodeWrapper(UI func, US last_func, UI cost, UI depth, bool xorable)
    {
        // Calculate hash based on the specified fields
        std::size_t seed = 0;
        hash_combine(seed, func);
        hash_combine(seed, last_func);
        hash_combine(seed, cost);
        hash_combine(seed, depth);
        hash_combine(seed, xorable);
        ULL hash = seed;

        // Check if a node with the same hash already exists in the map
        auto it = node_map_.find(hash);
        if (it != node_map_.end()) {
            // Node already exists in map, use existing node
            node_ = it->second;
        } else {
            // Node doesn't exist in map, create new node and add to map
            node_ = Node(func, last_func, cost, depth, xorable);
            node_map_[hash] = node_;
            Node & best_any = earliest_any[node_.func];
            Node & best_xor = earliest_xor[node_.func];
            // If the new node is better than any existing one
            if (std::tie(node_.lvl, node_.cost) < std::tie(best_any.lvl, best_any.cost))
            {
                earliest_any[node_.func] = node_;
                earliest_xor[node_.func] = node_;
            }
            // If the new node is better than any existing xorable one
            else if (std::tie(node_.lvl, node_.cost) < std::tie(best_xor.lvl, best_xor.cost))
            {
                earliest_xor[node_.func] = node_;
            }
        }
    }
    

    // ~NodeWrapper() {
    //     node_map_.erase(node_.hash);
    // }

    ULL hash() const {
        return node_.hash;
    }

    Node node() const {
        return node_;
    }

    US last_func() const {
        return node_.last_func;
    }

    size_t get_node_map_size()
    {
        return node_map_.size();
    }

    std::string to_str()
    {
        return node_.to_str();
    }

private:
    Node node_;
    std::unordered_map<ULL, Node>& node_map_ = GNM;
    std::array<Node, NUM_TT>& earliest_any = GEA;
    std::array<Node, NUM_TT>& earliest_xor = GEX;

    // Hash combiner
    template <typename T>
    static void hash_combine(std::size_t& seed, const T& val) {
        seed ^= std::hash<T>{}(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
};

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

std::vector<NodeWrapper> gen_pi(UI nvars)
{
    std::vector<NodeWrapper> pi;
    for (auto func : gen_pi_func(nvars))
    {
        std::cout << GNM.size() << std::endl;
        auto qqq = NodeWrapper(func, fPI, 0, 1, true);
        std::cout << GNM.size() << std::endl;

        pi.push_back( qqq );
    }
    return pi;
}

std::vector<ULL> select_lvl(US min_lvl, US max_lvl)
{
    std::vector<ULL> out;
    for (auto & [key, n] : GNM)
    {
        if (n.lvl >= min_lvl && n.lvl <= max_lvl){
            out.push_back(key);
        }
    }
    return out;
}

void cb_generation(US min_lvl, US max_lvl)
{
    std::vector<ULL> node_hashes = select_lvl(min_lvl, max_lvl);
    std::vector<ULL> new_hashes;
    UL num_nodes = node_hashes.size();
    UL N = num_nodes * (num_nodes - 1) / 2;
    UL i, j, k;
    for (UL k = 0u; k < N; k++)
    {
        i = N - 2 - floor(sqrt(-8*k + 4*N*(N-1)-7)/2.0 - 0.5);
        j = k + i + 1 - N*(N-1)/2 + (N-i)*((N-i)-1)/2;
        std::cout << i << " " << j << std::endl;
    }


    for (ULL i = 0; i < node_hashes.size(); i++)
        Node & ni = GNM[node_hashes[i]];
        // for (ULL j = i+1; j < node_hashes.size(); j++)
}


UI node_cost(const Node& n1, const Node& n2, UI gate_cost)
{
    std::cout << "\t Creating stack " << std::endl;
    std::stack<ULL> stack;
    std::cout << "\t Adding to stack " << std::endl;
    stack.push(n1.hash);
    stack.push(n2.hash);

    std::cout << "\t Creating counter " << std::endl;
    std::unordered_map<ULL, UI> ct_spl;
    std::unordered_set<ULL> non_splittable_nodes;

    while (!stack.empty())
    {
        ULL n_hash = stack.top();
        Node& n = GNM[n_hash];
        std::cout << "\t Got element " << n->func << std::endl;
        ct_spl[n_hash]++;
        std::cout << "\t Updated counter " << std::endl;
        if (ct_spl[n_hash] == 1)
        {
            gate_cost += COSTS[n.last_func];
            std::cout << "\t Added gate cost " << std::endl;
            std::cout << "\t Getting parents " << std::endl;
            for (const auto& parent : n.parents)
            {
                std::cout << "\t Pushing parent " << std::endl;
                std::cout << parent.to_str() << std::endl;
                stack.push(parent.hash);
            }
        }
        else
        {
            gate_cost += COSTS[fSPL];
        }
        if (n.last_func == fAND || n.last_func == fOR)
        {
            for (const auto& parent : n.parents)
            {
                non_splittable_nodes.emplace(parent.hash);
            }
        }
        stack.pop();
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


int main() {

    std::vector<NodeWrapper> pi = gen_pi(NUM_VARS);
    NodeWrapper CONST0 = NodeWrapper(   0, fPI, 0, 0, true);
    NodeWrapper CONST1 = NodeWrapper(ONES, fPI, 0, 0, true);

    // auto i = 0u;
    // for (auto & n : GNM)
    // {   
    //     std::cout << i++ << ": " << n.second.to_str() << std::endl;
    // }
    US lvl = 1;
    auto nodes = select_lvl(lvl - 1, lvl + 1);

    for (UL k = 0u, n = nodes.size(); k < nodes.size() * (nodes.size() - 1) / 2; k++)
    {
        UL i = n - 2 - UL(sqrt(4*n*(n-1) - 8*k - 7)/2.0 - 0.5);
        UL j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
        std::cout << k << " : " << i << " " << j << std::endl;
        ULL   hi = nodes[i];
        ULL   hj = nodes[j];
        Node& ni =  GNM[hi];
        Node& nj =  GNM[hj];
        TT func = ni.func | nj.func;
        bool xorable = func == (ni.func ^ nj.func);
        US last_func = fCB;
        US depth = ni.lvl * 3 + 1;
        // parents = (self, other)

        UI cost = node_cost(nodes[i], other, COSTS[last_func])
    }

    // for (auto i = 0u; i < NUM_TT; i++)
    // {
    //     std::cout << i << " A: " << GEA[i].to_str() << std::endl;
    //     std::cout << i << " X: " << GEX[i].to_str() << std::endl;
    // }
    // std::cout << GNM.size() << std::endl;
    return 0;
}
#include <kitty/static_truth_table.hpp>
#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <functional>
#include <stack>

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
constexpr UL INF = 99999;

#if false

class Node {
public:

    // static std::vector<Node*> empty;
    friend struct std::hash<Node>;

    UI node_cost(const Node& n2, UI gate_cost) const ;

    Node() : func(0), last_func(fNOFUNC), cost(INF), depth(INF), xorable(false), parents{}, hash(0), lvl(0)
    {
        calculateHash();
        calculateLvl();
    }
#pragma region init
    Node(const TT func, const bool xorable) : func(func), last_func(fNOFUNC), cost(INF), depth(INF), xorable(xorable), parents{}, hash(0), lvl(0)
    {
        calculateHash();
        calculateLvl();
    }
    Node(const TT func, const US last_func, const UI cost, const UI depth, const bool xorable, std::vector<std::shared_ptr<Node>> parents)
        : func(func), last_func(last_func), cost(cost), depth(depth), xorable(xorable), parents(parents), hash(0), lvl(0)
    {
        calculateHash();
        calculateLvl();
    }
    Node(const TT func, const US last_func, const UI cost, const UI depth, const bool xorable)
        : func(func), last_func(last_func), cost(cost), depth(depth), xorable(xorable), parents{}, hash(0), lvl(0)
    {
        calculateHash();
        calculateLvl();
    }

    Node& operator=(const Node& other) {
        if (this != &other) {   // check for self-assignment
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

#pragma endregion init

#pragma region get_functions
    ULL get_hash() const {
        return hash;
    }

    US get_lastfunc() const {
        return last_func;
    }

    US get_lvl() const {
        return lvl;
    }

    UI get_cost() const {
        return cost;
    }

    UI get_func() const {
        return func;
    }
    UI get_depth() const {
        return depth;
    }

    std::vector<std::shared_ptr<Node>> get_parents() const {
        return parents;
    }
#pragma endregion get_functions

#pragma region comparison
    bool operator==(const Node& other) const {
        return hash == other.hash;
    }

    bool operator!=(const Node& other) const {
        return hash != other.hash;
    }

    bool operator<(const Node& other) const {
        if ((*this).lvl < other.lvl)
        {
            return true;
        }
        else if ((*this).lvl == other.lvl)
        {
            if ((*this).cost < other.cost)
            {
                return true;
            }
            else if ((*this).cost == other.cost)
            {
                return (*this).func < other.func;
            }
        }
    }
#pragma endregion comparison

#pragma region operators
    Node operator + (const Node& other) const {
        assert( lvl == other.lvl );
        const TT new_func = func | other.func;
        const bool new_xorable = (new_func == (func ^ other.func));
        const UI new_cost = node_cost(other, COSTS[fCB]);
        const UI new_depth = lvl * 3 + 1;
        std::vector<std::shared_ptr<Node>> parents = { std::make_shared<Node>(*this), std::make_shared<Node>(other) };
        return Node(new_func, fCB, new_cost, new_depth, new_xorable, parents);
    }
    Node operator | (const Node& other) const {
        assert( depth == other.depth );
        assert( depth % 3 == 2 );
        const TT new_func = func | other.func;
        const UI new_cost = node_cost(other, COSTS[fOR]);
        std::vector<std::shared_ptr<Node>> parents = { std::make_shared<Node>(*this), std::make_shared<Node>(other) };
        return Node(new_func, fOR, new_cost, depth + 1, true, parents);
    }
    Node operator & (const Node& other) const {
        assert( depth == other.depth );
        assert( depth % 3 == 2 );
        const TT new_func = func & other.func;
        const UI new_cost = node_cost(other, COSTS[fAND]);
        std::vector<std::shared_ptr<Node>> parents = { std::make_shared<Node>(*this), std::make_shared<Node>(other) };
        return Node(new_func, fAND, new_cost, depth + 1, true, parents);
    }
    Node operator ^ (const Node& other) const {
        assert( lvl == other.lvl );
        const TT new_func = func ^ other.func;
        const UI new_cost = node_cost(other, COSTS[fXOR]);
        const UI new_depth = lvl * 3 + 2;
        std::vector<std::shared_ptr<Node>> parents = { std::make_shared<Node>(*this), std::make_shared<Node>(other) };
        return Node(new_func, fXOR, new_cost, new_depth, true, parents);
    }

    Node operator ~ () const {
        const TT new_func = func ^ ONES;
        const UI new_cost = cost + COSTS[fNOT];
        const UI new_depth = lvl * 3 + 2;
        std::vector<std::shared_ptr<Node>> parents = { std::make_shared<Node>(*this)};
        return Node(new_func, fXOR, new_cost, new_depth, true, parents);
    }

    Node operator + () const {
        const UI new_cost = cost + COSTS[fDFF];
        const UI new_depth = lvl * 3 + 2;
        std::vector<std::shared_ptr<Node>> parents = { std::make_shared<Node>(*this)};
        return Node(func, fXOR, new_cost, new_depth, true, parents);
    }
#pragma endregion operators

    std::string to_str()
    {
        return "Func: " + std::to_string(func) + "|Last: " + std::to_string(last_func) + "|Cost: " + std::to_string(cost) + "|Depth: " + std::to_string(depth) + "|X: " + std::to_string(xorable);
    }

private:
    TT func;
    US last_func;
    UI cost;
    UI depth;
    bool xorable;
    std::vector<std::shared_ptr<Node>> parents;
    ULL hash;
    UI lvl;

    void calculateHash() 
    {
        std::size_t seed = std::hash<TT>()(func) ^ std::hash<US>()(last_func) ^ std::hash<UI>()(cost) ^ std::hash<UI>()(depth) ^ std::hash<bool>()(xorable);
        if (!parents.empty()) {
            for (const auto & parent : parents) {
                seed ^= parent->hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
        }
        hash = seed;
    }

    void calculateLvl() 
    {
        lvl = (depth + 1) / 3;
    }
};

template<>
struct std::hash<Node> 
{
    std::size_t operator()(const Node& n) const 
    {
        return n.get_hash();
    }
};

class NodePtr {
public:
    NodePtr(const Node* node) : ptr( node) {}
    NodePtr(const std::shared_ptr<Node> sptr) : ptr(0) {
        auto a = *sptr;
        auto b = &a;
        ptr = b;
    }
    const Node* ptr;
    TT get_func() const
    {
        return (*ptr).get_func();
    }
    US get_lastfunc() const
    {
        return (*ptr).get_lastfunc();
    }
    UI get_cost() const
    {
        return (*ptr).get_cost();
    }
    UI get_depth() const
    {
        return (*ptr).get_depth();
    }
    UI get_hash() const
    {
        return (*ptr).get_hash();
    }
    std::vector<std::shared_ptr<Node>> get_parents() const
    {
        std::cout << "\t\t Entered " << std::endl;
        std::cout << ptr << std::endl;
        auto q = (*ptr);
        std::cout << "\t\t Dereferenced " << std::endl;
        auto parents = q.get_parents();
        std::cout << "\t\t Got parents " << std::endl;
        return parents;
    }
    bool operator == (const NodePtr& other) const
    {
        return get_hash() == other.get_hash();
    }
};

template<>
struct std::hash<NodePtr> 
{
    std::size_t operator()(const NodePtr& nptr) const 
    {
        return nptr.get_hash();
    }
};

UI Node::node_cost(const Node& n2, UI gate_cost) const
{
    // std::stack<std::shared_ptr<Node>> stack;
    // stack.push(std::make_shared(this));
    // stack.push(std::make_shared(&n2));
    std::cout << "\t Creating stack " << std::endl;
    std::stack<NodePtr> stack;
    std::cout << "\t Adding to stack " << std::endl;
    stack.push(NodePtr(this));
    stack.push(NodePtr(&n2));

    std::cout << "\t Creating counter " << std::endl;
    std::unordered_map<NodePtr, UI> ct_spl;
    std::unordered_set<NodePtr> non_splittable_nodes;

    while (!stack.empty())
    {
        auto n = stack.top();
        std::cout << "\t Got element " << n.get_func() << std::endl;
        ct_spl[n]++;
        std::cout << "\t Updated counter " << std::endl;
        if (ct_spl[n] == 1)
        {
            gate_cost += COSTS[n.get_lastfunc()];
            std::cout << "\t Added gate cost " << std::endl;
            auto qqq = n.get_parents();
            std::cout << "\t Getting parents " << std::endl;
            for (const auto& parent : qqq)
            {
                std::cout << "\t Pushing parent " << std::endl;
                std::cout << parent << std::endl;
                auto a = NodePtr(parent);
                stack.push(a);
            }
        }
        else
        {
            gate_cost += COSTS[fSPL];
        }
        if (n.get_lastfunc() == fAND || n.get_lastfunc() == fOR)
        {
            for (const auto& parent : n.get_parents())
            {
                non_splittable_nodes.emplace(parent);
            }
        }
        stack.pop();
    }

    for (const auto & n : non_splittable_nodes)
    {
        auto count = ct_spl[n];
        if (ct_spl[n] > 1)
        {
            auto lastfunc = n.ptr->get_lastfunc();
            gate_cost += COSTS[lastfunc] * (count - 1); // duplicate a gate
            if (lastfunc == fXOR)
            {
                gate_cost += (count - 1) * COSTS[fSPL]; // need to do twice more splittings for duplicating an XOR gate
            }
        }
    }
    return gate_cost;
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

std::vector<Node> gen_pi(UI nvars)
{
    std::vector<Node> pi;
    for (auto func : gen_pi_func(nvars))
    {
        pi.push_back(Node(func, fPI, 0, 0, true) );
    }
    return pi;
}

void update_func(std::vector<Node>& vec, const TT& func, const Node& new_n)
{
    vec.erase( vec.begin() + func );
    vec.insert( vec.begin() + func, new_n );
}

TT coverage(std::vector<Node*> vec)
{
    TT ct = 0;
    for (Node* n : vec)
    {
        ct += (*n).get_cost() < INF;
    }
    return ct;
}


int main()
{
    std::vector<Node> pi = gen_pi(NUM_VARS);
    Node CONST0 = Node(   0, fPI, 0, 0, true);
    Node CONST1 = Node(ONES, fPI, 0, 0, true);

    std::map<UI, std::set<Node>> all_nodes;
    std::vector<Node*> earliest_any;
    earliest_any.reserve(NUM_TT);
    std::vector<Node*> earliest_xor;
    earliest_xor.reserve(NUM_TT);

    for (auto i = 0u; i < NUM_TT; i++)
    {
        auto n_any = Node(i, false);
        auto n_xor = Node(i, true );
        earliest_any.push_back(&n_any);
        earliest_xor.push_back(&n_xor);
    }
    for (auto & n : pi)
    {
        auto func = n.get_func();
        all_nodes[1].emplace(n);
        earliest_any[func] = &n;
        earliest_xor[func] = &n;
    }
    earliest_any[   0] = &CONST0;
    earliest_xor[   0] = &CONST0;
    earliest_any[ONES] = &CONST1;
    earliest_xor[ONES] = &CONST1;

    UI lvl = 0;
    TT initial_num_nodes, final_sum_nodes;
    do
    {
        UI CB_tree_lvl = 0u;
        initial_num_nodes = coverage(earliest_xor);
        // UL old_len = all_nodes[lvl + 1].size();
        // UL new_len;
        bool update_flag; 
        do
        {
            CB_tree_lvl++;
            std::set<Node> nodes = all_nodes[lvl - 1];
            nodes.insert(all_nodes[lvl    ].begin(), all_nodes[lvl    ].end());
            nodes.insert(all_nodes[lvl + 1].begin(), all_nodes[lvl + 1].end());
            auto old_len = all_nodes[lvl + 1].size();
            update_flag = false;

            std::cout << "\t" << all_nodes[lvl + 1].size() << std::endl;
            std::cout << "\t" << nodes.size() << std::endl;

            for (auto i = nodes.begin(); i != nodes.end(); i = std::next(i))
            {
                auto ni = (*i);
                std::cout << "\t\t" << ni.to_str() << std::endl;
                for (auto j = std::next(i); j != nodes.end(); j = std::next(j))
                {
                    auto nj = (*j);
                    std::cout << "\t\t\t" << nj.to_str() << std::endl;
                    auto n = ni + nj;
                    auto result = all_nodes[lvl + 1].emplace(n);
                    update_flag |= result.second;
                    if (result.second)
                    {
                        std::cout << "\t\t\tAdded new node " << n.to_str() << std::endl;
                    }
                    else
                    {
                        std::cout << "\t\t\tNode " << n.to_str() << " already exists" << std::endl;
                    }
                }
            }
        } while (update_flag) ;
        final_sum_nodes = coverage(earliest_xor);
        std::cout << initial_num_nodes << " " << final_sum_nodes << std::endl;
    } while (initial_num_nodes != final_sum_nodes) ;
}

#endif
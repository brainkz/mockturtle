
#include <iostream>
#include <vector>
#include <functional>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <array>
#include <bitset>
#include <unordered_map>
#include <unordered_set>

#include "kitty/static_truth_table.hpp"

typedef unsigned short US;
typedef const unsigned short CUS;
typedef unsigned int UI;
typedef const unsigned int CUI;
typedef unsigned long UL;
typedef const unsigned long CUL;
typedef unsigned long long ULL;

constexpr US NUM_VARS = 4;

typedef kitty::static_truth_table<NUM_VARS> TT;

constexpr CUI INF = 9999;
constexpr CUI HEX_LEN = (1 << NUM_VARS); // 4
constexpr CUI ONES = (1 << (1 << NUM_VARS)) - 1;


const std::unordered_map<std::string, UI> COSTS = { {"MERGE", 8}, {"MERGE3", 11}, {"XOR", 7}, {"NOT", 9}, {"DFF", 7},  {"SPL", 7},  {"CB", 8},  {"AND", 8}, {"OR", 8}, {"PI", 0} };

class Node {
public:
    std::string last_func;

    Node() : hash(0), lvl(0), parents(nullptr) {}

    Node(TT func, std::string last_func, UI cost, UI depth, bool xorable, std::vector<Node>* parents)
        : func(func), last_func(last_func), cost(cost), depth(depth), xorable(xorable), parents(parents), hash(0), lvl(0)
    {
        calculateHash();
        calculateLvl();
    }

    ~Node() {
        if (parents != nullptr) {
            delete parents;
        }
    }

    ULL getHash() const {
        return hash;
    }

    // UI getCost() const {
    //     return COSTS[last_func];
    // }

    std::string getLastfunc() const {
        return last_func;
    }

    std::vector<Node>* getParents() const {
        return parents;
    }
    bool operator==(const Node& other) const {
        return this->hash == other.hash;
    }
    // bool operator|(const Node& other) const {

    //     func = func | other.func;
    //     last_func = "OR";
    //     // cost = node_cost(self, other, COSTS[last_func]);
    //     assert(self.depth == other.depth);
    //     assert(self.depth % 3 == 2);
    //     depth = self.depth + 1
    //     xorable = True
    //     parents = (self, other)
    //     return Node(func, last_func, cost, depth, xorable, parents)

    //     return this->hash == other.hash;
    // }

    UI node_cost(const Node& n2, UI gate_cost)
    {
        std::stack<Node> stack;
        stack.push(*this);
        stack.push(n2);
        std::unordered_map<Node, UI> ct_spl;
        std::unordered_set<Node> non_splittable_nodes;

        while (!stack.empty())
        {
            auto n = stack.top();
            ct_spl[n]++;
            if (ct_spl[n] == 1)
            {
                gate_cost += COSTS[n.last_func];
                for (const auto& parent : *n.getParents())
                {
                    stack.push(parent);
                }
            }
            else
            {
                gate_cost += COSTS["SPL"];
            }
            if (n.getLastfunc() == "AND" || n.getLastfunc() == "OR")
            {
                for (const auto& parent : *n.getParents())
                {
                    non_splittable_nodes.emplace(parent);
                }
            }
        }
    }

private:
    TT func;
    UI cost;
    UI depth;
    bool xorable;
    std::vector<Node>* parents;
    ULL hash;
    UI lvl;

    void calculateHash() {
        std::size_t seed = std::hash<TT>()(func._bits) ^ std::hash<std::string>()(last_func) ^ std::hash<UI>()(cost) ^ std::hash<UI>()(depth) ^ std::hash<bool>()(xorable);
        if (parents != nullptr) {
            for (const auto& parent : *parents) {
                seed ^= parent.hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
        }
        hash = seed;
    }

    void calculateLvl() {
        lvl = 0;
        if (parents != nullptr) {
            for (const auto& parent : *parents) {
                lvl = std::max(lvl, parent.lvl + 1);
            }
        }
    }
};

int main()
{
    return 0;
}


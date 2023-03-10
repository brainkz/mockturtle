
#if false

#pragma region imports
#include <iostream>
#include <fstream> //include the filestreamobject as the header files
#include <tuple>
#include <iterator>
#include <array>
#include <algorithm>
#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <bitset>
// #include <stdio.h>
#include <cstdio>
#include <fmt/format.h>
//  #include <kitty/kitty.hpp>
#pragma endregion imports

#pragma region constants

typedef unsigned short US;
typedef unsigned int   UI;
typedef unsigned long  UL;
typedef unsigned long long ULL;
typedef const unsigned short CUS;
typedef const unsigned int   CUI;
typedef const unsigned long  CUL;
typedef const unsigned long long  CULL;
typedef unsigned int  TT;

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

std::vector<UL> PIs = { 0xaaaa, 0xcccc, 0xf0f0, 0xff00 };


typedef const std::tuple<UI, US, UI, UI, bool, std::vector<ULL>, ULL, UI> Node;
/*
    TT func;
    US last_func;
    UI cost;
    UI depth;
    bool xorable;
    std::vector<std::shared_ptr<Node>> parents;
    ULL hash;
    UI lvl;
*/

struct node {
    TT func;
    US last_func;
    UI cost;
    UI depth;
    bool xorable;
    ULL hash;
    std::vector<ULL> parents;
    UI lvl;
};



// Address book. Keeps hashes to each Node object
std::unordered_map<ULL, node> hash_table;



std::vector<node> parents(node& n)
{
    std::vector<node> out;
    auto hashes = std::get<5>(n);
    for (auto hash : std::get<5>(n))
    {
        out.push_back(hash_table[hash]);
    }
}




UI node_cost(Node& n1, Node& n2, UI gate_cost)
{
    std::stack<Node> stack;
    stack.push(n1);
    stack.push(n2);
    std::unordered_map<Node, UI> ct_spl;
    std::unordered_set<Node> non_splittable_nodes;

    while (!stack.empty())
    {
        auto n = stack.top();
        ct_spl[n]++;
        if (ct_spl[n] == 1)
        {
            gate_cost += COSTS[(*n).last_func];
            for (const auto& parent : (*n).getParents())
            {
                stack.push(parent);
            }
        }
        else
        {
            gate_cost += COSTS[fSPL];
        }
        if (n.get_lastfunc() == fAND || n.get_lastfunc() == fOR)
        {
            for (const auto& parent : n.getParents())
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
            auto lastfunc = (*n).get_lastfunc();
            gate_cost += COSTS[lastfunc] * (count - 1); // duplicate a gate
            if (lastfunc == fXOR)
            {
                gate_cost += (count - 1) * COSTS[fSPL]; // need to do twice more splittings for duplicating an XOR gate
            }
        }
    }
    return gate_cost;
}

int main()
{   
    // Initialize

    // std::cout << kitty::to_binary(dummy) << " = " << kitty::to_hex(dummy) << std::endl;
}

#endif
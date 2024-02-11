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

constexpr unsigned short NUM_VARS = 4u;

constexpr unsigned short fDFF    = 0x0;
constexpr unsigned short fNOT    = 0x1;
constexpr unsigned short fOR     = 0x2;
constexpr unsigned short fAND    = 0x3;
constexpr unsigned short fMAJ    = 0x4;
constexpr unsigned short fXOR    = 0x5;
constexpr unsigned short fPI     = 0x6;
constexpr unsigned short fNOFUNC = 0x7;


// Hash combiner
template <typename T>
static void hash_combine(std::size_t& seed, const T& val) {
    seed ^= std::hash<T>{}(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

uint32_t calculate_hash(uint16_t func, uint8_t last_func, uint16_t cost, uint8_t depth, bool xorable, std::vector<uint32_t> parent_hashes = {})
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
    uint16_t func : 16  = 0u;
    uint8_t  last_func : 3  = fNOFUNC;
    uint16_t cost : 9 = 0b1'1111'1111;
    uint8_t  depth : 3= 0b111;
    bool xorable = false;
    uint32_t hash = 0u;
    std::vector<uint32_t> parent_hashes;

    Node() = default;
    Node(const Node& other) : func(other.func), last_func(other.last_func), cost(other.cost), depth(other.depth), xorable(other.xorable), hash(other.hash), parent_hashes(other.parent_hashes) {}

    Node& operator=(const Node& other) 
    {
        if (this != &other) 
        {
            func = other.func;
            last_func = other.last_func;
            cost = other.cost;
            depth = other.depth;
            xorable = other.xorable;
            hash = other.hash;
            parent_hashes = other.parent_hashes;
        }
        return *this;
    }

    Node(UI _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes(_parent_hashes)
    {
        // Calculate hash based on the hashes of parent_hashes and specified fields
        hash = calculate_hash(func, last_func, cost, depth, xorable, parent_hashes);
    }

    Node(UI _func, US _last_func, UI _cost, UI _depth, bool _xorable)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), parent_hashes{}
    {
        // It is assumed the node has no parents
        hash = calculate_hash(func, last_func, cost, depth, xorable);
    }

    Node(UI _func, US _last_func, UI _cost, UI _depth, bool _xorable, std::vector<ULL> _parent_hashes, ULL _hash)
        : func(_func), last_func(_last_func), cost(_cost), depth(_depth), xorable(_xorable), hash(_hash), parent_hashes(_parent_hashes)
    {
        // Hash is precalculated based on the hashes of parent_hashes and specified fields
    }

    ~Node() { }

    // Equality operator
    bool operator==(const Node& other) const 
    {
        return hash == other.hash;
    }

    ULL get_hash() 
    {
        return hash;
    }



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


int main()
{
    std::locale::global(std::locale("en_US.UTF-8"));
    std::vector<std::vector<uint8_t>> sets_of_levels { 
        {0,0,0,0},
        {0,0,0,1},
        {0,0,0,2}, //x
        {0,0,1,1},
        {0,0,1,2},
        {0,1,1,1},
        {0,0,1,3}, //x
        {0,0,2,2}, //x
        {0,1,1,2},
        {0,1,1,3}, //x
        {0,1,2,2},
        {0,1,2,3},
        {0,2,2,2}, //x
        {0,1,2,4}, //x
        {0,1,3,3}, //x
        {0,2,2,3}, //x
    };

    for (UI func : gen_pi_func(NUM_VARS))
    {
        create_node(func, fPI, 0, levels[i++]*3 + 1, true);
    }

}
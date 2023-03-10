#pragma region imports
#include <iostream>
#include <fstream> //include the filestreamobject as the header files
#include <tuple>
#include <iterator>
#include <array>
#include <algorithm>
#include <vector>
#include <set>
// #include <unordered_map>
#include <unordered_set>
#include <bitset>
// #include <stdio.h>
#include <cstdio>
#include <fmt/format.h>
 #include <kitty/kitty.hpp>
#pragma endregion imports

#pragma region constants

typedef unsigned short US;
typedef unsigned int   UI;
typedef unsigned long  UL;
typedef const unsigned short CUS;
typedef const unsigned int   CUI;
typedef const unsigned long  CUL;

constexpr CUS NVARS = 4u;
constexpr CUI BITLEN = 1 << NVARS;
constexpr CUL NTT = 1 << BITLEN;
constexpr CUI ONES = NTT - 1;
constexpr CUI INF = ONES;

unsigned int tt_levels[NTT];
// std::set<unsigned int> frontier { {0, 0xaa, 0xcc, 0xf0, 1} };
std::set<unsigned int> frontier { {0, 0xaaaa, 0xcccc, 0xf0f0, 0xff00, ONES} };

// std::set<kitty::static_truth_table<NVARS>> frontier {
//     kitty::static_truth_table<NVARS>,
//     kitty::static_truth_table<NVARS>,
//     kitty::static_truth_table<NVARS>,
//     kitty::static_truth_table<NVARS>,
//     kitty::static_truth_table<NVARS>,
//     kitty::static_truth_table<NVARS>,
// };
/*
template<typename TT>
UI p_repr(CUI tt_num, TT& tt)
{
    tt._bits[0] = tt_num;
    return std::get<0>( kitty::exact_p_canonization( tt ) )._bits[0];
}

std::set<unsigned int> filter_p(std::set<unsigned int>& frontier)
{   
    std::set<unsigned int> out;
    kitty::static_truth_table<NVARS> tt_obj;
    for (auto ii = frontier.begin(); ii != frontier.end(); ++ii) {
        out.emplace(p_repr(*ii, tt_obj));
    }
    return out;
}
*/
bool mergers()
{
    US i = 0;
    UI tt;
    bool status;
    do {
        status = false;
        std::cout << "Exploring depth " << ++i << std::endl;
        for (auto ii = frontier.begin(); ii != frontier.end(); ++ii) {
            for (auto jj = std::next(ii); jj != frontier.end(); ++jj) {
                tt = *ii | *jj;
                auto ret = frontier.emplace( tt );
                status |= ret.second;
            }
        }
    } while ( status );
    return status;
}

bool clocked()
{
    US i = 0;
    UI tt;
    bool status;
    std::set<unsigned int> new_tt;
    std::cout << "Exploring depth " << ++i << std::endl;
    for (auto ii = frontier.begin(); ii != frontier.end(); ++ii) {
        new_tt.emplace( ONES ^ *ii );
        for (auto jj = std::next(ii); jj != frontier.end(); ++jj) {
            new_tt.emplace( *ii ^ *jj );
        }
    }
    frontier.insert(new_tt.begin(), new_tt.end());
}

bool SA()
{
    US i = 0;
    UI tt;
    bool status;
    std::set<unsigned int> new_tt;

    std::cout << "Exploring depth " << ++i << std::endl;
    for (auto ii = frontier.begin(); ii != frontier.end(); ++ii) {
        for (auto jj = std::next(ii); jj != frontier.end(); ++jj) {
            new_tt.emplace( *ii & *jj );
        }
    }
    frontier.insert(new_tt.begin(), new_tt.end());
}

template <typename container>
void print(const container& frontier)
{
    auto i = 0u;
    for (auto &element : frontier)
    {
        // std::cout << "[" << i++ << "] Found: " << element << "\n";
        // std::cout << "[" << i++ << "] Found: " << element << "\n";
        if (NVARS == 3)
        {
            fmt::print( "[{0:03d}] Found: {1:08b} = {1:02x} = {1:d}\n", ++i, element);
        }
        else if (NVARS == 4)
        {
            fmt::print( "[{0:04d}] Found: {1:016b} = {1:04x} = {1:d}\n", ++i, element);
        }
    }
}

int main()
{   
    fmt::print("Initial TT-s: {}\n", frontier.size());
    print(filter_p(frontier));
    mergers();
    fmt::print("After initial AA phase: {}\n", frontier.size());
    print(filter_p(frontier));
    clocked();
    fmt::print("After AS phase: {}\n", frontier.size());
    print(filter_p(frontier));
    SA();
    fmt::print("After SA phase: {}\n", frontier.size());
    print(filter_p(frontier));
    mergers();
    fmt::print("After final AA phase: {}\n", frontier.size());
    print(filter_p(frontier));
}
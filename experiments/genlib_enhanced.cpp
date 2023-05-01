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

UI LOG_LVL = 1;

std::vector<std::vector<UI>> binned_indices(const std::vector<UI>& levels)
{
    std::vector<std::vector<UI>> binned_indices;

    for (UI i = 0; i < levels.size(); ++i) {
        // Check if this element is a duplicate
        if (i > 0 && levels[i] == levels[i-1]) {
            // If it is, add the index to the last inner vector
            binned_indices.back().push_back(i);
        } else {
            // If it's not, create a new inner vector
            binned_indices.push_back({i});
        }
    }
    return binned_indices;
}

std::vector<std::vector<UI>> get_permutations(std::vector<UI> elements)
{
    std::vector<std::vector<UI>> permutations;

    // Loop over all permutations and add them to the output vector
    do {
        permutations.push_back(elements);
    } while (std::next_permutation(elements.begin(), elements.end()));

    return permutations;
}

void cart_product(std::vector<std::vector<std::vector<UI>>>& perm_groups, UI perm_idx, std::vector<UI>& base_perm, std::vector<std::vector<UI>>& all_perms)
{
    if (perm_idx < perm_groups.size() - 1) //if the permutation group is not the last one
    {
        for (auto perm : perm_groups[perm_idx])
        {
            std::vector<UI> new_perm = base_perm;
            new_perm.insert(new_perm.end(), perm.begin(), perm.end());
            cart_product(perm_groups, perm_idx + 1, new_perm, all_perms);
        }
    }
    else
    {
        for (auto perm : perm_groups[perm_idx])
        {
            std::vector<UI> new_perm = base_perm;
            new_perm.insert(new_perm.end(), perm.begin(), perm.end());
            all_perms.push_back(new_perm);
        }
    }
}


std::vector<std::vector<UI>> pi_perms(const std::vector<UI> & levels) 
{
    std::vector<std::vector<UI>> idx_groups = binned_indices(levels);
    std::vector<std::vector<std::vector<UI>>> perm_groups(idx_groups.size());
    for (auto i = 0u; i < idx_groups.size(); ++i)
    {
        perm_groups[i] = get_permutations(idx_groups[i]);
    }
    std::vector<std::vector<UI>> all_perms;
    std::vector<UI> base_perm {};
    cart_product(perm_groups, 0, base_perm, all_perms);
    return all_perms;
}

std::vector<kitty::static_truth_table<NUM_VARS>> equivalent_tts(kitty::static_truth_table<NUM_VARS> Base_TT, std::vector<UI> levels)
{
    // TODO : USE KITTY SWAPS TO GENERATE EQUIVALENT TTs
    std::vector<kitty::static_truth_table<NUM_VARS>> eq_tt{ Base_TT };
    ULL old_size;
    do
    {
        old_size = eq_tt.size();
        for (auto & base : eq_tt)
        {
            for (auto i = 0u; i < levels.size(); ++i)
            {
                for (auto j = i + 1; j < levels.size(); ++j)
                {
                    if (levels[i] == levels[j])
                    {
                        auto new_tt = kitty::swap(base, i, j);
                        if (std::find(eq_tt.begin(), eq_tt.end(), new_tt) == eq_tt.end())
                        {
                            eq_tt.push_back(new_tt);
                        }
                    }
                }
            }
        }
    } while (eq_tt.size() > old_size);
    return eq_tt;
}

template <typename TT>
std::tuple<std::vector<uint8_t>, std::vector<UI>, std::vector<std::pair<uint8_t,uint8_t>>> min_base_perms( TT& tt , std::vector<UI> levels)
{
    std::vector<uint8_t> support;
    std::vector<std::pair<uint8_t,uint8_t>> swaps;

    if (LOG_LVL > 2){fmt::print("\t\t\t\t{:016b}\n", tt._bits);}
    
    auto k = 0u;
    for ( auto i = 0u; i < tt.num_vars(); ++i )
    {
    if ( !kitty::has_var( tt, i ) )
    {
        continue;
    }
    if ( k < i )
    {
        kitty::swap_inplace( tt, k, i );
        std::swap(levels[k], levels[i]);
        swaps.emplace_back(k,i);
        if (LOG_LVL > 2){fmt::print("\t\t\t\t{:016b}\n", tt._bits);}
    }
    support.push_back( i );
    ++k;
    }

    return std::make_tuple(support, levels, swaps);
}

union func_lvl 
{
    unsigned data : 32;
    struct
    {
        unsigned func : 16;
        unsigned support_size : 3;
        unsigned l1 : 3;
        unsigned l2 : 3;
        unsigned l3 : 3;
        unsigned l4 : 3;
    };
};

// template <typename TT>
// std::string to_genlib (UI orig_tt, TT & tt, func_lvl & signature)
// {

//     std::string eqn = GNM[orig_tt].genlib_eqn(nodemap, pis);
// }

#include <array>
#include <functional>

template <size_t N>
struct Delay 
{
    std::array<uint8_t, N> delays;
    Delay() = default;

    explicit Delay(const std::vector<uint8_t>& v) 
    {
        if (v.size() != N) 
        {
            throw std::invalid_argument("vector size does not match delay size");
        }
        std::copy(v.begin(), v.end(), delays.begin());
    }
    bool strictly_better_than(const Delay<N>& other) const 
    {
        bool strictly_smaller = false;
        for (size_t i = 0; i < N; ++i) 
        {
            if (delays[i] < other.delays[i]) 
            {
                strictly_smaller = true;
            } else if (delays[i] > other.delays[i]) 
            {
                return false;
            }
        }
        return strictly_smaller;
    }
    bool operator==(const Delay<N>& other) const 
    {
        return delays == other.delays;
    }
};

namespace std {
  template <size_t n>
  struct hash<Delay<n>> 
  {
    size_t operator()(const Delay<n>& d) const 
    {
      size_t h = 0;
      for (size_t i = 0; i < n; ++i) 
      {
        h = h * 31 + std::hash<uint8_t>()(d.delays[i]);
      }
      return h;
    }
  };
}

template <unsigned N>
void update_registry( Delay<N> D, uint32_t func, uint64_t hash,
        std::unordered_map <uint32_t, std::unordered_map <Delay<N>, uint64_t>> & functions,
        std::unordered_map<uint64_t, Node> & nodemap)
{
    bool add = true;
    std::vector<Delay<N>> to_delete;

    auto it = functions.find(func);
    if (it != functions.end()) 
    {
        std::unordered_map<Delay<N>, uint64_t>& inner_map = it->second;
        // Iterate over the entries in the inner map
        for (const auto& [other_D, other_hash] : inner_map) 
        {
            if (D.strictly_better_than(other_D))        // new delay pattern dominates the old delay pattern
            {
                to_delete.push_back(other_D); 
            }
            else if (other_D.strictly_better_than(D))   // old delay pattern dominates the new delay pattern
            {
                to_delete.push_back(D);
                add = false;
            }
            else if (D == other_D)                      // the delay patterns are equal, look at the cost
            {
                Node &       node = nodemap[      hash];
                Node & other_node = nodemap[other_hash];
                add = (node.cost < other_node.cost);
            }
        }
    }
    if (add)
    {
        // fmt::print("{}\n", node.cost);
        functions[func][D] = hash;
    }
    for (auto & bad_D : to_delete)
    {
        functions[func].erase(bad_D);
    }
}

template <unsigned N>
void force_update_registry( Delay<N> D, uint32_t func, uint64_t hash, std::unordered_map <uint32_t, std::unordered_map <Delay<N>, uint64_t>> & functions, std::unordered_map<uint64_t, Node> & nodemap)
{
    bool add = true;
    auto it = functions.find(func);
    if (it != functions.end()) 
    {
        std::unordered_map<Delay<N>, uint64_t>& inner_map = it->second;
        // Iterate over the entries in the inner map
        for (const auto& [other_D, other_hash] : inner_map) 
        {
            Node &       node = nodemap[      hash];
            Node & other_node = nodemap[other_hash];
            if (D.strictly_better_than(other_D))        // new delay pattern dominates the old delay pattern
            {
                // to_delete.push_back(other_D); 
            }
            else if (other_D.strictly_better_than(D))   // old delay pattern dominates the new delay pattern
            {
                // to_delete.push_back(D);
                // only add if the cost is better
                add = (node.cost < other_node.cost);
            }
            else if (D == other_D)                      // the delay patterns are equal, look at the cost
            {
                Node &       node = nodemap[      hash];
                Node & other_node = nodemap[other_hash];
                add = (node.cost < other_node.cost);
            }
        }
    }
    if (add)
    {
        functions[func][D] = hash; // fmt::print("{}\n", node.cost);
    }
}


template <unsigned N>
void filter_registry( std::unordered_map <uint32_t, std::unordered_map <Delay<N>, uint64_t>> & functions, std::unordered_map<uint64_t, Node> & nodemap)
{
    for (auto & [func, inner_map] : functions)
    {
        std::vector<Delay<N>> keys;
        std::vector<Delay<N>> to_delete;
        for (const auto & [D, hash] : inner_map) 
        {
            keys.push_back(D);
        }
        for (auto i = 0u; i < keys.size(); ++i)
        {
            auto Di = keys[i];
            auto hi = inner_map[Di];
            Node & ni = nodemap[hi];
            for (auto j = i + 1; j < keys.size(); ++j)
            {
                auto Dj = keys[j];
                auto hj = inner_map[Dj];
                Node & nj = nodemap[hj];
                if (Di.strictly_better_than(Dj) && ni.cost <= nj.cost)
                {
                    to_delete.push_back(Dj);
                }
                else if (Dj.strictly_better_than(Di) && nj.cost <= ni.cost)
                {
                    to_delete.push_back(Di);
                }
            }
        }
        for (auto D : to_delete)
        {
            inner_map.erase(D);
        }
        // TODO: delete the marked keys
    }
}

int main()
{
    std::ofstream outfile("LIBRARY_2023_05_01.genlib");

    std::vector<std::vector<UI>> sets_of_levels { { 
        {0,0,0,0},
        {0,0,0,1},
        {0,0,0,2},
        {0,0,1,1},
        {0,0,1,2},
        {0,1,1,1},
        {0,1,2,2},
        {0,1,1,3},
        {0,1,1,2},
        {0,1,2,3}
    } };
    std::reverse(sets_of_levels.begin(), sets_of_levels.end());


    std::unordered_map<ULL, Node> GNM_global;
    std::unordered_map <uint32_t, std::unordered_map <Delay<4>, uint64_t>> functions_4;
    std::unordered_map <uint32_t, std::unordered_map <Delay<3>, uint64_t>> functions_3;
    std::unordered_map <uint32_t, std::unordered_map <Delay<2>, uint64_t>> functions_2;
    std::unordered_map <uint32_t, std::unordered_map <Delay<1>, uint64_t>> functions_1;
    
    auto ctr = 0u;
    std::vector<func_lvl> seen_signatures;
    for (const std::vector<UI> & levels : sets_of_levels)
    {
        if (LOG_LVL > 0) {fmt::print("Processing patterns {}\n", fmt::join(levels, " "));}
        const std::vector<UI> PI_funcs {0x5555, 0x3333, 0x0F0F, 0x00FF};

        const std::string file_prefix = fmt::format("/Users/brainkz/Documents/GitHub/mockturtle/build/Golden_20230427/x3_{}_", fmt::join(levels, ""));

        std::unordered_map<ULL, Node> GNM = read_csv_gnm(fmt::format("{}gnm.csv", file_prefix));
        std::array<ULL, NUM_TT> GEX = read_csv_arr(fmt::format("{}gex.csv", file_prefix));

        GNM_global.insert(GNM.begin(), GNM.end());
        // for (auto [hash, node] : GNM)
        // {
        //     GNM_global[hash] = node;
        // }

        // const std::vector<std::vector<UI>> perms = pi_perms(levels);
        // // std::vector<kitty::static_truth_table<NUM_VARS>> skip_tt;
        // std::vector<kitty::static_truth_table<NUM_VARS>> skip_tt;

        for (auto & [hash, node] : GNM)
        {
            if (!node.xorable)
            {
                continue;
            }

            auto [status, support_size, func, DP, pi_map] = process_node(hash, GNM);
            if (node.func == 0xd597 || node.func == 0xd957 || node.func == 0xe557 || node.func == 0xe337 || node.func == 0xb397 ||
                node.func == 0xb937 || node.func == 0xad1f || node.func == 0xcb1f || node.func == 0xc1bf || node.func == 0xa1df ||
                node.func == 0x89f7 || node.func == 0x8f97 || node.func == 0xb397)
            {
                fmt::print("Found {}\n", node.func);
                fmt::print("\tStack: {}\n", node.to_stack(GNM, pi_map));
                fmt::print("\tCharacteristics: {}\n", node.to_str());
                fmt::print("\tDelays: {}\n", fmt::join(DP, ", "));
            }
            if (status)
            {
                if (support_size == 4)
                {
                    force_update_registry<4>( Delay<4>(DP), func, hash, functions_4, GNM_global);
                }
                else if (support_size == 3)
                {
                    force_update_registry<3>( Delay<3>(DP), func, hash, functions_3, GNM_global);
                }
                else if (support_size == 2)
                {
                    force_update_registry<2>( Delay<2>(DP), func, hash, functions_2, GNM_global);
                }
                else if (support_size == 1)
                {
                    force_update_registry<1>( Delay<1>(DP), func, hash, functions_1, GNM_global);
                }
            }
        }
    }

    filter_registry<1>(functions_1, GNM_global);
    filter_registry<2>(functions_2, GNM_global);
    filter_registry<3>(functions_3, GNM_global);
    filter_registry<4>(functions_4, GNM_global);

    auto out_ctr = 0u;
    fmt::print("1-input functions\n");
    for (auto & [func, inner_map] : functions_1)
    {
        fmt::print("{0}: Function: {1:01x}={1:02b}\n", out_ctr++, func); 
        for (auto & [D, hash] : inner_map)
        {
            Node & node = GNM_global[hash];
            fmt::print("{}\n", node.to_stack(GNM_global));
            fmt::print("\tStack: {}\n", node.to_stack(GNM_global));
            fmt::print("\tCost: {}\n", node.cost);
            fmt::print("\tCharacteristics: {}\n", node.to_str());
            fmt::print("\tDelays: {}\n\n", fmt::join(D.delays, ","));
        }
    }
    fmt::print("2-input functions\n");
    for (auto & [func, inner_map] : functions_2)
    {
        fmt::print("{0}: Function: {1:01x}={1:04b}\n", out_ctr++, func); 
        for (auto & [D, hash] : inner_map)
        {
            Node & node = GNM_global[hash];
            fmt::print("{}\n", node.to_stack(GNM_global));
            fmt::print("\tStack: {}\n", node.to_stack(GNM_global));
            fmt::print("\tCost: {}\n", node.cost);
            fmt::print("\tCharacteristics: {}\n", node.to_str());
            fmt::print("\tDelays: {}\n\n", fmt::join(D.delays, ","));
        }
    }
    fmt::print("3-input functions\n");
    for (auto & [func, inner_map] : functions_3)
    {
        fmt::print("{0}: Function: {1:02x}={1:08b}\n", out_ctr++, func); 
        for (auto & [D, hash] : inner_map)
        {
            Node & node = GNM_global[hash];
            fmt::print("{}\n", node.to_stack(GNM_global));
            fmt::print("\tStack: {}\n", node.to_stack(GNM_global));
            fmt::print("\tCost: {}\n", node.cost);
            fmt::print("\tCharacteristics: {}\n", node.to_str());
            fmt::print("\tDelays: {}\n\n", fmt::join(D.delays, ","));
        }
    }
    fmt::print("4-input functions\n");
    for (auto & [func, inner_map] : functions_4)
    {
        fmt::print("{0}: Function: {1:04x}={1:016b}\n", out_ctr++, func); 
        for (auto & [D, hash] : inner_map)
        {
            Node & node = GNM_global[hash];
            fmt::print("{}\n", node.to_stack(GNM_global));
            fmt::print("\tStack: {}\n", node.to_stack(GNM_global));
            fmt::print("\tCost: {}\n", node.cost);
            fmt::print("\tCharacteristics: {}\n", node.to_str());
            fmt::print("\tDelays: {}\n\n", fmt::join(D.delays, ","));
        }
    }
}

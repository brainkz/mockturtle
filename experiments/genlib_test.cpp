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


#if false 
    int main()
    {
        std::string gnm_file = "/Users/brainkz/Documents/GitHub/mockturtle/build/20230320_vec/x3_0123_gnm.csv";
        std::string gex_file = "/Users/brainkz/Documents/GitHub/mockturtle/build/20230320_vec/x3_0123_gex.csv";

        std::unordered_map<ULL, Node> GNM = read_csv_gnm(gnm_file);
        // std::array<ULL, NUM_TT> GEA = read_csv_arr(fmt::format("{}gea.csv", file_prefix));
        std::array<ULL, NUM_TT> GEX = read_csv_arr(gex_file);
        UI i = 0;
        for (auto hash : GEX)
        {
            Node & n = GNM.at(hash);
            // if (LOG_LVL > 2){fmt::print("{:>5d}:\t{}\t{}\n", i++, hash, n.to_str());}
            if (n.lvl == 1) if (LOG_LVL > 2){fmt::print("{:>5d}:\t{}\n", i++, n.to_str());}
        }
    }
#else
    int main()
    {
        std::ofstream outfile("LIBRARY_2023_03_24.genlib");

        std::vector<std::vector<UI>> sets_of_levels { { {0,1,2,3}, {0,0,0,0}, {0,0,0,1}, {0,0,0,2}, {0,0,0,3}, {0,0,0,4}, {0,0,1,1}, {0,0,1,2}, {0,0,1,3}, {0,0,1,4}, {0,0,2,2}, {0,0,2,3}, {0,0,2,4}, {0,0,3,3}, {0,0,3,4}, {0,0,4,4}, {0,1,1,1}, {0,1,1,2}, {0,1,1,3}, {0,1,1,4}, {0,1,2,2}, {0,1,2,3}, {0,1,2,4}, {0,1,3,3}, {0,1,3,4}, {0,1,4,4}, {0,2,2,2}, {0,2,2,3}, {0,2,2,4}, {0,2,3,3}, {0,2,3,4}, {0,2,4,4}, {0,3,3,3}, {0,3,3,4}, {0,3,4,4}, {0,4,4,4} } };
        
        auto ctr = 0u;
        std::vector<func_lvl> seen_signatures;
        for (const std::vector<UI> levels : sets_of_levels)
        {
            if (LOG_LVL > 0) {fmt::print("Processing patterns {}\n", fmt::join(levels, " "));}
            const std::vector<UI> PI_funcs {0x5555, 0x3333, 0x0F0F, 0x00FF};

            std::string file_prefix = fmt::format("/Users/brainkz/Documents/GitHub/mockturtle/build/20230320_vec/x3_{}_", fmt::join(levels, ""));

            std::unordered_map<ULL, Node> GNM = read_csv_gnm(fmt::format("{}gnm.csv", file_prefix));
            // std::array<ULL, NUM_TT> GEA = read_csv_arr(fmt::format("{}gea.csv", file_prefix));
            std::array<ULL, NUM_TT> GEX = read_csv_arr(fmt::format("{}gex.csv", file_prefix));

            // const std::vector<std::vector<UI>> perms = pi_perms(levels);
            // // std::vector<kitty::static_truth_table<NUM_VARS>> skip_tt;
            // std::vector<kitty::static_truth_table<NUM_VARS>> skip_tt;
            for (auto nhash : GEX)
            {
                Node & n = GNM[nhash];
                auto [is_ok, sup_size] = n.redundancy_check(GNM);
                if (LOG_LVL > 2){fmt::print("\tOK: {} | support_size: {}\n\n", is_ok, (is_ok?fmt::format("{}",sup_size):"N/A"));}
                if (!is_ok || sup_size == 0) continue;
                kitty::static_truth_table<NUM_VARS> Base_TT;
                std::vector<UI> words {n.func};
                kitty::create_from_words(Base_TT, words.begin(), words.end());
                if (LOG_LVL > 2){fmt::print("Base_TT: 0x{:04x} \n", Base_TT._bits);}

                // IMPORTANT: here the min_base is determined
                auto [support_idx, perm_levels, swaps] = min_base_perms(Base_TT, levels);
                auto hexlen = (sup_size > 1 ? ( NUM_VARS >> (NUM_VARS - sup_size) ) : 1);
                if (LOG_LVL > 2){fmt::print("Reduced TT: 0x{:0{}x} \n", Base_TT._bits  & ((1 << (1 << sup_size)) - 1), hexlen);}
                if (LOG_LVL > 2){fmt::print("Permuted levels: {} \n", fmt::join(perm_levels, " "));}
                bool is_lifted = false;    
                for (auto idx : support_idx)
                {
                    if (levels[idx] == 0)
                    {
                        is_lifted = true;
                        break;
                    }
                }
                if (!is_lifted)
                {
                    if (LOG_LVL > 2){fmt::print("Permuted levels do not contain 0. Skipping\n");}
                    continue;
                }

                assert(sup_size == support_idx.size());
                func_lvl signature;
                signature.data = 0u;
                signature.func = Base_TT._bits & ((1 << (1 << sup_size)) - 1);
                signature.support_size = sup_size;
                if (sup_size > 0) 
                {
                    signature.l1 = perm_levels[0];
                    if (sup_size > 1)
                    {
                        signature.l2 = perm_levels[1];
                        if (sup_size > 2)
                        {
                            signature.l3 = perm_levels[2];
                            if (sup_size > 3)
                            {
                                signature.l4 = perm_levels[3];
                            }
                        }
                    }
                }

                if (std::find_if(seen_signatures.begin(), seen_signatures.end(), [signature](func_lvl a) {return (a.data == signature.data);}) == seen_signatures.end())
                {
                    auto eq_tt = equivalent_tts(Base_TT, levels);
                    seen_signatures.push_back(signature);

                    if (LOG_LVL > 0) {
                        fmt::print("Writing TT #{:>5d}: {}\n", ++ctr, n.to_str());
                        fmt::print("\t{}\t{}\n", sup_size, fmt::join(perm_levels.begin(), perm_levels.begin() + sup_size, " "));
                        fmt::print("\t{}\n", n.to_stack(GNM));
                    }
                    outfile << fmt::format("{}\n\n", n.to_genlib(GNM, levels, PI_funcs));
                }

                // fmt::format("{:032b}", signature.func);
                // outfile << fmt::format("{}\n\n", n.to_genlib(GNM, levels, PI_funcs));
            }
        }
    }
#endif
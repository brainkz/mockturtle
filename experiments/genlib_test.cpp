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

#define LOG_PRINT(lvl, ...) \
  do { \
      if ((lvl) <= LOG_LVL) { \
          fmt::print(__VA_ARGS__); \
      } \
  } while (false)

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
struct spec
{
  // std::string stack;
  // std::string eqn;
  ULL hash;
  UI func;
  UI cost;
  UI sup_size; 
  std::vector<UI> delays; 
  std::vector<UI> indices; 
  std::vector<std::string> symbol; 
  bool has_better_lvl_than(spec& other) 
  {
    for (std::size_t idx = 0; idx < delays.size(); ++idx) 
    {
      if (delays[idx] > other.delays[idx]) 
      {
        return false;
      }
    }
    return true;
  }
  bool has_equal_level(spec& other) 
  {
    for (std::size_t idx = 0; idx < delays.size(); ++idx) 
    {
      if (delays[idx] != other.delays[idx]) 
      {
        return false;
      }
    }
    return true;
  }
  bool has_not_equal_level(spec& other) 
  {
    for (std::size_t idx = 0; idx < delays.size(); ++idx) 
    {
      if (delays[idx] != other.delays[idx]) 
      {
        return true;
      }
    }
    return false;
  }
};

#if false
  template<typename TT>
  std::tuple<TT, std::vector<UI>, std::vector<uint8_t>> p_canonize( const TT& tt , std::vector<UI> levels)
  {
      auto t1 = tt;
      auto tmin = t1;
      auto tminbase = t1;

      const auto& swaps = kitty::detail::swaps[NUM_VARS - 2u];

      int best_swap = -1;
      std::vector<uint8_t> support;
      for ( std::size_t i = 0; i < swaps.size(); ++i )
      {
          const auto pos = swaps[i];
          swap_adjacent_inplace( t1, pos );
          support = kitty::min_base_inplace(tminbase);

          if ( t1 < tmin && tminbase == t1)
          {
              best_swap = static_cast<int>( i );
              tmin = t1;
          }
      }
      for ( auto i = 0; i <= best_swap; ++i )
      {
          std::swap(levels[swaps[i]], levels[swaps[i] + 1]);
      }

      return std::make_tuple( tmin, levels, support );
  }

  spec canonize( Node & n, std::vector<UI> levels )
  {
      UI func = n.func;
      kitty::static_truth_table<NUM_VARS> t1;
      kitty::create_from_words(t1, &func, &func + 1);

      std::vector<UI> indices(NUM_VARS);
      std::vector<bool> support(NUM_VARS, false);
      std::vector<UI> delays(NUM_VARS, INF);
      LOG_PRINT(0, "levels: ({})\n", fmt::join(levels, ","));
      for ( auto i = 0; i < NUM_VARS; ++i )
      {
        indices[i] = i;
        if (kitty::has_var(t1, i))
        {
          support[i] = true;
          delays[i] = ( n.depth + 2 - (3 * levels[i])) / 3;
        }
        fmt::print("{}: delay: {}\tdepth: {}\t3xlvl: {}\n", i, delays[i], n.depth, (3 * levels[i]));
      }
      LOG_PRINT(0, "Started pattern {0:04x}={0:016b} | support: ({1}) | delay:({2}) | indices:({3})\n", func, fmt::join(support, ","), fmt::join(delays, ","), fmt::join(indices, ","));

      const auto& swaps = kitty::detail::swaps[NUM_VARS - 2u];

      auto tmin = t1;
      int best_swap = -1;
      for ( std::size_t i = 0; i < swaps.size(); ++i )
      {
          const auto pos = swaps[i];
          swap_adjacent_inplace( t1, pos );
          kitty::static_truth_table<NUM_VARS> tminbase;
          kitty::min_base_inplace(tminbase);

          if ( t1 < tmin && tminbase == t1)
          {
              best_swap = static_cast<int>( i );
              tmin = t1;
          }
      }

      // fmt::print("Func: {:016b}\nTmin: {:016b}\n ", func, tmin._bits);

      // const UI true_lvl = (n.depth + 3) / 3;/
      for ( auto i = 0; i <= best_swap; ++i )
      {
          auto idx = swaps[i];
          std::swap(indices[idx], indices[idx + 1]);
          std::swap(delays[idx], delays[idx + 1]);
          std::swap(support[idx], support[idx + 1]);
      }

      spec S;
      S.func = tmin._bits;
      S.hash = n.hash;
      S.support = support;
      S.delays = delays;
      S.indices = indices;
      return S;
  }
  inline bool is_lifted(const std::vector<UI> & levels, const std::vector<uint8_t> & support)
  {
      for (const auto idx : support)
      {
          if (levels[idx] == 0)
          {
              return false;
          }
      }
      return true;
  }

#endif
std::vector<std::string> permute(const std::vector<std::string>& a, const std::vector<UI>& indices) {
    std::vector<std::string> result(a.size());
    for (int i = 0; i < indices.size(); i++) {
        result[i] = a[indices[i]];
    }
    return result;
}

template <typename TT>
void canonize( spec & S , TT & t1)
{
  UI func = S.func;
  kitty::create_from_words(t1, &func, &func + 1);

  // LOG_PRINT(0, "Started  pattern 0x{0:04x}=0b{0:016b} | delay: ({1}) | indices:({2})\n", S.func, fmt::join(S.delays, ","), fmt::join(S.indices, ","));

  const auto& swaps = kitty::detail::swaps[t1.num_vars() - 2u];

  TT tmin = t1;
  int best_swap = -1;
  for ( std::size_t i = 0; i < swaps.size(); ++i )
  {
      const auto pos = swaps[i];
      swap_adjacent_inplace( t1, pos );

      if ( t1 < tmin )
      {
          best_swap = static_cast<int>( i );
          tmin = t1;
      }
  }
  for ( auto i = 0; i <= best_swap; ++i )
  {
      auto idx = swaps[i];
      std::swap(S.indices[idx], S.indices[idx + 1]);
      std::swap(S.delays[idx], S.delays[idx + 1]);
      std::swap(S.symbol[idx], S.symbol[idx + 1]);
  }
  S.func = tmin._bits;
  S.symbol = permute(S.symbol, S.indices);
  // LOG_PRINT(0, "Finished pattern 0x{0:04x}=0b{0:016b} | delay: ({1}) | indices:({2})\n", S.func, fmt::join(S.delays, ","), fmt::join(S.indices, ","));
}

spec min_base_perms(const Node & n, std::vector<UI> levels, std::unordered_map<ULL, Node>& nodemap )
{
    UI func = n.func;
    kitty::static_truth_table<NUM_VARS> tt;
    kitty::create_from_words(tt, &func, &func + 1);

    std::vector<UI> indices(NUM_VARS);
    std::vector<UI> delays(NUM_VARS, INF);
    // std::vector<UI> pi = { 0x5555, 0x3333, 0x0F0F, 0x00FF };
    // std::vector<std::string> symbol = { "a", "b", "c", "d" };
    std::vector<std::string> symbol = { "d", "c", "b", "a" };
    // LOG_PRINT(0, "levels: ({})\n", fmt::join(levels, ","));
    for ( auto i = 0; i < NUM_VARS; ++i )
    {
      indices[i] = i;
      delays[i] = ( n.depth + 2 - (3 * levels[i])) / 3;
      // fmt::print("{}: delay: {}\tdepth: {}\t3xlvl: {}\n", i, delays[i], n.depth, (3 * levels[i]));
    }
    // LOG_PRINT(0, "Started pattern {0:04x}={0:016b} | delay: ({1}) | indices:({2})\n", func, fmt::join(delays, ","), fmt::join(indices, ","));

    auto sup_size = 0u;
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
        std::swap(delays[k], delays[i]);
        std::swap(indices[k], indices[i]);
        std::swap(symbol[k], symbol[i]);
        LOG_PRINT(2, "\t\t\t\t{:016b}\n", tt._bits);
      }
      ++k;
    }
    spec S;
    if (k == 4)
    {
      S.func = func;
      S.hash = n.hash;
      S.cost = n.cost;
      // S.stack = n.to_stack(nodemap, symbol);
      // S.eqn = n.genlib_eqn(nodemap, symbol);
      S.sup_size = k;
      S.delays.insert(S.delays.end(), delays.begin(), delays.begin() + k);
      S.indices.insert(S.indices.end(), indices.begin(), indices.begin() + k);
      S.symbol.insert(S.symbol.end(), symbol.begin(), symbol.begin() + k);
    }
    else
    {
      S.func = tt._bits & ((1 << (1 << k)) - 1);
      S.hash = n.hash;
      S.cost = n.cost;
      // S.stack = n.to_stack(nodemap, symbol);
      // S.eqn = n.genlib_eqn(nodemap, symbol);
      S.sup_size = k;
      S.delays.insert(S.delays.end(), delays.begin(), delays.begin() + k);
      S.indices.insert(S.indices.end(), indices.begin(), indices.begin() + k);
      S.symbol.insert(S.symbol.end(), symbol.begin(), symbol.begin() + k);
    }
    if (S.sup_size == 4)
    {
      kitty::static_truth_table<4> tt;
      canonize(S, tt);
      if (S.indices[0] == 0 && S.indices[1] == 1 && S.indices[2] == 2 && S.indices[3] == 3) 
      {
        eqn_dict[S.func] = nodemap[S.hash].genlib_eqn(nodemap, S.symbol, S.indices);
      }
    }
    else if (S.sup_size == 3)
    {
      kitty::static_truth_table<3> tt;
      canonize(S, tt);
    }
    else if (S.sup_size == 2)
    {
      kitty::static_truth_table<2> tt;
      canonize(S, tt);
    }
    return S;
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

#define mode 1234

#if (mode == 1234)

  int main()
  {
    std::ofstream outfile("LIBRARY_2023_03_28.genlib");
    std::vector<UI> PI = {0x5555, 0x3333, 0x0F0F, 0x00FF};
    std::vector<std::vector<UI>> sets_of_levels { 
      {0,1,2,3},
      {0,0,0,0},
      {0,0,0,1},
      {0,0,0,2},
      {0,0,0,3},
      {0,0,0,4},
      {0,0,1,1},
      {0,0,1,2},
      {0,0,1,3},
      {0,0,1,4},
      {0,0,2,2},
      {0,0,2,3},
      {0,0,2,4},
      {0,0,3,3},
      {0,0,3,4},
      {0,0,4,4},
      {0,1,1,1},
      {0,1,1,2},
      {0,1,1,3},
      {0,1,1,4},
      {0,1,2,2},
      {0,1,2,3},
      {0,1,2,4},
      {0,1,3,3},
      {0,1,3,4},
      {0,1,4,4},
      {0,2,2,2},
      {0,2,2,3},
      {0,2,2,4},
      {0,2,3,3},
      {0,2,3,4},
      {0,2,4,4},
      {0,3,3,3},
      {0,3,3,4},
      {0,3,4,4},
      {0,4,4,4} 
    };

    std::unordered_map<ULL, Node> GNM;
    std::vector<std::array<ULL, NUM_TT>> GEX_global;
    // std::vector<std::vector<func_lvl>> delays_global(NUM_TT, std::vector<func_lvl>());
    // std::unordered_map<UI, std::vector<spec>> best_specs;
    // std::unordered_map<UI, std::vector<std::pair<std::string, spec>>> best_specs_2;
    // std::unordered_map<UI, std::vector<std::pair<std::string, spec>>> best_specs_3;
    // std::unordered_map<UI, std::vector<std::pair<std::string, spec>>> best_specs_4;
    std::vector<spec> best_specs_2;
    std::vector<spec> best_specs_3;
    std::vector<spec> best_specs_4;
    auto lvl_idx = 0u;
    for (std::vector<UI> & levels : sets_of_levels)
    {
      // if (lvl_idx > 3) break;
      std::string file_prefix = fmt::format("/Users/brainkz/Documents/GitHub/mockturtle/build/20230320_vec/x3_{}_", fmt::join(levels, ""));
      auto NM = read_csv_gnm(fmt::format("{}gnm.csv", file_prefix));
      GEX_global.push_back(read_csv_arr(fmt::format("{}gex.csv", file_prefix)));
      UI func = 10;
      ULL nhash = GEX_global.back()[func];
      Node & n = NM[nhash]; 
      // fmt::print(" nhash: {}\nn.hash: {}\n", nhash, n.hash);
      // fmt::print("  func: {}\nn.func: {}\n",  func, n.func);
      // fmt::print("\t{}\n",  n.to_str());
      GNM.merge(NM);
      Node & q = GNM[nhash]; 
      fmt::print("{}: GNM SIZE IS {}\n", fmt::join(levels, ""), GNM.size());
      // fmt::print(" nhash: {}\nq.hash: {}\n", nhash, q.hash);
      // fmt::print("  func: {}\nq.func: {}\n",  func, q.func);
      // fmt::print("\t{}\n",  q.to_str());
      assert(GNM[GEX_global.back()[10]].func == 10);
      lvl_idx++;
    }
    fmt::print("DONE: GNM SIZE IS {}\n", GNM.size());
    lvl_idx = 0u;
    for (std::vector<UI> & levels : sets_of_levels)
    {
      // if (lvl_idx > 3) break;
      LOG_PRINT(0, "Processing patterns {}\n", fmt::join(levels, " "));
      // UI func = 0u;
      for (auto nhash : GEX_global[lvl_idx])
      {
        Node & n = GNM[nhash];
        if (n.func == 0 || n.func == ONES || n.func == 0x5555 || n.func == 0x3333 || n.func == 0x0F0F || n.func == 0x00FF) continue;
        // fmt::print(" nhash: {}\nn.hash: {}\n", nhash, n.hash);
        // fmt::print("  func: {}\nn.func: {}\n",  func, n.func);
        spec S = min_base_perms( n, levels, GNM );
        if (std::find(S.delays.begin(), S.delays.end(), 0) != S.delays.end())
        {
          continue;
        }
        if (S.sup_size == 4)
        {
          best_specs_4.emplace_back( S );
        }
        assert(S.hash == n.hash && nhash == n.hash);
        // best_specs[ S.func ].push_back(S);
        // LOG_PRINT(0, "Emplaced pattern {0:04x}={0:016b} | support: ({1}) | delay:({2}) | indices:({3})\n", S.func, S.sup_size, fmt::join(S.delays, ","), fmt::join(S.indices, ","));
        // fmt::print("{}\n\n", n.to_str());
        // fmt::print("{}\n\n", n.to_genlib(GNM, levels, PI));
      }
      lvl_idx++;
    }

    auto spec_ctr = 0u;
    for (UI func = 0u; func < NUM_TT; ++func )
    {
      std::vector<spec> specs;
      for (auto & S : best_specs_4)
      {
        if (S.func == func)
        {
          specs.push_back(S);
        }
      }
      if (specs.empty()) continue;
      spec_ctr++;
      std::unordered_set<UI> to_remove;
      for (int i = 0u; i < specs.size(); ++i)
      {
        spec Si = specs[i];
        for (int j = i + 1; j < specs.size(); ++j)
        {
          spec Sj = specs[j];
          if (!Si.has_equal_level(Sj))
          {
            if (Si.has_better_lvl_than(Sj))
            {   
              to_remove.emplace(j);
            }
            if (Sj.has_better_lvl_than(Si))
            {   
              to_remove.emplace(i);
            }
          }
          else //if the levels are equal
          {
            if (Si.cost < Sj.cost)
            {
              to_remove.emplace(j);
            }
            if (Sj.cost < Si.cost)
            {
              to_remove.emplace(i);
            }
          }
        }
      }
      LOG_PRINT(0, "Removing {} | {} out of {}\n", fmt::join(to_remove, ","), to_remove.size(), specs.size());

      for (auto i = 0u; i < specs.size(); ++i)
      {
        // std::array<ULL, NUM_TT> GEX = GEX_global[i];
        if (std::find(to_remove.begin(), to_remove.end(), i) == to_remove.end())
        {
          // std::unordered_map<ULL, Node> GNM = read_csv_gnm(fmt::format("/Users/brainkz/Documents/GitHub/mockturtle/build/20230320_vec/x3_{}_", fmt::join(levels, "")));
          spec S = specs[i];
          Node & n = GNM[S.hash];
          assert(n.hash == S.hash);
          LOG_PRINT(0, "{4:>4}: Writing pattern {0:04x}={0:016b} | support: ({1}) | delay:({2}) | indices:({3})\n", S.func, S.sup_size, fmt::join(S.delays, ","), fmt::join(S.indices, ","), spec_ctr);
          fmt::print("{}\n", n.to_stack(GNM, S.symbol));
          fmt::print("{}\n", n.to_str());
          // std::unordered_map<UI, std::string> pi_map;
          // std::vector<std::string> pi_str {"a", "b", "c", "d"};
          // for (auto i = 0u; i < NUM_VARS; ++i)
          // {
          //   pi_map[i] = pi_str[S.indices[i]];
          // }
          fmt::print("GATE 0x{:04x}_{} {} O={};\n", func, fmt::join(S.delays, ""), S.cost, eqn_dict[S.func]);
          outfile << fmt::format("GATE 0x{:04x}_{} {} O={};\n", func, fmt::join(S.delays, ""), S.cost, eqn_dict[S.func]);
          std::vector<std::string> _symbol = { "d", "c", "b", "a" };
          for (auto i=0u; i < S.delays.size(); ++i)
          {
            outfile << fmt::format("\tPIN {} {} {} {} {:d} {:0.3f} {:d} {:0.3f}\n", 
            _symbol[i], GENLIB_PHASE, GENLIB_INPUT_LOAD, GENLIB_MAX_LOAD, 
            S.delays[i], GENLIB_RISE_FANOUT_DELAY, S.delays[i], GENLIB_FALL_FANOUT_DELAY);
          }
        }
      }
    }
  }

#elif (mode == 999)
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
#elif (mode == 923)
    int main()
    {
        std::ofstream outfile("LIBRARY_2023_03_24.genlib");

        std::vector<std::vector<UI>> sets_of_levels { { {0,1,2,3}, {0,0,0,0}, {0,0,0,1}, {0,0,0,2}, {0,0,0,3}, {0,0,0,4}, {0,0,1,1}, {0,0,1,2}, {0,0,1,3}, {0,0,1,4}, {0,0,2,2}, {0,0,2,3}, {0,0,2,4}, {0,0,3,3}, {0,0,3,4}, {0,0,4,4}, {0,1,1,1}, {0,1,1,2}, {0,1,1,3}, {0,1,1,4}, {0,1,2,2}, {0,1,2,3}, {0,1,2,4}, {0,1,3,3}, {0,1,3,4}, {0,1,4,4}, {0,2,2,2}, {0,2,2,3}, {0,2,2,4}, {0,2,3,3}, {0,2,3,4}, {0,2,4,4}, {0,3,3,3}, {0,3,3,4}, {0,3,4,4}, {0,4,4,4} } };
        
        auto ctr = 0u;
        std::array<UI, NUM_TT> ctr_arr;
        std::fill(std::begin(ctr_arr), std::end(ctr_arr), 0);
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
                const auto [min_TT, perm_levels, support_idx] = p_canonize(Base_TT, levels);
                // auto [support_idx, perm_levels] = p_canonize(Base_TT, levels);
                auto hexlen = (sup_size > 1 ? ( NUM_VARS >> (NUM_VARS - sup_size) ) : 1);
                if (LOG_LVL > 2){fmt::print("Reduced TT: 0x{:0{}x} \n", min_TT._bits  & ((1 << (1 << sup_size)) - 1), hexlen);}
                if (LOG_LVL > 2){fmt::print("Permuted levels: {} \n", fmt::join(perm_levels, " "));}
                
                if (is_lifted(perm_levels, support_idx))
                {
                    if (LOG_LVL > 2){fmt::print("Permuted levels do not contain 0. Skipping\n");}
                    continue;
                }

                assert(sup_size == support_idx.size());
                func_lvl signature;
                signature.data = 0u;
                signature.func = min_TT._bits & ((1 << (1 << sup_size)) - 1);
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

                    if (LOG_LVL > 0) 
                    {
                        // fmt::print("Writing TT #{:>5d}: {}\n", ++ctr, n.to_str());
                        UI data = signature.data;
                        fmt::print("Writing TT #{:>5d} (#{:d}): {}, {:032b}\n", ++ctr, ++ctr_arr[min_TT._bits], min_TT._bits, data); // n.to_str()
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
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

constexpr uint8_t LOG_LVL = 1;

#define continue_if(condition) if (condition) continue
#define break_if(condition) if (condition) break
#define LOG_PRINT(lvl, ...) \
  do { \
      if ((lvl) <= LOG_LVL) { \
          fmt::print(__VA_ARGS__); \
      } \
  } while (false)

std::unordered_map<uint32_t, std::string> eqn_dict;

struct spec
{
  uint64_t hash;
  TT func;
  uint32_t cost;
  uint8_t sup_size; 
  std::vector<uint8_t> delays; 
  std::vector<uint8_t> indices; 
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

std::vector<std::string> permute(const std::vector<std::string>& a, const std::vector<uint8_t>& indices) {
    std::vector<std::string> result(a.size());
    for (uint64_t i = 0; i < indices.size(); i++) {
        result[i] = a[indices[i]];
    }
    return result;
}

void canonize( spec & S )
{
  TT t1 = S.func;
  TT tmin = t1;
  const auto& swaps = kitty::detail::swaps[t1.num_vars() - 2u];

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
  S.func = tmin;
}

spec min_base_perms(const Node & n, std::vector<uint8_t> levels, std::unordered_map<uint64_t, Node>& nodemap )
{
    TT tt = n.func;

    std::vector<uint8_t> indices(NUM_VARS);
    std::vector<uint8_t> delays(NUM_VARS, INF8);
    std::vector<std::string> symbol = { "d", "c", "b", "a" };
    for ( auto i = 0; i < NUM_VARS; ++i )
    {
      indices[i] = i;
      delays[i] = ( n.depth + 2 - (3 * levels[i])) / 3;
    }

    uint8_t k = 0u;
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
      S.func = tt;
      S.hash = n.hash;
      S.cost = n.cost;
      S.sup_size = k;
      S.delays.insert(S.delays.end(), delays.begin(), delays.begin() + k);
      S.indices.insert(S.indices.end(), indices.begin(), indices.begin() + k);
      S.symbol.insert(S.symbol.end(), symbol.begin(), symbol.begin() + k);
    }
    else
    {
      // TODO: adapt for 3-input TT
      S.func = tt; //  & ((1 << (1 << k)) - 1)
      S.hash = n.hash;
      S.cost = n.cost;
      S.sup_size = k;
      S.delays.insert(S.delays.end(), delays.begin(), delays.begin() + k);
      S.indices.insert(S.indices.end(), indices.begin(), indices.begin() + k);
      S.symbol.insert(S.symbol.end(), symbol.begin(), symbol.begin() + k);
    }
    canonize(S);
    if (S.indices[0] == 0 && S.indices[1] == 1 && S.indices[2] == 2 && S.indices[3] == 3) 
    {
      eqn_dict[S.func._bits] = nodemap[S.hash].genlib_eqn(nodemap, S.symbol, S.indices);
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

int main()
{
  std::ofstream outfile("LIBRARY_2023_03_28.genlib");
  std::vector<std::vector<uint8_t>> sets_of_levels { 
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

  std::unordered_map<uint64_t, Node> GNM;
  std::vector<std::array<uint64_t, NUM_TT>> GEX_global;
  std::vector<spec> best_specs;
  auto lvl_idx = 0u;
  for (std::vector<uint8_t> & levels : sets_of_levels)
  {
    std::string file_prefix = fmt::format("/Users/brainkz/Documents/GitHub/mockturtle/build/20230320_vec/x3_{}_", fmt::join(levels, ""));
    auto NM = read_csv_gnm(fmt::format("{}gnm.csv", file_prefix));
    GEX_global.push_back(read_csv_arr(fmt::format("{}gex.csv", file_prefix)));
    /*
      TT func = num2tt(10);
      uint64_t nhash = GEX_global.back()[func._bits];
      Node & n = NM[nhash]; 
      fmt::print(" nhash: {}\nn.hash: {}\n", nhash, n.hash);
      fmt::print("  func: {}\nn.func._bits: {}\n",  func, n.func._bits);
      fmt::print("\t{}\n",  n.to_str());*/
    GNM.merge(NM);
    fmt::print("{}: GNM SIZE IS {}\n", fmt::join(levels, ""), GNM.size());
    /* 
      fmt::print(" nhash: {}\nq.hash: {}\n", nhash, q.hash);
      fmt::print("  func: {}\nq.func._bits: {}\n",  func, q.func._bits);
      fmt::print("\t{}\n",  q.to_str());*/
    assert(GNM[GEX_global.back()[10]].func._bits == 10);
    lvl_idx++;
  }
  fmt::print("DONE: GNM SIZE IS {}\n", GNM.size());
  lvl_idx = 0u;
  for (std::vector<uint8_t> & levels : sets_of_levels)
  {
    LOG_PRINT(0, "Processing patterns {}\n", fmt::join(levels, " "));
    for (auto nhash : GEX_global[lvl_idx])
    {
      Node & n = GNM[nhash];
      continue_if (is_pi(n.func));
      spec S = min_base_perms( n, levels, GNM );
      continue_if (std::find(S.delays.begin(), S.delays.end(), 0) != S.delays.end());
      best_specs.push_back( S );
      assert(S.hash == n.hash && nhash == n.hash);
    }
    lvl_idx++;
  }

  auto spec_ctr = 0u;
  for (uint32_t func_idx = 0u; func_idx < NUM_TT; ++func_idx )
  {
    std::vector<spec> specs;
    for (spec & S : best_specs)
    {
      if (S.func._bits == func_idx)
      {
        specs.push_back(S);
      }
    }
    if (specs.empty()) continue;
    spec_ctr++;
    std::unordered_set<uint32_t> to_remove;
    for (uint64_t i = 0u; i < specs.size(); ++i)
    {
      spec Si = specs[i];
      for (uint64_t j = i + 1; j < specs.size(); ++j)
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

    for (uint64_t i = 0u; i < specs.size(); ++i)
    {
      if (std::find(to_remove.begin(), to_remove.end(), i) == to_remove.end())
      {
        spec S = specs[i];
        Node & n = GNM[S.hash];
        assert(n.hash == S.hash);
        LOG_PRINT(0, "{4:>4}: Writing pattern {0:04x}={0:016b} | support: ({1}) | delay:({2}) | indices:({3})\n", S.func._bits, S.sup_size, fmt::join(S.delays, ","), fmt::join(S.indices, ","), spec_ctr);
        fmt::print("{}\n", n.to_stack(GNM, S.symbol));
        fmt::print("{}\n", n.to_str());
        fmt::print("GATE 0x{:04x}_{} {} O={};\n", S.func._bits, fmt::join(S.delays, ""), S.cost, eqn_dict[S.func._bits]);
        outfile << fmt::format("GATE 0x{:04x}_{} {} O={};\n", S.func._bits, fmt::join(S.delays, ""), S.cost, eqn_dict[S.func._bits]);
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
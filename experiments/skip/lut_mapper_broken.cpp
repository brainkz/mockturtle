/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2019  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include <kitty/kitty.hpp>

#include <fmt/format.h>
#include <lorina/aiger.hpp>

#include <mockturtle/generators/arithmetic.hpp>

#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/io/write_bench.hpp>
// #include <mockturtle/io/write_blif.hpp>

// #include <mockturtle/utils/four_input_tt_maj.hpp>
// #include <mockturtle/utils/four_input_tt_maj_and_or.hpp>
// #include <mockturtle/utils/arrays_with_inv.hpp>
// #include <mockturtle/utils/arrays_OR_only_xorables.hpp>
// #include <mockturtle/utils/tt4_with_xorables.hpp> //ONLY mergers 
// #include <mockturtle/utils/tt4_first_2_3_with_xorables.hpp>
#include <mockturtle/utils/tt4_first_2_3_with_xorables.hpp> //AND3/OR3/MAJ3 + mergers
const std::string suffix = "first_2_3_with_xorables";
const unsigned short CLA_NBITS = 8;
// #include <mockturtle/utils/tt4_first_2_with_xorables.hpp> //AND2/OR2 + mergers


// #include <mockturtle/utils/arrays.hpp>
#include <experiments.hpp>

/*
struct lut_sfq_cost
{
  std::pair<uint32_t, uint32_t> operator()( uint32_t num_leaves ) const
  {
    if ( num_leaves < 2u )
      return { 0u, 0u };
    return { 1u, 1u }; // area, delay
  }

  std::pair<uint64_t, uint64_t> operator()( kitty::dynamic_truth_table const& tt ) const
  {
    
    uint64_t area = area_delay[tt._bits[0]][0];
    uint64_t delay = area_delay[tt._bits[0]][1];
    return { area , delay };
    // std::cout << "CALLED " << delay << std::endl;.
    // if ( delay > 1u)
    // {
    //   return { UINT16_MAX, UINT16_MAX };
    // }
    // if ( tt.num_vars() < 2u )
    //   return { 0u, 0u };
    // return { 1u, 1u }; // area, delay
  }
};
*/

template<typename Ntk, typename AdderFn>
Ntk create_adder( uint32_t width, AdderFn&& adder )
{
  Ntk ntk;

  std::vector<typename Ntk::signal> a( width ), b( width );
  std::generate( a.begin(), a.end(), [&ntk]() { return ntk.create_pi(); } );
  std::generate( b.begin(), b.end(), [&ntk]() { return ntk.create_pi(); } );
  auto carry = ntk.get_constant( false );

  adder( ntk, a, b, carry );

  std::for_each( a.begin(), a.end(), [&]( auto f ) { ntk.create_po( f ); } );
  ntk.create_po( carry );

  assert( ntk.num_pis() == 2 * width );
  assert( ntk.num_pos() == width + 1 );

  return ntk;
}


mockturtle::klut_network run_mapping(mockturtle::aig_network& aig, const mockturtle::lut_map_params& ps, mockturtle::lut_map_stats& st, const std::string benchmark, experiments::experiment<std::string, uint32_t, uint32_t, uint32_t, double, bool>& exp, const bool cec_test = false)
{
  mockturtle::mapping_view<mockturtle::aig_network, true> mapped_aig{ aig };
  mockturtle::lut_map<decltype( mapped_aig ), true, lut_sfq_cost>( mapped_aig, ps, &st );
  const auto klut = *mockturtle::collapse_mapped_network<mockturtle::klut_network>( mapped_aig );

  mockturtle::depth_view<mockturtle::klut_network> klut_d{ klut };

  auto const cec = (cec_test) ? ((benchmark == "hyp") ? true : experiments::abc_cec( klut, benchmark )) : true;
  fmt::print( "Check 1 \n" );
  exp( benchmark, klut.num_gates(), klut_d.depth(), st.edges, mockturtle::to_seconds( st.time_total ), cec );
  fmt::print( "Check 2 \n" );
  return klut_d;
}

mockturtle::klut_network analyze_adder(experiments::experiment<std::string, uint32_t, uint32_t, uint32_t, double, bool>& exp, const mockturtle::lut_map_params& ps, mockturtle::lut_map_stats& st, const std::string name, const std::string suffix, const bool write)
{
    // Routine to create CLA 
    const unsigned short NBITS = 8;
    fmt::print( "[i] processing {}-bit CLA\n", NBITS );
    mockturtle::aig_network aig = create_adder<mockturtle::aig_network>( NBITS, mockturtle::carry_ripple_adder_inplace<mockturtle::aig_network> );

    auto klut_d = run_mapping(aig, ps, st, name, exp, true);

    if (write)
    {
      mockturtle::write_bench(klut_d, name + "_" + suffix + ".bench");
    }
    return klut_d;
}

mockturtle::klut_network analyze_benchmark(experiments::experiment<std::string, uint32_t, uint32_t, uint32_t, double, bool>& exp, mockturtle::aig_network& aig, const mockturtle::lut_map_params& ps, mockturtle::lut_map_stats& st, const std::string& benchmark, const std::string suffix, const bool write )
{
    fmt::print( "[i] processing {}\n", benchmark );

    auto klut_d = run_mapping(aig, ps, st, benchmark, exp, true);
    fmt::print( "Check 3 \n" );

    if (write)
    {
      mockturtle::write_bench(klut_d, benchmark + "_" + suffix + ".bench");
    }
    return klut_d;
}

template<typename TT>
std::vector<int> redundant_variables(TT tt)
{
  std::vector<int> redundant_vars;
  for ( auto i = 0u; i < tt.num_vars(); ++i )
  {
    auto copy = kitty::flip(tt, i);
    if (copy == tt)
    {
      redundant_vars.push_back(i);
    }
  }
  return redundant_vars;
}

void simplify_benchmark(mockturtle::klut_network& klut, const uint32_t area_delay[NTT][12])
{
  // TODO: For each complex node within network:
  // TODO:  if so, replace the node with its preds + func
  // 1. create pred and func nodes
  // 2. replace old node with func node.
  return;
}

int main()
{
  using namespace experiments;
  using namespace mockturtle;
  // auto qwer = kitty::create_nth_var();
  experiment<std::string, uint32_t, uint32_t, uint32_t, double, bool> exp( "lut_mapper", "benchmark", "luts", "lut_depth", "edges", "runtime", "equivalent" );

  lut_map_params ps;
  ps.cut_enumeration_ps.cut_size = 4u;
  ps.cut_enumeration_ps.cut_limit = 10u;
  ps.recompute_cuts = false;
  ps.area_oriented_mapping = true; //false;
  ps.cut_expansion = true;
  ps.edge_optimization = true; //can try set to true
  ps.remove_dominated_cuts = false; //can try set to true
  ps.cost_cache_vars = 0u; 
  ps.verbose = true;
  lut_map_stats st;

  const bool write = false;

  #pragma region CLA8bit
  fmt::print( "[i] processing {}-bit CLA\n", CLA_NBITS );
  aig_network create_adder<aig_network>( CLA_NBITS, carry_ripple_adder_inplace<aig_network> );

  mapping_view<aig_network, false> mapped_aig{ aig };
  lut_map<decltype( mapped_aig ), false>( mapped_aig, ps, &st );
  const auto klut = *collapse_mapped_network<klut_network>( mapped_aig );

  depth_view<klut_network> klut_d{ klut };
  exp( "CLA", klut.num_gates(), klut_d.depth(), st.edges, to_seconds( st.time_total ), true );
  #pragma endregion CLA8bit

  if (write)
  {
    mockturtle::write_bench(klut_d, name + "_" + suffix + ".bench");
  }

  #pragma region epfl
  for ( auto const& benchmark : epfl_benchmarks() )
    {
      fmt::print( "[i] processing {}\n", benchmark );
      aig_network aig;
      if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
      {
        continue;
      }

    depth_view<klut_network> klut_d{ klut };

    auto const cec = benchmark == "hyp" ? true : abc_cec( klut, benchmark );

    exp( benchmark, klut.num_gates(), klut_d.depth(), st.edges, to_seconds( st.time_total ), cec );
    #pragma endregion epfl
  }
  // fmt::print( "[i] processing \n" );
  
  exp.save();
  exp.table();

  return 0;
}

/*
int main_unused()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, double, bool> exp( "lut_mapper", "benchmark", "luts", "lut_depth", "edges", "runtime", "equivalent" );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    lut_map_params ps;
    ps.cut_enumeration_ps.cut_size = 6u;
    ps.cut_enumeration_ps.cut_limit = 8u;
    ps.recompute_cuts = true;
    ps.area_oriented_mapping = false;
    ps.cut_expansion = true;
    lut_map_stats st;
    mapping_view<aig_network, false> mapped_aig{ aig };
    lut_map<decltype( mapped_aig ), false>( mapped_aig, ps, &st );
    const auto klut = *collapse_mapped_network<klut_network>( mapped_aig );

    depth_view<klut_network> klut_d{ klut };

    auto const cec = benchmark == "hyp" ? true : abc_cec( klut, benchmark );

    exp( benchmark, klut.num_gates(), klut_d.depth(), st.edges, to_seconds( st.time_total ), cec );
  }

  exp.save();
  exp.table();

  return 0;
}
*/
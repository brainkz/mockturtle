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

#include <fmt/format.h>
#include <lorina/aiger.hpp>

#include <mockturtle/generators/arithmetic.hpp>

#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
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

template<typename num>
num maj(num a, num b, num c) 
{
  return (a & (b | c)) | (b & c);
}

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

void analyze_adder(experiments::experiment<std::string, uint32_t, uint32_t, uint32_t, double, bool>& exp, const mockturtle::lut_map_params& ps, mockturtle::lut_map_stats& st, const std::string name, const std::string suffix, const bool write)
{
    // Routine to create CLA 
    const unsigned short NBITS = 8;
    fmt::print( "[i] processing {}-bit CLA\n", NBITS );
    mockturtle::aig_network aig = create_adder<mockturtle::aig_network>( NBITS, mockturtle::carry_ripple_adder_inplace<mockturtle::aig_network> );

    mockturtle::mapping_view<mockturtle::aig_network, true> mapped_aig{ aig };

    mockturtle::lut_map<decltype( mapped_aig ), true, lut_sfq_cost>( mapped_aig, ps, &st );
    const auto klut = *mockturtle::collapse_mapped_network<mockturtle::klut_network>( mapped_aig );

    mockturtle::depth_view<mockturtle::klut_network> klut_d{ klut };

    exp( name, klut.num_gates(), klut_d.depth(), st.edges, mockturtle::to_seconds( st.time_total ), true );

    // mockturtle::write_bench(klut_d, benchmark + "_OR_only_xorables.bench");
    if (write)
    {
      mockturtle::write_bench(klut_d, name + "_" + suffix + ".bench");
    }
}

void analyze_benchmark(experiments::experiment<std::string, uint32_t, uint32_t, uint32_t, double, bool>& exp, const mockturtle::lut_map_params& ps, mockturtle::lut_map_stats& st, const std::string& benchmark, const std::string suffix, const bool write )
{
    fmt::print( "[i] processing {}\n", benchmark );
    mockturtle::aig_network aig;
    if ( lorina::read_aiger( experiments::benchmark_path( benchmark ), mockturtle::aiger_reader( aig ) ) != lorina::return_code::success )
    {
      return;
    }

    mockturtle::mapping_view<mockturtle::aig_network, true> mapped_aig{ aig };

    mockturtle::lut_map<decltype( mapped_aig ), true, lut_sfq_cost>( mapped_aig, ps, &st );
    const auto klut = *mockturtle::collapse_mapped_network<mockturtle::klut_network>( mapped_aig );

    mockturtle::depth_view<mockturtle::klut_network> klut_d{ klut };

    auto const cec = (benchmark == "hyp") ? true : experiments::abc_cec( klut, benchmark );

    exp( benchmark, klut.num_gates(), klut_d.depth(), st.edges, mockturtle::to_seconds( st.time_total ), cec );

    // mockturtle::write_bench(klut_d, benchmark + "_OR_only_xorables.bench");
    if (write)
    {
      mockturtle::write_bench(klut_d, benchmark + "_" + suffix + ".bench");
    }
}



auto func( mockturtle::klut_network& klut, mockturtle::klut_network::signal node)
// uint32_t& area_delay[NTT][12],
// mockturtle::klut_network::signal n 
{ 
  fmt::print( "[debug] Processing node {}\n", node);
  auto ct = 0u;
  auto tt = klut.node_function(node);
  fmt::print( "[debug] Function : {}\n", kitty::to_hex(tt));
  auto num_in = klut.fanin_size(node);
  fmt::print( "[debug] num_in : {}\n", num_in);
  uint32_t tt_idx = tt._bits[0];
  // Determine tt entry
  // fmt::print( "[debug] num_in : {}\n", num_in);
  // fmt::print( "[debug] Initial tt: {}\n", kitty::to_binary(tt));
  // fmt::print( "[debug] Initial tt_idx: {:016b}\n", tt_idx);

  for (auto width = (1 << num_in); width < 16; width *= 2) // convert to 4-LUT if fewer variables
  {
    tt_idx = (tt_idx << width) + tt_idx;
    // fmt::print( "[debug] new tt_idx: {:016b}\n", tt_idx);
  }


  // The redundant variables are set to the latter variables (lower frequency)

  auto p1 = area_delay[tt_idx][3];
  auto p2 = area_delay[tt_idx][4];
  auto p3 = area_delay[tt_idx][5];
    
  mockturtle::klut_network::signal p;
  kitty::dynamic_truth_table tt1, tt2, tt3;

  tt1 = kitty::dynamic_truth_table( num_in );
  tt1._bits[0] = p1;
  // auto n1 = klut.create_node(node.children, tt1);
  std::vector<mockturtle::klut_network::signal> children;
  klut.foreach_fanin(
    node, [&] (auto n) { children.push_back(n); }
  );


  fmt::print( "[debug] the preds of {} are 0x{:04x}, 0x{:04x}, 0x{:04x}\n", kitty::to_hex(tt), area_delay[tt_idx][3], area_delay[tt_idx][4], area_delay[tt_idx][5]);
  if (p2 == INF) // 1-input function, should be an inverter
  {
    assert(tt == ~tt1);
    fmt::print( "[debug] detected inverter\n");
    fmt::print( "[debug] creating node for {}\n", kitty::to_hex(tt1));
    auto n1 = klut.create_node(children, tt1);
    p = klut.create_not( n1 ); 
  }
  else if (p3 == INF) // 2-input function
  {
    fmt::print( "[debug] detected binary function\n");
    tt2 = kitty::dynamic_truth_table( num_in );
    tt2._bits[0] = p2;
    fmt::print( "[debug] creating node for {}\n", kitty::to_hex(tt1));
    auto n1 = klut.create_node(children, tt1);
    fmt::print( "[debug] creating node for {}\n", kitty::to_hex(tt2));
    auto n2 = klut.create_node(children, tt2);
    if      (tt == (tt1 | tt2)) 
    {
      fmt::print( "[debug] detected OR2\n");
      p = klut.create_or( n1, n2 ); 
    }
    else if (tt == (tt1 & tt2)) 
    { 
      fmt::print( "[debug] detected AND2\n");
      p = klut.create_and( n1, n2 ); 
    }
    else if (tt == (tt1 ^ tt2)) 
    {
      fmt::print( "[debug] detected XOR2\n");
      p = klut.create_xor( n1, n2 ); 
    }
    else
    {
      throw("The 2-input TT could not be recognized from predecessors");
    }
  }
  else  // 3-input function
  {
    tt2 = kitty::dynamic_truth_table( num_in );
    tt2._bits[0] = p2;
    fmt::print( "[debug] detected ternary function\n");
    tt3 = kitty::dynamic_truth_table( num_in );
    tt3._bits[0] = p3;
    fmt::print( "[debug] creating node for {}\n", kitty::to_hex(tt1));
    auto n1 = klut.create_node(children, tt1);
    fmt::print( "[debug] creating node for {}\n", kitty::to_hex(tt2));
    auto n2 = klut.create_node(children, tt2);
    fmt::print( "[debug] creating node for {}\n", kitty::to_hex(tt3));
    auto n3 = klut.create_node(children, tt3);
    if      (tt == (tt1 | tt2 | tt3)) 
    {
      fmt::print( "[debug] detected OR3\n");
      p = klut.create_nary_or( {n1, n2, n3} ); 
    }
    else if (tt == (tt1 & tt2 & tt3)) 
    {
      fmt::print( "[debug] detected AND3\n");
      p = klut.create_nary_and( {n1, n2, n3} );
    }
    else if (tt == maj(tt1, tt2, tt3)) 
    {
      fmt::print( "[debug] detected MAJ3\n");
      p = klut.create_maj( n1, n2, n3 ); 
    }
    else {
      throw("The 3-input TT could not be recognized from predecessors");
    }
  }
  fmt::print( "[debug] Check {}\n", ++ct);
  klut.substitute_node(node, p);
  fmt::print( "[debug] Fanout of an old node is {}\n", klut.fanout_size(node));
  
}

void substitute_nodes(mockturtle::klut_network& klut)
{
  // klut._storage
  klut.foreach_gate(
      [&]( auto n ) { func(klut, n); }
  );
};

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, double, bool> exp( "lut_mapper", "benchmark", "luts", "lut_depth", "edges", "runtime", "equivalent" );

  lut_map_params ps;
  ps.cut_enumeration_ps.cut_size = 4u;
  ps.cut_enumeration_ps.cut_limit = 10u;
  ps.recompute_cuts = false;
  ps.area_oriented_mapping = false; //true;
  ps.cut_expansion = true;
  ps.edge_optimization = true; //can try set to true
  ps.remove_dominated_cuts = false; //can try set to true
  ps.cost_cache_vars = 0u; 
  ps.verbose = true;
  lut_map_stats st;

  // fmt::print( "[i] processing CLA\n" );

  bool write = false;

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    // Routine to run benchmakrs
    // fmt::print( "[i] processing {}\n", benchmark );
    // aig_network aig;
    // if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    // {
    //   continue;
    // }

    // if (benchmark != "bar") continue; 


    const unsigned short NBITS = 4;
    fmt::print( "[i] processing {}-bit CLA\n", NBITS );
    mockturtle::aig_network aig = create_adder<mockturtle::aig_network>( NBITS, mockturtle::carry_ripple_adder_inplace<mockturtle::aig_network> );
    // aig_network aig;
    // lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) );
    // fmt::print( "[i] processing {}\n", benchmark );

    mapping_view<aig_network, true> mapped_aig{ aig };

    lut_map<decltype( mapped_aig ), true, lut_sfq_cost>( mapped_aig, ps, &st );
    auto klut = *collapse_mapped_network<klut_network>( mapped_aig );

    depth_view<klut_network> klut_d{ klut };

    // auto const cec = (benchmark == "hyp") ? true : experiments::abc_cec( klut, benchmark );

    exp( benchmark, klut.num_gates(), klut_d.depth(), st.edges, to_seconds( st.time_total ), true );

    // mockturtle::write_bench(klut_d, benchmark + "_OR_only_xorables.bench");
    // if (write) write_bench(klut_d, benchmark + "_" + suffix + ".bench");
    
    write_bench(klut_d, benchmark + "_reduced_delay_raw_" + suffix + ".bench");
    fmt::print( "[i] Completed mapping {}\n", benchmark );

    fmt::print( "[debug] 1st pass\n" );
    klut.foreach_gate(
      [&]( auto n ) { 
        func(klut, n); 
      }
    ); 
    fmt::print( "[debug] Cleaning up dangling\n" );
    klut = mockturtle::cleanup_dangling(klut, true, true);
    fmt::print( "[debug] Cleaning up LUTs\n" );
    klut = mockturtle::cleanup_luts(klut);

    fmt::print( "[debug] 2nd pass\n" );
    klut.foreach_gate(
      [&]( auto n ) { 
        func(klut, n); 
      }
    ); 
    fmt::print( "[debug] Cleaning up dangling\n" );
    klut = mockturtle::cleanup_dangling(klut, true, true);
    fmt::print( "[debug] Cleaning up LUTs\n" );
    klut = mockturtle::cleanup_luts(klut);

    write_bench(klut, benchmark + "_reduced_delay_" + suffix + ".bench");

    break;
  }
  
  exp.save();
  exp.table();

  return 0;
}

  // ps.area_oriented_mapping = false; //false;
  // ps.edge_optimization = false; //can try set to true
  // ps.remove_dominated_cuts = false; //can try set to true
  // analyze_adder(exp, ps, st, "CLA0", suffix, write );
  // ps.area_oriented_mapping = false; //false;
  // ps.edge_optimization = false; //can try set to true
  // ps.remove_dominated_cuts = true; //can try set to true
  // analyze_adder(exp, ps, st, "CLA1", suffix, write );
  // ps.area_oriented_mapping = false; //false;
  // ps.edge_optimization = true; //can try set to true
  // ps.remove_dominated_cuts = false; //can try set to true
  // analyze_adder(exp, ps, st, "CLA2", suffix, write );
  // ps.area_oriented_mapping = false; //false;
  // ps.edge_optimization = true; //can try set to true
  // ps.remove_dominated_cuts = true; //can try set to true
  // analyze_adder(exp, ps, st, "CLA3", suffix, write );
  // ps.area_oriented_mapping = true; //false;
  // ps.edge_optimization = false; //can try set to true
  // ps.remove_dominated_cuts = false; //can try set to true
  // analyze_adder(exp, ps, st, "CLA4", suffix, write );
  // ps.area_oriented_mapping = true; //false;
  // ps.edge_optimization = false; //can try set to true
  // ps.remove_dominated_cuts = true; //can try set to true
  // analyze_adder(exp, ps, st, "CLA5", suffix, write );
  // ps.area_oriented_mapping = true; //false;
  // ps.edge_optimization = true; //can try set to true
  // ps.remove_dominated_cuts = false; //can try set to true
  // analyze_adder(exp, ps, st, "CLA6", suffix, write );
  // ps.area_oriented_mapping = true; //false;
  // ps.edge_optimization = true; //can try set to true
  // ps.remove_dominated_cuts = true; //can try set to true
  // analyze_adder(exp, ps, st, "CLA7", suffix, write );

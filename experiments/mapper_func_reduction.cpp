/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/rsfq/rsfq_network_conversion.hpp>
#include <mockturtle/algorithms/rsfq/rsfq_path_balancing.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/nodes.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/algorithms/retiming.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/generators/arithmetic.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/rsfq_view.hpp>
#include <mockturtle/algorithms/functional_reduction.hpp>

#include <mockturtle/io/auxiliary_genlib.hpp>

// #include <mockturtle/utils/GNM_global.hpp> // GNM global is stored here

#include <experiments.hpp>

// constexpr int fDFF   = 0;
// // constexpr int fNOT   = 1;
// // constexpr int fMERGE = 2;
// // constexpr int fOR    = 3;
// // constexpr int fAND   = 4;
// // constexpr int fXOR   = 5;
// // constexpr int fOR3   = 6;
// // constexpr int fAND3  = 7;
// // constexpr int fMAJ3  = 8;
// // constexpr int fCB    = 9;
// constexpr int fSPL   = 10;
// // constexpr int fPI    = 11;
// // constexpr int fNOFUNC= 99;


typedef mockturtle::klut_network klut;
typedef mockturtle::xag_network   xag;
typedef mockturtle::xmg_network   xmg;
typedef mockturtle::mig_network   mig;
typedef mockturtle::aig_network   aig;

// constexpr std::array<int,12> COSTS = {7, 9, 8, 8, 8, 7, 11, 11, 11, 8, 7, 0};
// constexpr std::array<int,12> COSTS_CONNECT = {6, 10, 7, 7, 7, 11, 999, 999, 999, 7, 3, 0};

// Removed input buffers in AND/OR gates
constexpr std::array<int,12> COSTS_CONNECT = {6, 9, 7, 3, 3, 11, 999, 999, 999, 7, 3, 0};


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

template <typename Ntk>
std::tuple<mockturtle::binding_view<klut>, mockturtle::map_stats, double, double, bool> map_with_pb 
( 
  const std::string & benchmark, 
  const Ntk & input_ntk, 
  const mockturtle::tech_library<4u, mockturtle::classification_type::p_configurations> & tech_lib, 
  std::unordered_map<std::string, int> & nDFF_global, 
  bool area_oriented = false 
)
{
  mockturtle::map_params ps;
  ps.cut_enumeration_ps.minimize_truth_table = true;
  ps.cut_enumeration_ps.cut_limit = 24;
  ps.buffer_pis = false;
  if (area_oriented)
  {
      ps.skip_delay_round = true;
      ps.required_time = std::numeric_limits<float>::max();
  }
  mockturtle::map_stats st;
  fmt::print("Started mapping of {}\n", benchmark);
  mockturtle::binding_view<klut> res = map( input_ntk, tech_lib, ps, &st );
  fmt::print("Finished mapping of {}\n", benchmark);
  mockturtle::depth_view<mockturtle::binding_view<klut>> dv { res };

  std::map<klut::node, int> dff_count;
  std::map<klut::node, int> fanout_count;

  /* RSFQ path balancing */
  fmt::print("Started RSFQ path balancing of {}\n", benchmark);
  auto balanced_res = mockturtle::rsfq_path_balancing( res );
  fmt::print("Finished RSFQ path balancing of {}\n", benchmark);

  /* retiming */
  mockturtle::retime_params rps;
  rps.iterations = 100;
  rps.verbose = true;
  
  mockturtle::retime_stats rst;
  fmt::print("Started rsfq_generic_network_create_from_mapped of {}->net\n", benchmark);
  auto net = mockturtle::rsfq_generic_network_create_from_mapped( balanced_res );
  fmt::print("Finished rsfq_generic_network_create_from_mapped of {}->net\n", benchmark);
  fmt::print("Started retime of {}\n", benchmark);
  mockturtle::retime( net, rps, &rst );
  fmt::print("Finished retime of {}\n", benchmark);
  fmt::print("Started rsfq_generic_network_create_from_mapped of net->{}\n", benchmark);
  auto retime_res = mockturtle::rsfq_mapped_create_from_generic_network( net );
  fmt::print("Finished rsfq_generic_network_create_from_mapped of net->{}\n", benchmark);

  uint32_t num_ext_dffs = retime_res.num_dffs();
  
  uint32_t num_int_dffs = 0;

  retime_res.foreach_node( 
    [&]( auto const& n ) 
    {
      if ( !retime_res.has_binding( n ) )
        return;
      auto const& g = retime_res.get_binding( n );
      num_int_dffs += nDFF_global[g.name];
      // fmt::print("Node {}\tGate {}\tnDFF {}\n", n, g.name, nDFF_global.at(g.name));
    } 
  );

  /* RSFQ splitter insertion */
  uint32_t num_splitters = 0;
  retime_res.foreach_node( [&]( auto const& n ) {
    if ( !retime_res.is_constant( n ) )
      num_splitters += retime_res.fanout_size( n ) - 1;
  } );

  fmt::print("Started rsfq_check_buffering of {}\n", benchmark);
  bool cec = rsfq_check_buffering( retime_res );
  fmt::print("Finished rsfq_check_buffering of {}\n", benchmark);
  fmt::print("Started abc_cec of {}\n", benchmark);
  cec &= benchmark == "hyp" ? true : experiments::abc_cec( retime_res, benchmark );
  fmt::print("Finished abc_cec of {}\n", benchmark);

  // Internal DFF area is already counted in the library
  // External DFF area is already counted after retiming
  double total_ndff = num_int_dffs + num_ext_dffs;
  double total_area = st.area + COSTS_CONNECT[fSPL] * num_splitters;
  //  +  COSTS_CONNECT[fDFF] * num_ext_dffs;
  fmt::print("\t{} : Int: {}, Ext: {}, ratio: {}\n", benchmark, num_int_dffs, num_ext_dffs, (float)num_int_dffs / (num_int_dffs + num_ext_dffs) );
  return std::make_tuple( res, st, total_ndff, total_area, cec );
}

template <typename Ntk>
std::tuple<mockturtle::binding_view<klut>, mockturtle::map_stats, int, int, bool> map_wo_pb 
( 
  const std::string & benchmark, 
  const Ntk & input_ntk, 
  const mockturtle::tech_library<4u, mockturtle::classification_type::p_configurations> & tech_lib, 
  std::unordered_map<std::string, int> & nDFF_global, 
  bool area_oriented = false 
)
{
  mockturtle::map_params ps;
  ps.cut_enumeration_ps.minimize_truth_table = true;
  ps.cut_enumeration_ps.cut_limit = 24;
  ps.buffer_pis = false;
  if (area_oriented)
  {
      ps.skip_delay_round = true;
      ps.required_time = std::numeric_limits<float>::max();
  }
  mockturtle::map_stats st;
  mockturtle::binding_view<klut> res = map( input_ntk, tech_lib, ps, &st );
  mockturtle::depth_view<mockturtle::binding_view<klut>> dv { res };

  std::map<klut::node, int> dff_count;

  uint32_t num_int_dffs = 0;
  res.foreach_node( [&]( auto const& n ) {
    if ( !res.has_binding( n ) )
      return;
    auto const& g = res.get_binding( n );
    num_int_dffs += nDFF_global[res.get_binding( n ).name];
    // fmt::print("Node {}\tGate {}\tnDFF {}\n", n, g.name, nDFF_global.at(g.name));
  } );

  /* RSFQ splitter insertion */
  uint32_t num_splitters = 0;
  res.foreach_node( [&]( auto const& n ) {
    if ( !res.is_constant( n ) )
      num_splitters += res.fanout_size( n ) - 1;
  } );

  bool cec = rsfq_check_buffering( res );
  cec &= benchmark == "hyp" ? true : experiments::abc_cec( res, benchmark );

  // int internal_dff = 0;

  // dv.foreach_node( [&]( auto const& node ) 
  // {
  //     if ( dv.has_binding( node ) )
  //     {
  //       auto const& g = dv.get_binding( node );
        // fmt::print("Node {}\tGate {}\tnDFF {}\n", node, g.name, nDFF_global.at(g.name));
        // internal_dff += nDFF_global[g.name];//[];
        // dv.foreach_fanin(node, [&](auto pred, auto i) {
        //   int expected_lvl = dv.level(node) - g.pins[i].rise_block_delay;//g.num_vars-1-
        //   int actual_lvl = dv.level(pred);
        //   int extra_delay = expected_lvl - actual_lvl;

        //   dff_count[pred] = std::max(dff_count[pred], extra_delay);
          // fanout_count[pred]++;
        // });
  //     }
  // } );

  // int PB_DFF = std::accumulate(dff_count.begin(), dff_count.end(), 0, [](int a, const auto& p){ return a + p.second; });
  // PB_DFF += internal_dff;
  // int ext_SPL = std::accumulate(fanout_count.begin(), fanout_count.end(), 0, [](int a, const auto& p){ 
  //     return (p.second > 1)?(a + p.second - 1):(a); 
  //     });
  
  // Internal DFF area is already counted in the library
  const int num_ext_dffs = 0;
  int total_area = st.area + COSTS_CONNECT[fSPL] * num_splitters;
  // +  COSTS_CONNECT[fDFF] * num_ext_dffs;
  fmt::print("\t{} : Int: {}, Ext: {}, ratio: {}\n", benchmark, num_int_dffs, num_ext_dffs, (float)num_int_dffs / (num_int_dffs + num_ext_dffs) );
  return std::make_tuple( res, st, num_int_dffs + num_ext_dffs, total_area, cec );
}

// mockturtle::binding_view<klut> recompose_cells(

// )

mockturtle::binding_view<klut> decompose_cells(
  mockturtle::binding_view<klut> ntk, 
  std::unordered_map<ULL, Node> nodemap,
  std::unordered_map<std::string, LibEntry> entries,
  int max_replacements
  )
{

  int replacements = 0;
  ntk.foreach_node( [&]( const klut::node & n ) 
  {
    if ( !ntk.has_binding( n ) || replacements >= max_replacements)
    {
      return;
    }
    std::vector<klut::node> klut_fanins;
    ntk.foreach_fanin(n, [&klut_fanins] (const klut::node & parent) 
    {
      klut_fanins.push_back(parent);
    } );  

    if (klut_fanins.size() == 2)
    {
      replacements++;
    }
    else
    {
      return;
    }
    auto const& g = ntk.get_binding( n );

    // fmt::print( "Analyzing node {0} using cell {1}\n", n, g.name );

    LibEntry entry = entries.at(g.name);
    // fmt::print( "\tLibEntry hash {0}\n", entry.hash );
    Node & root = nodemap[entry.hash];
    // fmt::print( "\tNode {0}\n", root.to_str() );
    // std::unordered_map<uint16_t, uint8_t> pi_func_map = {
    //   {0x5555, entry.chars[0] - 'a'},
    //   {0x3333, entry.chars[1] - 'a'},
    //   {0x0F0F, entry.chars[2] - 'a'},
    //   {0x00FF, entry.chars[3] - 'a'}
    // };
    fmt::print("Replacing node {0} ({1}) with node 0x{2:x} ({3}) fcode={4} \n", n, g.name, g.function._bits[0], root.to_str(), root.last_func);
    // Now, get the list of gates in topological order
    std::vector<ULL> topo_order;
    root.topo_sort(nodemap, topo_order, true);

    std::unordered_map<uint64_t, klut::node> node2klut;
    for (auto & hash : topo_order)
    {
      Node & node = nodemap[hash];
      klut::node klut_node;
      switch (node.last_func)
      {
        case fPI:
        {
          auto pi_idx = entry.chars[PIFUNC2IDX.at(node.func)] - 'a';
          klut_node = klut_fanins[pi_idx];
          break;
        }
        case fDFF:
        {
          klut::node input_klut_node = node2klut.at(node.parent_hashes.front());
          klut_node = ntk.create_buf(input_klut_node);
          break;
        }
        case fNOT:
        {
          klut::node input_klut_node = node2klut.at(node.parent_hashes.front());
          klut_node = ntk.create_not(input_klut_node);
          break;
        }
        case fAND:
        {
          klut::node input_klut_node_1 = node2klut.at(node.parent_hashes.front());
          klut::node input_klut_node_2 = node2klut.at(node.parent_hashes.back() );
          klut_node = ntk.create_and(input_klut_node_1, input_klut_node_2);
          break;
        }
        case fOR:
        // case fMERGE: //likely unused
        {
          klut::node input_klut_node_1 = node2klut.at(node.parent_hashes.front());
          klut::node input_klut_node_2 = node2klut.at(node.parent_hashes.back() );
          klut_node = ntk.create_or(input_klut_node_1, input_klut_node_2);
          break;
        }
        case fCB:
        {
          klut::node input_klut_node_1 = node2klut.at(node.parent_hashes.front());
          klut::node input_klut_node_2 = node2klut.at(node.parent_hashes.back() );
          klut_node = ntk.create_or(input_klut_node_1, input_klut_node_2);
          break;
        }
        case fXOR:
        {
          klut::node input_klut_node_1 = node2klut.at(node.parent_hashes.front());
          klut::node input_klut_node_2 = node2klut.at(node.parent_hashes.back() );
          klut_node = ntk.create_xor(input_klut_node_1, input_klut_node_2);
          break;
        }
        default:
        {
          throw "Unsupported function!";
          break;
        }
      }
      node2klut.emplace(hash, klut_node);
    }
    uint64_t root_hash = topo_order.back();
    klut::node klut_root = node2klut[root_hash];
    ntk.substitute_node(n, klut_root);
  } );
    // LibEntry & entry = 

    // fmt::print("Node {}\tGate {}\tnDFF {}\n", n, g.name, nDFF_global.at(g.name));
  return ntk;
}

// union NodeParams
// {
//   uint64_t value;
//   struct
//   {
//     unsigned node : 32; // id of the underlying node in the network
//     unsigned func : 16; // 4-input function
//     unsigned type : 2; // AA - 0, AS - 1, SA - 2
//     unsigned fanin: 3; //maximum of four inputs - i.e. 3 bits are sufficient
//   };
template <typename Ntk>
struct NodeParams
{
  using ntk_signal = typename Ntk::signal;
  ntk_signal signal; // id of the underlying node in the network
  uint16_t func; // 4-input function
  uint8_t type; // AA - 0, AS - 1, SA - 2
  uint8_t fanin; //maximum of four inputs - i.e. 3 bits are sufficient
  NodeParams(ntk_signal _signal, uint16_t _func, uint8_t _type, uint8_t _fanin) : 
                 signal(_signal),    func(_func),   type(_type),  fanin(_fanin) {};
};



std::tuple<klut, std::unordered_map<klut::node, NodeParams<klut>>> decompose_to_klut(
  mockturtle::binding_view<klut> src, 
  std::unordered_map<ULL, Node> nodemap,
  std::unordered_map<std::string, LibEntry> entries,
  bool add_buffers = false
  )
{
  std::unordered_map<klut::node, NodeParams<klut>> tgt_node_params;

  std::unordered_map<klut::signal, klut::signal> src2tgt;

  mockturtle::rsfq_view<klut> tgt;
  src.foreach_pi( [&]( const klut::node & src_pi ) 
  {
    klut::signal tgt_pi = tgt.create_pi();
    src2tgt.emplace(src_pi, tgt_pi);
  } );
  // src.foreach_po( [&]( const klut::node & src_po ) 
  // {
  //   klut::signal tgt_po = tgt.create_po();
  //   src2tgt.emplace(src_po, tgt_po);
  // } );

  /* The bindings will be replaced in forward topological order.  */

  std::vector<klut::node> src_po;
  src_po.reserve(src.num_pos());
  src.foreach_po( [&] (auto src_n)
  {
    src_po.push_back(src_n);
  } );

  // std::vector<klut::node> src_top_order;
  // src_top_order.reserve( src.size() );
  mockturtle::topo_view<klut>( src ).foreach_node( [&]( auto src_node ) 
  {
    // fmt::print(
    //   "\nSRC  function : 0x{0:x} with {1} fanins\n", 
    //   src.node_function(src_node)._bits[0], 
    //   src.fanin_size(src_node)
    // );
    if ( !src.has_binding( src_node ) )
    {
      return;
    }

    auto const& g = src.get_binding( src_node );

    if (g.name == "one" && g.expression == "CONST1")
    {
      klut::node tgt_node = 1;
      src2tgt.emplace(src_node, tgt_node);
      return;
    }
    else if (g.name == "zero" && g.expression == "CONST0")
    {
      klut::node tgt_node = 0;
      src2tgt.emplace(src_node, tgt_node);
      return;
    }

    LibEntry entry = entries.at(g.name);
    Node & root = nodemap[entry.hash];


    // fmt::print("Binding : {0}\n", entry.to_str());
    // fmt::print("NODE function : 0x{0:{2}x}, {1}\n", root.func, root.to_str(), std::max(1 << (src.fanin_size(src_node) - 2), 1) );

    // kitty::dynamic_truth_table func(src.fanin_size(src_node));
    // func._bits = root.func;
    // tgt._create_node(parent_tgt_1, parent_tgt_2);

    /* Get the topological order of internal nodes (i.e., nodes inside the cell) */
    std::vector<ULL> topo_order;
    root.topo_sort(nodemap, topo_order, true);

    std::vector<klut::node> tgt_fanins;
    src.foreach_fanin(src_node, [&](const auto & src_fanin)
    {
      tgt_fanins.push_back(src2tgt.at(src_fanin));
      // fmt::print("Recorded fanin of src_node {0}:\n\tsrc_fanin: {1}\n\ttgt_fanin: {2}\n", src_node, src_fanin, tgt_fanins.back());
    } );

    std::unordered_map<uint64_t, klut::node> node2tgt;
    for (ULL hash : topo_order)
    {
      Node & node = nodemap.at(hash);
      if (node.last_func == fPI)
      {
        /*determine which index of the PI is needed */
        // auto pi_idx = entry.chars[PIFUNC2IDX.at(node.func)] - 'a';
        char letter;
        if      (node.func == 0x5555) { letter = 'a'; }
        else if (node.func == 0x3333) { letter = 'b'; }
        else if (node.func == 0x0F0F) { letter = 'c'; }
        else if (node.func == 0x00FF) { letter = 'd'; }

        auto pi_idx = std::find(entry.chars.begin(), entry.chars.end(), letter) - entry.chars.begin();
        node2tgt.emplace(hash, tgt_fanins[pi_idx]);
      }
      else if (node.last_func == fDFF)
      {
        klut::node parent_tgt = node2tgt.at(node.parent_hashes.front());
        // klut::node tgt_node = tgt.create_buf(parent_tgt);
        klut::node tgt_node;
        if (add_buffers)
        {
          kitty::dynamic_truth_table _DFF_TT( 1u ); 
          _DFF_TT._bits[0] = 0b10;
          tgt_node = tgt.create_node({ parent_tgt }, _DFF_TT);
          // fmt::print("parent node: {0}, new_node: {1}\n", parent_tgt, tgt_node);
        }
        else
        {
          tgt_node = tgt.create_buf( parent_tgt );
        }
        fmt::print("Created node n{0} = BUF({1})\n", tgt_node, parent_tgt);
        node2tgt.emplace(hash, tgt_node);

        tgt_node_params.emplace(tgt_node, NodeParams<klut>(tgt_node, 0xAAAAu, 1u, 1u)); //0xAAAA s simply a DFF TT
      }
      else if (node.last_func == fNOT)
      {
        klut::node parent_tgt = node2tgt.at(node.parent_hashes.front());
        klut::node tgt_node = tgt.create_not(parent_tgt);
        // fmt::print("Created node n{0} = NOT({1})\n", tgt_node, parent_tgt);
        node2tgt.emplace(hash, tgt_node);
        tgt_node_params.emplace(tgt_node, NodeParams<klut>(tgt_node, 0x5555u, 1u, 1u));
      }
      else if (node.last_func == fAND)
      {
        klut::node parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::node parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::node tgt_node = tgt.create_and(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = AND({1}, {2})\n", tgt_node, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_node);
        tgt_node_params.emplace(tgt_node, NodeParams<klut>(tgt_node, 0x8888u, 2u, 2u));
      }
      else if (node.last_func == fOR)
      {
        klut::node parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::node parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::node tgt_node = tgt.create_or(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = OR({1}, {2})\n", tgt_node, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_node);
        tgt_node_params.emplace(tgt_node, NodeParams<klut>(tgt_node, 0xEEEEu, 2u, 2u));
      }
      else if (node.last_func == fCB)
      {
        klut::node parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::node parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::node tgt_node = tgt.create_or(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = CB({1}, {2})\n", tgt_node, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_node);
        tgt_node_params.emplace(tgt_node, NodeParams<klut>(tgt_node, 0xEEEEu, 0u, 2u));
      }
      else if (node.last_func == fXOR)
      {
        klut::node parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::node parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::node tgt_node = tgt.create_xor(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = XOR({1}, {2})\n", tgt_node, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_node);
        tgt_node_params.emplace(tgt_node, NodeParams<klut>(tgt_node, 0x6666u, 1u, 2u));
      }
      // else if 
      else
      {
        throw "Unsupported function";
      }
    }
    uint64_t root_hash = topo_order.back();
    klut::node tgt_node = node2tgt.at(root_hash);
    src2tgt.emplace(src_node, tgt_node);
  } );

  src.foreach_po([&](auto const & src_po)
  {
    tgt.create_po(src2tgt.at(src_po));
  } );

  return std::make_pair(tgt, tgt_node_params);
}

// buffer(A) is expressed as XOR3(false, <invert?>, A) 
xmg::signal create_buffer(xmg & tgt, xmg::signal & fanin, bool invert)
{
  return tgt.create_xor3( tgt.get_constant(false), tgt.get_constant(invert), fanin );
  // fmt::print("Created buffer n{0} = XOR3( 0, {2:d}, {1})\n", tgt.get_node(fanout), tgt.get_node(fanin), invert );
}

// buffer(A) in XAG is expressed as XOR(false, <invert?>, A) 
xag::signal create_buffer(xag & tgt, xag::signal & fanin, bool invert)
{
  return tgt.create_xor( tgt.get_constant(invert), fanin );
  // fmt::print("Created buffer n{0} = XOR( 0, {2:d}, {1})\n", tgt.get_node(fanout), tgt.get_node(fanin), invert );
}

aig::signal create_buffer(aig & tgt, aig::signal & fanin, bool invert)
{
  if (!invert)
  {
    return tgt.create_and( tgt.get_constant( true ), fanin );
  }
  else
  {
    return !tgt.create_and( tgt.get_constant( true ), !fanin );
  }
  // fmt::print("Created buffer n{0} = XOR3( 0, {2:d}, {1})\n", tgt.get_node(fanout), tgt.get_node(fanin), invert );
}


template <typename Ntk>
std::tuple<Ntk, std::unordered_map<typename Ntk::signal, NodeParams<Ntk>>> decompose_to_ntk(
  mockturtle::binding_view<klut> src, 
  std::unordered_map<ULL, Node> nodemap,
  std::unordered_map<std::string, LibEntry> entries,
  bool add_buffers = false,
  bool add_inverters = true
  )
{
  using signal = typename Ntk::signal;
  std::unordered_map<signal, NodeParams<Ntk>> tgt_sig_params;

  std::unordered_map<klut::signal, signal> src2tgt;
  
  Ntk tgt;
  src.foreach_pi( [&]( const klut::signal & src_pi ) 
  {
    signal tgt_pi = tgt.create_pi();
    src2tgt.emplace(src_pi, tgt_pi);
  } );

  /* The bindings will be replaced in forward topological order.  */
  std::vector<klut::signal> src_po;
  src_po.reserve(src.num_pos());
  src.foreach_po( [&] (auto src_n)
  {
    src_po.push_back(src_n);
  } );

  // std::vector<klut::node> src_top_order;
  // src_top_order.reserve( src.size() );
  mockturtle::topo_view<klut>( src ).foreach_node( [&]( auto src_node ) 
  {
    // fmt::print(
    //   "\nSRC  function : 0x{0:x} with {1} fanins\n", 
    //   src.node_function(src_node)._bits[0], 
    //   src.fanin_size(src_node)
    // );
    if ( !src.has_binding( src_node ) )
    {
      return;
    }

    auto const& g = src.get_binding( src_node );

    /* handing constants */
    if (g.name == "one" && g.expression == "CONST1")
    {
      signal tgt_sig = tgt.get_constant( true );
      src2tgt.emplace(src_node, tgt_sig);
      return;
    }
    else if (g.name == "zero" && g.expression == "CONST0")
    {
      signal tgt_sig = tgt.get_constant( false );
      src2tgt.emplace(src_node, tgt_sig);
      return;
    }
    
    LibEntry entry = entries.at(g.name);
    Node & root = nodemap[entry.hash];

    // fmt::print("Binding : {0}\n", entry.to_str());
    // fmt::print("NODE function : 0x{0:{2}x}, {1}\n", root.func, root.to_str(), std::max(1 << (src.fanin_size(src_node) - 2), 1) );

    // kitty::dynamic_truth_table func(src.fanin_size(src_node));
    // func._bits = root.func;
    // tgt._create_node(parent_tgt_1, parent_tgt_2);

    /* Get the topological order of internal nodes (i.e., nodes inside the cell) */
    std::vector<ULL> topo_order;
    root.topo_sort(nodemap, topo_order, true);

    std::vector<signal> tgt_fanins;
    src.foreach_fanin(src_node, [&](const auto & src_fanin)
    {
      tgt_fanins.push_back(src2tgt.at(src_fanin));
      // fmt::print("\tRecorded fanin of src_node {0}:\n\t\tsrc_fanin: {1}\n\t\ttgt_fanin: {2}\n", src_node, src_fanin, tgt_fanins.back().data);
    } );

    std::unordered_map<ULL, signal> node2tgt;
    for (ULL hash : topo_order)
    {
      Node & node = nodemap.at(hash);
      if (node.last_func == fPI)
      {
        /*determine which index of the PI is needed */
        char letter = node.pi_letter();

        auto pi_idx = std::find(entry.chars.begin(), entry.chars.end(), letter) - entry.chars.begin();
        node2tgt.emplace(hash, tgt_fanins[pi_idx]);
      }
      else if (node.last_func == fDFF)
      {
        signal parent_tgt = node2tgt.at(node.parent_hashes.front());
        signal tgt_sig = add_buffers ? create_buffer(tgt, parent_tgt, false) : tgt.create_buf( parent_tgt );
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, NodeParams<Ntk>(tgt_sig, 0xAAAAu, 1u, 1u)); //0xAAAA s simply a DFF TT
      }
      else if (node.last_func == fNOT)
      {
        signal parent_tgt = node2tgt.at(node.parent_hashes.front());
        signal tgt_sig = add_inverters ? create_buffer(tgt, parent_tgt, true ) : tgt.create_not( parent_tgt );
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, NodeParams<Ntk>(tgt_sig, 0x5555u, 1u, 1u));
      }
      else if (node.last_func == fAND)
      {
        signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        signal tgt_sig = tgt.create_and(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = AND({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, NodeParams<Ntk>(tgt_sig, 0x8888u, 2u, 2u));
      }
      else if (node.last_func == fOR)
      {
        signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        signal tgt_sig = tgt.create_or(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = OR({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, NodeParams<Ntk>(tgt_sig, 0xEEEEu, 2u, 2u));
      }
      else if (node.last_func == fCB)
      {
        signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        signal tgt_sig = tgt.create_or(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = CB({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, NodeParams<Ntk>(tgt_sig, 0xEEEEu, 0u, 2u));
      }
      else if (node.last_func == fXOR)
      {
        signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        signal tgt_sig = tgt.create_xor(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = XOR({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, NodeParams<Ntk>(tgt_sig, 0x6666u, 1u, 2u));
      }
      // else if 
      else
      {
        throw "Unsupported function";
      }
    }
    ULL root_hash = topo_order.back();
    signal tgt_sig = node2tgt.at(root_hash);
    src2tgt.emplace(src_node, tgt_sig);
  } );

  src.foreach_po([&](auto const & src_po)
  {
    tgt.create_po(src2tgt.at(src_po));
  } );

  return std::make_pair(tgt, tgt_sig_params);
}
  // TODO: 
  // TODO : Here I need to replace each binding with the appropriate primitives.
  // TODO : record the timing constraints for the ILP
  // TODO : record the mapping from src to tgt network
  // TODO : record the data pertaining to each node :
  //        - whether the element is AA, AS, or SA
  //          - AA elements are placed at the phase of the latest input
  //          - SA elements tie the preceding AS elements to itself to ensure simultaneous arrival of pulses  
  //        - anything else???

  // /* start replacing cells with regular nodes */
  // src.foreach_node( [&]( const klut::node & n ) 
  // {
  //   if ( !src.has_binding( n ) )
  //   {
  //     return;
  //   }
  //   /* Get fanins of the node */
  //   std::vector<klut::node> src_klut_fanins;
  //   src.foreach_fanin(n, [&src_klut_fanins] (const klut::node & parent) 
  //   {
  //     src_klut_fanins.push_back(parent);
  //   } );  ;

  //   /* get library cell information */
  //   auto const& g = src.get_binding( n ); // fmt::print( "Analyzing node {0} using cell {1}\n", n, g.name );
  //   LibEntry entry = entries.at(g.name);  // fmt::print( "\tLibEntry hash {0}\n", entry.hash );
  //   Node & root = nodemap[entry.hash];    // fmt::print( "\tNode {0}\n", root.to_str() );
  // //   fmt::print("Replacing node {0} ({1}) with node 0x{2:x} ({3}) fcode={4} \n", n, g.name, g.function._bits[0], root.to_str(), root.last_func);

  //   /* Get the topological order of internal nodes (i.e., nodes inside the cell) */
  //   std::vector<ULL> topo_order;
  //   root.topo_sort(nodemap, topo_order, true);

  //   /* mapping from Node to tgt klut*/
  //   std::unordered_map<uint64_t, klut::node> node_to_tgt_klut;
  //   for (auto & hash : topo_order)
  //   {
  //     Node & node = nodemap[hash];
  //     klut::node tgt_node;
  //     switch (node.last_func)
  //     {
  //       case fPI:
  //       {
  //         auto pi_idx = entry.chars[PIFUNC2IDX.at(node.func)] - 'a';
  //         tgt_node = src_klut_fanins[pi_idx];
  //         break;
  //       }
  //       case fDFF:
  //       {
  //         klut::node input_tgt_node = node_to_tgt_klut.at(node.parent_hashes.front());
  //         tgt_node = tgt.create_buf(input_tgt_node);
  //         break;
  //       }
  //       case fNOT:
  //       {
  //         klut::node input_tgt_node = node_to_tgt_klut.at(node.parent_hashes.front());
  //         tgt_node = src.create_not(input_tgt_node);
  //         break;
  //       }
  //       case fAND:
  //       {
  //         klut::node input_tgt_node_1 = node_to_tgt_klut.at(node.parent_hashes.front());
  //         klut::node input_tgt_node_2 = node_to_tgt_klut.at(node.parent_hashes.back() );
  //         tgt_node = src.create_and(input_tgt_node_1, input_tgt_node_2);
  //         break;
  //       }
  //       case fOR:
  //       // case fMERGE: //likely unused
  //       {
  //         klut::node input_tgt_node_1 = node_to_tgt_klut.at(node.parent_hashes.front());
  //         klut::node input_tgt_node_2 = node_to_tgt_klut.at(node.parent_hashes.back() );
  //         tgt_node = src.create_or(input_tgt_node_1, input_tgt_node_2);
  //         break;
  //       }
  //       case fCB:
  //       {
  //         klut::node input_tgt_node_1 = node_to_tgt_klut.at(node.parent_hashes.front());
  //         klut::node input_tgt_node_2 = node_to_tgt_klut.at(node.parent_hashes.back() );
  //         tgt_node = src.create_or(input_tgt_node_1, input_tgt_node_2);
  //         break;
  //       }
  //       case fXOR:
  //       {
  //         klut::node input_tgt_node_1 = node_to_tgt_klut.at(node.parent_hashes.front());
  //         klut::node input_tgt_node_2 = node_to_tgt_klut.at(node.parent_hashes.back() );
  //         tgt_node = src.create_xor(input_tgt_node_1, input_tgt_node_2);
  //         break;
  //       }
  //       default:
  //       {
  //         throw "Unsupported function!";
  //         break;
  //       }
  //     }
  //     node_to_tgt_klut.emplace(hash, tgt_node);
  //   }
  //   uint64_t root_hash = topo_order.back();
  //   klut::node klut_root = node_to_tgt_klut[root_hash];
  //   src.substitute_node(n, klut_root);
  // } );
  //   // LibEntry & entry = 

  //   // fmt::print("Node {}\tGate {}\tnDFF {}\n", n, g.name, nDFF_global.at(g.name));
  // return src;
// }

// Function to read unordered_map from CSV file
std::unordered_map<std::string, int> readCSV(const std::string& filename) 
{
    std::ifstream infile(filename);             // Open the input file stream
    std::unordered_map<std::string, int> map;   // Create the unordered_map
    
    std::string line;
    std::getline(infile, line);                 // Ignore the header row

    // Read each subsequent row and add the key-value pair to the unordered_map
    while (std::getline(infile, line)) 
    {
        std::stringstream ss(line);
        std::string key;
        int value;
        std::getline(ss, key, ',');
        ss >> value;
        map[key] = value;
    }
    infile.close(); // Close the input file stream
    return map;
}


/* Gate costs are based on CONNECT library (from Japan) */
// const std::string DATABASE_PATH { "LIBRARY_2023_05_19_CONNECT.genlib" } ;
/* CONNECT library (from Japan) */
// const std::string DATABASE_PATH { "/Users/brainkz/Documents/GitHub/mockturtle_alessandro/build/LIBRARY_VANILLA_CONNECT.genlib" } ;
// LIBRARY_2023_05_19_CONNECT.genlib
// const std::string DATABASE_PATH { "/Users/brainkz/Documents/GitHub/mockturtle_alessandro/build/LIBRARY_2023_06_26_CONNECT_1111.genlib" } ;
// const std::string DATABASE_PATH { "/Users/brainkz/Documents/GitHub/mockturtle_alessandro/build/LIBRARY_2023_05_19_CONNECT.genlib" } ;
const std::string DATABASE_PATH { "/Users/brainkz/Documents/GitHub/mockturtle/build/LIBRARY_2023_06_27_CONNECT_CONSERVATIVE.genlib" } ;
/*The number of internal DFFs within each cell. 
Some of them are necessary not only for path balancing but also 
for synchronizing the pulses for AND gates. I include them 
in total DFF count */
// const std::string NDFF_PATH { "/Users/brainkz/Documents/GitHub/mockturtle_alessandro/build/nDFF_2023_05_08_CONNECT.csv" } ; 
const std::string NDFF_PATH { "/Users/brainkz/Documents/GitHub/mockturtle_alessandro/build/nDFF_2023_06_27_CONNECT_CONSERVATIVE.csv" } ; 

void GNM_to_hpp()
{
  /* Import the existing data */
  const std::string nodemap_prefix = "/Users/brainkz/Documents/GitHub/mockturtle/build/Golden_20230427/x3";
  const std::vector<std::vector<UI>> sets_of_levels { { {0,0,0,0}, {0,0,0,1}, {0,0,0,2}, {0,0,1,1}, {0,0,1,2}, {0,1,1,1}, {0,1,1,2}, {0,1,2,2}, {0,1,2,3} } }; // {0,1,1,3}, 
  std::unordered_map<ULL, Node> GNM_global = read_global_gnm( sets_of_levels, nodemap_prefix );
  std::ofstream gnm_file ("/Users/brainkz/Documents/GitHub/mockturtle/include/mockturtle/utils/GNM_global.hpp");

  gnm_file << fmt::format(
    "#pragma once\n"
    "#include <iostream>\n"
    "#include <unordered_map>\n"
    "#include <vector>\n"
    "#include <mockturtle/algorithms/nodes.hpp>\n"
    "std::unordered_map<uint64_t, Node> get_map()\n{{\n\tstd::unordered_map<uint64_t, Node> GNM_global;\n");

  for (auto & [hash, node] : GNM_global)
  {
    std::vector<std::string> phashes(node.parent_hashes.size());
    for (auto i = 0u; i < node.parent_hashes.size(); ++i)
    {
      phashes[i] = fmt::format("0x{0:016x}", node.parent_hashes[i]);
    }
    gnm_file << fmt::format(
      "\tGNM_global.emplace(0x{0:016x}, Node(0x{1:04x}, {2:>2d}, {3:>3d}, {4:>2d}, {5:>5}, {{ {6} }}, 0x{0:016x} ) );\n",
      // "GNM_global[0x{0:016x}] = Node(0x{1:04x}, {2:>2d}, {3:>3d}, {4:>2d}, {5:>5}, {{ {6} }}, 0x{0:016x} );\n",
      hash,             // 0
      node.func,        // 1
      node.last_func,   // 2
      node.cost,        // 3
      node.depth,       // 4
      node.xorable,     // 5
      fmt::join(phashes, ", ") // 6
    );
  }
  gnm_file << fmt::format("\treturn GNM_global;\n}};");
  return;
}



int main()  //int argc, char* argv[]
{
  // // Check if the required number of command-line arguments is provided
  // if (argc < 2) 
  // {
  //     std::cout << "Usage: ./program_name max_replacements" << std::endl;
  //     return 1;
  // }

  // // Parse the command-line argument into an integer
  // int max_replacements = std::stoi(argv[1]);

  using namespace experiments;
  using namespace mockturtle;
  
  
    std::unordered_map<std::string, std::tuple<double,double,double>> PBMAP;
    PBMAP["int2float"] = std::make_tuple(  270,   6432,  16);
    PBMAP["priority"]  = std::make_tuple( 9064, 102085, 127);
    PBMAP["sin"]       = std::make_tuple(13666, 215318, 182);
    PBMAP["cavlc"]     = std::make_tuple(  522,  16339,  17);
    PBMAP["dec"]       = std::make_tuple(    8,   5469,   4);
    PBMAP["c499"]      = std::make_tuple(  476,   7758,  13);
    PBMAP["c880"]      = std::make_tuple(  774,  12909,  22);
    PBMAP["c1908"]     = std::make_tuple(  696,  12013,  20);
    PBMAP["c3540"]     = std::make_tuple( 1159,  28300,  31);
    PBMAP["c5315"]     = std::make_tuple( 2908,  52033,  23);
    PBMAP["c7552"]     = std::make_tuple( 2429,  48482,  19);
    

  // experiment<std::string, 
  //           double, double, double, 
  //           double, double, double, 
  //           double, double, double, 
  //           float> exp(
  //     "mapper", "benchmark", 
  //     "#DFF (base)", "#DFF (our)", "#DFF (ratio)", 
  //     "area (base)", "area (our)", "area (ratio)", 
  //     "delay (base)", "delay (our)", "delay (ratio)", 
  //     "time");

  experiment<std::string, 
            double, double, double, 
            double, double, double, 
            double, double, double, 
            double, double, double> exp(
      "mapper", "benchmark", 
      "#DFF (AIG)", "#DFF (XMG0)", "#DFF (XMG_FR)", 
      "area (AIG)", "area (XMG0)", "area (XMG_FR)", 
      "delay (AIG)", "delay (XMG0)", "delay (XMG_FR)", 
      "time (AIG)", "time (XMG0)", "time (XMG_FR)" 
      );

  fmt::print( "[i] processing technology library\n" );

  /* library to map to technology */
  std::vector<gate> gates;
  std::ifstream inputFile( DATABASE_PATH );
  // std::unordered_map<std::string, int> nDFF_global = readCSV( NDFF_PATH );
  std::unordered_map<std::string, int> nDFF_global;

  // Prints # DFF for each cell in the library
  /*
    int ctr = 0;
    for (auto & [key, value]: nDFF_global)
    {
      fmt::print("{}:\t{}\t", key, value);
      if (++ctr % 8 == 0) 
      {
        std::cout << std::endl;
      }
    }
  */
  
  if ( lorina::read_genlib( inputFile, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }


  // const std::vector<std::vector<UI>> sets_of_levels { { {0,1,2,3} } };
  const std::string LibEntry_file = "/Users/brainkz/Documents/GitHub/mockturtle/build/LibEntry_2023_06_27_CONNECT_CONSERVATIVE.csv";


  mockturtle::tech_library_params tps; // tps.verbose = true;
  tech_library<4, mockturtle::classification_type::p_configurations> tech_lib( gates, tps );

  // // Where to find the MCNC benchmarks?
  // auto benchmarks2 = epfl_benchmarks( experiments::sin | experiments::cavlc | experiments::int2float | experiments::priority | experiments::i2c | experiments::voter | experiments::dec );
  // //   auto benchmarks1 = epfl_benchmarks( experiments::epfl & ~experiments::div & ~experiments::hyp & ~experiments::log2 & ~experiments::sqrt );
  // auto benchmarks1 = iscas_benchmarks( experiments::c499 | experiments::c432 | experiments::c880 | experiments::c1908 | experiments::c3540 | experiments::c5315 | experiments::c7552 );
  // benchmarks1.insert(benchmarks1.end(), benchmarks2.begin(), benchmarks2.end());
  
  auto benchmarks1 = epfl_benchmarks(); //0xFFFFFFFFFFFFFFFF ^ experiments::div
  //   auto benchmarks1 = epfl_benchmarks( experiments::epfl & ~experiments::div & ~experiments::hyp & ~experiments::log2 & ~experiments::sqrt );
  auto benchmarks2 = iscas_benchmarks();
  benchmarks1.insert(benchmarks1.end(), benchmarks2.begin(), benchmarks2.end());

  // TODO :  Debug experiments::c5315 with GATE 0x9696_112_0 not showing up in GNM . 
  // TODO :  Perhaps genlib_hierarchical.cpp is broken, need to debug
  
  const std::string nodemap_prefix = "/Users/brainkz/Documents/GitHub/mockturtle/build/Golden_20230427/x3";
  const std::vector<std::vector<UI>> sets_of_levels { { {0,0,0,0}, {0,0,0,1}, {0,0,0,2}, {0,0,1,1}, {0,0,1,2}, {0,1,1,1}, {0,1,1,2}, {0,1,2,2}, {0,1,2,3} } }; //  {0,1,1,3},
  std::unordered_map<ULL, Node> GNM_global = read_global_gnm( sets_of_levels, nodemap_prefix );
  std::ofstream gnm_file ("/Users/brainkz/Documents/GitHub/mockturtle/include/mockturtle/utils/GNM_global.hpp");
  std::unordered_map<std::string, LibEntry> entries = read_LibEntry_map(LibEntry_file);

  for ( auto const& benchmark : benchmarks1 )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig aig_original;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig_original ) ) != lorina::return_code::success )
    {
      continue;
    }
    xmg xmg_original;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xmg_original ) ) != lorina::return_code::success )
    {
      continue;
    }

    // auto [res_area , st_area , PB_DFF_area , nJJ_area , cec_area ] = map_with_pb(benchmark, aig_original, tech_lib, nDFF_global, true );
    auto [res_no_pb, st_no_pb, PB_DFF_no_pb, nJJ_no_pb, cec_no_pb] = map_wo_pb(benchmark, aig_original, tech_lib, nDFF_global, false);

    if (abc_cec(res_no_pb, benchmark))
    {
      fmt::print("RES_no_pb is equivalent to the original\n");
    }
    else
    {
      fmt::print("ERROR: RES_no_pb is not equivalent\n");
    }

    auto [xmg_decomposed, tgt_node_params] = decompose_to_ntk<xmg>(res_no_pb, GNM_global, entries, false, false);

    if (abc_cec_cmp(xmg_decomposed, res_no_pb))
    {
      fmt::print("SUCCESS: Decomposed network is equivalent to the original\n");
    }
    else
    {
      fmt::print("ERROR: Not equivalent\n");
    }

    fmt::print("Benchmark \"{0}\"\n", benchmark);
    fmt::print("\tXMG size before: \"{0}\"\n", xmg_decomposed.size());
    functional_reduction_params fr_ps;
    fr_ps.verbose = true;
    // ps.max_TFI_nodes = 10000;
    // ps.skip_fanout_limit = 1000;
    // ps.conflict_limit = 1000;
    // ps.max_clauses = 10000;
    // ps.num_patterns = 1024; //4x default
    // ps.max_patterns = 4096; //4x default
    functional_reduction_stats fr_st;
    fmt::print("Started functional reduction of {}\n", benchmark);
    functional_reduction(xmg_decomposed, fr_ps, &fr_st);
    fmt::print("Finished functional reduction of {}\n", benchmark);

    fmt::print("Started cleanup dangling of {}\n", benchmark);
    xmg_decomposed = cleanup_dangling( xmg_decomposed );
    fmt::print("Finished cleanup dangling of {}\n", benchmark);

    fmt::print("\tXMG size after: \"{0}\"\n", xmg_decomposed.size());

    auto [res_aig_direct, st_aig_direct, PB_DFF_aig_direct, nJJ_aig_direct, cec_aig_direct] = map_with_pb(benchmark, aig_original, tech_lib, nDFF_global, false);

    auto [res_xmg_direct, st_xmg_direct, PB_DFF_xmg_direct, nJJ_xmg_direct, cec_xmg_direct] = map_with_pb(benchmark, xmg_original, tech_lib, nDFF_global, false);

    auto [res_after, st_after, PB_DFF_after, nJJ_after, cec_after] = map_with_pb(benchmark, xmg_decomposed, tech_lib, nDFF_global, false);

    auto t_base = to_seconds( st_aig_direct.time_total );

    exp( benchmark, PB_DFF_aig_direct, 
                    PB_DFF_xmg_direct/PB_DFF_aig_direct, 
                    PB_DFF_after /PB_DFF_aig_direct, 

                    nJJ_aig_direct, 
                    nJJ_xmg_direct/nJJ_aig_direct, 
                    nJJ_after /nJJ_aig_direct, 

                    st_aig_direct.delay, 
                    st_xmg_direct.delay/st_aig_direct.delay, 
                    st_after.delay/st_aig_direct.delay,

                    t_base * 1000 , 
                    to_seconds( st_xmg_direct.time_total ) / t_base, 
                    1 + (to_seconds( st_after.time_total ) + to_seconds( fr_st.time_total ))  / t_base  
                    );

    // fmt::print("\tNTK DFF BEFORE: \"{0}\"\n", PB_DFF_original);
    // fmt::print("\tNTK DFF AFTER: \"{0}\"\n\n", PB_DFF_after);
    // fmt::print("\tNTK AREA BEFORE: \"{0}\"\n", nJJ_original);
    // fmt::print("\tNTK AREA AFTER: \"{0}\"\n\n", nJJ_after);
    // fmt::print("\tNTK DELAY BEFORE: \"{0}\"\n", st_original.delay);
    // fmt::print("\tNTK DELAY AFTER: \"{0}\"\n\n", st_after.delay);
    // fmt::print("\tNTK CEC BEFORE: \"{0}\"\n\n", cec_original);
    // fmt::print("\tNTK CEC AFTER: \"{0}\"\n\n", cec_after);

//     if (PBMAP.find(benchmark) != PBMAP.end())
//     {
//         const double ndff  = std::get<0>(PBMAP[benchmark]);
//         const double njj   = std::get<1>(PBMAP[benchmark]);
//         const double delay = std::get<2>(PBMAP[benchmark]);
        
// //         exp( benchmark + " area" , ndff, PB_DFF_area,     PB_DFF_area/ndff, 
// //                                    njj,   nJJ_area ,        nJJ_area/njj,
// //                                    delay, st_area.delay, st_area.delay/delay,
// //                                    to_seconds( st_area.time_total  )  );
        
//         exp( benchmark + " delay", ndff, PB_DFF_no_pb,     PB_DFF_no_pb/ndff, 
//                                   njj,   nJJ_no_pb ,        nJJ_no_pb/njj,
//                                   delay, st_no_pb.delay, st_no_pb.delay/delay,
//                                   to_seconds( st_no_pb.time_total * 1000 )  );
//     }    
//     else
//     {        
// //         exp( benchmark + " area" , 0xFFFFFFFF, PB_DFF_area,    0xFFFFFFFF, 
// //                                    0xFFFFFFFF,   nJJ_area ,    0xFFFFFFFF,
// //                                    0xFFFFFFFF, st_area.delay,  0xFFFFFFFF,
// //                                    to_seconds( st_area.time_total  )  );
        
//         exp( benchmark + " delay", 0xFFFFFFFF, PB_DFF_no_pb,   0xFFFFFFFF, 
//                                   0xFFFFFFFF,   nJJ_no_pb ,   0xFFFFFFFF,
//                                   0xFFFFFFFF, st_no_pb.delay, 0xFFFFFFFF,  
//                                   to_seconds( st_no_pb.time_total * 1000 ) ); //
//     }    
    // Write .bench file for inspection
    // mockturtle::write_bench(res_area , benchmark + "_connect_area"  + ".bench");
    // mockturtle::write_bench(res_no_pb, benchmark + "_connect_full" + ".bench");
    mockturtle::write_bench(res_after, benchmark + "_xag_one_round" + ".bench");
    exp.save();
    exp.table();
    // break;
  }
  // exp.table();

  return 0;
}

/* results with full PB
|       benchmark |   #DFF (base) | #DFF (our) |  #DFF (ratio) |   area (base) | area (our) |  area (ratio) |  delay (base) | delay (our) | delay (ratio) |   time |
|       sin delay |      13666.00 |    9299.00 |          0.68 |     215318.00 |   84754.00 |          0.39 |        182.00 |       92.00 |          0.51 | 546.45 |
|     cavlc delay |        522.00 |     490.00 |          0.94 |      16339.00 |   12716.00 |          0.78 |         17.00 |       15.00 |          0.88 |  20.12 |
|       dec delay |          8.00 |       8.00 |          1.00 |       5469.00 |    6276.00 |          1.15 |          4.00 |        4.00 |          1.00 |  10.56 |
|       i2c delay | 4294967296.00 |    2348.00 | 4294967296.00 | 4294967296.00 |   24117.00 | 4294967296.00 | 4294967296.00 |       16.00 | 4294967296.00 |  30.90 |
| int2float delay |        270.00 |     211.00 |          0.78 |       6432.00 |    4683.00 |          0.73 |         16.00 |       15.00 |          0.94 |   9.17 |
|  priority delay |       9064.00 |    7323.00 |          0.81 |     102085.00 |   24185.00 |          0.24 |        127.00 |       84.00 |          0.66 |  31.67 |
|     voter delay | 4294967296.00 |    3655.00 | 4294967296.00 | 4294967296.00 |  164539.00 | 4294967296.00 | 4294967296.00 |       41.00 | 4294967296.00 | 924.69 |
|      c432 delay | 4294967296.00 |     596.00 | 4294967296.00 | 4294967296.00 |    3389.00 | 4294967296.00 | 4294967296.00 |       21.00 | 4294967296.00 |   8.60 |
|      c499 delay |        476.00 |     256.00 |          0.54 |       7758.00 |    4065.00 |          0.52 |         13.00 |        8.00 |          0.62 |  39.59 |
|      c880 delay |        774.00 |     586.00 |          0.76 |      12909.00 |    5125.00 |          0.40 |         22.00 |       15.00 |          0.68 |  15.88 |
|     c1908 delay |        696.00 |     401.00 |          0.58 |      12013.00 |    3141.00 |          0.26 |         20.00 |       13.00 |          0.65 |  28.96 |
|     c3540 delay |       1159.00 |     765.00 |          0.66 |      28300.00 |   16293.00 |          0.58 |         31.00 |       23.00 |          0.74 |  52.46 |
|     c5315 delay |       2908.00 |    1817.00 |          0.62 |      52033.00 |   22308.00 |          0.43 |         23.00 |       15.00 |          0.65 | 113.08 |
|     c7552 delay |       2429.00 |    2234.00 |          0.92 |      48482.00 |   17692.00 |          0.36 |         19.00 |       19.00 |          1.00 | 134.29 |



*/
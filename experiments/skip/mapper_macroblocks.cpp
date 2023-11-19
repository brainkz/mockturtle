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

typedef mockturtle::klut_network klut;
typedef mockturtle::xag_network   xag;
typedef mockturtle::xmg_network   xmg;
typedef mockturtle::mig_network   mig;
typedef mockturtle::aig_network   aig;

// constexpr std::array<int,12> COSTS = {7, 9, 8, 8, 8, 7, 11, 11, 11, 8, 7, 0};
// constexpr std::array<int,12> COSTS_CONNECT = {6, 10, 7, 7, 7, 11, 999, 999, 999, 7, 3, 0};

// Removed input buffers in AND/OR gates
constexpr std::array<int,12> COSTS_CONNECT = {6, 9, 7, 3, 3, 11, 999, 999, 999, 7, 3, 0};

template <typename Ntk>
std::tuple<mockturtle::binding_view<klut>, mockturtle::map_stats> map_wo_pb 
( 
  const Ntk & input_ntk, 
  const mockturtle::tech_library<4u, mockturtle::classification_type::p_configurations> & tech_lib, 
  const bool area_oriented = false,
  const bool verbose = false
)
{
  mockturtle::map_params ps;
  ps.verbose = verbose;
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
  return std::make_tuple( res, st );
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
  fmt::print("Started mapping of {}\n", benchmark);
  auto [res, st] = map_wo_pb(benchmark, input_ntk, tech_lib, area_oriented);
  fmt::print("Finished mapping of {}\n", benchmark);

  std::map<klut::node, int> dff_count;
  std::map<klut::node, int> fanout_count;

  /* RSFQ path balancing */
  fmt::print("Started RSFQ path balancing of {}\n", benchmark);
  auto balanced_res = mockturtle::rsfq_path_balancing( res );
  fmt::print("Finished RSFQ path balancing of {}\n", benchmark);

  mockturtle::retime_params rps;
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


// mockturtle::binding_view<klut> recompose_cells(

// )

template <typename Ntk>
struct Primitive
{
  using ntk_signal = typename Ntk::signal;
  ntk_signal sig; // id of the underlying node in the network
  uint16_t  func; // 4-input function
  uint8_t   type; // AA - 0, AS - 1, SA - 2
  std::vector<ntk_signal> fanins;
  Primitive(ntk_signal _sig, uint16_t _func, uint8_t _type, std::vector<ntk_signal> _fanins) : 
                   sig(_sig),    func(_func),   type(_type),                 fanins(_fanins) {};
};

template <typename Ntk>
struct ASBlk
{
  using ntk_signal = typename Ntk::signal;
  ntk_signal  fanout;                // fanout node in the network
  std::vector<ntk_signal> internals; // internal nodes
  std::vector<ntk_signal> fanins;    // fanout nodes
  ASBlk(ntk_signal _fanout, std::vector<ntk_signal> _internals, std::vector<ntk_signal> _fanins) : 
        fanout(_fanout), internals(_internals), fanins(_fanins) {};
};


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

// #if false
std::pair<klut, std::unordered_map<klut::signal, Primitive<klut>>> decompose_to_klut(mockturtle::binding_view<klut> src, std::unordered_map<ULL, Node> nodemap,std::unordered_map<std::string, LibEntry> entries)
{
  std::unordered_map<klut::signal, Primitive<klut>> tgt_sig_params;
  std::unordered_map<klut::signal, klut::signal> src2tgt;
  
  klut tgt;
  src.foreach_pi( [&]( const klut::signal & src_pi ) 
  {
    klut::signal tgt_pi = tgt.create_pi();
    src2tgt.emplace(src_pi, tgt_pi);
  } );

  // The bindings are replaced in forward topological order.
  std::vector<klut::signal> src_po;
  src_po.reserve(src.num_pos());
  src.foreach_po( [&] (auto src_n)
  {
    src_po.push_back(src_n);
  } );

  std::unordered_set<klut::signal> OR_replacement_candidates;
  std::vector<klut::signal> XOR_gates;

  auto num_node_src = src.size();
  auto ctr = 0u;
  mockturtle::topo_view<klut>( src ).foreach_node( [&]( auto src_node ) 
  {
    fmt::print("Processing node {0} ({1} out of {2})\r", src_node, ++ctr, num_node_src);
    if (ctr == num_node_src)
    {
      fmt::print("\n");
    }
    if ( !src.has_binding( src_node ) )
    {
      return;
    }

    auto const& g = src.get_binding( src_node );

    // handing constants
    if (g.name == "one" && g.expression == "CONST1")
    {
      klut::signal tgt_sig = tgt.get_constant( true );
      src2tgt.emplace(src_node, tgt_sig);
      return;
    }
    else if (g.name == "zero" && g.expression == "CONST0")
    {
      klut::signal tgt_sig = tgt.get_constant( false );
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

    // Get the topological order of internal nodes (i.e., nodes inside the cell)
    std::vector<ULL> topo_order;
    root.topo_sort(nodemap, topo_order, true);

    std::vector<klut::signal> tgt_fanins;
    src.foreach_fanin(src_node, [&](const auto & src_fanin)
    {
      tgt_fanins.push_back(src2tgt.at(src_fanin));
      // fmt::print("\tRecorded fanin of src_node {0}:\n\t\tsrc_fanin: {1}\n\t\ttgt_fanin: {2}\n", src_node, src_fanin, tgt_fanins.back().data);
    } );

    std::unordered_map<ULL, klut::signal> node2tgt;

    for (ULL hash : topo_order)
    {
      Node & node = nodemap.at(hash);
      if (node.last_func == fPI)
      {
        // determine which index of the PI is needed 
        char letter = node.pi_letter();

        auto pi_idx = std::find(entry.chars.begin(), entry.chars.end(), letter) - entry.chars.begin();
        node2tgt.emplace(hash, tgt_fanins[pi_idx]);
      }
      else if (node.last_func == fDFF)
      {
        klut::signal parent_tgt = node2tgt.at(node.parent_hashes.front());
        node2tgt.emplace(hash, parent_tgt);
        // tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0xAAAAu, 1u, { parent_tgt }));
        //  nothing else to do, no new constraints, no contribution to objective function
      }
      else if (node.last_func == fNOT)
      {
        klut::signal parent_tgt = node2tgt.at(node.parent_hashes.front());
        klut::signal tgt_sig = tgt.create_not( parent_tgt );
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0x5555u, 1u, { parent_tgt }));
      }
      else if (node.last_func == fAND)
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_and(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = AND({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0x8888u, 2u, { parent_tgt_1, parent_tgt_2 }));
      }
      else if (node.last_func == fOR) 
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_or(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = OR({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0xEEEEu, 2u, { parent_tgt_1, parent_tgt_2 }));
        OR_replacement_candidates.emplace(tgt_sig);
      }
      else if (node.last_func == fCB)
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_or(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = CB({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0xEEEEu, 0u, { parent_tgt_1, parent_tgt_2 }));
      }
      else if (node.last_func == fXOR)
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_xor(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = XOR({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0x6666u, 1u, { parent_tgt_1, parent_tgt_2 }));
        XOR_gates.push_back(tgt_sig);
      }
      // else if 
      else
      {
        throw "Unsupported function";
      }
    }
    ULL root_hash = topo_order.back();
    klut::signal tgt_sig = node2tgt.at(root_hash);
    src2tgt.emplace(src_node, tgt_sig);
  } );

  for (klut::signal const & xor_sig : XOR_gates)
  {
    std::vector<klut::signal> stack;
    std::vector<klut::signal> visited;
    Primitive<klut> & xor_prim = tgt_sig_params.at(xor_sig);

    for (klut::signal const & fanin : xor_prim.fanins)
    {
      // Primitive<klut> & fanin_prim = tgt_sig_params.at(fanin);
      // if the OR gate 
      auto it = std::find(OR_replacement_candidates.begin(), OR_replacement_candidates.end(), fanin);
      if (it != OR_replacement_candidates.end())
      {
        OR_replacement_candidates.erase(it);
      }
    }
  }

  for (klut::signal const & sig : OR_replacement_candidates)
  {
    Primitive<klut> & prim = tgt_sig_params.at(sig);
    prim.type = 0u;
  }

  src.foreach_po([&](auto const & src_po)
  {
    tgt.create_po(src2tgt.at(src_po));
  } );
  
  return std::make_pair(tgt, tgt_sig_params);
}
// #endif 

// std::vector<ASBlk> merge_blk(klut ntk, std::unordered_map<klut::signal, Primitive<klut>> sig_params)
// {
//   std::vector<klut::signal> stack;
//   stack.reserve(ntk.num_pos());
//   ntk.foreach_po([&] (auto node) 
//   {
//     stack.push_back(node);
//   });
//   // build macroblocks with DFS

//   while (!stack.empty())
//   {
    
//   }

//   // std::vector<klut::signal> reverse_topo( ntk.size() );

//   // size_t ctr = ntk.size();
//   // mockturtle::topo_view<klut>( ntk ).foreach_node( [&]( auto ntk_node ) 
//   // {
//   //   reverse_topo[--ctr] = ntk_node;
//   // });

//   // for (klut::signal & node : reverse_topo)
//   // {

//   // }

// }

template <typename Ntk>
std::tuple<Ntk, std::unordered_map<typename Ntk::signal, Primitive<Ntk>>> decompose_to_ntk(
  mockturtle::binding_view<klut> src, 
  std::unordered_map<ULL, Node> nodemap,
  std::unordered_map<std::string, LibEntry> entries,
  bool add_buffers = false,
  bool add_inverters = true
  )
{
  using signal = typename Ntk::signal;
  std::unordered_map<signal, Primitive<Ntk>> tgt_sig_params;

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
        tgt_sig_params.emplace(tgt_sig, Primitive<Ntk>(tgt_sig, 0xAAAAu, 1u, 1u)); //0xAAAA s simply a DFF TT
      }
      else if (node.last_func == fNOT)
      {
        signal parent_tgt = node2tgt.at(node.parent_hashes.front());
        signal tgt_sig = add_inverters ? create_buffer(tgt, parent_tgt, true ) : tgt.create_not( parent_tgt );
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<Ntk>(tgt_sig, 0x5555u, 1u, 1u));
      }
      else if (node.last_func == fAND)
      {
        signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        signal tgt_sig = tgt.create_and(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = AND({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<Ntk>(tgt_sig, 0x8888u, 2u, 2u));
      }
      else if (node.last_func == fOR)
      {
        signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        signal tgt_sig = tgt.create_or(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = OR({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<Ntk>(tgt_sig, 0xEEEEu, 2u, 2u));
      }
      else if (node.last_func == fCB)
      {
        signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        signal tgt_sig = tgt.create_or(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = CB({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<Ntk>(tgt_sig, 0xEEEEu, 0u, 2u));
      }
      else if (node.last_func == fXOR)
      {
        signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        signal tgt_sig = tgt.create_xor(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = XOR({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<Ntk>(tgt_sig, 0x6666u, 1u, 2u));
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
  
  tgt.foreach_node([&](auto const & node) 
  {
    tgt.foreach_fanin(node, [&] (auto const & sig )
    {
      if (tgt.is_complemented(sig) && !tgt.is_constant(tgt.get_node( sig )))
      {
        fmt::print("Signal: {0} is complemented. [fanin of {1}]\n", sig.data, node);
      }
    });
  });

  return std::make_pair(tgt, tgt_sig_params);
}

void write_klut_specs(const klut & ntk, const std::unordered_map<klut::signal, Primitive<klut>> & sig_params, const std::string & filename)
{
  std::ofstream spec_file (filename);
  ntk.foreach_pi([&] (const auto & sig) {
    spec_file << fmt::format("PI {0}\n", sig);
  } );
  for (const auto & [sig, prim] : sig_params)
  {
    spec_file << fmt::format("{0},{1},{2},{3}\n", sig, prim.func, prim.type, fmt::join(prim.fanins, "|"));
  }
}

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

int main()  //int argc, char* argv[]
{
  using namespace experiments;
  using namespace mockturtle;

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
  if ( lorina::read_genlib( inputFile, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  // std::unordered_map<std::string, int> nDFF_global = readCSV( NDFF_PATH );
  std::unordered_map<std::string, int> nDFF_global;


  // const std::vector<std::vector<UI>> sets_of_levels { { {0,1,2,3} } };
  const std::string LibEntry_file = "/Users/brainkz/Documents/GitHub/mockturtle/build/LibEntry_2023_06_27_CONNECT_CONSERVATIVE.csv";


  mockturtle::tech_library_params tps;
  // tps.verbose = true;
  tech_library<4, mockturtle::classification_type::p_configurations> tech_lib( gates, tps );

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

    auto [res_no_pb, st_no_pb] = map_wo_pb(aig_original, tech_lib, false, true);

    // if (abc_cec(res_no_pb, benchmark))
    // {
    //   fmt::print("RES_no_pb is equivalent to the original\n");
    // }
    // else
    // {
    //   fmt::print("ERROR: RES_no_pb is not equivalent\n");
    // }

    auto [klut_decomposed, klut_prim_params] = decompose_to_klut(res_no_pb, GNM_global, entries);

    // Here, I need to decompose the network into the macroblocks



    fmt::print("Decomposition complete\n");
    // if (abc_cec_cmp(klut_decomposed, res_no_pb))
    // {
    //   fmt::print("SUCCESS: Decomposed network is equivalent to the original\n");
    // }
    // else
    // {
    //   fmt::print("ERROR: Not equivalent\n");
    // }
    write_klut_specs(klut_decomposed, klut_prim_params, fmt::format("{0}_specs.csv", benchmark));

    // fmt::print("Benchmark \"{0}\"\n", benchmark);
    // fmt::print("\tXMG size before: \"{0}\"\n", xmg_decomposed.size());
    // functional_reduction_params ps;
    // ps.verbose = true;
    // ps.max_TFI_nodes = 10000;
    // ps.skip_fanout_limit = 1000;
    // ps.conflict_limit = 1000;
    // ps.max_clauses = 10000;
    // ps.num_patterns = 1024; //4x default
    // ps.max_patterns = 4096; //4x default
    // fmt::print("Started functional reduction of {}\n", benchmark);
    // functional_reduction(xmg_decomposed, ps);
    // fmt::print("Finished functional reduction of {}\n", benchmark);

    // fmt::print("Started cleanup dangling of {}\n", benchmark);
    // xmg_decomposed = cleanup_dangling( xmg_decomposed );
    // fmt::print("Finished cleanup dangling of {}\n", benchmark);


    // fmt::print("\tXMG size after: \"{0}\"\n", xmg_decomposed.size());

    // auto [res_aig_direct, st_aig_direct, PB_DFF_aig_direct, nJJ_aig_direct, cec_aig_direct] = map_with_pb(benchmark, aig_original, tech_lib, nDFF_global, false);

    // auto [res_xmg_direct, st_xmg_direct, PB_DFF_xmg_direct, nJJ_xmg_direct, cec_xmg_direct] = map_with_pb(benchmark, xmg_original, tech_lib, nDFF_global, false);

    // auto [res_after, st_after, PB_DFF_after, nJJ_after, cec_after] = map_with_pb(benchmark, klut_decomposed, tech_lib, nDFF_global, false);

    // auto t_base = to_seconds( st_aig_direct.time_total );

    // exp( benchmark, 0, //PB_DFF_aig_direct, 
    //                 0,//PB_DFF_xmg_direct/PB_DFF_aig_direct, 
    //                 PB_DFF_after,// /PB_DFF_aig_direct, 

    //                 0,//nJJ_aig_direct, 
    //                 0,//nJJ_xmg_direct/nJJ_aig_direct, 
    //                 nJJ_after,// /nJJ_aig_direct, 

    //                 0,//st_aig_direct.delay, 
    //                 0,//st_xmg_direct.delay/st_aig_direct.delay, 
    //                 st_after.delay,// /st_aig_direct.delay,

    //                 0,//t_base * 1000 , 
    //                 0,//to_seconds( st_xmg_direct.time_total ) / t_base, 
    //                 to_seconds( st_after.time_total )//  / t_base  
    //                 );

    // // mockturtle::write_bench(res_after, benchmark + "_xag_one_round" + ".bench");
    // exp.save();
    // exp.table();
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
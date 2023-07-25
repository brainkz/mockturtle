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
#include <set>
#include <cstdio>

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
  // ps.cut_enumeration_ps.very_verbose = true;
  ps.cut_enumeration_ps.verbose = true;
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
  auto [res, st] = map_wo_pb(input_ntk, tech_lib, area_oriented);
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

constexpr uint8_t PI_GATE = 0u;
constexpr uint8_t PO_GATE = 1u;
constexpr uint8_t AA_GATE = 2u;
constexpr uint8_t AS_GATE = 3u;
constexpr uint8_t SA_GATE = 4u;

const std::vector<std::string> GATE_TYPE { "PI", "PO", "AA", "AS", "SA" }; 

typedef uint32_t glob_phase_t;

template <typename Ntk>
struct Primitive
{
  using ntk_signal = typename Ntk::signal;
  ntk_signal sig; // id of the underlying node in the network
  uint16_t  func; // 4-input function
  uint8_t   type; // AA - 0, AS - 1, SA - 2
  std::vector<ntk_signal> fanins;
  bool is_spl = false;

  Primitive(ntk_signal _sig, uint16_t _func, uint8_t _type, std::vector<ntk_signal> _fanins) : 
                   sig(_sig),    func(_func),   type(_type), fanins(_fanins), is_spl(false) {};

  explicit Primitive(ntk_signal _sig, uint16_t _func, uint8_t _type, std::vector<ntk_signal> _fanins, bool _is_spl) : 
                   sig(_sig),    func(_func),   type(_type), fanins(_fanins), is_spl(_is_spl) {};
  Primitive( const Primitive & other )
  {
    sig = other.sig; // id of the underlying node in the network
    func = other.func; // 4-input function
    type = other.type; 
    fanins = other.fanins;
    is_spl = other.is_spl;
  }
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
    tgt_sig_params.emplace( tgt_pi, Primitive<klut>(tgt_pi, 0xAAAA, AS_GATE, {}));
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
        
        // tgt_sig_params.emplace();
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
        tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0x5555u, AS_GATE, { parent_tgt }));
      }
      else if (node.last_func == fAND)
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_and(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = AND({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0x8888u, SA_GATE, { parent_tgt_1, parent_tgt_2 }));
      }
      else if (node.last_func == fOR) 
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_or(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = OR({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0xEEEEu, SA_GATE, { parent_tgt_1, parent_tgt_2 }));
        OR_replacement_candidates.emplace(tgt_sig);
      }
      else if (node.last_func == fCB)
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_or(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = CB({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0xEEEEu, AA_GATE, { parent_tgt_1, parent_tgt_2 }));
      }
      else if (node.last_func == fXOR)
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_xor(parent_tgt_1, parent_tgt_2);
        // fmt::print("Created node n{0} = XOR({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0x6666u, AS_GATE, { parent_tgt_1, parent_tgt_2 }));
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
    prim.type = AA_GATE;
  }

  src.foreach_po([&](auto const & src_po)
  {
    tgt.create_po(src2tgt.at(src_po));
  } );
  
  return std::make_pair(tgt, tgt_sig_params);
}

// void insert_splitters(klut ntk, const std::unordered_map<klut::signal, Primitive<klut>> & sig_params, const uint8_t n_phases, const bool verbose = false)
// {

// }

std::unordered_map<klut::signal, glob_phase_t> greedy_assemble(const klut & ntk, const std::unordered_map<klut::signal, Primitive<klut>> & sig_params, const uint8_t n_phases, const bool verbose = false)
{
  mockturtle::topo_view<klut> ntk_topo ( ntk );
  std::unordered_map<klut::signal, glob_phase_t> glob_phase;
  glob_phase.reserve(ntk_topo.size());

  ntk_topo.foreach_node([&] ( const klut::signal & node ) 
  {
    if ( ntk_topo.is_constant( node ) )
    {
      if (verbose)
      {
        fmt::print("Skipping constant {}\n", node);
      }
      return;
    }
    if ( ntk_topo.is_pi( node ) )
    {
      glob_phase[node] = 0u;
      if (verbose)
      {
        fmt::print("PI {} placed at ɸ=0 [S=0, φ=0]\n", node);
      }
      return;
    }
    const Primitive<klut> & node_params = sig_params.at( node );
    if (verbose)
    {
      fmt::print("{} GATE {}:\n", GATE_TYPE.at(node_params.type), node);
    }

    if ( node_params.type == AA_GATE )
    {
      // place at the same phase as the latest fanin
      glob_phase_t phase = 0u;
      ntk_topo.foreach_fanin(node, [&] ( const klut::signal & parent )
      {
        if ( ntk_topo.is_constant( parent ) )
        {
          if (verbose)
          {
            fmt::print("\t\tfanin {} is constant, skipping...\n", parent);
          }
          return;
        }
        phase = std::max(phase, glob_phase.at( parent ));
        if (verbose)
        {
          glob_phase_t gp = glob_phase.at( parent );
          fmt::print("\t\tfanin {} ɸ={} [S={}, φ={}]\n", parent, gp, gp/n_phases, gp%n_phases);
        }

      });
      glob_phase[node] = phase;
      if (verbose)
      {
        fmt::print("\tAA GATE {} placed at ɸ={} [S={}, φ={}]\n", node, phase, phase/n_phases, phase%n_phases);
      }
    }
    else if ( node_params.type == SA_GATE )
    {
      // Place at the earliest feasible phase
      glob_phase_t phase = 0u;
      ntk_topo.foreach_fanin(node, [&] (const klut::signal & parent)
      {
        if ( ntk_topo.is_constant( parent ) )
        {
          if (verbose)
          {
            fmt::print("\t\tfanin {} is constant, skipping...\n", parent);
          }
          return;
        }
        // if d==false, SA gate can be directly connected to fanin
        //cannot directly connect PI
        //cannot directly connect split signal 
        //cannot directly connect SA/AA gates to SA gates
        bool d = (ntk_topo.is_pi(parent)) || (ntk_topo.fanout_size(parent) > 1) || (sig_params.at(parent).type != AS_GATE);                     

        phase = std::max(phase, glob_phase.at( parent ) + static_cast<int>(d) );
        if (verbose)
        {
          glob_phase_t gp = glob_phase.at( parent );
          std::string gt;
          fmt::print("\t\tfanin: {4} {0} ɸ={1} [S={2}, φ={3}], fanout={5} \n", 
          parent, 
          gp, 
          gp/n_phases, 
          gp%n_phases, 
          ntk_topo.is_pi(parent) ? "PI" : GATE_TYPE.at( sig_params.at(parent).type ),
          ntk_topo.fanout_size(parent));
        }
      });
      glob_phase[node] = phase;
      if (verbose)
      {
        fmt::print("\tSA GATE {} placed at ɸ={} [S={}, φ={}]\n", node, phase, phase/n_phases, phase%n_phases);
      }
    }
    else if ( node_params.type == AS_GATE )
    {
      // Place at the earliest feasible phase
      glob_phase_t phase = 0u;
      ntk_topo.foreach_fanin(node, [&] (const klut::signal & parent)
      {
        if ( ntk_topo.is_constant( parent ) )
        {
          if (verbose)
          {
            fmt::print("\t\tfanin {} is constant, skipping...\n", parent);
          }
          return;
        }
        phase = std::max(phase, glob_phase.at( parent ) );
        if (verbose)
        {
          glob_phase_t gp = glob_phase.at( parent );
          fmt::print("\t\tfanin {} ɸ={} [S={}, φ={}]\n", parent, gp, gp/n_phases, gp%n_phases);
        }
      });
      phase++;
      glob_phase[node] = phase;
      if (verbose)
      {
        fmt::print("\tAS GATE {} placed at ɸ={} [S={}, φ={}]\n", node, phase, phase/n_phases, phase%n_phases);
      }
    }
    else 
    {
      fmt::print("\t GATE {} : fanins[{}], type[{}], func[{}]\n", node, fmt::join(node_params.fanins, ","), node_params.type, node_params.func);
      throw "Unsupported gate type!";
    }
  });
  return glob_phase;
}

struct Splitter
{
  // the object has an ID value that is a nonexistent node in an ntk
  klut::signal id;
  klut::signal fanin;
  std::vector<klut::signal> fanouts;
  Splitter(const klut::signal & _id, const klut::signal & _fanin, const std::vector<klut::signal>& _fanouts)
    : id(_id), fanin(_fanin), fanouts(_fanouts) {}
};

struct Path
{
  std::set<klut::signal> sources;   // AS/SA gates
  std::set<klut::signal> internals; // AA    gates
  std::set<klut::signal> targets;   // AS/SA gates
  Path(const std::set<klut::signal> & _sources, const std::set<klut::signal>& _internals, const std::set<klut::signal>& _targets)
    : sources(_sources), internals(_internals), targets(_targets) {}
  Path() : sources({}), internals({}), targets({}) {}
  
  void absorb(Path & other)
  {
    sources.insert(other.sources.begin(), other.sources.end());
    internals.insert(other.internals.begin(), other.internals.end());
    targets.insert(other.targets.begin(), other.targets.end());
  }

  void print() const
  {
    fmt::print("Path from [{}]\n\tvia [{}]\n\tto [{}]\n", fmt::join(sources, ","), fmt::join(internals, ","), fmt::join(targets, ","));
  }

  std::vector<klut::signal> preds(const klut::signal & sig, const klut & ntk) const
  {
    if ( sources.count(sig) != 0)
    {
      return {};
    }
    std::vector<klut::signal> predecessors;
    ntk.foreach_fanin(sig, [&] (const klut::signal & parent) {
      if ( internals.count(parent) != 0 || sources.count(parent) != 0 )
      {
        predecessors.push_back( parent );
      }
    });
    return predecessors;
  }

  std::vector<klut::signal> src_int() const
  {
    std::vector<klut::signal> out;
    out.insert(out.end(), sources.begin(), sources.end());
    out.insert(out.end(), internals.begin(), internals.end());
    return out;
  }
  std::vector<klut::signal> int_tgt() const
  {
    std::vector<klut::signal> out;
    out.insert(out.end(), internals.begin(), internals.end());
    out.insert(out.end(), targets.begin(), targets.end());
    return out;
  }

  void dff_cost(klut ntk, std::unordered_map<klut::signal, glob_phase_t> glob_phase)
  {
    std::vector<Splitter> splitters;
    std::vector<std::pair<klut::signal, klut::signal>> leq_constr;
    klut::signal spl_idx; 

    std::unordered_map<klut::signal, std::set<klut::signal>> ntk_fanouts;
    ntk.foreach_node([&](const klut::signal & node)
    {
      spl_idx = std::max(spl_idx, node);
      ntk.foreach_fanin(node, [&](const klut::signal & parent)
      {
        ntk_fanouts[parent].insert(node);
      });
    });
    spl_idx++;


    std::set<klut::signal> src_int = sources; 
    src_int.insert(internals.begin(), internals.end());
    for (auto & sig : src_int)
    {
      std::set<klut::signal> fanouts = ntk_fanouts.at( sig );
      if (fanouts.size() == 0)
      {
        continue;
      }

      if (fanouts.size() == 1)
      {
        // No need to insert splitters, add constraint directly
        leq_constr.emplace_back( sig, *(fanouts.begin()) );
      }

      if (fanouts.size() > 1)
      {
        // Need to insert splitters
        // first, sort fanouts based on glob_phase
        std::vector<klut::signal> fanouts_vec (fanouts.begin(), fanouts.end());
        std::sort(fanouts_vec.begin(), fanouts_vec.end(), 
        [&] (const klut::signal & sig1, const klut::signal & sig2)
        {
          return glob_phase.at( sig1 ) < glob_phase.at( sig2 );
        });
        
        // the splitters are placed between each pair of the fanouts

        klut::signal last_node = sig;
        uint16_t fo_size = fanouts_vec.size();

        auto nspl = fanouts.size() - 1;

        // Names of splitters
        std::vector<klut::signal> spl_ids ( nspl );
        std::iota(spl_ids.begin(), spl_ids.end(), spl_idx);

        // Fanins are the first node and the splitters up to pre-last  
        std::vector<klut::signal> spl_fi;
        spl_fi.reserve( nspl );
        spl_fi.push_back( sig );
        spl_fi.insert(spl_fi.end(), spl_ids.begin(), spl_ids.end()-1);

        // Fanout 1 are the splitters from second to last, and the last fanout
        std::vector<klut::signal> spl_fo_1;
        spl_fo_1.reserve( nspl );
        spl_fo_1.insert(spl_fo_1.end(), spl_ids.begin()+1, spl_ids.end());
        spl_fo_1.push_back( fanouts_vec.back() );
        
        // Fanout 2 are the fanouts from first to pre-last
        std::vector<klut::signal> spl_fo_2;
        spl_fo_2.reserve( nspl );
        spl_fo_2.insert(spl_fo_2.end(), fanouts_vec.begin(), fanouts_vec.end()-1);
        
        leq_constr.emplace_back( sig , spl_ids.front() );

        for (auto i = 0u; i < nspl; ++i)
        {
          Splitter new_spl { spl_ids[i], spl_fi[i], { spl_fo_1[i], spl_fo_2[i] } };
          splitters.emplace_back( new_spl );
          leq_constr.emplace_back( spl_ids[i] , spl_fo_1[i] );
          leq_constr.emplace_back( spl_ids[i] , spl_fo_2[i] );
        }
      }
    }
    
  }
};

template <typename sig>
bool haveCommonElements(const std::set<sig>& set1, const std::set<sig>& set2) 
{
    for (const auto& element : set1) 
    {
        if (set2.count(element) > 0) 
        {
            return true;  // Found a common element
        }
    }
    return false;  // No common elements found
}

std::vector<Path> extract_paths(const klut & ntk, const std::unordered_map<klut::signal, Primitive<klut>> & sig_params, bool verbose = false)
{
  // std::unordered_map<klut::signal, std::vector<klut::signal>> ancestors;
  std::vector<Path> paths;
  ntk.foreach_gate([&] (const klut::signal & node )
  {
    if (ntk.is_constant( node ) || ntk.is_pi( node ) )
    {
      return;
    }

    const Primitive<klut> & params = sig_params.at( node );
    if (params.type != AS_GATE && params.type != SA_GATE)
    {
      return;
    }
    // at this point, the node should be a AS/SA gate

    // Create a separate path for each fanin of the node
    std::vector<Path> aa_paths;
    aa_paths.reserve( ntk.fanin_size(node) );
    ntk.foreach_fanin( node, [&](const klut::signal & parent) 
    { 
      // Create a path object with only a target
      Path node_path;
      node_path.targets.emplace( node );
      
      std::vector<klut::signal> stack;
      stack.push_back( parent );
      
      fmt::print("Node {}, traversal for parent {}\n", node, parent);

      std::set<klut::signal> seen;
      while (!stack.empty())
      {
        klut::signal & n = stack.back();
        stack.pop_back();

        fmt::print("\t[Parent {}]: Analyzing node {}\n", parent, n);
        // A constant does not have any effect on the DFF placement, we can skip it
        if ( ntk.is_constant( n ) )
        {
          continue;
        }
        // Found a source of the path, add to sources, do not continue traversal
        else if ( ntk.is_pi( n ) || sig_params.at( n ).type == AS_GATE || sig_params.at( n ).type == SA_GATE )
        {
          fmt::print("\t\t{} is a source node \n", n);
          node_path.sources.emplace( n );
        }
        // Found AA gate, add to internal nodes, add parents for further traversal
        else if ( sig_params.at( n ).type == AA_GATE )
        {
          fmt::print("\t\t{} is INTERNAL with function {}, adding fanins \n", n, kitty::to_hex(ntk.node_function(n)) );
          node_path.internals.emplace( n );
          ntk.foreach_fanin( n , [&](const klut::signal & sig )
          {
            stack.push_back( sig );
            // if (std::find( seen.begin(), seen.end(), sig ) != seen.end() )
            // {
            //   fmt::print("\t\t\tAdded {} to stack \n", sig);
            // }
          });
        }
        else
        {
          fmt::print("Signal {}: {} is not recognized \n", n, GATE_TYPE.at( sig_params.at( n ).type ));
          throw "Unsupported case";
        }
        seen.emplace( n );
      }
      // aa_paths.push_back(node_path);
      // Identify overlapping paths
      std::vector<size_t> to_merge;
      for (size_t i = 0u; i < paths.size(); ++i)
      {
        Path & known_paths = paths[i];
        // merge if there are sources in common
        if( haveCommonElements( known_paths.sources, node_path.sources) )
        {
          to_merge.push_back(i);
        }
      }
      // Merge overlapping paths into the path object and remove the path object
      // Iterating in reverse order to ensure seamless deletion
      for (auto it = to_merge.rbegin(); it != to_merge.rend(); ++it) 
      {
        auto idx = *it;
        // fmt::print("Before absorption\n");
        // node_path.print();
        node_path.absorb(paths[idx]);
        // fmt::print("After absorption\n");
        // node_path.print();
        paths.erase(paths.begin() + idx);
      }
      paths.push_back(node_path);
    });
  });
  return paths;
}

uint64_t calculate_hash(glob_phase_t phase, klut::signal fanout, klut::signal fanin)
{
    std::size_t seed = 0;
    hash_combine(seed, phase);
    hash_combine(seed, fanout);
    hash_combine(seed, fanin);
    return seed;
};

struct Chain
{
  glob_phase_t min_phase;
  glob_phase_t max_phase;
  klut::signal fanout; // successor gate (not DFF!, can even be a splitter)
  klut::signal fanin;  // predecessor gate (not DFF!, can even be a splitter)
  uint64_t hash;
  std::vector<uint64_t> parent_hashes;
  Chain(glob_phase_t _min_phase, glob_phase_t _max_phase, klut::signal _fanout, klut::signal _fanin, uint64_t _hash, const std::vector<uint64_t>& _parent_hashes)
    : min_phase(_min_phase), max_phase(_max_phase), fanout(_fanout), fanin(_fanin), hash(_hash), parent_hashes(_parent_hashes)
  {}
  Chain(glob_phase_t _min_phase, glob_phase_t _max_phase, klut::signal _fanout, klut::signal _fanin)
    : min_phase(_min_phase), max_phase(_max_phase), fanout(_fanout), fanin(_fanin), parent_hashes({})
  {
    std::size_t seed = 0;
    hash_combine(seed, _fanout);
    hash_combine(seed, _fanin);
    hash = seed;
  }
  void print() const
  {
    fmt::print("DFF chain between {} and {} at phases from {} to {}\n", fanin, fanout, min_phase, max_phase);
  }
};

struct PB_DFF
{
  glob_phase_t phase;
  klut::signal fanout; // successor gate (not DFF!, can even be a splitter)
  klut::signal fanin;  // predecessor gate (not DFF!, can even be a splitter)
  uint64_t hash;
  std::vector<std::tuple<glob_phase_t, klut::signal, klut::signal>> parent_dffs;

  explicit PB_DFF (glob_phase_t _phase, klut::signal _fanout, klut::signal _fanin) 
  : phase(_phase), fanout(_fanout), fanin(_fanin), parent_dffs({}) {
    hash = calculate_hash(_phase, _fanout, _fanin);
  };
  void print() const
  {
    fmt::print("DFF location between {} and {} at phase {}\n", fanin, fanout, phase);
  }
};

struct Chain_REGISTRY
{
  std::unordered_map<uint64_t, Chain> reg;
  // std::unordered_map<uint64_t, uint64_t> adjacent_hashes;

  const Chain & at( klut::signal _fanout, klut::signal _fanin ) const
  {
    std::size_t hash = 0;
    hash_combine(hash, _fanout);
    hash_combine(hash, _fanin);
    return reg.at( hash );
  }
  const Chain & at( uint64_t hash ) const
  {
    return reg.at( hash );
  }
  uint64_t add( glob_phase_t _min_phase, glob_phase_t _max_phase, klut::signal _fanout, klut::signal _fanin )
  {
    Chain chain { _min_phase, _max_phase, _fanout, _fanin };
    reg.emplace( chain.hash, chain );
    return chain.hash;
  }
  void print() const
  {
    for (const auto & [hash, chain] : reg)
    {
      fmt::print("{} : ", hash);
      chain.print();
    }
  }
};

void insert_splitters( klut & ntk, klut::node const& node, std::unordered_map<klut::signal, glob_phase_t> & glob_phase, std::unordered_map<mockturtle::klut_network::signal, Primitive<klut>> & klut_prim_params, std::vector<klut::signal> & fanouts)
{
  auto buf_func = kitty::dynamic_truth_table( 1u );
  buf_func._bits[0] = 0x2;

  unsigned int fo_size = fanouts.size();
  int fanout_diff = ntk.fanout_size(node) - fo_size;
  fmt::print("\t\tReported fanout is greater by {} \n", fanout_diff);

  fmt::print("\t\tFound a total of {} fanouts \n", fo_size);

  // Create splitter nodes
  std::vector<klut::signal> splitters;
  splitters.reserve( fo_size );
  klut::signal last_spl_fanin = node;
  for (auto i = 0u; i < fo_size-1; ++i)
  {
    klut::signal spl = ntk.create_node({ last_spl_fanin }, buf_func);
    assert( (spl != node) && (spl != last_spl_fanin) );
    fmt::print("\t\tConnected splitter #{} to {}\n", spl, last_spl_fanin);
    splitters.push_back( spl );
    last_spl_fanin = spl;

    auto obj = Primitive<klut>(spl, 0xAAAA, AA_GATE, {last_spl_fanin});
    obj.is_spl = true;
    // (spl, 0xAAAA, AA_GATE, { last_spl_fanin }, true );
    klut_prim_params.emplace(spl, obj);
  }
  // add the last splitter twice to connect the two last fanouts
  splitters.push_back( last_spl_fanin );

  std::sort(fanouts.begin(), fanouts.end(), [&](const klut::signal & a, const klut::signal & b)
  {
    if (ntk.is_constant( a ))
    {
      return false;
    }
    else if (ntk.is_constant( b ))
    {
      return true;
    }

    return (glob_phase.at( a ) - (klut_prim_params.at( a ).type == SA_GATE)) < (glob_phase.at( b ) - (klut_prim_params.at( b ).type == SA_GATE));
  });
  fmt::print("\t\tSorting fanouts -> {}\n", fmt::join(fanouts, ","));

  // redirect the fanouts to splitters
  for ( auto i = 0u; i < fo_size; ++i )
  {
    klut::signal fo = fanouts[i];
    auto& fo_cfg = ntk._storage->nodes[fo];
    fmt::print("\t\t\tFanout: {}\n", fo);

    klut::signal spl = splitters[i];
    auto& spl_cfg = ntk._storage->nodes[spl];
    fmt::print("\t\t\tSplitter: {}\n", spl);

    for ( auto& child : fo_cfg.children )
    {
      if ( child == node )
      {
        fmt::print("BEFORE:\tnode {} fo={}\tspl {} fo={}\n", node, static_cast<int>(ntk._storage->nodes[node].data[0].h1), spl, static_cast<int>(ntk._storage->nodes[spl].data[0].h1));
        
        fmt::print("\t\t\tAdded {} as a fanin to {}\n", spl, fo);
        child = spl;
        // increment fan-out of splitter
        ntk._storage->nodes[spl].data[0].h1++;
        // decrement fan-out of old node
        ntk._storage->nodes[node].data[0].h1--;
        fmt::print("AFTER :\tnode {} fo={}\tspl {} fo={}\n", node, static_cast<int>(ntk._storage->nodes[node].data[0].h1), spl, static_cast<int>(ntk._storage->nodes[spl].data[0].h1));
      }
    }
  }
  for ( auto i = 0u; i < fo_size-1; ++i )
  {
    klut::signal fo = fanouts[i];
    klut::signal spl = splitters[i];
    // assign splitter's phase to be equal to the fanout's phase
    fmt::print("\t\t\tAssigned phase {} to splitter {}\n", glob_phase.at( fo ), spl);
    glob_phase[spl] = glob_phase.at( fo );
  }

  fmt::print("\t\t\tFinal node {} fanout: {}\n", node, ntk.fanout_size( node ));
  assert((ntk.fanout_size( node ) - fanout_diff) == 1);
  for ( const auto & spl : splitters )
  {
    fmt::print("\t\t\tFinal spl {} fanout: {}\n", spl, ntk.fanout_size( spl ));
    assert(ntk.fanout_size( spl ) == 2);
  }
  // return std::make_tuple( ntk, glob_phase, klut_prim_params );
}

/// @brief Fills DFFs along the path
/// @param ntk - decomposed klut network
/// @param path - Path object 
/// @param glob_phase - contains information about all phases in network
/// @param klut_prim_params - has gate parameters for each gate 
/// @param n_phases - number of phases in the network
Chain_REGISTRY fill_dffs( const klut & ntk, const Path & path, const std::unordered_map<klut::signal, glob_phase_t> & glob_phase, const std::unordered_map<klut::signal, Primitive<klut>> & klut_prim_params)
{
  Chain_REGISTRY REG;

  std::vector<std::pair<klut::signal, klut::signal>> stack;
  for (const klut::signal & fanout: path.targets)
  {
    for (const klut::signal & fanin : path.preds(fanout, ntk))
    {
      stack.emplace_back(fanout, fanin);
    }
  }
  std::set<std::pair<klut::signal, klut::signal>> seen;

  while (!stack.empty())
  {
    auto [fanout, fanin] = stack.back();
    stack.pop_back();

    const bool fanin_is_AA = (klut_prim_params.at( fanin ).type == AA_GATE);
    const glob_phase_t fanin_phase = glob_phase.at( fanin ); //  + fanin_is_AA;

    const bool fanout_is_AS = ( klut_prim_params.at( fanout ).type == AS_GATE );
    const glob_phase_t fanout_phase = glob_phase.at( fanout ) - fanout_is_AS;

    REG.add(fanin_phase, fanout_phase, fanout, fanin);

    if ( fanin_is_AA && seen.count(std::make_pair(fanout, fanin)) == 0 ) 
    {
      for (const klut::signal & predecessor : path.preds(fanin, ntk))
      {
        stack.emplace_back(fanin, predecessor);
      }
    }
    seen.emplace(fanout, fanin);
  }
  // std::vector<klut::signal> stack;
  // stack.insert(stack.end(), path.targets.begin(), path.targets.end());
  // std::set<klut::signal> seen;
  // while (!stack.empty())
  // {
  //   klut::signal fanout = stack.back();
  //   stack.pop_back();
  //   bool fanout_is_AS = klut_prim_params.at( fanout ).type == AS_GATE;
  //   glob_phase_t fanout_phase = glob_phase.at( fanout ) - fanout_is_AS;

  //   for (const klut::signal & fanin : path.preds(fanout, ntk))
  //   {
  //     bool fanin_is_AA = (klut_prim_params.at( fanin ).type == AA_GATE);
  //     bool fanout_is_AS = (klut_prim_params.at( fanout ).type == AS_GATE);
  //     glob_phase_t fanin_phase = glob_phase.at( fanin ) + fanin_is_AA;
  //     glob_phase_t fanout_phase = glob_phase.at( fanout ) - fanout_is_AS;

  //     REG.add(fanin_phase, fanout_phase, fanout, fanin);

  //     // if there are unseen chains in the back
  //     if ( klut_prim_params.at( fanin ).type == AA_GATE && seen.count(fanin) == 0 ) 
  //     {
  //       stack.push_back(fanin);
  //     }
  //   }
  //   seen.emplace(fanout);
  // }


  // // at this point, all dffs are registered and internal connections are made
  // // now, tail nodes need to be connected to preceding chains
  // for (const auto & [phase, fanout, mid] : tails)
  // {
  //   PB_DFF & tail_dff = REG.at(phase, fanout, mid);
  //   for (const klut::signal & fanin : path.preds(mid, ntk))
  //   {
  //     // PB_DFF & head_dff = REG.at( phase, mid, fanin );
  //     tail_dff.parent_dffs.emplace_back( phase, mid, fanin );
  //   }
  // }
  REG.print();
  return REG;
  // TODO : extract chains and formulate SAT (make sure to consider special cases, like a phase with AA->AS connection)
}
#if false
struct DFF_var 
{
  klut::signal fanin;
  klut::signal fanout;
  glob_phase_t phase;
  bool negated = false;

  DFF_var(klut::signal _fanin, klut::signal _fanout, glob_phase_t _phase)
      : fanin(_fanin), fanout(_fanout), phase(_phase) {}

  DFF_var(klut::signal _fanin, klut::signal _fanout, glob_phase_t _phase, bool _negated)
      : fanin(_fanin), fanout(_fanout), phase(_phase), negated(_negated) {}

  DFF_var(const DFF_var& other)
      : fanin(other.fanin), fanout(other.fanout), phase(other.phase), negated(other.negated) {}

  std::string str()
  {
    if (negated)
    {
      return fmt::format( "var_{}_{}_{}.Not()", fanin, fanout, phase );
    }
    return fmt::format( "var_{}_{}_{}", fanin, fanout, phase );
  }
};
#endif
struct DFF_var 
{
  klut::signal fanin;
  klut::signal fanout;
  glob_phase_t phase;
  std::vector<uint64_t> parent_hashes;
  uint64_t hash;

  DFF_var(klut::signal _fanin, klut::signal _fanout, glob_phase_t _phase)
      : fanin(_fanin), fanout(_fanout), phase(_phase), parent_hashes({}) 
      {
        hash = calculate_hash(phase, fanout, fanin);
      }

  DFF_var(klut::signal _fanin, klut::signal _fanout, glob_phase_t _phase, std::vector<uint64_t> _parent_hashes, uint64_t _hash)
      : fanin(_fanin), fanout(_fanout), phase(_phase), parent_hashes(_parent_hashes), hash(_hash) {}

  DFF_var(const DFF_var& other)
      : fanin(other.fanin), fanout(other.fanout), phase(other.phase),parent_hashes(other.parent_hashes), hash(other.hash) {}

  std::string str(bool negated = false)
  {
    if (negated)
    {
      return fmt::format( "var_{}_{}_{}.Not()", fanin, fanout, phase );
    }
    return fmt::format( "var_{}_{}_{}", fanin, fanout, phase );
  }
};

struct DFF_registry
{
  std::unordered_map<uint64_t, DFF_var> variables;

  DFF_var & at(klut::signal _fanin, klut::signal _fanout, glob_phase_t _phase)
  {
    return variables.at( calculate_hash(_phase, _fanout, _fanin) );
  } 
  DFF_var & at(uint64_t _hash)
  {
    return variables.at( _hash );
  } 
  uint64_t add(klut::signal _fanin, klut::signal _fanout, glob_phase_t _phase, std::vector<uint64_t> _parent_hashes = {})
  {
    uint64_t _hash = calculate_hash(_phase, _fanout, _fanin);
    DFF_var temp { _fanin, _fanout, _phase, _parent_hashes, _hash };
    // variables.emplace(_hash, std::move(temp));
    variables.insert(std::make_pair(_hash, temp));
    return _hash;
  } 
  std::vector<uint64_t> from_chain(const Chain & chain, const std::unordered_map<klut::signal, glob_phase_t> & glob_phase, const std::unordered_map<mockturtle::klut_network::signal, Primitive<klut>> & klut_prim_params, uint64_t & precalc_ndff)
  {
    const glob_phase_t fanin_phase = glob_phase.at( chain.fanin );
    const uint8_t fanin_type = klut_prim_params.at( chain.fanin ).type;

    const glob_phase_t fanout_phase = glob_phase.at( chain.fanout );
    const uint8_t fanout_type = klut_prim_params.at( chain.fanout ).type;

    // check if the chain is straight - #DFF is just floor(delta-phase)
    if (fanout_type != AA_GATE && fanin_type != AA_GATE)
    {
      if (fanout_phase == fanin_phase)
      {
        // special case when an AS gate feeds directly into SA gate
        assert(fanin_type == AS_GATE && fanout_type == SA_GATE);
        // do nothing, no additional DFFs needed
      }
      else
      {
        // straight chain, just floor the difference!
        precalc_ndff += (glob_phase.at( chain.fanout ) - glob_phase.at( chain.fanin ))/3 + (fanout_type == SA_GATE); //extra DFF before SA gate
      }
      return {};
    }
    // create the earliest DFF
    glob_phase_t earliest_phase = fanin_phase + (fanin_type != AA_GATE);
    glob_phase_t latest_phase = fanout_phase - (fanout_type == AS_GATE);
    std::vector<uint64_t> out_hashes;
    uint64_t last_hash = add(chain.fanin, chain.fanout, earliest_phase);
    out_hashes.push_back(last_hash);
    for (glob_phase_t phase = fanin_phase+1; phase <= latest_phase; ++phase)
    {
      last_hash = add(chain.fanin, chain.fanout, phase, {last_hash});
      out_hashes.push_back(last_hash);
    }
    return out_hashes;
  }
  std::vector<uint64_t> preds_at_phase(uint64_t hash, glob_phase_t phase)
  {
    std::vector<uint64_t> out;
    preds_at_phase(hash, phase, out);
    return out;
  }

  void preds_at_phase(uint64_t hash, glob_phase_t phase, std::vector<uint64_t> & out)
  {
    DFF_var & dff = variables.at( hash );
    if (dff.phase < phase)
    {
      return;
    }
    for (const uint64_t & phash: dff.parent_hashes)
    {
      const DFF_var & p = variables.at( phash );
      if (p.phase == phase)
      {
        out.push_back(p.hash);
        preds_at_phase(p.hash, phase, out);
      }
    }
  }
  std::string str(uint64_t hash, bool negated = false)
  {
    DFF_var & dff = variables.at(hash);
    if (negated)
    {
      return fmt::format( "var_{}_{}_{}.Not()", dff.fanin, dff.fanout, dff.phase );
    }
    return fmt::format( "var_{}_{}_{}", dff.fanin, dff.fanout, dff.phase );
  }
};

std::pair<DFF_registry, uint64_t> generate_dff_vars(const klut & ntk, const Chain_REGISTRY & REG, const Path & path, const std::unordered_map<klut::signal, glob_phase_t> & glob_phase, const std::unordered_map<mockturtle::klut_network::signal, Primitive<klut>> & klut_prim_params)
{
  std::unordered_map<klut::signal, std::vector<uint64_t>> node2chains;
  for (const auto & [hash, chain] : REG.reg)
  {
    node2chains[chain.fanout].push_back(hash);
  }
  
  DFF_registry DFF_REG;

  uint64_t precalc_ndff = 0; //precalculated number of DFFs for trivial cases

  //map from fanin 
  std::unordered_map<klut::signal, std::unordered_set<uint64_t>> earliest_dffs; 
  std::unordered_map<klut::signal, std::unordered_set<uint64_t>> latest_dffs;
  std::vector<uint64_t> stack;
  std::unordered_set<uint64_t> seen;
  for (const klut::signal & fanout : path.targets)
  {
    for (const uint64_t chain_hash : node2chains[ fanout ]) //changed to brackets
    {
      stack.push_back(chain_hash);
    }
  }
  while (!stack.empty())
  {
    uint64_t chain_hash = stack.back();
    stack.pop_back();
    seen.emplace(chain_hash);

    const Chain & chain = REG.at( chain_hash );
    //create local dffs
    fmt::print("[generate_dff_vars]: Adding hashes from chain ");
    chain.print();
    fmt::print("[generate_dff_vars]: Old DFF_REG.size = {}\n", DFF_REG.variables.size());
    std::vector<uint64_t> new_hashes = DFF_REG.from_chain(chain, glob_phase, klut_prim_params, precalc_ndff);
    fmt::print("[generate_dff_vars]: Added {} new hashes\n", new_hashes.size());
    fmt::print("[generate_dff_vars]: NEW DFF_REG.size = {}\n", DFF_REG.variables.size());
    if (new_hashes.empty()) 
    {
      continue;
    }
    earliest_dffs[chain.fanin].insert( new_hashes.front() );
    latest_dffs[chain.fanout].insert( new_hashes.back() );
    for (const uint64_t parent_chain_hash : node2chains[chain.fanin])
    {
      // add the new chain to stack for further processing 
      if (seen.count(parent_chain_hash) == 0)
      {
        stack.push_back(parent_chain_hash);
      }
    }
  }

  for (const klut::signal & mid : path.internals)
  {
    for ( const uint64_t & hash_after_mid : earliest_dffs[ mid ] )
    {
      DFF_var & dff_after_mid = DFF_REG.at( hash_after_mid );
      for (const uint64_t & hash_before_mid : latest_dffs[ mid ] )
      {
        // DFF_var & dff_before_mid = DFF_REG.at( hash_before_mid );
        dff_after_mid.parent_hashes.push_back( hash_before_mid );
      }
    }
  }
  for (const auto & [hash, dff]: DFF_REG.variables)
  {
    if (dff.parent_hashes.empty())
    {
      auto it = std::find(path.sources.begin(), path.sources.end(), dff.fanin);
      assert(it != path.sources.end());
    }
  }

  return std::make_pair(DFF_REG, precalc_ndff);
}

std::vector<std::vector<uint64_t>> worm_snake_2(DFF_registry & DFF_REG, const std::unordered_map<klut::signal, glob_phase_t> & glob_phase, const Path & path, const uint8_t n_phases)
{
  fmt::print("Starting extraction of worms \n");
  std::vector<std::vector<uint64_t>> worms;
  std::vector<std::vector<uint64_t>> stack;
  std::vector<std::vector<uint64_t>> seen;

  // Add variables near the targets
  for (const auto & [hash, dff]: DFF_REG.variables)
  {
    auto it = std::find(path.targets.begin(), path.targets.end(), dff.fanout);
    if (it != path.targets.end() && (dff.phase == glob_phase.at( dff.fanout )))
    {
      stack.push_back( { { hash } } );
    }
  }

}


std::vector<std::vector<uint64_t>> worm_snake(DFF_registry & DFF_REG, const std::unordered_map<klut::signal, glob_phase_t> & glob_phase, const Path & path, const uint8_t n_phases, const std::string cfg_name)
{
  fmt::print("[{}]: Starting extraction of worms \n", cfg_name);
  std::vector<std::vector<uint64_t>> worms;
  std::vector<std::vector<uint64_t>> stack;
  std::vector<uint64_t> pi_dffs;
  std::vector<std::vector<uint64_t>> seen;
  for (const auto & [hash, dff]: DFF_REG.variables)
  {
    auto it = std::find(path.targets.begin(), path.targets.end(), dff.fanout);
    if (it != path.targets.end() && (dff.phase == glob_phase.at( dff.fanout )))
    {
      stack.push_back( { { hash } } );
    }
  }
  fmt::print("Stack size is {} \n", stack.size());
  while (!stack.empty())
  {
    std::vector<uint64_t> worm = stack.back();
    stack.pop_back();
    seen.push_back(worm);

    // DFF_var & worm_tail = DFF_REG.at( worm.front() );
    DFF_var & worm_head = DFF_REG.at( worm.back() );
    glob_phase_t old_head_phase = worm_head.phase;
    DFF_var & worm_tail = DFF_REG.at( worm.front() );
    glob_phase_t old_tail_phase = worm_tail.phase;
    bool is_full_diff = ((old_tail_phase - old_head_phase) == (n_phases - 1));

    fmt::print("\tCurrent worm size {}, between phases {} and {} \n", worm.size(), DFF_REG.str(worm.front()), DFF_REG.str(worm.back()));

    bool the_worm_is_good = worm_head.parent_hashes.empty() && is_full_diff;

    for (const uint64_t phash : worm_head.parent_hashes)
    {
      DFF_var & new_head = DFF_REG.at( phash );

      glob_phase_t head_phase = new_head.phase;
      glob_phase_t tail_phase = head_phase + n_phases - 1;
      
      the_worm_is_good |= (head_phase < old_head_phase && is_full_diff );

      std::vector<uint64_t> new_worm;
      for (auto it = worm.begin(); it < worm.end(); ++it)
      {
        DFF_var & dff = DFF_REG.at( *it );
        if (dff.phase <= tail_phase)
        {
          new_worm.insert(new_worm.end(), it, worm.end());
          fmt::print("\tAdding new worm with size {}, between phases {} and {} \n", worm.size(), dff.phase, head_phase);
          break;
        }
      }
      new_worm.push_back( phash );
      
      auto it = std::find(seen.begin(), seen.end(), new_worm);
      if (it == seen.end())
      {
        stack.push_back( new_worm );
      }
    }
    if (the_worm_is_good)
    {
      fmt::print("\tAdding the worm with size {}, between phases {} and {} to the output vector\n", worm.size(), old_tail_phase, old_head_phase);
      worms.push_back(worm);
    }
  }
  if (worms.size() > 0)
  {
    std::ofstream spec_file (cfg_name);

    std::vector<std::vector<uint64_t>> single_phase_conflicts;
    for (const auto & worm : worms)
    {
      std::map<glob_phase_t, std::vector<uint64_t>> buckets;
      for (auto it = worm.begin(); it < worm.end(); ++it)
      {
        DFF_var & dff = DFF_REG.at(*it);
        buckets[dff.phase].push_back(*it);
      }

      std::vector<std::string> vars_bucket;
      for (const auto & [phase, bucket] : buckets)
      {
        std::vector<std::string> vars;
        for (uint64_t hash : bucket)
        {
          vars.push_back(DFF_REG.str(hash));
        }
        vars_bucket.push_back(fmt::format(vars.size()>1?"({})":"{}", fmt::join(vars, "+")));

        auto it = std::find(single_phase_conflicts.begin(), single_phase_conflicts.end(), bucket);
        if (it == single_phase_conflicts.end() && vars.size() > 1)
        {
          single_phase_conflicts.push_back(bucket);
          fmt::print("New single phase conflict : {}≤1\n", fmt::join(vars, "+"));
          spec_file << fmt::format("PHASE,{}\n", fmt::join(vars, ","));
        }
      }
      std::reverse(vars_bucket.begin(), vars_bucket.end());
      fmt::print("New buffer requirement : ~({})\n", fmt::join(vars_bucket, "|"));
      spec_file << fmt::format("BUFFER,{}\n", fmt::join(vars_bucket, ","));
    }
  }

  return worms;
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


  // TODO: 
  // TODO : Here I need to replace each binding with the appropriate primitives.
  // TODO : record the timing constraints for the ILP
  // TODO : record the mapping from src to tgt network
  // TODO : record the data pertaining to each node :
  //        - whether the element is AA, AS, or SA
  //          - AA elements are placed at the phase of the latest input
  //          - SA elements tie the preceding AS elements to itself to ensure simultaneous arrival of pulses  
  //        - anything else???

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

int cpsat_ortools(const std::string & cfg_name) 
{
  std::string command = fmt::format("/Users/brainkz/anaconda3/bin/python /Users/brainkz/Documents/GitHub/ortools_python/config_solver.py {}", cfg_name);
  std::string pattern = "Objective value: (\\d+)";

  // Run the command and capture its output
  FILE* pipe = popen(command.c_str(), "r");
  if (!pipe) 
  {
    std::cerr << "Error running the command." << std::endl;
    return -1;
  }

  char buffer[128];
  std::string output;
  while (fgets(buffer, sizeof(buffer), pipe) != nullptr) 
  {
    output += buffer;
  }
  fmt::print(output);

  int result = pclose(pipe);
  if (result == -1) 
  {
    std::cerr << "Error closing the command pipe." << std::endl;
    return -1;
  }

  // Use regex to find the objective value in the output
  std::regex regex(pattern);
  std::smatch match;
  if (std::regex_search(output, match, regex) && match.size() > 1) 
  {
    std::string value_str = match[1];
    return std::stoi(value_str);
  } 
  else 
  {
    std::cerr << "Objective value not found in the output." << std::endl;
    return -1;
  }
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

constexpr uint8_t n_phases = 4u;

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

  // auto benchmarks1 = epfl_benchmarks( experiments::adder | experiments::sin | experiments::cavlc | experiments::int2float | experiments::priority | experiments::i2c | experiments::voter | experiments::dec );
  // experiments::adder | experiments::sin | experiments::cavlc | | experiments::dec | experiments::i2c | experiments::adder | experiments::int2float | | experiments::voter | experiments::c432 | experiments::c499 | | experiments::priority
  auto benchmarks1 = epfl_benchmarks(  experiments::dec  );
  // //   auto benchmarks1 = epfl_benchmarks( experiments::epfl & ~experiments::div & ~experiments::hyp & ~experiments::log2 & ~experiments::sqrt );
  // experiments::c880 | experiments::c1908 | experiments::c3540 |
  auto benchmarks2 = iscas_benchmarks(   experiments::c5315 | experiments::c7552 );

  benchmarks1.insert(benchmarks1.end(), benchmarks2.begin(), benchmarks2.end());



  fmt::print("{}\n", fmt::join(benchmarks1, "\n"));

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
      fmt::print("Failed to read {}\n", benchmark);
      continue;
    }
    fmt::print("Started mapping of {}\n", benchmark);
    // auto [res_no_pb, st_no_pb] = 
    auto [res_w_pb, st_w_pb, total_ndff_w_pb, total_area_w_pb, cec_w_pb ] = map_with_pb(benchmark, aig_original, tech_lib, nDFF_global, false); //, true
    fmt::print("Finished mapping of {}\n", benchmark);

    auto [klut_decomposed, klut_prim_params] = decompose_to_klut(res_w_pb, GNM_global, entries);
    fmt::print("Decomposition complete\n");

    std::unordered_map<klut::signal, glob_phase_t> glob_phase = greedy_assemble(klut_decomposed, klut_prim_params, n_phases, true);

    fmt::print("[i] FINISHED PHASE ASSIGNMENT\n");

    klut spl_klut = klut_decomposed;
    std::unordered_map<klut::signal, glob_phase_t> spl_phase = glob_phase;
    std::unordered_map<mockturtle::klut_network::signal, Primitive<klut>> spl_params = klut_prim_params;

    std::map<klut::signal, unsigned int> multifanouts;
    spl_klut.foreach_node([&] (const klut::signal & n)
    {
      if ( !(spl_klut.is_constant( n )) && !( spl_params.at( n ).is_spl ) && spl_klut.fanout_size( n ) > 1)
      {
        multifanouts.emplace(n, spl_klut.fanout_size( n ));
      }
    });

    auto size_before = spl_klut.size();

    // Find fanouts of each node
    std::map<klut::signal, std::vector<klut::signal>> fanout_map;
    spl_klut.foreach_node( [&](const klut::signal & n)
    {
      if (spl_klut.is_constant(n))
      {
        fmt::print("\t\tNot considering fanouts of constants {}\n", n);
        return;
      }
      spl_klut.foreach_fanin(n, [&](const klut::signal & p)
      {
        fmt::print("\t\tAdding {} to fanouts of {} \n", n, p);
        fanout_map[p].push_back( n );
      });
    });

    auto ctr = 0u;
    for (auto & [n, fo_size] : multifanouts)
    {
      fmt::print("Splitter insertion #{} out of {}\n", ++ctr, multifanouts.size());
      insert_splitters(spl_klut, n, spl_phase, spl_params, fanout_map[n]);
    }

    std::map<klut::signal, unsigned int> multifanouts_after;
    spl_klut.foreach_node([&] (const klut::signal & n)
    {
      if ( !(spl_klut.is_constant( n )) && !( spl_params.at( n ).is_spl ) && spl_klut.fanout_size( n ) > 1)
      {
        multifanouts_after.emplace(n, spl_klut.fanout_size( n ));
      }
    });

    bool spl_cec = experiments::abc_cec( spl_klut, benchmark );
    if (!(spl_cec))
    {
      throw "Networks are not equivalent\n";
    }

    if (multifanouts_after.size() > 0)
    {
      for (auto & [k,v]:multifanouts_after)
      {
        fmt::print("Warning - reported fanout of node {} is {}\n", k, v);
      }
      // throw "Splitters not inserted\n";
    }

    std::vector<Path> paths = extract_paths(spl_klut, spl_params);

    uint64_t total_num_dff = 0u;
    auto file_ctr = 0u;
    for (const Path & path : paths)
    {
      path.print();
      Chain_REGISTRY path_reg = fill_dffs( spl_klut, path, spl_phase, spl_params );

      auto [DFF_REG, precalc_ndff] = generate_dff_vars(spl_klut, path_reg, path, spl_phase, spl_params);
      fmt::print("Precalculated {} DFFs!\n", precalc_ndff);

      std::string cfg_name = fmt::format("/Users/brainkz/Documents/GitHub/ortools_python/{}_cfg_{}.csv", benchmark, file_ctr++);

      std::vector<std::vector<uint64_t>> worms = worm_snake(DFF_REG, spl_phase, path, n_phases, cfg_name);
      
      int num_dff = 0u;
      if (!worms.empty())
      {
        num_dff = cpsat_ortools(cfg_name);
        // std::string command_str = fmt::format("/Users/brainkz/anaconda3/bin/python /Users/brainkz/Documents/GitHub/ortools_python/config_solver.py {}", cfg_name);
        // const char* command = command_str.c_str();

        // std::string pattern = "Objective value: (\\d+)";

        // num_dff = std::system(command);
        fmt::print("OR Tools optimized to {} DFF\n", num_dff);
      }

      total_num_dff += precalc_ndff + num_dff;
    }
    fmt::print("#DFF for {} is {}\n", benchmark, total_num_dff);
  }
  // TODO : now, count #DFFs in extracted paths. Perhaps, the function can be written within the Path object.
  return 0;
}

// Adder 18542 vs 24384

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

typedef mockturtle::klut_network klut;
typedef mockturtle::xag_network   xag;
typedef mockturtle::xmg_network   xmg;
typedef mockturtle::mig_network   mig;
typedef mockturtle::aig_network   aig;
typedef uint64_t node_t;

template <typename elem>
bool haveCommonElements(const std::set<elem>& set1, const std::set<elem>& set2) 
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
  uint8_t   type; 
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

// #if false
std::pair<klut, std::unordered_map<klut::signal, Primitive<klut>>> decompose_to_klut(mockturtle::binding_view<klut> src, std::unordered_map<ULL, Node> nodemap,std::unordered_map<std::string, LibEntry> entries, bool verbose = false)
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
    if (verbose) fmt::print("Processing node {0} ({1} out of {2})\r", src_node, ++ctr, num_node_src);
    if (ctr == num_node_src)
    {
      if (verbose) fmt::print("\n");
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

struct NtkNode
{
  node_t              id;
  std::vector<node_t> fanins;
  std::vector<node_t> fanouts;
  glob_phase_t        phase;
  uint8_t             type; 
  klut::signal        ntk_sig;
  bool                is_spl;
  bool                is_pi;
  bool                is_const;
  bool                is_valid;

  NtkNode(const klut::signal & _id, const std::vector<klut::signal>& _fanins, const std::vector<klut::signal>& _fanouts, glob_phase_t _phase, uint8_t _type, klut::signal _ntk_sig, bool _is_spl, bool _is_pi, bool _is_const)
    : id(_id), fanins(_fanins), fanouts(_fanouts), phase(_phase), type(_type), ntk_sig(_ntk_sig), is_spl(_is_spl), is_pi(_is_pi), is_const(_is_const), is_valid(true) {}
  NtkNode(const klut::signal & _id)
    : id(_id), fanins({}), fanouts({}), phase(0), type(0), ntk_sig(0), is_spl(false), is_pi(false), is_const(false), is_valid(false) {}
};

std::unordered_map<node_t, NtkNode> greedy_assemble(const klut & ntk, const std::unordered_map<klut::signal, Primitive<klut>> & sig_params, const uint8_t n_phases, const bool verbose = false)
{
  mockturtle::topo_view<klut> ntk_topo ( ntk );
  std::unordered_map<node_t, NtkNode> NR( ntk_topo.size() );

  ntk_topo.foreach_node([&] ( const klut::signal & node ) 
  {
    NR.emplace(node, NtkNode(node));
    if ( ntk_topo.is_constant( node ) )
    {
      if (verbose)
      {
        if (verbose) fmt::print("Skipping constant {}\n", node);
      }
      NR.at(node).id = node;
      NR.at(node).ntk_sig = node;
      NR.at(node).is_const = true;
      NR.at(node).is_valid = true;
      return;
    }
    if ( ntk_topo.is_pi( node ) )
    {
      NR.at(node).id = node;
      NR.at(node).phase = 0;
      NR.at(node).type = AS_GATE;
      NR.at(node).ntk_sig = node;
      NR.at(node).is_pi = true;
      NR.at(node).is_valid = true;

      if (verbose)
      {
        if (verbose) fmt::print("PI {} placed at ɸ=0 [S=0, φ=0]\n", node);
      }
      return;
    }
    const Primitive<klut> & node_params = sig_params.at( node );
    if (verbose)
    {
      if (verbose) fmt::print("{} GATE {}:\n", GATE_TYPE.at(node_params.type), node);
    }

    if ( node_params.type == AA_GATE )
    {
      // place at the same phase as the latest fanin

      glob_phase_t phase = 0u;
      ntk_topo.foreach_fanin(node, [&] ( const klut::signal & parent )
      {
        NR.at(node).fanins.push_back(NR.at(parent).id);
        NR.at(parent).fanouts.push_back(NR.at(node).id);
        if ( ntk_topo.is_constant( parent ) )
        {
          if (verbose) fmt::print("\t\tfanin {} is constant, skipping...\n", parent);
          return;
        }
        phase = std::max( phase, NR.at( parent ).phase );
        if (verbose)
        {
          glob_phase_t gp = NR.at( parent ).phase;
          fmt::print("\t\tfanin {} ɸ={} [S={}, φ={}]\n", parent, gp, gp/n_phases, gp%n_phases);
        }

      });
      NR.at(node).id = node;
      NR.at(node).phase = phase;
      NR.at(node).type = AA_GATE;
      NR.at(node).ntk_sig = node;
      NR.at(node).is_valid = true;
      fmt::print("\tAA GATE {} placed at ɸ={} [S={}, φ={}]\n", node, phase, phase/n_phases, phase%n_phases);
    }
    else if ( node_params.type == SA_GATE )
    {
      // Place at the earliest feasible phase
      glob_phase_t phase = 0u;
      ntk_topo.foreach_fanin(node, [&] (const klut::signal & parent)
      {
        NR.at(node).fanins.push_back(NR.at(parent).id);
        NR.at(parent).fanouts.push_back(NR.at(node).id);
        if ( ntk_topo.is_constant( parent ) )
        {
          if (verbose) fmt::print("\t\tfanin {} is constant, skipping...\n", parent);
          return;
        }
        // if d==false, SA gate can be directly connected to fanin
        //cannot directly connect PI
        //cannot directly connect split signal 
        //cannot directly connect SA/AA gates to SA gates
        bool d = (ntk_topo.is_pi(parent)) || (ntk_topo.fanout_size(parent) > 1) || (sig_params.at(parent).type != AS_GATE);                     

        phase = std::max(phase, NR.at( parent ).phase + static_cast<int>(d) );
        if (verbose)
        {
          glob_phase_t gp = NR.at( parent ).phase;
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
      NR.at(node).id = node;
      NR.at(node).phase = phase;
      NR.at(node).type = SA_GATE;
      NR.at(node).ntk_sig = node;
      NR.at(node).is_valid = true;
      if (verbose) fmt::print("\tSA GATE {} placed at ɸ={} [S={}, φ={}]\n", node, phase, phase/n_phases, phase%n_phases);
    }
    else if ( node_params.type == AS_GATE )
    {
      // Place at the earliest feasible phase
      glob_phase_t phase = 0u;
      ntk_topo.foreach_fanin(node, [&] (const klut::signal & parent)
      {
        NR.at(node).fanins.push_back(NR.at(parent).id);
        NR.at(parent).fanouts.push_back(NR.at(node).id);
        if ( ntk_topo.is_constant( parent ) )
        {
          if (verbose) fmt::print("\t\tfanin {} is constant, skipping...\n", parent);
          return;
        }
        phase = std::max(phase, NR.at( parent ).phase );
        if (verbose)
        {
          glob_phase_t gp = NR.at( parent ).phase;
          fmt::print("\t\tfanin {} ɸ={} [S={}, φ={}]\n", parent, gp, gp/n_phases, gp%n_phases);
        }
      });
      phase++;
      NR.at(node).id = node;
      NR.at(node).phase = phase;
      NR.at(node).type = AS_GATE;
      NR.at(node).ntk_sig = node;
      NR.at(node).is_valid = true;
      if (verbose) fmt::print("\tAS GATE {} placed at ɸ={} [S={}, φ={}]\n", node, phase, phase/n_phases, phase%n_phases);
    }
    else 
    {
      fmt::print("\t GATE {} : fanins[{}], type[{}], func[{}]\n", node, fmt::join(node_params.fanins, ","), node_params.type, node_params.func);
      throw "Unsupported gate type!";
    }
  });
  return NR;
}

bool phase_comparison(const node_t & a, const node_t & b, const std::unordered_map<node_t, NtkNode> & NR)
{
  const NtkNode & a_node = NR.at(a);
  const NtkNode & b_node = NR.at(b);
  if (a_node.is_const)
  {
    return true;
  }
  else if (b_node.is_const)
  {
    return false;
  }
  return a_node.phase < b_node.phase;
}

void splitter_insertion( std::unordered_map<node_t, NtkNode> & NR , bool verbose = false)
{
  auto phase_comp = [&](const node_t & a, const node_t & b){ return phase_comparison(a, b, NR); };

  std::vector<node_t> stack;
  std::vector<node_t> seen;
  // record the PO-s of the network
  for (auto & [id, ntk_node] : NR)
  {
    if (ntk_node.is_const)
    {
      continue;
    }
    std::sort(ntk_node.fanins.begin(), ntk_node.fanins.end(), phase_comp);
    std::sort(ntk_node.fanouts.begin(), ntk_node.fanouts.end(), phase_comp);
    if (ntk_node.fanouts.size() > 1)
    {
      stack.push_back(id);
    }
  }

  node_t spl_id = NR.size();
  while (!stack.empty())
  {
    if (verbose) fmt::print("Stack size: {}\n", stack.size());
    // std::sort(stack.begin(), stack.end(), phase_comp);
    node_t parent_id = stack.back();
    NtkNode & parent = NR.at( parent_id );

    if (verbose) fmt::print("\t{} fanouts: [{}]\n", parent.fanouts.size(), fmt::join(parent.fanouts, ","));
    // create a splitter object, increment splitter counter
    
    if (verbose) fmt::print("\t\tAdding splitter {} to NR (size before = {}, ", spl_id, NR.size());
    NR.emplace(spl_id, spl_id);
    if (verbose) fmt::print("size after = {})\n", NR.size());
    NtkNode & spl = NR.at(spl_id);
    spl_id++;

    if (verbose) fmt::print("\tAdding fanouts to splitter {} -> [{}]\n", spl.id, fmt::join(parent.fanouts.end()-2, parent.fanouts.end(), ","));
    // splitter ----> fanouts 
    spl.fanouts.insert( spl.fanouts.end(), parent.fanouts.end()-2, parent.fanouts.end() );


    if (verbose) fmt::print("\tRemoving connection   parent {} <---- fanouts [{}]\n", parent.id, fmt::join(spl.fanouts, ","));
    if (verbose) fmt::print("\tCreating connection splitter {} <---- fanouts [{}]\n",    spl.id, fmt::join(spl.fanouts, ","));
    // parent  <-/-/- fanouts
    // splitter <---- fanouts
    for (auto it = parent.fanouts.end()-2; it < parent.fanouts.end(); ++it)
    {
      node_t fo_id = *it;
      NtkNode & fo = NR.at( fo_id );
      // remove parent from the fanins of the ntk_node
      // replace with the splitter parent from the fanins of the ntk_node
      if (verbose) fmt::print("\t\t\t{} fanins before replacement: [{}]\n", fo.id, fmt::join(fo.fanins, ","));
      auto parent_it = std::find(fo.fanins.begin(), fo.fanins.end(), parent_id);
      *parent_it = spl.id; // fo.fanins.erase(parent_it);
      if (verbose) fmt::print("\t\t\t{} fanins after  replacement: [{}]\n", fo.id, fmt::join(fo.fanins, ","));
    }
    // parent  -/-/-> fanouts
    if (verbose) fmt::print("\tRemoving connection   parent {} ----> fanouts [{}]\n", parent.id, fmt::join(spl.fanouts, ","));
    if (verbose) fmt::print("\t\t\tParent {} fanouts before replacement: [{}]\n", parent.id, fmt::join(parent.fanouts, ","));
    parent.fanouts.erase(parent.fanouts.end()-2, parent.fanouts.end());
    if (verbose) fmt::print("\t\t\tParent {} fanouts after  replacement: [{}]\n", parent.id, fmt::join(parent.fanouts, ","));
    
    if (verbose) fmt::print("\tConnecting parent {} <---> splitter {}\n", parent.id, spl.id);
    // parent <----> splitter
    if (verbose) fmt::print("\t\t\tSplitter {} fanins before connection: [{}]\n", spl.id, fmt::join(spl.fanins, ","));
    spl.fanins.push_back(parent_id);  // connect spl to parent 
    if (verbose) fmt::print("\t\t\tSplitter {} fanins after  connection: [{}]\n", spl.id, fmt::join(spl.fanins, ","));
    
    if (verbose) fmt::print("\t\t\tParent {} fanouts before adding splitter {}: [{}]\n", parent.id, spl.id, fmt::join(parent.fanouts, ","));
    parent.fanouts.push_back(spl.id); // connect parent to spl
    if (verbose) fmt::print("\t\t\tParent {} fanouts after  adding splitter {}: [{}]\n", parent.id, spl.id, fmt::join(parent.fanouts, ","));
    // get maximum phase
    std::vector<glob_phase_t> phases;
    for (const node_t & fo_id : spl.fanouts)
    {
      phases.push_back(NR.at( fo_id ).phase);
    }
    glob_phase_t phase = std::min(phases);

    // Set splitter parameters
    if (verbose) fmt::print("\tAssigning phase {}\n", phase);
    spl.phase    = phase;
    spl.type     = AA_GATE;
    spl.ntk_sig  = spl.id;
    spl.is_spl   = true;
    spl.is_valid = true;

    // sort fanins/fanouts
    std::sort(spl.fanins.begin(), spl.fanins.end(), phase_comp);
    std::sort(spl.fanouts.begin(), spl.fanouts.end(), phase_comp);
    std::sort(parent.fanouts.begin(), parent.fanouts.end(), phase_comp);
    for (const node_t & fo_id : spl.fanouts)
    {
      auto & fo = NR.at( fo_id );
      std::sort(fo.fanins.begin(), fo.fanins.end(), phase_comp);
    }
    if (parent.fanouts.size() == 1)
    {
      stack.pop_back();
    }
    // stack.push_back( spl.id );
  }
}

uint64_t calculate_hash(glob_phase_t phase, klut::signal fanout, klut::signal fanin)
{
    std::size_t seed = 0;
    hash_combine(seed, phase);
    hash_combine(seed, fanout);
    hash_combine(seed, fanin);
    return seed;
};

struct DFF_var 
{
  node_t fanin;
  node_t fanout;
  glob_phase_t phase;
  std::unordered_set<uint64_t> parent_hashes;
  uint64_t hash;

  DFF_var(node_t _fanin, node_t _fanout, glob_phase_t _phase)
      : fanin(_fanin), fanout(_fanout), phase(_phase), parent_hashes({}) 
      {
        hash = calculate_hash(phase, fanout, fanin);
      }

  DFF_var(node_t _fanin, node_t _fanout, glob_phase_t _phase, std::unordered_set<uint64_t> _parent_hashes, uint64_t _hash)
      : fanin(_fanin), fanout(_fanout), phase(_phase), parent_hashes(_parent_hashes), hash(_hash) {}

  DFF_var(const DFF_var& other)
      : fanin(other.fanin), fanout(other.fanout), phase(other.phase), parent_hashes(other.parent_hashes), hash(other.hash) {}

  std::string str()
  {
    return fmt::format( "var_{}_{}_{}", fanin, fanout, phase );
  }
};

uint64_t dff_hash(node_t _fanin, node_t _fanout, glob_phase_t _phase)
{
  return ( (uint64_t)_fanin << 40 ) | ( (uint64_t)_fanout << 16 ) | _phase;
}

struct DFF_registry
{
  std::unordered_map<uint64_t, DFF_var> variables;

  DFF_var & at(node_t _fanin, node_t _fanout, glob_phase_t _phase)
  {
    return variables.at( dff_hash(_phase, _fanout, _fanin) );
  } 
  DFF_var & at(uint64_t _hash)
  {
    return variables.at( _hash );
  } 
  uint64_t add(node_t _fanin, node_t _fanout, glob_phase_t _phase, std::unordered_set<uint64_t> _parent_hashes = {})
  {
    uint64_t _hash = dff_hash(_phase, _fanout, _fanin);
    DFF_var temp { _fanin, _fanout, _phase, _parent_hashes, _hash };
    variables.emplace(_hash, temp);
    // variables.insert(std::make_pair(_hash, temp));
    return _hash;
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


struct Path
{
  std::set<node_t> sources;   // AS/SA gates
  std::set<node_t> internals; // AA    gates
  std::set<node_t> targets;   // AS/SA gates
  Path(const std::set<node_t> & _sources, const std::set<node_t>& _internals, const std::set<node_t>& _targets)
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
  std::string format() const
  {
    return fmt::format("Path from [{}]\n\tvia [{}]\n\tto [{}]\n", fmt::join(sources, ","), fmt::join(internals, ","), fmt::join(targets, ","));
  }

  std::vector<node_t> preds(const node_t & sig_id, const std::unordered_map<node_t, NtkNode> & NR) const
  {
    if ( sources.count(sig_id) != 0)
    {
      return {};
    }
    const NtkNode & sig = NR.at( sig_id );
    std::vector<node_t> predecessors;
    for (const node_t & parent : sig.fanins)
    {
      if ( internals.count(parent) != 0 || sources.count(parent) != 0 )
      {
        predecessors.push_back( parent );
      }
    }
    return predecessors;
  }

  std::vector<node_t> src_int() const
  {
    std::vector<node_t> out;
    out.insert(out.end(), sources.begin(), sources.end());
    out.insert(out.end(), internals.begin(), internals.end());
    return out;
  }
  std::vector<node_t> int_tgt() const
  {
    std::vector<node_t> out;
    out.insert(out.end(), internals.begin(), internals.end());
    out.insert(out.end(), targets.begin(), targets.end());
    return out;
  }
};

std::vector<Path> extract_paths(const std::unordered_map<node_t, NtkNode> & NR, bool verbose = false)
{
  // std::unordered_map<klut::signal, std::vector<klut::signal>> ancestors;
  std::vector<Path> paths;

  for (const auto & [fo_id, fo_node] : NR)
  {
    if (fo_node.is_const || fo_node.is_pi || (fo_node.type != AS_GATE && fo_node.type != SA_GATE) )
    {
      continue;
    }
    // at this point, the node should be AS/SA

    // Create a separate path for each fanin of the node
    std::vector<Path> aa_paths; 
    aa_paths.reserve( fo_node.fanins.size() );

    for (const node_t & fi_id : fo_node.fanins)
    {
      // Create a path object with only a target
      Path node_path;
      node_path.targets.emplace( fo_id );
      
      std::vector<node_t> stack;
      stack.push_back( fi_id );
      
      if ( verbose ) fmt::print("Node {}, traversal for parent {}\n", fo_id, fi_id);

      std::set<node_t> seen;
      while (!stack.empty())
      {
        node_t & node_id = stack.back();
        const NtkNode & n = NR.at( node_id );
        stack.pop_back();

        if ( verbose ) fmt::print("\t[Parent {}]: Analyzing node {}\n", fi_id, node_id);
        // A constant does not have any effect on the DFF placement, we can skip it
        if ( n.is_const )
        {
          continue;
        }
        // Found a source of the path, add to sources, do not continue traversal
        else if ( n.is_pi || n.type == AS_GATE || n.type == SA_GATE )
        {
          if ( verbose ) fmt::print("\t\t{} is a source node \n", node_id);
          node_path.sources.emplace( node_id );
        }
        // Found AA gate, add to internal nodes, add parents for further traversal
        else if ( n.type == AA_GATE )
        {
          if ( verbose ) fmt::print("\t\t{} is INTERNAL adding fanins \n", node_id);
          node_path.internals.emplace( node_id );
          for (const node_t & sig_id : n.fanins)
          {
            stack.push_back( sig_id );
          }
        }
        else
        {
          if ( verbose ) fmt::print("Signal {}: {} is not recognized \n", node_id, GATE_TYPE.at( n.type ));
          throw "Unsupported case";
        }
        seen.emplace( node_id );
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
      // Iterating in reverse order to preserve order
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
    }
  }
  return paths;
}

union Edge
{
  struct 
  {
    uint32_t fanin;
    uint32_t fanout;
  };
  uint64_t value;

  Edge(uint64_t _fanin, uint64_t _fanout) : fanin(_fanin), fanout(_fanout) {};
  Edge(uint64_t _value) : value(_value) {};

  bool operator==(const Edge & other) const 
  {
    return value == other.value;
  }
};


#if false
  std::pair<DFF_registry, uint64_t> dff_vars(const std::unordered_map<node_t, NtkNode> & NR, const std::vector<Path> & paths, const uint8_t n_phases)
  {
    DFF_registry DFF_REG;
    uint64_t precalc_ndff = 0; //precalculated number of DFFs for trivial cases

    for (const Path & path : paths)
    {
      std::vector<std::tuple<node_t, uint64_t>> stack;
      for (const node_t & tgt_id : path.targets)
      {
        stack.emplace_back(tgt_id, 0);
      }
      std::vector<Edge> seen;
      
      while (!stack.empty())
      {
        const auto & [ fo_id, earliest_child_hash ] = stack.back();
        stack.pop_back();

        const NtkNode & fo_node = NR.at( fo_id );
        glob_phase_t latest_phase = fo_node.phase - (fo_node.type == AS_GATE);
        for (const node_t & fi_id : fo_node.fanins)
        {
          const NtkNode & fi_node = NR.at( fi_id );
          glob_phase_t earliest_phase = fi_node.phase + (fi_node.type != AA_GATE);

          Edge edge { fi_id, fo_id }; 
          
          // check if the chain is straight - #DFF is just floor(delta-phase), no need to create the dff vars
          if (fo_node.type != AA_GATE && fi_node.type != AA_GATE)
          {
            // special case when an AS gate feeds directly into SA gate
            if (fo_node.phase == fi_node.phase)
            {
              // do nothing, no additional DFFs needed
              assert(fo_node.type == SA_GATE && fi_node.type == AS_GATE && fi_node.fanouts.size() == 1);
            }
            else
            {
              // straight chain, just floor the difference!
              precalc_ndff += (fo_node.phase - fi_node.phase)/n_phases + (fo_node.type == SA_GATE); //extra DFF before SA gate
            }
          }
          
          // std::vector<uint64_t> out_hashes;
          uint64_t prev_hash = DFF_REG.add(fi_id, fo_id, earliest_phase);
          uint64_t earliest_hash = prev_hash;
          // out_hashes.push_back(prev_hash);
          for (glob_phase_t phase = earliest_phase + 1; phase <= latest_phase; ++phase)
          {
            prev_hash = DFF_REG.add(fi_id, fo_id, phase, { prev_hash });
            // out_hashes.push_back(last_hash);
          }

          // if the node is internal, connect with the fanout phase
          if (fo_node.type == AA_GATE)
          {
            DFF_var & child_dff = DFF_REG.at(earliest_child_hash);
            child_dff.parent_hashes.emplace( prev_hash );
          }
          auto it = std::find(seen.begin(), seen.end(), edge);
          if (it == seen.end())
          {
            stack.emplace_back( fi_id, earliest_hash );
            seen.push_back( edge );
          }
        }
      }
    }
    return std::make_pair(DFF_REG, precalc_ndff);
  }
#endif

// std::pair<std::vector<DFF_var>, uint64_t> 
std::tuple<DFF_registry, uint64_t, std::vector<uint64_t>> dff_vars_single(const Path & path, const std::unordered_map<node_t, NtkNode> & NR, const uint8_t n_phases, bool verbose = false)
{
  DFF_registry DFF_REG;
  std::vector<uint64_t> required_SA_DFFs;

  std::vector<std::tuple<node_t, uint64_t>> stack;
  for (const node_t & tgt_id : path.targets)
  {
    stack.emplace_back(tgt_id, 0);
  }
  if (verbose) fmt::print("[DFF] Target nodes: {}\n", fmt::join(path.targets, ","));

  std::vector<Edge> seen;

  auto precalc_ndff = 0u;
  
  while (!stack.empty())
  {
    if (verbose) fmt::print("STACK :\n");
    for (const auto & [ fo_id, earliest_child_hash ] : stack)
    {
      const NtkNode & fo_node = NR.at( fo_id );
      if (earliest_child_hash != 0)
      {
        if (verbose) fmt::print("\t{}({})[{}], {}\n", GATE_TYPE.at(fo_node.type), fo_id, fo_node.phase, DFF_REG.at(earliest_child_hash).str());
      }
      else
      {
        if (verbose) fmt::print("\t{}({})[{}]\n", GATE_TYPE.at(fo_node.type), fo_id, fo_node.phase);
      }
    }
    const auto & [ fo_id, earliest_child_hash ] = stack.back();
    stack.pop_back();

    // DFF_var & child_dff = DFF_REG.at( earliest_child_hash );
    // if (verbose) fmt::print("Prior node is {}[{}]\n", child_dff.str(), child_dff.phase); 
    // assert(child_dff.fanin == fo_id);

    const NtkNode & fo_node = NR.at( fo_id );
    glob_phase_t latest_phase = fo_node.phase - (fo_node.type == AS_GATE);
    if (verbose) fmt::print("[DFF] Analyzing child: {}({})[{}]\n", GATE_TYPE.at(fo_node.type), fo_id, fo_node.phase);
    for (const node_t & fi_id : fo_node.fanins)
    {
      const NtkNode & fi_node = NR.at( fi_id );
      glob_phase_t earliest_phase = fi_node.phase + (fi_node.type != AA_GATE);

      if (verbose) fmt::print("\t[DFF] Analyzing parent: {}({})[{}]\n", GATE_TYPE.at(fi_node.type), fi_id, fi_node.phase);
      Edge edge { fi_id, fo_id }; 
      
      // check if the chain is straight - #DFF is just floor(delta-phase), no need to create the dff vars
      if (fo_node.type != AA_GATE && fi_node.type != AA_GATE)
      {
        // special case when an AS gate feeds directly into SA gate
        if (fo_node.phase == fi_node.phase)
        {
          if (verbose) fmt::print("\t[DFF] Straight chain: AS{} -> SA{}\n", fi_id, fo_id);
          // do nothing, no additional DFFs needed
          assert(fo_node.type == SA_GATE && fi_node.type == AS_GATE && fi_node.fanouts.size() == 1);
        }
        else
        {
          if (verbose) fmt::print("\t[DFF] Straight chain: {}[{}] -> {}[{}]\n", GATE_TYPE.at(fi_node.type), fi_node.phase, GATE_TYPE.at(fo_node.type), fo_node.phase);
          // straight chain, just floor the difference!
          precalc_ndff += (fo_node.phase - fi_node.phase)/n_phases + (fo_node.type == SA_GATE); //extra DFF before SA gate
        }
        continue;
      }
      

      if (verbose) fmt::print("\t[DFF] Non-straight chain: {}[{}] -> {}[{}]\n", GATE_TYPE.at(fi_node.type), fi_node.phase, GATE_TYPE.at(fo_node.type), fo_node.phase);
      std::vector<uint64_t> out_hashes;
      if (verbose) fmt::print("\tAdding new DFFs [reg size = {}]\n", DFF_REG.variables.size());
      for (glob_phase_t phase = earliest_phase; phase <= latest_phase; ++phase)
      {
        uint64_t new_hash = DFF_REG.add(fi_id, fo_id, phase);
        out_hashes.push_back(new_hash);
        if (verbose) fmt::print("\tAdded new DFFs at phase {} [reg size = {}]\n", phase, DFF_REG.variables.size());
      }

      if (verbose) fmt::print("\tConnecting new DFFs\n");
      for (auto i = 1u; i < out_hashes.size(); ++i)
      {
        DFF_var & dff = DFF_REG.at( out_hashes[i] );
        dff.parent_hashes.emplace(out_hashes[i-1]);
      }

      if (fo_node.type == SA_GATE)
      {
        assert( !out_hashes.empty() );
        required_SA_DFFs.push_back(out_hashes.back());
      }

      uint64_t earliest_hash = (out_hashes.empty()) ? earliest_child_hash : out_hashes.front();
      // if the node is internal, connect with the fanout phase
      if (fo_node.type == AA_GATE && !out_hashes.empty() && earliest_hash != 0 && earliest_child_hash != 0)
      {
        DFF_var & child_dff = DFF_REG.at( earliest_child_hash );
        if (verbose) fmt::print("\tPrior node is {}[{}]\n", child_dff.str(), child_dff.phase); 
        // assert(child_dff.fanin == fo_id);
        child_dff.parent_hashes.emplace( out_hashes.back() );
      }
      if (fi_node.type == AA_GATE)
      {
        stack.emplace_back( fi_id, earliest_hash );
        if (verbose) fmt::print("\tEmplacing {}({})[{}], {}\n", GATE_TYPE.at(fi_node.type), fi_id, fi_node.phase, (earliest_hash!=0)?DFF_REG.at(earliest_hash).str():"");
      }

      // auto it = std::find(seen.begin(), seen.end(), edge);
      // if (it == seen.end())
      // {
      //   seen.push_back( edge );
      // }
    }
  }
  return std::make_tuple(DFF_REG, precalc_ndff, required_SA_DFFs);
}

struct Snake
{
  std::deque<std::vector<uint64_t>> sections;

  Snake(): sections({}) {}
  Snake( const uint64_t head ): sections({ { head } }) {}
  Snake( const std::deque<std::vector<uint64_t>> _sections ): sections(_sections) {}
  Snake( const Snake & _other ): sections(_other.sections) {}

  bool append(const uint64_t dff_hash, DFF_registry &DFF_REG, const uint8_t n_phases)
  {
    DFF_var & dff = DFF_REG.at( dff_hash );

    std::vector<uint64_t> & head_section = sections.back();
    uint64_t & head_hash = head_section.back();
    DFF_var & head_dff = DFF_REG.at( head_hash );
    if (dff.phase == head_dff.phase) // add to the same section
    {
      head_section.push_back( dff_hash );
      return false;
    }
    else
    {
      // assert( head_dff.phase - dff.phase == 1 );
      sections.push_back( { dff_hash } );
      if (sections.size() > n_phases)
      {
        sections.pop_front();
      }
      return true;
    }
  }

};

void write_snakes(const std::vector<Snake> & snakes, const std::unordered_map<node_t, NtkNode> & NR, DFF_registry & DFF_REG, const std::vector<uint64_t> & required_SA_DFFs, const std::string cfg_name, uint8_t n_phases, bool verbose = false)
{
  std::ofstream spec_file (cfg_name);

  for (const Snake & snake : snakes)
  {
    std::vector<std::string> vars_bucket;
    for (const std::vector<uint64_t> & section : snake.sections)
    {
      std::vector<std::string> vars;
      for (uint64_t hash : section)
      {
        vars.push_back(DFF_REG.str( hash ));
      }
      if (verbose) fmt::print("New single phase conflict : {}≤1\n", fmt::join(vars, "+"));
      vars_bucket.push_back(fmt::format(vars.size()>1?"({})":"{}", fmt::join(vars, "+")));
      if (vars.size() > 1)
      {
        spec_file << fmt::format("PHASE,{}\n", fmt::join(vars, ","));
      }
    }
    std::reverse(vars_bucket.begin(), vars_bucket.end());
    if (verbose) fmt::print("New buffer requirement : ({})\n", fmt::join(vars_bucket, "|"));
    if (vars_bucket.size() == n_phases)
    {
      spec_file << fmt::format("BUFFER,{}\n", fmt::join(vars_bucket, ","));
    }
  }

  for (const uint64_t & hash : required_SA_DFFs)
  {
    if (verbose) fmt::print("New SA_REQUIRED : {}\n", DFF_REG.str( hash ));
    spec_file << fmt::format("SA_REQUIRED,{}\n", DFF_REG.str( hash ));
  }
}

std::vector<Snake> sectional_snake(const Path & path, std::unordered_map<node_t, NtkNode> &NR,  DFF_registry & DFF_REG, uint8_t n_phases, bool verbose = false)
{
  std::vector<Snake> out_snakes; 
  std::vector<Snake> stack;
  
  if (verbose) fmt::print("[i]: Starting extraction of worms \n");
  // get all DFFs 
  for (const auto & [hash, dff]: DFF_REG.variables)
  {
    auto fanout_phase = NR.at( dff.fanout ).phase - ( NR.at( dff.fanout ).type == AS_GATE );
    auto it = std::find(path.targets.begin(), path.targets.end(), dff.fanout);
    if (it != path.targets.end() && (dff.phase >= fanout_phase ))
    {
      stack.emplace_back( hash );
    }
  }
  
  while (!stack.empty())
  {
    if (verbose) fmt::print("[i] Stack size is {} \n", stack.size());
    Snake snake = stack.back();
    stack.pop_back();

    if (verbose) fmt::print("\t[i] The snake has {} sections\n", snake.sections.size());
    uint64_t hash = snake.sections.back().back();
    DFF_var & dff = DFF_REG.at( hash );

    // fmt::print("\tCurrent worm size {}, between phases {} and {} \n", worm.size(), DFF_REG.str(worm.front()), DFF_REG.str(worm.back()));


    if (verbose) fmt::print("\t\t[i] The DFF {} has {} parents\n", DFF_REG.at( hash ).str(),  dff.parent_hashes.size() );

    bool returned_current_snake = false;
    for (const uint64_t parent_hash : dff.parent_hashes)
    {
      Snake snake_copy = snake; 
      if (verbose) fmt::print("\t\t[i] Advancing towards fanin {}\n", DFF_REG.at( parent_hash ).str() );
      bool status = snake_copy.append(parent_hash, DFF_REG, n_phases);
      if (verbose) fmt::print((status) ? "\t\t\tAdded new section!\n" :"\t\t\tExtended existing section!\n"  );
      if (verbose) fmt::print("\t\t\tThe new length is {}\n", snake_copy.sections.size() );
      
      stack.push_back( snake_copy );
      if (status && !returned_current_snake && snake_copy.sections.size() == n_phases)
      {
        if (verbose) fmt::print("\t\tAdding the snake to the output\n");
        out_snakes.push_back(snake);
        returned_current_snake = true;
      }
    }
  }
  return out_snakes;
}

std::vector<std::pair<std::vector<std::vector<uint64_t>>,std::string>> snake(DFF_registry & DFF_REG, const std::unordered_map<node_t, NtkNode> & NR, const std::vector<Path> & paths, const uint8_t n_phases, const std::string & cfg_prefix)
{
  std::vector< std::pair< std::vector< std::vector<uint64_t> >, std::string > >  global_worms;

  auto file_ctr = 0u;
  for (const Path & path : paths)
  {
    path.print(); 

    fmt::print("[i]: Starting extraction of worms \n");
    std::vector<std::vector<uint64_t>> worms;
    std::vector<std::vector<uint64_t>> stack;
    std::vector<std::vector<uint64_t>> seen;
    
    for (const auto & [hash, dff]: DFF_REG.variables)
    {
      auto it = std::find(path.targets.begin(), path.targets.end(), dff.fanout);
      if (it != path.targets.end() && (dff.phase == NR.at( dff.fanout ).phase ))
      {
        stack.push_back( { hash } );
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
    std::string cfg_name;
    if (worms.size() > 0)
    {
      cfg_name = fmt::format("{}{}.csv", cfg_prefix, file_ctr++);
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
    global_worms.emplace_back(worms, cfg_name);
  }
  return global_worms;
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

constexpr uint8_t N_PHASES = 5u;

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
  // experiments::adder | experiments::sin | experiments::cavlc | | experiments::dec | experiments::i2c | experiments::adder | experiments::int2float | | experiments::voter | experiments::c432 | experiments::c499 | | |experiments::dec |
  auto benchmarks1 = epfl_benchmarks(   experiments::sin | experiments::priority | experiments::voter | experiments::int2float | experiments::i2c  );
  // //   auto benchmarks1 = epfl_benchmarks( experiments::epfl & ~experiments::div & ~experiments::hyp & ~experiments::log2 & ~experiments::sqrt );
  // experiments::c880 | experiments::c1908 | experiments::c3540 |
  auto benchmarks2 = iscas_benchmarks( experiments::c432 | experiments::c499 | experiments::c1908 | experiments::c3540 |  experiments::c5315 | experiments::c7552 );

  benchmarks1.insert(benchmarks1.end(), benchmarks2.begin(), benchmarks2.end());



  fmt::print("{}\n", fmt::join(benchmarks1, "\n"));

  // TODO :  Debug experiments::c5315 with GATE 0x9696_112_0 not showing up in GNM . 
  // TODO :  Perhaps genlib_hierarchical.cpp is broken, need to debug
  
  const std::string nodemap_prefix = "/Users/brainkz/Documents/GitHub/mockturtle/build/Golden_20230427/x3";
  const std::vector<std::vector<UI>> sets_of_levels { { {0,0,0,0}, {0,0,0,1}, {0,0,0,2}, {0,0,1,1}, {0,0,1,2}, {0,1,1,1}, {0,1,1,2}, {0,1,2,2}, {0,1,2,3} } }; //  {0,1,1,3},
  // const std::vector<std::vector<UI>> sets_of_levels { { {0,1,2,3} } }; //  {0,1,1,3},
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
    auto [res_w_pb, st_w_pb] = map_wo_pb(aig_original, tech_lib, false); //benchmark, true, nDFF_global, total_ndff_w_pb, total_area_w_pb, cec_w_pb 
    fmt::print("Finished mapping of {}\n", benchmark);

    auto [klut_decomposed, klut_prim_params] = decompose_to_klut(res_w_pb, GNM_global, entries);
    fmt::print("Decomposition complete\n");

    std::unordered_map<node_t, NtkNode> NR = greedy_assemble(klut_decomposed, klut_prim_params, N_PHASES, true);

    splitter_insertion( NR );

    for (auto & [id, ntk_node] : NR)
    {
      if (ntk_node.is_spl)
      {
        assert(ntk_node.fanouts.size() == 2);
      }
      else
      {
        assert(ntk_node.fanouts.size() <= 1);
      }
    }
    fmt::print("[i] FINISHED PHASE ASSIGNMENT\n");

    fmt::print("[i] EXTRACTING PATHS\n");
    std::vector<Path> paths = extract_paths(NR, true);
    // auto [DFF_REG, precalc_ndff] = dff_vars(NR, paths, N_PHASES);

    auto total_num_dff = 0u;
    auto file_ctr = 0u;
    for (const Path & path : paths)
    {
      auto [DFF_REG, precalc_ndff, required_SA_DFFs] = dff_vars_single(path, NR, N_PHASES);
      total_num_dff += precalc_ndff;
      fmt::print("[i]: Precalculated {} DFFs\n", precalc_ndff);

      std::vector<Snake> snakes = sectional_snake(path, NR, DFF_REG, N_PHASES);
      fmt::print("\tCreated {} snakes\n", snakes.size());
      if (!snakes.empty())
      {
        std::string cfg_file = fmt::format("/Users/brainkz/Documents/GitHub/ortools_python/{}_cfgNR_{}.csv", benchmark, file_ctr++);
        write_snakes(snakes, NR, DFF_REG, required_SA_DFFs, cfg_file, N_PHASES);
        auto num_dff = cpsat_ortools(cfg_file);
        fmt::print("OR Tools optimized to {} DFF\n", num_dff);
        total_num_dff += num_dff;
      }
    }
    fmt::print("#DFF for {} is {}\n", benchmark, total_num_dff);

    // std::string cfg_prefix = fmt::format("/Users/brainkz/Documents/GitHub/ortools_python/{}_cfgNR_", benchmark);

    // std::vector<std::pair<std::vector<std::vector<uint64_t>>, std::string>> global_worms = snake(DFF_REG, NR, paths, N_PHASES, cfg_prefix);

    // auto total_num_dff = precalc_ndff;
    // for (const auto & [worms, cfg_name] : global_worms)
    // {
    //   if (cfg_name.size() != 0)
    //   {
    //     auto num_dff = cpsat_ortools(cfg_name);
    //     fmt::print("OR Tools optimized to {} DFF\n", num_dff);
    //     total_num_dff += precalc_ndff + num_dff;
    //   }
    // }
    // fmt::print("#DFF for {} is {}\n", benchmark, total_num_dff);

    // for (const Path & path : paths)
    // {
    //   path.print();
    //   Chain_REGISTRY path_reg = fill_dffs( spl_klut, path, spl_phase, spl_params );

    //   auto [DFF_REG, precalc_ndff] = generate_dff_vars(spl_klut, path_reg, path, spl_phase, spl_params);
    //   fmt::print("Precalculated {} DFFs!\n", precalc_ndff);

    //   std::string cfg_name = fmt::format("/Users/brainkz/Documents/GitHub/ortools_python/{}_cfg_{}.csv", benchmark, file_ctr++);

    //   std::vector<std::vector<uint64_t>> worms = worm_snake(DFF_REG, spl_phase, path, n_phases, cfg_name);
      
    //   int num_dff = 0u;
    //   if (!worms.empty())
    //   {
    //     num_dff = cpsat_ortools(cfg_name);
    //     // std::string command_str = fmt::format("/Users/brainkz/anaconda3/bin/python /Users/brainkz/Documents/GitHub/ortools_python/config_solver.py {}", cfg_name);
    //     // const char* command = command_str.c_str();

    //     // std::string pattern = "Objective value: (\\d+)";

    //     // num_dff = std::system(command);
    //     fmt::print("OR Tools optimized to {} DFF\n", num_dff);
    //   }

    //   total_num_dff += precalc_ndff + num_dff;
    // }
    // fmt::print("#DFF for {} is {}\n", benchmark, total_num_dff);
  }
  // TODO : now, count #DFFs in extracted paths. Perhaps, the function can be written within the Path object.
  return 0;
}

// Adder 18542 vs 24384

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
#include <filesystem>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/blif.hpp>
#include <lorina/genlib.hpp>
#include <lorina/bench.hpp>
#include <mockturtle/algorithms/rsfq/rsfq_network_conversion.hpp>
#include <mockturtle/algorithms/rsfq/rsfq_path_balancing.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/nodes.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/algorithms/retiming.hpp>
// #include <mockturtle/algorithms/refactoring.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/bench_reader.hpp>
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
#include <mockturtle/algorithms/klut_to_graph.hpp>

// mockturtle/algorithms/mig_algebraic_rewriting.hpp

#include <mockturtle/io/auxiliary_genlib.hpp>

// #include <mockturtle/utils/GNM_global.hpp> // GNM global is stored here

#include <mockturtle/utils/misc.hpp>

#include <experiments.hpp>

#define DEBUG_PRINT(format, ...) if (verbose) fmt::print(format, ##__VA_ARGS__)

typedef mockturtle::klut_network klut;
typedef mockturtle::xag_network   xag;
typedef mockturtle::xmg_network   xmg;
typedef mockturtle::mig_network   mig;
typedef mockturtle::aig_network   aig;
typedef uint64_t node_t;

// Sunmagnetics Technology Library
constexpr std::array<int,12> COSTS_MAP = {7, 9, 8, 8, 12, 8, 999, 999, 999, 8, 3, 0};

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
  double total_area = st.area + COSTS_MAP[fSPL] * num_splitters;
  //  +  COSTS_MAP[fDFF] * num_ext_dffs;
  fmt::print("\t{} : Int: {}, Ext: {}, ratio: {}\n", benchmark, num_int_dffs, num_ext_dffs, (float)num_int_dffs / (num_int_dffs + num_ext_dffs) );
  return std::make_tuple( res, st, total_ndff, total_area, cec );
}


constexpr uint8_t PI_GATE = 0u;
constexpr uint8_t AA_GATE = 1u;
constexpr uint8_t AS_GATE = 2u;
constexpr uint8_t SA_GATE = 3u;
constexpr uint8_t FA_GATE = 4u;

union NodeData 
{
  struct 
  {
    unsigned int sigma : 29;
    unsigned int type : 3;
  };
  uint32_t value;
  NodeData() : value(0) {}
  NodeData(uint32_t _sigma, uint8_t _type) : sigma(_sigma), type(_type) {}
  NodeData(uint32_t _value) : value(_value) {}
};

const std::vector<std::string> GATE_TYPE { "PI", "AA", "AS", "SA", "FA" }; //, "PO"

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

template <typename Ntk>
class fmt::formatter<Primitive<Ntk>> {
public:
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

  template <typename Context>
  constexpr auto format(const Primitive<Ntk>& primitive, Context& ctx) const {
    return format_to(ctx.out(), "Primitive{{sig={}, func={}, type={}, fanins={}, is_spl={}}}",
      primitive.sig, primitive.func, static_cast<int>(primitive.type),
      fmt::join(primitive.fanins, ", "), primitive.is_spl);
  }
};

// #if false
std::tuple<klut, std::unordered_map<klut::signal, Primitive<klut>>, int64_t> decompose_to_klut(mockturtle::binding_view<klut> src, std::unordered_map<ULL, Node> nodemap, std::unordered_map<std::string, LibEntry> entries, bool verbose = false)
{
  /* record the information of each signal in tgt, including: signal index (why?), */
  /* truth table ( represented using an uint16 ), type ( PI/AA/AS/AS ), fan-ins,   */
  /* and if it is a splitter                                                       */
  std::unordered_map<klut::signal, Primitive<klut>> tgt_sig_params;
  /* record the correspondence between nodes in src and node in tgt                */
  std::unordered_map<klut::signal, klut::signal> src2tgt;
  int64_t area = 0;

  /* src: a 4-LUT with gates in the target technology libarary binded to each node */
  /* tgt: a 2-LUT ( essentially 4-LUT ), each node is associated with a NodeData,  */
  /* which contains the assigned stage and the gate type ( PI/SA/AS/AA ), denoted  */
  /* using 'value' since src is converted into tgt by decomposing each 4-LUT into  */
  /* 2-LUTs, for each node in src, there is a node in tgt implementing the same    */
  /* function                                                                      */
  klut tgt;
  src.foreach_pi( [&]( const klut::signal & src_pi ) 
  {
    klut::signal tgt_pi = tgt.create_pi();
    tgt.set_value(tgt_pi, NodeData(0, PI_GATE).value);
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
    DEBUG_PRINT("Processing node {0} ({1} out of {2})\r", src_node, ++ctr, num_node_src);
    if (ctr == num_node_src)
    {
      DEBUG_PRINT("\n");
    }
    if ( !src.has_binding( src_node ) )
    {
      return;
    }

    /* get gate type in the technology library */
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
    
    /* for each cell in the technology library, there is a primitive gates-based     */
    /* implementation, prepared in advance and ready for a look-up                   */
    LibEntry entry = entries.at(g.name);
    Node & root = nodemap[entry.hash];

    // fmt::print("Binding : {0}\n", entry.to_str());
    // fmt::print("NODE function : 0x{0:{2}x}, {1}\n", root.func, root.to_str(), std::max(1 << (src.fanin_size(src_node) - 2), 1) );

    // kitty::dynamic_truth_table func(src.fanin_size(src_node));
    // func._bits = root.func;
    // tgt._create_node(parent_tgt_1, parent_tgt_2);

    // _create_node( { a }, 2 )

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

    /* for each cell/compound gate in src, implement it using primitive gates in tgt */
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
        tgt.set_value(tgt_sig, NodeData(0, AS_GATE).value);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0x5555u, AS_GATE, { parent_tgt }));
        area += COSTS_MAP[fNOT];
        DEBUG_PRINT("ADDED NOT = {}\n", area);
      }
      else if (node.last_func == fAND)
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_and(parent_tgt_1, parent_tgt_2);
        tgt.set_value(tgt_sig, NodeData(0, SA_GATE).value);
        // fmt::print("Created node n{0} = AND({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0x8888u, SA_GATE, { parent_tgt_1, parent_tgt_2 }));
        area += COSTS_MAP[fAND];
        DEBUG_PRINT("ADDED AND = {}\n", area);
      }
      else if (node.last_func == fOR) 
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_or(parent_tgt_1, parent_tgt_2);
        tgt.set_value(tgt_sig, NodeData(0, SA_GATE).value);
        // fmt::print("Created node n{0} = OR({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0xEEEEu, SA_GATE, { parent_tgt_1, parent_tgt_2 }));
        OR_replacement_candidates.emplace(tgt_sig);
        area += COSTS_MAP[fOR];
        DEBUG_PRINT("ADDED OR  = {}\n", area);
      }
      else if (node.last_func == fCB)
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_or(parent_tgt_1, parent_tgt_2);
        tgt.set_value(tgt_sig, NodeData(0, AA_GATE).value);
        // fmt::print("Created node n{0} = CB({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0xEEEEu, AA_GATE, { parent_tgt_1, parent_tgt_2 }));
        area += COSTS_MAP[fCB];
        DEBUG_PRINT("ADDED CB  = {}\n", area);
      }
      else if (node.last_func == fXOR)
      {
        klut::signal parent_tgt_1 = node2tgt.at(node.parent_hashes.front());
        klut::signal parent_tgt_2 = node2tgt.at(node.parent_hashes.back());
        klut::signal tgt_sig = tgt.create_xor(parent_tgt_1, parent_tgt_2);
        tgt.set_value(tgt_sig, NodeData(0, AS_GATE).value);
        // fmt::print("Created node n{0} = XOR({1}, {2})\n", tgt_sig, parent_tgt_1, parent_tgt_2);
        node2tgt.emplace(hash, tgt_sig);
        tgt_sig_params.emplace(tgt_sig, Primitive<klut>(tgt_sig, 0x6666u, AS_GATE, { parent_tgt_1, parent_tgt_2 }));
        XOR_gates.push_back(tgt_sig);
        area += COSTS_MAP[fXOR];
        DEBUG_PRINT("ADDED XOR = {}\n", area);
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
    //std::vector<klut::signal> stack;
    //std::vector<klut::signal> visited;
    Primitive<klut> & xor_prim = tgt_sig_params.at(xor_sig);

    for (klut::signal const & fanin : xor_prim.fanins)
    {
      // Primitive<klut> & fanin_prim = tgt_sig_params.at(fanin);
      /* hazzards potentially produced by a CB can lead to incorrect behavior of its   */
      /* following XOR, such a CB should therefore be replaced by an OR */
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
    tgt.set_value(sig, NodeData(0, AA_GATE).value);
    prim.type = AA_GATE;
    area += COSTS_MAP[fCB];
    area -= COSTS_MAP[fOR];
    DEBUG_PRINT("REPLACED OR WITH CB = {}\n", area);
  }

  src.foreach_po([&](auto const & src_po)
  {
    tgt.create_po(src2tgt.at(src_po));
  } );
  
  return std::make_tuple(tgt, tgt_sig_params, area);
}

#if false
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
    bool                is_po;
    bool                is_const;
    bool                is_valid;

    NtkNode(const klut::signal & _id, const std::vector<klut::signal>& _fanins, const std::vector<klut::signal>& _fanouts, glob_phase_t _phase, uint8_t _type, klut::signal _ntk_sig, bool _is_spl, bool _is_pi, bool _is_const)
      : id(_id), fanins(_fanins), fanouts(_fanouts), phase(_phase), type(_type), ntk_sig(_ntk_sig), is_spl(_is_spl), is_pi(_is_pi), is_const(_is_const), is_valid(true) {}
    NtkNode(const klut::signal & _id)
      : id(_id), fanins({}), fanouts({}), phase(0), type(0), ntk_sig(0), is_spl(false), is_pi(false), is_const(false), is_valid(false) {}
  };
#endif

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
  std::string format() const
  {
    return fmt::format("Path from [{}]\n\tvia [{}]\n\tto [{}]\n", fmt::join(sources, ","), fmt::join(internals, ","), fmt::join(targets, ","));
  }

  std::vector<klut::signal> preds(const klut::signal & sig, const klut & ntk) const
  {
    if ( sources.count(sig) != 0)
    {
      return {};
    }

    std::vector<klut::signal> predecessors;

    ntk.foreach_fanin(sig, [&](const klut::signal & parent){
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
};

struct Standard_T1_CELL
{
  // 3 bits to store input negations
  unsigned in_phase : 3;
  /* truth tables of the 3 outputs under the current input phase */
  uint8_t sum_truth_table{};
  uint8_t carry_truth_table{};
  uint8_t cbar_truth_table{};

  Standard_T1_CELL( const unsigned in_phase_, const uint8_t sum_truth_table_, const uint8_t carry_truth_table_, const uint8_t cbar_truth_table_ )
    : in_phase( in_phase_ ), sum_truth_table( sum_truth_table_ ), carry_truth_table( carry_truth_table_ ), cbar_truth_table( cbar_truth_table_ ) {}
};

static std::array<Standard_T1_CELL, 8> standard_T1_cells = 
{
  Standard_T1_CELL( 0, 0x96, 0xe8, 0xfe ), 
  Standard_T1_CELL( 1, 0x69, 0xd4, 0xfd ), 
  Standard_T1_CELL( 2, 0x69, 0xb2, 0xfb ), 
  Standard_T1_CELL( 3, 0x96, 0x71, 0xf7 ), 
  Standard_T1_CELL( 4, 0x69, 0x8e, 0xef ), 
  Standard_T1_CELL( 5, 0x96, 0x4d, 0xdf ), 
  Standard_T1_CELL( 6, 0x69, 0xb2, 0xbf ), 
  Standard_T1_CELL( 7, 0x69, 0x17, 0x7f )
};

struct T1_OUTPUTS
{
  // 3 bits to store input negations
  unsigned in_phase          : 3;
  // whether carry (maj3) output is used with inverter, feel free to rename if other naming is more convenient for you
  bool     has_carry_inverted: 1;
  // whether cbar (or3) output is used with inverter, feel free to rename if other naming is more convenient for you
  bool     has_cbar_inverted : 1;
  bool     has_sum           : 1;
  // whether carry (maj3) output is used with DFF
  bool     has_carry         : 1;
  // whether cbar (or3) output is used with DFF
  bool     has_cbar          : 1;
  /* TTs of the 3 outputs */
  uint8_t   sum_truth_table{};
  uint8_t carry_truth_table{};
  uint8_t  cbar_truth_table{};
  /* usage of the 5 output ports */
  klut::node       sum_to;
  klut::node     carry_to;
  klut::node inv_carry_to;
  klut::node      cbar_to;
  klut::node  inv_cbar_to;

  // TODO : calculate cost, add mask for specific outputs
  T1_OUTPUTS() : in_phase(0u), has_carry_inverted(false), has_cbar_inverted(false), has_carry(false), has_cbar(false) {}

  // Explicit Constructor
  T1_OUTPUTS(const uint8_t in, const bool p_c, const bool p_cb, const bool h_c, const bool h_cb)
    : in_phase(in), has_carry_inverted(p_c), has_cbar_inverted(p_cb), has_carry(h_c), has_cbar(h_cb) {}

  /* Explicit constructor with output functions provided */
  T1_OUTPUTS( const uint8_t in, const bool p_c, const bool p_cb, const bool h_s, const bool h_c, const bool h_cb, const klut::node sum_node, 
              const klut::node carry_node, const klut::node inv_carry_node, const klut::node cbar_node, const klut::node inv_cbar_node )
    : in_phase( in ), has_carry_inverted( p_c ), has_cbar_inverted( p_cb ), has_sum( h_s ), has_carry( h_c ), has_cbar( h_cb ), 
      sum_to( sum_node ), carry_to( carry_node ), inv_carry_to( inv_carry_node ), cbar_to( cbar_node ), inv_cbar_to( inv_cbar_node ) {}
    
  // returns the TT of the sum output of the T1 cell
  /* Q: does tt information really matter? */
  // uint8_t sum_tt() const
  // {
  //   return (__builtin_popcount(in_phase) & 1) ? 0x69 : 0x96;
  // }
  // // returns the TT of the carry output of the T1 cell
  // uint8_t carry_tt(const bool output_phase = false) const
  // {
  //   TT3 tt_in;
  //   tt_in._bits = MAJ3;
  //   const uint32_t phase = in_phase | (output_phase << CUT_SIZE);
    
  //   TT3 tt_out = kitty::create_from_npn_config( std::make_tuple( tt_in, phase, perm ) );
  //   return static_cast<uint8_t>(tt_out._bits);
  // }
  // // returns the TT of the cbar output of the T1 cell
  // uint8_t cbar_tt(const bool output_phase = false) const
  // {
  //   TT3 tt_in;
  //   tt_in._bits = OR3;
  //   const uint32_t phase = in_phase | (output_phase << CUT_SIZE);
    
  //   TT3 tt_out = kitty::create_from_npn_config( std::make_tuple( tt_in, phase, perm ) );
  //   return static_cast<uint8_t>(tt_out._bits);
  // }
  // uint8_t cost() const
  // {
  //   return T1_COST // base cost (2x merger + T1 cell)
  //     + __builtin_popcount(in_phase) * (COSTS_MAP[fNOT] - COSTS_MAP[fDFF]) // cost of negating the inputs
  //     + COSTS_MAP[fDFF] * (has_carry + has_cbar)                          // cost of MAJ3 and OR3
  //     + COSTS_MAP[fNOT] * (has_carry_inverted + has_cbar_inverted)        // cost of ~(MAJ3) and ~(OR3)
  //     + COSTS_MAP[fSPL] * ((has_carry & has_carry_inverted) + (has_cbar & has_cbar_inverted)); //cost of splitters, if needed 
  // }
  void report() const
  {
    fmt::print("\t\tInput phases : {{ {}0,{}1,{}2 }}\n", ( ( in_phase & 1 ) ? "~" : " " ), ( ( in_phase >> 1 & 1 ) ? "~" : " " ), ( ( in_phase >> 2 & 1 ) ? "~" : " " ) );
    fmt::print("\t\t      Sum ( Node {} ): {0:08b} ( 0x{0:02x} )\n",       sum_to,    sum_truth_table );
    fmt::print("\t\t    Carry ( Node {} ): {0:08b} ( 0x{0:02x} )\n",     carry_to,  carry_truth_table );
    fmt::print("\t\tInv carry ( Node {} ): {0:08b} ( 0x{0:02x} )\n", inv_carry_to, ~carry_truth_table );
    fmt::print("\t\t     Cbar ( Node {} ): {0:08b} ( 0x{0:02x} )\n",      cbar_to,   cbar_truth_table );
    fmt::print("\t\tInv carry ( Node {} ): {0:08b} ( 0x{0:02x} )\n", inv_carry_to,  ~cbar_truth_table );

  }
};

template <typename Ntk>
glob_phase_t latest_fanin_phase(const Ntk & ntk, const typename Ntk::signal & node, const uint8_t n_phases, const uint8_t type, const bool verbose)
{
  bool valid = false;
  uint32_t phase = 0u;
  ntk.foreach_fanin(node, [&] ( const typename Ntk::signal & parent )
  {
    if ( ntk.is_constant( parent ) )
    {
      DEBUG_PRINT("\t\tfanin {} is constant, skipping...\n", parent);
      return;
    }
    valid = true;
    NodeData parent_data = ntk.value( parent );

    // if d==false, SA gate can be directly connected to fanin
    //cannot directly connect PI (convention)
    //cannot directly connect split signal 
    //cannot directly connect SA/AA gates to SA gates
    int d = (type == SA_GATE) && ( (ntk.is_pi(parent)) || (ntk.fanout_size(parent) > 1) || (parent_data.type != AS_GATE) );

    phase = generic_max( phase, parent_data.sigma + d );

    if (verbose)
    {
      unsigned int gp = static_cast<unsigned int>(parent_data.sigma);
      fmt::print("\t\tfanin {} ɸ={} [S={}, φ={}]\n", parent, gp, gp/n_phases, gp%n_phases);
    }
  });
  assert(valid);

  // DEBUG_PRINT("\t{} GATE {} placed at ɸ={} [S={}, φ={}]\n",  node, phase, phase/n_phases, phase%n_phases);

  return phase;
}

void greedy_ntk_assign(const klut & ntk, const std::unordered_map<klut::signal, Primitive<klut>> & sig_params, const uint8_t n_phases, const std::unordered_map<unsigned int, unsigned int> & phase_assignment, const bool verbose = false)
{
  mockturtle::topo_view<klut> ntk_topo ( ntk );

  // ntk.set_value(sig, node_data.value);

  ntk_topo.foreach_node([&] ( const klut::signal & node ) 
  {
    if ( ntk_topo.is_constant( node ) )
    {
      return;
    }

    if ( ntk_topo.is_pi( node ) )
    {
      NodeData node_data;
      node_data.type = AS_GATE;
      auto ct = phase_assignment.count(node);
      node_data.sigma = (ct != 0) ? phase_assignment.at( node ) : 0;
      ntk.set_value(node, node_data.value);

      DEBUG_PRINT("PI {} placed at ɸ=0 [S=0, φ=0]\n", node);
      return;
    }

    // const Primitive<klut> & node_params = sig_params.at( node );
    // if (verbose) fmt::print("{} GATE {}:\n", GATE_TYPE.at(node_params.type), node);

    uint8_t node_type = NodeData( ntk.value( node ) ).type;
    // const Primitive<klut> & prim = sig_params.at( node );
    DEBUG_PRINT("{} GATE {}:\n", GATE_TYPE.at(node_type), node);

    NodeData node_data;
    node_data.type = node_type;

    auto ct = phase_assignment.count(node);
    if (ct != 0) // there is a precalculated phase assignment
    {
      node_data.sigma = phase_assignment.at( node );
      ntk.set_value(node, node_data.value);
    }
    else
    {
      if ( node_type == AA_GATE )
      {
        node_data.sigma = latest_fanin_phase(ntk_topo, node, n_phases, AA_GATE, verbose);
        uint32_t sigma = static_cast<uint32_t>(node_data.sigma);
        DEBUG_PRINT("\tAA GATE {} placed at ɸ={} [S={}, φ={}]\n", node, sigma, sigma/n_phases, sigma%n_phases);
      }
      else if ( node_type == SA_GATE )
      {
        node_data.sigma = latest_fanin_phase(ntk_topo, node, n_phases, SA_GATE, verbose);
        uint32_t sigma = static_cast<uint32_t>(node_data.sigma);
        DEBUG_PRINT("\tSA GATE {} placed at ɸ={} [S={}, φ={}]\n", node, sigma, sigma/n_phases, sigma%n_phases);
      }
      else if ( node_type == AS_GATE )
      {
        node_data.sigma = latest_fanin_phase(ntk_topo, node, n_phases, SA_GATE, verbose) + 1;
        uint32_t sigma = static_cast<uint32_t>(node_data.sigma);
        DEBUG_PRINT("\tAS GATE {} placed at ɸ={} [S={}, φ={}]\n", node, sigma, sigma/n_phases, sigma%n_phases);
      }
      else 
      {
        std::vector<klut::signal> fanins;
        ntk.foreach_fanin(node, [&](const klut::signal & fanin){ fanins.push_back(fanin); });
        fmt::print("\t GATE {0} : fanins[{1}], type[{2}], func[{3}=0x{3:x}=0b{3:b}]\n", node, fmt::join(fanins, ","), node_type, ntk.node_function(node)._bits[0]);
        throw "Unsupported gate type!";
      }
    }
  });
}


std::unordered_map<klut::signal, glob_phase_t> greedy_assign(const klut & ntk, const std::unordered_map<klut::signal, Primitive<klut>> & sig_params,  const bool verbose = false) //const uint8_t n_phases,
{
  mockturtle::topo_view<klut> ntk_topo ( ntk );
  std::unordered_map<klut::signal, glob_phase_t> glob_phase;
  glob_phase.reserve(ntk_topo.size());

  printUnorderedMap(sig_params);

  ntk_topo.foreach_node([&] ( const klut::signal & node ) 
  {
    if ( ntk_topo.is_constant( node ) )
    {
      DEBUG_PRINT("Skipping constant {}\n", node);
      return;
    }
    if ( ntk_topo.is_pi( node ) )
    {
      glob_phase[node] = 0u;
      DEBUG_PRINT("PI {} placed at ɸ=0 [S=0, φ=0]\n", node);
      return;
    }
    const Primitive<klut> & node_params = sig_params.at( node );
    DEBUG_PRINT("{} GATE {}:\n", GATE_TYPE.at(node_params.type), node);

    if ( node_params.type == AA_GATE )
    {
      // place at the same phase as the latest fanin
      glob_phase_t phase = 0u;
      ntk_topo.foreach_fanin(node, [&] ( const klut::signal & parent )
      {
        if ( ntk_topo.is_constant( parent ) )
        {
        DEBUG_PRINT("\t\tfanin {} is constant, skipping...\n", parent);
        return;
        }
        phase = std::max(phase, glob_phase.at( parent ));
        DEBUG_PRINT("\t\tfanin {} ɸ={}\n", parent, glob_phase.at( parent ));
        // fmt::print("\t\tfanin {} ɸ={} [S={}, φ={}]\n", parent, gp, gp/n_phases, gp%n_phases);
      });
      glob_phase[node] = phase;
    if (verbose)
    {
      fmt::print("\tAA GATE {} placed at ɸ={}\n", node, phase);
      // fmt::print("\tAA GATE {} placed at ɸ={} [S={}, φ={}]\n", node, phase, phase/n_phases, phase%n_phases);
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
        // fmt::print("\t\tfanin: {4} {0} ɸ={1} [S={2}, φ={3}], fanout={5} \n",
        fmt::print("\t\tfanin: {2} {0} ɸ={1}, fanout={3} \n", 
            parent, 
            gp, 
            // gp/n_phases, 
            // gp%n_phases, 
            ntk_topo.is_pi(parent) ? "PI" : GATE_TYPE.at( sig_params.at(parent).type ),
            ntk_topo.fanout_size(parent)
        );
        }
    });
    glob_phase[node] = phase;
    if (verbose)
    {
        // fmt::print("\tSA GATE {} placed at ɸ={} [S={}, φ={}]\n", node, phase, phase/n_phases, phase%n_phases);
        fmt::print("\tSA GATE {} placed at ɸ={}\n", node, phase);
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
        // fmt::print("\t\tfanin {} ɸ={} [S={}, φ={}]\n", parent, gp, gp/n_phases, gp%n_phases);
        fmt::print("\t\tfanin {} ɸ={}\n", parent, gp);
        }
    });
    phase++;
    glob_phase[node] = phase;
    if (verbose)
    {
        fmt::print("\tAS GATE {} placed at ɸ={}\n", node, phase);
        // fmt::print("\tAS GATE {} placed at ɸ={} [S={}, φ={}]\n", node, phase, phase/n_phases, phase%n_phases);
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

bool phase_ntk_comparison(const klut::signal & a, const klut::signal & b, const klut & ntk )
{
  if (ntk.is_constant(a))
  {
    return true;
  }
  else if (ntk.is_constant(b))
  {
    return false;
  }
  NodeData a_data = ntk.value(a);
  NodeData b_data = ntk.value(b);

  return a_data.sigma < b_data.sigma;
}

// Function to insert splitter nodes in a KLUT-based network (klut = k-input lookup table network).
void splitter_ntk_insertion(klut & ntk, const bool verbose = false)
{
  // Lambda function for comparing the phases of two signals in the network.
  auto phase_comp = [&](const klut::signal & a, const klut::signal & b)
  {
    return phase_ntk_comparison(a, b, ntk);
  };

  // Create a view of the network that provides access to fanout information.
  auto ntk_fo = mockturtle::fanout_view<klut>(ntk);

  // For each node in the fanout view:
  ntk_fo.foreach_node([&](const klut::signal & node)
  {
    // Get the number of fanouts for the current node.
    uint32_t fo_size = ntk_fo.fanout_size(node);

    DEBUG_PRINT("\t[NODE {}] FANOUT SIZE = {}\n", node, fo_size);
    // If the current node is a constant or it has fanout ≤ 1, skip to the next node.
    if (ntk_fo.is_constant(node) || fo_size <= 1)
    {
      return;
    }

    // Populate the fanouts vector.
    std::vector<klut::signal> fanouts;
    fanouts.reserve(fo_size);
    ntk_fo.foreach_fanout(node, [&](const klut::signal & fo_node)
    {
      fanouts.push_back(fo_node);
      DEBUG_PRINT("\t\t[NODE {}] ADDING FANOUT\n", node, fo_node);
    });

    // Fix the fanout count (bugged fanouts_size()?)
    if ( fanouts.size() != fo_size )
    {
      ntk._storage->nodes[node].data[0].h1 = fanouts.size();
    }

    // Sort the fanouts using the phase comparison function.
    std::sort(fanouts.begin(), fanouts.end(), phase_comp);
    DEBUG_PRINT("\t[NODE {}] SORTED FANOUTS:\n", node);
    printVector(fanouts, 2);

    // Create [fo_size - 1] splitter nodes.
    klut::signal last_spl = node;
    std::vector<klut::signal> splitters;
    splitters.reserve( fanouts.size() - 1 );
    for (auto it = fanouts.begin(); it < fanouts.end() - 1; it++)
    {
      DEBUG_PRINT("\t\t[NODE {}] LAST SPL: {}\n", node, last_spl);
      // Copy sigma and type data from the first fanout for the splitter node data.
      NodeData spl_data = ntk.value(*it);
      // Change gate type to AA for the splitter node.
      spl_data.type = AA_GATE;
      // Create a new splitter node connected to 'last_spl'.

      DEBUG_PRINT("\t\t[NODE {}] CREATING SPL FOR {}\n", node, *it);
      DEBUG_PRINT("\t\t[NODE {}] LAST_SPL {} FANOUT BEFORE: {}\n", node, last_spl, ntk.fanout_size(last_spl));
      const klut::signal spl = ntk._create_node({last_spl}, 2, last_spl);
      ntk.set_value(spl, spl_data.value);
      DEBUG_PRINT("\t\t[NODE {}] CREATED SPL {}\n", node, spl);
      DEBUG_PRINT("\t\t[NODE {}] LAST_SPL {} FANOUT AFTER: {}\n", node, last_spl, ntk.fanout_size(last_spl));
      DEBUG_PRINT("\t\t[NODE {}] SPL {} FANIN: {}\n", node, spl, ntk.fanin_size(spl));

      // auto n = ntk._storage->nodes[*it];
      // auto& preds = n.children;

      // Update the connections to reflect the splitter's presence.
      DEBUG_PRINT("\t\t[NODE {}, LAST_SPL {}] UPDATING CONNECTIONS\n", node, last_spl);
      // for (auto& pred : preds)
      for (auto pred_it = ntk._storage->nodes[*it].children.begin();
                pred_it < ntk._storage->nodes[*it].children.end(); pred_it++ )
      {
        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PRED {}\n", node, last_spl, pred_it->data);
        if (pred_it->data == node)
        {
          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] RECORDING OLD PREDS\n", node, last_spl);
          // Store the previous connections.
          std::vector<klut::signal> old_preds(ntk._storage->nodes[*it].children.size());
          std::transform(ntk._storage->nodes[*it].children.begin(), ntk._storage->nodes[*it].children.end(), old_preds.begin(), [](auto c) { return c.index; });
          printVector(old_preds, 4);

          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS BEFORE\n", node, last_spl);
          for (const auto& entry : ntk._storage->nodes[*it].children)  { fmt::print("\t\t\t\t{}\n", entry.data); }
          pred_it->data = spl;                             // Replace the connection with the splitter.
          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS AFTER\n", node, last_spl);
          for (const auto& entry : ntk._storage->nodes[*it].children)  { fmt::print("\t\t\t\t{}\n", entry.data); }

          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] SPL {} FANOUT BEFORE: {}\n", node, last_spl, spl, static_cast<int>(ntk._storage->nodes[spl].data[0].h1));
          DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(spl));
          ntk._storage->nodes[spl].data[0].h1++;  // Increment fan-out of the splitter.
          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] SPL {} FANOUT  AFTER: {}\n", node, last_spl, spl, static_cast<int>(ntk._storage->nodes[spl].data[0].h1));
          DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(spl));

          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] NODE {} FANOUT BEFORE: {}\n", node, last_spl, node, static_cast<int>(ntk._storage->nodes[node].data[0].h1));
          DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(node));
          ntk._storage->nodes[node].data[0].h1--; // Decrement fan-out of the current node.
          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] NODE {} FANOUT  AFTER: {}\n", node, last_spl, node, static_cast<int>(ntk._storage->nodes[node].data[0].h1));
          DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(node));

          // // Notify listeners of the modification.
          // for (auto const& fn : ntk._events->on_modified)
          // {
          //   (*fn)(*it, old_preds);
          // }
        }
      }
      last_spl = spl;
      splitters.push_back(spl);
    }

    // Process the last fanout.
    // auto& preds = ntk._storage->nodes[fanouts.back()].children;
    DEBUG_PRINT("\t\t[NODE {}, LAST_SPL {}] UPDATING CONNECTIONS TO LAST FANOUT {}\n", node, last_spl, fanouts.back());
    // for (auto& pred : preds)
    for (auto pred_it = ntk._storage->nodes[fanouts.back()].children.begin();
              pred_it < ntk._storage->nodes[fanouts.back()].children.end(); pred_it++ )
    {
      if (pred_it->data == node)
      {
        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] RECORDING OLD PREDS\n", node, last_spl);
        // Store the previous connections.
        std::vector<klut::signal> old_preds(ntk._storage->nodes[fanouts.back()].children.size());
        std::transform(ntk._storage->nodes[fanouts.back()].children.begin(), ntk._storage->nodes[fanouts.back()].children.end(), old_preds.begin(), [](auto c) { return c.index; });
        printVector(old_preds, 4);

        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS BEFORE\n", node, last_spl);
        for (const auto& entry : ntk._storage->nodes[fanouts.back()].children)  { fmt::print("\t\t\t\t{}\n", entry.data); }
        pred_it->data = last_spl;                            // Replace the connection with the last splitter.
        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS AFTER\n", node, last_spl);
        for (const auto& entry : ntk._storage->nodes[fanouts.back()].children)  { fmt::print("\t\t\t\t{}\n", entry.data); }     


        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] SPL {} FANOUT BEFORE: {}\n", node, last_spl, last_spl, static_cast<int>(ntk._storage->nodes[last_spl].data[0].h1));
        DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(last_spl));
        ntk._storage->nodes[last_spl].data[0].h1++; // Increment fan-out of the last splitter.
        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] SPL {} FANOUT  AFTER: {}\n", node, last_spl, last_spl, static_cast<int>(ntk._storage->nodes[last_spl].data[0].h1));
        DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(last_spl));

        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] NODE {} FANOUT BEFORE: {}\n", node, last_spl, node, static_cast<int>(ntk._storage->nodes[node].data[0].h1));
        DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(node));
        ntk._storage->nodes[node].data[0].h1--;     // Decrement fan-out of the current node.
        DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] NODE {} FANOUT  AFTER: {}\n", node, last_spl, node, static_cast<int>(ntk._storage->nodes[node].data[0].h1));
        DEBUG_PRINT("\t\t\t\t Call: {}\n", ntk.fanout_size(node));

        // Notify listeners of the modification.
        for (auto const& fn : ntk._events->on_modified)
        {
          (*fn)(fanouts.back(), old_preds);
        }
      }
    }

    auto fo_ntk { mockturtle::fanout_view( ntk ) };
    for (const auto & node : fanouts)
    {
      DEBUG_PRINT("\t\t\t node: {}\n", node);
      
      fo_ntk.foreach_fanin( node, [&](const klut::signal & fi_node)
      {
        DEBUG_PRINT("\t\t\t\t fanin : {}\n", fi_node);
      });
      fo_ntk.foreach_fanout( node, [&](const klut::signal & fo_node)
      {
        DEBUG_PRINT("\t\t\t\t fanout: {}\n", fo_node);
      });
    }
    for (const auto & node : splitters)
    {
      DEBUG_PRINT("\t\t\t spl : {}\n", node);
      
      fo_ntk.foreach_fanin( node, [&](const klut::signal & fi_node)
      {
        DEBUG_PRINT("\t\t\t\t fanin : {}\n", fi_node);
      });
      fo_ntk.foreach_fanout( node, [&](const klut::signal & fo_node)
      {
        DEBUG_PRINT("\t\t\t\t fanout: {}\n", fo_node);
      });
    }
    // Ensure that the current node's fan-out count is now 1 (since all other fanouts have been replaced by splitters).
    assert(ntk._storage->nodes[node].data[0].h1 == 1);
  });
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
  klut::signal fanin;
  klut::signal fanout;
  uint32_t sigma;
  std::unordered_set<uint64_t> parent_hashes;
  uint64_t hash;

  DFF_var(klut::signal _fanin, klut::signal _fanout, uint32_t _sigma)
      : fanin(_fanin), fanout(_fanout), sigma(_sigma), parent_hashes({}) 
      {
        hash = calculate_hash(_sigma, _fanout, _fanin);
      }

  DFF_var(klut::signal _fanin, klut::signal _fanout, uint32_t _sigma, std::unordered_set<uint64_t> _parent_hashes, uint64_t _hash)
      : fanin(_fanin), fanout(_fanout), sigma(_sigma), parent_hashes(_parent_hashes), hash(_hash) {}

  DFF_var(const DFF_var& other)
      : fanin(other.fanin), fanout(other.fanout), sigma(other.sigma), parent_hashes(other.parent_hashes), hash(other.hash) {}

  std::string str() 
  {
      return fmt::format("var_{}_{}_{}", fanin, fanout, sigma);
  }
};

uint64_t dff_hash(klut::signal _fanin, klut::signal _fanout, uint32_t _sigma)
{
  return ( (uint64_t)_fanin << 40 ) | ( (uint64_t)_fanout << 16 ) | _sigma;
}

struct DFF_registry
{
  std::unordered_map<uint64_t, DFF_var> variables;

  DFF_var & at(node_t _fanin, node_t _fanout, glob_phase_t _sigma)
  {
    return variables.at( dff_hash(_fanin, _fanout, _sigma) );
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
      return fmt::format( "var_{}_{}_{}.Not()", dff.fanin, dff.fanout, dff.sigma );
    }
    return fmt::format( "var_{}_{}_{}", dff.fanin, dff.fanout, dff.sigma );
  }
};


std::vector<Path> extract_paths(const klut & ntk, bool verbose = false)
{
  DEBUG_PRINT("\t[i] ENTERED FUNCTION extract_paths\n");
  std::vector<Path> paths;

  ntk.foreach_node([&](const klut::signal & fo_node)
  {
    DEBUG_PRINT("\t\t[i] PROCESSING NODE {}\n", fo_node);
    if (ntk.is_constant(fo_node) || ntk.is_pi(fo_node))
    {
      DEBUG_PRINT("\t\t\t[NODE {}] the node is IS CONSTANT\n", fo_node);
      return;
    }
    NodeData fo_node_data = ntk.value(fo_node);
    if (fo_node_data.type != AS_GATE && fo_node_data.type != SA_GATE) 
    {
      DEBUG_PRINT("\t\t\t[NODE {}] the node is AA, skipping\n", fo_node);
      return;
    }
    // at this point, the node should be AS/SA
    DEBUG_PRINT("\t\t[NODE {}] the node is AS/SA, continuing...\n", fo_node);

    // Create a separate path for each fanin of the node
    std::vector<Path> aa_paths; 
    aa_paths.reserve( ntk.fanin_size(fo_node) );

    ntk.foreach_fanin(fo_node, [&](const klut::signal & fi_node)
    {
      DEBUG_PRINT("\t\t\t[NODE {}] processing fanin {}\n", fo_node, fi_node);
      // Create a path object with only a target
      Path node_path;
      node_path.targets.emplace( fo_node );

      std::vector<klut::signal> stack { fi_node };
      
      DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}] created stack\n", fo_node, fi_node);
      
      std::set<klut::signal> seen;
      while (!stack.empty())
      {
        DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}] stack contents:\n", fo_node, fi_node);
        printVector(stack, 4);

        const klut::signal & n = stack.back();
        stack.pop_back();

        DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}]: Analyzing node {}\n", fo_node, fi_node, n);

        // A constant does not have any effect on the DFF placement, we can skip it
        if ( ntk.is_constant( n ) )
        {
          continue;
        }        
        const NodeData n_data { ntk.value(n) };
        // Found a source of the path, add to sources, do not continue traversal
        if ( ntk.is_pi(n) || n_data.type == AS_GATE || n_data.type == SA_GATE )
        {
          DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}]: node {} is a source \n", fo_node, fi_node, n);
          node_path.sources.emplace( n );
        }
        // Found AA gate, add to internal nodes, add parents for further traversal
        else if ( n_data.type == AA_GATE )
        {
          DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}]: node is INTERNAL adding fanins \n", fo_node, fi_node, n);
          node_path.internals.emplace( n );

          ntk.foreach_fanin(n, [&](const klut::signal & sig){
            stack.push_back( sig );
          });
        }
        else
        {
          DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}]: Signal {}: {} is not recognized \n", fo_node, fi_node, n, GATE_TYPE.at( n_data.type ));
          throw "Unsupported case";
        }
        seen.emplace( n );
      }

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
    });
  });
  return paths;
}

/// @brief Create binary variables for DFF placement in a given path
/// @param path - a path object to insert DFFs into
/// @param NR - unordered_map of NtkNode objects 
/// @param n_phases - # of phases
/// @param verbose - prints debug messages if set to *true*
/// @return 
std::tuple<DFF_registry, uint64_t, std::vector<uint64_t>> dff_vars_single_paths(const Path & path, const klut & ntk, const uint8_t n_phases, bool verbose = false)
{
  DFF_registry DFF_REG;
  std::vector<uint64_t> required_SA_DFFs;

  std::vector<std::tuple<klut::signal, uint64_t>> stack;
  for (const klut::signal & tgt : path.targets)
  {
    stack.emplace_back(tgt, 0);
  }
  DEBUG_PRINT("[DFF] Target nodes: {}\n", fmt::join(path.targets, ","));

  auto precalc_ndff = 0u;
  
  while (!stack.empty())
  { 
    if (verbose)
    {
      DEBUG_PRINT("STACK :\n");
      for (const auto & [ fo_node, earliest_child_hash ] : stack)
      {
        NodeData fo_data { ntk.value( fo_node ) };
        if (earliest_child_hash != 0)
        {
          DEBUG_PRINT("\t{}({})[{}], {}\n", GATE_TYPE.at((int)fo_data.type), fo_node, (int)fo_data.sigma, DFF_REG.at(earliest_child_hash).str());
        }
        else
        {
          DEBUG_PRINT("\t{}({})[{}]\n", GATE_TYPE.at((int)fo_data.type), fo_node, (int)fo_data.sigma);
        }
      }
    }
    // fixing this stupid clang bug with structured bindings
    auto node_tuple = stack.back();
    const klut::signal fo_node = std::get<0>(node_tuple);
    const uint64_t earliest_child_hash = std::get<1>(node_tuple);
    // const auto & [ fo_node, earliest_child_hash ] = node_tuple ;
    stack.pop_back();
    NodeData fo_data { ntk.value( fo_node ) };

    uint32_t latest_sigma = fo_data.sigma - (fo_data.type == AS_GATE);
    DEBUG_PRINT("[DFF] Analyzing child: {}({})[{}]\n", GATE_TYPE.at(fo_data.type), fo_node, (int)fo_data.sigma);

    ntk.foreach_fanin(fo_node, [&](const klut::signal & fi_node){
      NodeData fi_data { ntk.value( fi_node ) };
      uint32_t earliest_sigma = fi_data.sigma + (fi_data.type != AA_GATE);

      DEBUG_PRINT("\t[DFF] Analyzing parent: {}({})[{}]\n", GATE_TYPE.at(fi_data.type), fi_node, (int)fi_data.sigma);

      // check if the chain is straight - #DFF is just floor(delta-phase), no need to create the dff vars
      if (fo_data.type != AA_GATE && fi_data.type != AA_GATE)
      {
        // special case when an AS gate feeds directly into SA gate
        if (fo_data.sigma == fi_data.sigma)
        {
          DEBUG_PRINT("\t[DFF] Straight chain: AS{} -> SA{}\n", fi_node, fo_node);
          // do nothing, no additional DFFs needed
          assert(fo_data.type == SA_GATE && fi_data.type == AS_GATE && ntk.fanout_size(fi_node) == 1);
        }
        else
        {
          DEBUG_PRINT("\t[DFF] Straight chain: {}[{}] -> {}[{}]\n", GATE_TYPE.at(fi_data.type), (int)fi_data.sigma, GATE_TYPE.at(fo_data.type), (int)fo_data.sigma);
          // straight chain, just floor the difference!
          precalc_ndff += (fo_data.sigma - fi_data.sigma)/n_phases + (fo_data.type == SA_GATE); //extra DFF before SA gate
        }
        return;
      }

      DEBUG_PRINT("\t[DFF] Non-straight chain: {}[{}] -> {}[{}]\n", GATE_TYPE.at(fi_data.type), (int)fi_data.sigma, GATE_TYPE.at(fo_data.type), (int)fo_data.sigma);
      std::vector<uint64_t> out_hashes;
      DEBUG_PRINT("\tAdding new DFFs [reg size = {}]\n", DFF_REG.variables.size());

      for (glob_phase_t sigma = earliest_sigma; sigma <= latest_sigma; ++sigma)
      {
        uint64_t new_hash = DFF_REG.add(fi_node, fo_node, sigma);
        out_hashes.push_back(new_hash);
        DEBUG_PRINT("\tAdded new DFFs at phase {} [reg size = {}]\n", sigma, DFF_REG.variables.size());
      }
      DEBUG_PRINT("\tConnecting new DFFs\n");
      for (auto i = 1u; i < out_hashes.size(); ++i)
      {
        DFF_var & dff = DFF_REG.at( out_hashes[i] );
        dff.parent_hashes.emplace(out_hashes[i-1]);
      }
      if (fo_data.type == SA_GATE)
      {
        assert( !out_hashes.empty() );
        required_SA_DFFs.push_back(out_hashes.back());
      }

      uint64_t earliest_hash = (out_hashes.empty()) ? earliest_child_hash : out_hashes.front();
      // if the node is internal, connect with the fanout phase
      if (fo_data.type == AA_GATE && !out_hashes.empty() && earliest_hash != 0 && earliest_child_hash != 0)
      {
        DFF_var & child_dff = DFF_REG.at( earliest_child_hash );
        DEBUG_PRINT("\tPrior node is {}[{}]\n", child_dff.str(), (int)child_dff.sigma); 
        // assert(child_dff.fanin == fo_id);
        child_dff.parent_hashes.emplace( out_hashes.back() );
      }
      if (fi_data.type == AA_GATE)
      {
        stack.emplace_back( fi_node, earliest_hash );
        DEBUG_PRINT("\tEmplacing {}({})[{}], {}\n", GATE_TYPE.at(fi_data.type), fi_node, (int)fi_data.sigma, (earliest_hash!=0)?DFF_REG.at(earliest_hash).str():"");
      }
    });
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
    if (dff.sigma == head_dff.sigma) // add to the same section
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

void write_snakes(const std::vector<Snake> & snakes, DFF_registry & DFF_REG, const std::vector<uint64_t> & required_SA_DFFs, const std::string cfg_name, uint8_t n_phases, bool verbose = false)
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
      DEBUG_PRINT("New single phase conflict : {}≤1\n", fmt::join(vars, "+"));
      vars_bucket.push_back(fmt::format(vars.size()>1?"({})":"{}", fmt::join(vars, "+")));
      if (vars.size() > 1)
      {
        spec_file << fmt::format("PHASE,{}\n", fmt::join(vars, ","));
      }
    }
    std::reverse(vars_bucket.begin(), vars_bucket.end());
    DEBUG_PRINT("New buffer requirement : ({})\n", fmt::join(vars_bucket, "|"));
    if (vars_bucket.size() == n_phases)
    {
      spec_file << fmt::format("BUFFER,{}\n", fmt::join(vars_bucket, ","));
    }
  }

  for (const uint64_t & hash : required_SA_DFFs)
  {
    DEBUG_PRINT("New SA_REQUIRED : {}\n", DFF_REG.str( hash ));
    spec_file << fmt::format("SA_REQUIRED,{}\n", DFF_REG.str( hash ));
  }
}

/// @brief 
/// @param path 
/// @param NR 
/// @param DFF_REG 
/// @param n_phases 
/// @param verbose 
/// @return 
std::vector<Snake> sectional_snake(const Path & path, klut & ntk,  DFF_registry & DFF_REG, uint8_t n_phases, bool verbose = false)
{
  std::vector<Snake> out_snakes; 
  std::vector<Snake> stack;
  
  DEBUG_PRINT("[i]: Starting extraction of worms \n");
  // get all DFFs 
  for (const auto & [hash, dff]: DFF_REG.variables)
  {
    NodeData fo_data { ntk.value( dff.fanout ) };
    auto fanout_sigma = fo_data.sigma - ( fo_data.type == AS_GATE );
    auto it = std::find(path.targets.begin(), path.targets.end(), dff.fanout);
    if (it != path.targets.end() && (dff.sigma >= fanout_sigma ))
    {
      stack.emplace_back( hash );
    }
  }
  
  while (!stack.empty())
  {
    DEBUG_PRINT("[i] Stack size is {} \n", stack.size());
    Snake snake = stack.back();
    stack.pop_back();

    DEBUG_PRINT("\t[i] The snake has {} sections\n", snake.sections.size());
    uint64_t hash = snake.sections.back().back();
    DFF_var & dff = DFF_REG.at( hash );

    // fmt::print("\tCurrent worm size {}, between phases {} and {} \n", worm.size(), DFF_REG.str(worm.front()), DFF_REG.str(worm.back()));


    DEBUG_PRINT("\t\t[i] The DFF {} has {} parents\n", DFF_REG.at( hash ).str(),  dff.parent_hashes.size() );

    bool returned_current_snake = false;
    for (const uint64_t parent_hash : dff.parent_hashes)
    {
      Snake snake_copy = snake; 
      DEBUG_PRINT("\t\t[i] Advancing towards fanin {}\n", DFF_REG.at( parent_hash ).str() );
      bool status = snake_copy.append(parent_hash, DFF_REG, n_phases);
      DEBUG_PRINT((status) ? "\t\t\tAdded new section!\n" :"\t\t\tExtended existing section!\n"  );
      DEBUG_PRINT("\t\t\tThe new length is {}\n", snake_copy.sections.size() );
      
      stack.push_back( snake_copy );
      if (status && !returned_current_snake && snake_copy.sections.size() == n_phases)
      {
        DEBUG_PRINT("\t\tAdding the snake to the output\n");
        out_snakes.push_back(snake);
        returned_current_snake = true;
      }
    }
  }
  return out_snakes;
}

void write_klut_specs(const klut & ntk, const std::unordered_map<klut::signal, Primitive<klut>> & sig_params, const std::unordered_map<klut::signal, glob_phase_t> &glob_phase, const std::string & filename)
{
  auto ntk_fo = mockturtle::fanout_view<klut, false>( ntk );
  std::ofstream spec_file (filename);
  // ntk_fo.foreach_pi([&] (const auto & sig) 
  // {
  //   std::vector<klut::node> fanouts;
  //   ntk_fo.foreach_fanout(sig, [&](const klut::node fo)
  //   {
  //     fanouts.push_back(fo);
  //   });
  //   spec_file << fmt::format("PI {0},{1}\n", sig, fmt::join(fanouts, "|"));
  // } );
  for (const auto & [sig, prim] : sig_params)
  {
    std::vector<klut::node> fanouts;
    ntk_fo.foreach_fanout(sig, [&](const klut::node fo)
    {
      fanouts.push_back(fo);
    });

    uint8_t prim_type = prim.fanins.empty() ? PI_GATE : prim.type;
    spec_file << fmt::format("{0},{1},{2},{3},{4},{5}\n", sig, prim.func, prim_type, fmt::join(prim.fanins, "|"), glob_phase.at(sig), fmt::join(fanouts, "|"));
  }
}

std::tuple<int, std::unordered_map<unsigned int, unsigned int>, std::string>  cpsat_macro_opt(const std::string & cfg_name, uint8_t n_phases) 
{
  std::string command = fmt::format("{} {} {} {}", PYTHON_EXECUTABLE, PYTHON_PHASE_ASSIGNMENT, n_phases, cfg_name);
  
  std::string pattern = "Objective value: (\\d+)";

  // Run the command and capture its output
  FILE* pipe = popen(command.c_str(), "r");
  if (!pipe) 
  {
    std::cerr << "Error running the command." << std::endl;
    throw;
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
    throw;
  }

  // Parse output
  std::istringstream iss(output);

  std::string line;
  std::string solve_status;
  int objective_value;
  std::unordered_map<unsigned int, unsigned int> key_value_pairs;

  while (std::getline(iss, line))
  {
    if (line.find("Solve status:") != std::string::npos) 
    {
      if (line.find("OPTIMAL") || line.find("FEASIBLE"))
      {
        solve_status = "SUCCESS";
        break;
      }
      else
      {
        solve_status = "UNKNOWN";
        return {0, {}, ""};
      }
    }
  }

  // Parse the second line (Objective value)
  std::getline(iss, line);
  if (line.find("Objective value: ") != std::string::npos) 
  {
      objective_value = std::stoi(line.substr(17));
  } 
  else 
  {
      // Handle missing or incorrect format for objective value
      std::cerr << "Error: Objective value not found or invalid format." << std::endl;
      return {0, {}, ""};
  }

  // Parse the key-value pairs in subsequent lines
  while (std::getline(iss, line)) 
  {
      std::istringstream line_stream(line);
      unsigned int key, value;
      char colon;
      if (line_stream >> key >> colon >> value && colon == ':') 
      {
          key_value_pairs[key] = value;
      } 
      else 
      {
          // Handle incorrect format for key-value pairs
          std::cerr << "Error: Invalid format for key-value pairs." << std::endl;
          return {0, {}, ""};
      }
  }

  return {objective_value, key_value_pairs, solve_status};
}

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
  std::string command = fmt::format("{} {} {}", PYTHON_EXECUTABLE, PYTHON_DFF_PLACEMENT, cfg_name);
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

void config_t1_ports_connection( klut const& ntk, const bool output_phase, 
                                 const uint32_t target_node, 
                                 klut::node& port, klut::node& inv_port )
{
  if ( !output_phase )
  {
    /* connect the node to the output port */
    port = ntk.index_to_node( target_node );
    /* if (1) the target node is an inverter, and    */
    /* (2) its fanin node has more than one fanout,  */
    /* the fanin node of the target node shall be    */
    /* connected to the output port of inverted port */
    if ( ntk.fanin_size( port ) == 1 )
    {
      /* TO CHECK: a valid standard to dectect INV?  */
      ntk.foreach_fanin( port, [&inv_port, &ntk]( auto const& ni ) {
        if ( ntk.fanout_size( ni ) > 1 )
        {
          inv_port = ni;
        }
      } );
    }
  }
  else
  {
    inv_port = ntk.index_to_node( target_node );
    if ( ntk.fanin_size( inv_port ) == 1 )
    {
      ntk.foreach_fanin( inv_port, [&port, &ntk]( auto const& ni ) {
        if ( ntk.fanout_size( ni ) > 1 )
        {
          port = ni;
        }
      } );
    }
  }
}

uint32_t get_node_cost( const uint32_t gate_type )
{
  switch( gate_type )
  {
  case 3:
    return COSTS_MAP[fNOT];
  case 4:
    return COSTS_MAP[fAND];
  case 6:
    return COSTS_MAP[fOR];
  case 12:
    return COSTS_MAP[fXOR];
  default:
    std::cerr << "[e] detected unsupported gate type in the decomposed 2-LUT\n";
    return 0u;
  }
}

uint32_t deref_node( klut& ntk, const std::array<uint32_t, 3> leaves, const uint32_t n )
{
  /* return the number of nodes in the original cut that can be removed, */
  /* if the function of current cut is implemented using an T1 cell,     */
  /* i.e., the MFFC size                                                 */
  if ( auto it_find = std::find( leaves.begin(), leaves.end(), n ); it_find != leaves.end() )
  {
    return 0u;
  }

  /* TODO: is there a better way to access to the node function?         */
  uint32_t gain = get_node_cost( ntk.index_to_node( n ).data[1].h1 );
  ntk.foreach_fanin( n, [&]( auto const& ni ) {
    if ( ntk.decr_fanout_size( ni ) == 0 )
    {
      gain += deref_node( ntk, leaves, ntk.node_to_index( ni ) );
    }
    else
    {
      /* even if a node is not in the MFFC, reduce its fanout size by 1  */
      /* can save a splitter                                             */
      gain += COSTS_MAP[fSPL];
    }
  } );

  return gain;
}

void reref_node( klut& ntk, const std::array<uint32_t, 3> leaves, const uint32_t n )
{
  /* the inverse process of 'deref_root_node', invoked if the T1 cell    */
  /* based implementation turns out to be more expensive                 */
  if ( auto it_find = std::find( leaves.begin(), leaves.end(), n ); it_find != leaves.end() )
  {
    return;
  }

  ntk.foreach_fanin( n, [&]( auto const& ni ) {
    reref_node( ntk, leaves, ntk.node_to_index( ni ) );
  } );
}

bool t1_usage_sanity_check( klut& ntk, std::tuple<std::array<uint32_t, 3>, T1_OUTPUTS>& t1_candidate )
{
  int32_t gain{ 0 };
  const std::array<uint32_t, 3> leaves = t1_candidate.first;
  T1_OUTPUTS& t1_outputs = t1_candidate.second;
  std::array<uint32_t, 3> roots;
  roots[0] = ntk.node_to_index( t1_outputs.sum_to );
  roots[1] = std::max( ntk.node_to_index( t1_outputs.carry_to ), ntk.node_to_index( t1_outputs.inv_carry_to ) );
  roots[2] = std::max( ntk.node_to_index(  t1_outputs.cbar_to ), ntk.node_to_index(  t1_outputs.inv_cbar_to ) );

  for ( const uint32_t root : roots )
  {
    if ( root != ntk.node_to_index( ntk.get_constant( false ) ) )
    {
      gain += deref_root_node( ntk, leaves, root );
    }
  }
  gain += __builtin_popcount( t1_outputs.in_phase ) * ( COSTS_MAP[fNOT] - COSTS_MAP[fDFF] );

  if ( gain -= COSTS_MAP[fFA]; gain < 0 )
  {
    for ( const uint32_t root : roots )
    {
      reref_root_node( ntk, leaves, root );
    }

    return false;
  }

  /* update gate type for the committed T1 cells */
  if ( t1_outputs.has_sum ) { t1_outputs.sum_to.type = FA_GATE; }
  if ( t1_outputs.has_carry ) { t1_outputs.carry_to.type = FA_GATE; }
  if ( t1_outputs.has_carry_inverted ) { t1_outputs.inv_carry_to.type = FA_GATE; }
  if ( t1_outputs.has_cbar ) { t1_outputs.cbar_to.type = FA_GATE; }
  if ( t1_outputs.has_cbar_inverted ) { t1_outputs.inv_cbar_to.type = FA_GATE; }

  return true;
}

void write_klut_specs_supporting_t1( klut const& ntk, phmap::flat_hash_map<std::array<uint32_t,3>, T1_OUTPUTS, ArrayHash<3>> const& t1s, std::string const& filename )
{
  std::ofstream spec_file( filename );

  ntk.foreach_gate( [&]( auto const& n ) {
    if ( ntk.fanout_size( n ) == 0 && !ntk.is_po( n ) )
    {
      /* a dangling node due to the usage of T1 cells */
      return true;
    }

    if ( n.type == FA_GATE )
    {
      /* T1 cells would be handled together later     */
      return true;
    }

    std::vector<uint32_t> n_fanins;
    ntk.foreach_fanin( n, [&ntk, &n_fanins]( auto const& ni ) {
      /* notice that 'signal' and 'node' are equal    */
      /* in kluts                                     */
      n_fanins.push_back( ntk.node_to_index( ni ) );
    } );

    spec_file << fmt::format( "{0},{1},{2}\n", ntk.node_to_index( n ), n.type, fmt::join( n_fanins, "|" ) );
  } );

  /* write information of the commited T1 cells into  */
  /* the csv file                                     */
  for ( auto const& t1 : t1s )
  {
    /* the 5 output ports are in the order of: sum,   */
    /* carry, inverted carry, cbar, and inverted cbar */
    std::vector<uint32_t> output_ports( 5, 0u );
    output_ports[0] = ntk.node_to_index(       t1.second.sum_to );
    output_ports[1] = ntk.node_to_index(     t1.second.carry_to );
    output_ports[2] = ntk.node_to_index( t1.second.inv_carry_to );
    output_ports[3] = ntk.node_to_index(      t1.second.cbar_to );
    output_ports[4] = ntk.node_to_index(  t1.second.inv_cbar_to );

    /* combine input phases with input ports          */
    std::array<std::string, 3> input_ports;
    input_ports[0] = ( ( t1.second.in_phase >> 0 & 1 ) ? "~" : "" ) + std::to_string( t1.first[0] );
    input_ports[1] = ( ( t1.second.in_phase >> 1 & 1 ) ? "~" : "" ) + std::to_string( t1.first[1] );
    input_ports[2] = ( ( t1.second.in_phase >> 2 & 1 ) ? "~" : "" ) + std::to_string( t1.first[2] );

    spec_file << fmt::format( "{0},{1},{2}\n", fmt::join( output_ports, "|" ), 4u, fmt::join( input_ports, "|" ) );
  }
}

int main(int argc, char* argv[])  //
{
  using namespace experiments;
  using namespace mockturtle;

  TT3 XOR3, MAJ3, OR3;
  XOR3._bits = 0x96;
  MAJ3._bits = 0xe8;
  OR3._bits  = 0xfe;

  std::vector<TT3> xor3_tts;
  std::vector<TT3> maj3_tts;
  std::vector<TT3>  or3_tts;
  auto add_xor3 = [&xor3_tts](const TT3 & tt){xor3_tts.push_back(tt);};
  auto add_maj3 = [&maj3_tts](const TT3 & tt){maj3_tts.push_back(tt);};
  auto add_or3  =  [&or3_tts](const TT3 & tt){ or3_tts.push_back(tt);};

  kitty::exact_npn_canonization(XOR3, add_xor3);
  kitty::exact_npn_canonization(MAJ3, add_maj3);
  kitty::exact_npn_canonization( OR3, add_or3 );

  // fmt::print("Compatible TTs:\n");
  // for (auto i = 0u; i < xor3_tts.size(); ++i)
  // {
  //   fmt::print("\t[{}]:\n", i);
  //   fmt::print("\t\tXOR3: {0:08b}={0:02x}={0:d}\n", xor3_tts[i]._bits);
  //   fmt::print("\t\tMAJ3: {0:08b}={0:02x}={0:d}\n", maj3_tts[i]._bits);
  //   fmt::print("\t\t OR3: {0:08b}={0:02x}={0:d}\n",  or3_tts[i]._bits);
  // }
  // return 0;

  experiment<std::string, double, double, double, double> exp( "mapper", "benchmark", "N_PHASES", "#DFF", "area", "delay");

  // uint8_t MIN_N_PHASES = std::stoi(argv[1]);
  // uint8_t MAX_N_PHASES = std::stoi(argv[2]);

  std::vector<uint8_t> PHASES;
  for (auto i = 1; i < argc; ++i)
  {
    PHASES.push_back(std::stoi(argv[i]));
  }
  fmt::print("Phases to analyze: [{}]\n", fmt::join(PHASES, ", "));

  fmt::print( "[i] processing technology library\n" );

  // library to map to technology
  std::vector<gate> gates;
  std::ifstream inputFile( DATABASE_PATH );
  if ( lorina::read_genlib( inputFile, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  // std::unordered_map<std::string, int> nDFF_global = readCSV( NDFF_PATH );
  std::unordered_map<std::string, int> nDFF_global;

  mockturtle::tech_library_params tps; // tps.verbose = true;
  tech_library<NUM_VARS, mockturtle::classification_type::p_configurations> tech_lib( gates, tps );

  #pragma region benchmark_parsing
    // *** BENCHMARKS OF INTEREST ***
    auto benchmarks1 = epfl_benchmarks( experiments::int2float | experiments::priority | experiments::voter);
    auto benchmarks2 = iscas_benchmarks( experiments::c432 | experiments::c880 | experiments::c1908 | experiments::c1355 | experiments::c3540 );
    benchmarks1.insert(benchmarks1.end(), benchmarks2.begin(), benchmarks2.end());

    // *** OPENCORES BENCHMARKS (DO NOT LOOK GOOD) ***
    const std::vector<std::string> BEEREL_BENCHMARKS 
    {
      "simple_spi-gates",
      "des_area-gates",
      "pci_bridge32-gates",
      "spi-gates",
      "mem_ctrl-gates"
    };

    // *** ISCAS89 SEQUENTIAL BENCHMARKS (DO NOT LOOK GOOD) ***
    const std::vector<std::string> ISCAS89_BENCHMARKS {"s382.aig", "s5378.aig", "s13207.aig"};
    benchmarks1.insert(benchmarks1.end(), ISCAS89_BENCHMARKS.begin(), ISCAS89_BENCHMARKS.end());
    // std::reverse(benchmarks1.begin(), benchmarks1.end());

    // *** LIST ALL CONSIDERED BENCHMARKS ***
    fmt::print("Benchmarks:\n\t{}\n", fmt::join(benchmarks1, "\n\t"));
    
    // *** READ COMPOUND GATE LIBRARY ***
    const std::vector<std::vector<UI>> sets_of_levels { { {0,0,0,0}, {0,0,0,1}, {0,0,0,2}, {0,0,1,1}, {0,0,1,2}, {0,1,1,1}, {0,1,1,2}, {0,1,2,2}, {0,1,2,3} } }; //  {0,1,1,3},

    std::unordered_map<ULL, Node> GNM_global = read_global_gnm( sets_of_levels, NODEMAP_PREFIX );
    std::unordered_map<std::string, LibEntry> entries = read_LibEntry_map(LibEntry_file);

  #pragma endregion benchmark_parsing

  // *** START PROCESSING BECNHMARKS ***
  for ( auto const& benchmark : benchmarks1 )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    // *** LOAD NETWORK INTO MIG ***
    mig ntk_original;
    if (benchmark.find("-gates") != std::string::npos) 
    {
      fmt::print("USING THE BLIF READER\n");

      std::string abc_command = fmt::format("abc -c \"read_blif {}{}.blif\" -c strash -c \"write_aiger temp.aig\" ", OPENCORES_FOLDER, benchmark);

      std::system(abc_command.c_str());

      klut temp_klut;
      if ( lorina::read_aiger( "temp.aig", aiger_reader( ntk_original ) ) != lorina::return_code::success )
      {
        fmt::print("Failed to read {}\n", benchmark);
        continue;
      }
    }
    else if ( benchmark.find(".aig") != std::string::npos ) // ISCAS89 benchmark
    {
      fmt::print("USING THE BENCH READER\n");
      std::string path = fmt::format("{}{}", ISCAS89_FOLDER, benchmark);
      if ( lorina::read_aiger( path, aiger_reader( ntk_original ) ) != lorina::return_code::success )
      {
        fmt::print("Failed to read {}\n", benchmark);
        continue;
      }
      // convert_klut_to_graph<mig>(ntk_original, temp_klut);
    }
    else // regular benchmark
    {
      fmt::print( "USING THE AIGER READER\n" );
      if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( ntk_original ) ) != lorina::return_code::success )
      {
        fmt::print("Failed to read {}\n", benchmark);
        continue;
      }
      // convert_klut_to_graph<mig>(ntk_original, temp_klut);
    }

    // *** MAP, NO NEED FOR RETIMING/PATH BALANCING ***
    fmt::print("Started mapping {}\n", benchmark);
    auto [res_wo_pb, st_wo_pb] = map_wo_pb(ntk_original, tech_lib, false); //benchmark, true, nDFF_global, total_ndff_w_pb, total_area_w_pb, cec_w_pb 
    fmt::print("Finished mapping {}\n", benchmark);


    // *** DECOMPOSE COMPOUND GATES INTO PRIMITIVES, REMOVE DFFS, REPLACE OR GATES WITH CB WHERE POSSIBLE ***
    auto _result = decompose_to_klut(res_wo_pb, GNM_global, entries);
    auto klut_decomposed = std::get<0>(_result);
    auto klut_prim_params = std::get<1>(_result);
    auto raw_area = std::get<2>(_result);
    fmt::print("Decomposition complete\n");

    cut_enumeration_params ce_params; 
    ce_params.cut_size = 3u;
  
    auto cuts = mockturtle::cut_enumeration<klut, true>( klut_decomposed, ce_params );
    klut_decomposed.foreach_node( [&]( auto node ) {
      auto idx = klut_decomposed.node_to_index( node );
      auto & node_cuts = cuts.cuts( idx );
      std::cout << node_cuts << "\n";
    } );

    phmap::flat_hash_map<std::array<uint32_t,3>, std::tuple<kitty::dynamic_truth_table, uint32_t>, ArrayHash<3>> maj3_cuts;
    phmap::flat_hash_map<std::array<uint32_t,3>, std::tuple<kitty::dynamic_truth_table, uint32_t>, ArrayHash<3>> xor3_cuts;
    phmap::flat_hash_map<std::array<uint32_t,3>, std::tuple<kitty::dynamic_truth_table, uint32_t>, ArrayHash<3>>  or3_cuts;
    // std::unordered_set<std::array<uint32_t,4>> maj3_cuts;
    // std::unordered_set<std::array<uint32_t,4>> xor3_cuts;
    // phmap::flat_hash_set<cut_type<true, empty_cut_data>> maj3_cuts;
    // phmap::flat_hash_set<cut_type<true, empty_cut_data>> xor3_cuts;
    klut_decomposed.foreach_node( [&]( const klut::signal & node ) { // 
      const auto idx = klut_decomposed.node_to_index( node );
      const auto & node_cut_set = cuts.cuts( idx );
      // fmt::print("[{}] Extracted {} cuts for node {} (idx={})\n", benchmark, node_cut_set.size(), node, idx);
      // if (node_cut_set.size() == 0) { return; }

      for (const auto & cut_entry : node_cut_set)
      {
        //auto qq = *cut_entry;
        const auto tt = cuts.truth_table( *cut_entry );
        /* check if the function of the cut is NPN-equivalent to XOR3 */
        if (std::find_if(xor3_tts.begin(), xor3_tts.end(), [&](const TT3 & xor3_tt){return tt._bits.front() == xor3_tt._bits;}) != xor3_tts.end())
        {
          std::array<uint32_t,3> item;
          std::copy(cut_entry->begin(), cut_entry->end(), item.begin());
          xor3_cuts.emplace( item, std::make_tuple(tt, idx) );

          std::ostringstream oss;
          oss << *cut_entry;
          fmt::print("[{}]\t XOR3 cut {} with tt 0b{}=0x{}\n", benchmark, oss.str(), kitty::to_binary( tt ), kitty::to_hex( tt ));
        }
        if (std::find_if(maj3_tts.begin(), maj3_tts.end(), [&](const TT3 & maj3_tt){return tt._bits.front() == maj3_tt._bits;}) != maj3_tts.end())
        {
          std::array<uint32_t,3> item;
          std::copy(cut_entry->begin(), cut_entry->end(), item.begin());
          maj3_cuts.emplace( item, std::make_tuple(tt, idx) );

          std::ostringstream oss;
          oss << *cut_entry;
          fmt::print("[{}]\t MAJ3 cut {} with tt 0b{}=0x{}\n", benchmark, oss.str(), kitty::to_binary( tt ), kitty::to_hex( tt ));
        }
        if ( std::find_if( or3_tts.begin(), or3_tts.end(), [&]( TT3 const& or3_tt ) { return tt._bits.front() == or3_tt._bits; } ) != or3_tts.end() )
        {
          std::array<uint32_t, 3> item;
          std::copy( cut_entry->begin(), cut_entry->end(), item.begin() );
          or3_cuts.emplace( item, std::make_tuple( tt, idx ) );
        }
        // std::ostringstream oss;
        // oss << *cut_entry;
        // fmt::print("[{}] \tCut {} with tt 0b{}=0x{}\n", benchmark, oss.str(), kitty::to_binary( tt ), kitty::to_hex( tt ));
      }
      // flat_hash_set<std::array<uint32_t, 3>> arraySet;
      // auto tt = cuts.truth_table( node_cut_set );
      // std::cout << node_cut_set << " -> " <<  << "\n";
      // cuts.truth_table( node_cut_set );
    } );

    /* TODO: adopt the assumption that, it is a good deal if an implementation can make use of more than 2 out of the 3 outputs of T1,  */
    /* then we would need four rounds of matching: (1) 3 leaves using all 3 outputs; (2) ...using XOR and MAJ; (3) ...using MAJ and OR; */
    /* (4) ...using XOR and OR. For each 3 leaves found, check validity by deciding which of the 8 T1s to choose                        */
    phmap::flat_hash_map<std::array<uint32_t,3>, T1_OUTPUTS, ArrayHash<3>> t1_candidates;
    // std::vector<std::array<uint32_t, 3>> checked_leaves;

    for ( auto const& xor3_cut : xor3_cuts )
    {
      auto target_leaves = xor3_cut.first;
      if ( maj3_cuts.contains( target_leaves ) )
      {
        if ( or3_cuts.contains( target_leaves ) )
        {
          /* find a candidate potentially using all 3 outputs of an T1 cell */
          /* check whether the combination of truth table is valid          */
          bool maj_phase{ false }, or_phase{ false };
          for ( auto const& t1_cell : standard_T1_cells )
          {
            maj_phase = or_phase = false;
            if ( std::get<0>( xor3_cut.second )._bits[0] != t1_cell.sum_truth_table )
            {
              continue;
            }

            if ( std::get<0>( maj3_cuts[target_leaves] )._bits[0] == ~t1_cell.carry_truth_table )
            {
              maj_phase = true;
            }
            else if ( std::get<0>( maj3_cuts[target_leaves] )._bits[0] != t1_cell.carry_truth_table )
            {
              continue;
            }

            if ( std::get<0>( or3_cuts[target_leaves] )._bits[0] == ~t1_cell.cbar_truth_table )
            {
              or_phase = true;
            }
            else if ( std::get<0>( or3_cuts[target_leaves] )._bits[0] == t1_cell.cbar_truth_table )
            {
              continue;
            }

            /* managed to find a T1 cell-based implementation      */
            /* figure out the nodes connecting to the output ports */
            klut::node sum_to = klut_decomposed.index_to_node( std::get<1>( xor3_cut.second ) );
            klut::node carry_to = klut_decomposed.get_constant( false );
            klut::node inv_carry_to = klut_decomposed.get_constant( false );
            klut::node cbar_to = klut_decomposed.get_constant( false );
            klut::node inv_cbar_to = klut_decomposed.get_constant( false );
            config_t1_ports_connection( klut_decomposed, maj_phase, std::get<1>( maj3_cuts[target_leaves] ), carry_to, inv_carry_to );
            config_t1_ports_connection( klut_decomposed,  or_phase, std::get<1>(  or3_cuts[target_leaves] ),  cbar_to,  inv_cbar_to );

            /* instantiate an T1 cell */
            t1_candidates.emplace( target_leaves, T1_OUTPUTS( t1_cell.in_phase, inv_carry_to == klut_decomposed.get_constant( false ), inv_cbar_to == klut_decomposed.get_constant( false ), 
                                                              sum_to == klut_decomposed.get_constant( false ), carry_to == klut_decomposed.get_constant( false ), cbar_to == klut_decomposed.get_constant( false ), 
                                                              sum_to, carry_to, inv_carry_to, cbar_to, inv_cbar_to ) );
            break;
            /* TODO: think about if it is possible that a design can be implemented under */
            /* different input phases of a T1 cell; if yes, it is better to enumerate all */
            /* the potential input phases and choose the one requires the minimum cost    */
          }
        }
      }
    }

    for ( auto const& xor3_cut : xor3_cuts )
    {
      auto target_leaves = xor3_cut.first;
      if ( maj3_cuts.contains( target_leaves ) )
      {
        if ( !t1_candidates.contains( target_leaves ) )
        {
          /* find a candidate potentially using 2 outputs of an T1 cell  */
          /* ( XOR3, MAJ3 ), check whether the combnation of tt is valid */
          bool maj_phase{ false };
          for ( auto const& t1_cell : standard_T1_cells )
          {
            if ( std::get<0>( xor3_cut.second )._bits[0] != t1_cell.sum_truth_table )
            {
              continue;
            }

            if ( std::get<0>( maj3_cuts[target_leaves] )._bits[0] == ~t1_cell.carry_truth_table )
            {
              maj_phase = true;
            }
            else if ( std::get<0>( maj3_cuts[target_leaves] )._bits[0] != t1_cell.carry_truth_table )
            {
              continue;
            }

            /* managed to find a T1 cell-based implementation      */
            /* figure out the nodes connecting to the output ports */
            klut::node sum_to = klut_decomposed.index_to_node( std::get<1>( xor3_cut.second ) );
            klut::node carry_to = klut_decomposed.get_constant( false );
            klut::node inv_carry_to = klut_decomposed.get_constant( false );
            config_t1_ports_connection( klut_decomposed, maj_phase, std::get<1>( maj3_cuts[target_leaves] ), carry_to, inv_carry_to );

            /* instantiate an T1 cell */
            t1_candidates.emplace( target_leaves, T1_OUTPUTS( t1_cell.in_phase, inv_carry_to == klut_decomposed.get_constant( false ), false, 
                                                              carry_to == klut_decomposed.get_constant( false ), false, 
                                                              sum_to, carry_to, inv_carry_to, klut_decomposed.get_constant( false ), klut_decomposed.get_constant( false ) ) );
            break;
          }
        }
      }
    }

    for ( auto const& xor3_cut : xor3_cuts )
    {
      auto target_leaves = xor3_cut.first;
      if ( or3_cuts.contains( target_leaves ) )
      {
        if ( !t1_candidates.contains( target_leaves ) )
        {
          /* find a candidate potentially using 2 outputs of an T1 cell  */
          /* ( XOR3,  OR3 ), check whether the combnation of tt is valid */
          bool or_phase{ false };
          for ( auto const& t1_cell : standard_T1_cells )
          {
            if ( std::get<0>( xor3_cut.second )._bits[0] != t1_cell.sum_truth_table )
            {
              continue;
            }

            if ( std::get<0>( or3_cuts[target_leaves] )._bits[0] == ~t1_cell.cbar_truth_table )
            {
              or_phase = true;
            }
            else if ( std::get<0>( or3_cuts[target_leaves] )._bits[0] != t1_cell.cbar_truth_table )
            {
              continue;
            }

            /* managed to find a T1 cell-based implementation      */
            /* figure out the nodes connecting to the output ports */
            klut::node sum_to = klut_decomposed.index_to_node( std::get<1>( xor3_cut.second ) );
            klut::node cbar_to = klut_decomposed.get_constant( false );
            klut::node inv_cbar_to = klut_decomposed.get_constant( false );
            config_t1_ports_connection( klut_decomposed, or_phase, std::get<1>( or3_cuts[target_leaves] ), cbar_to, inv_cbar_to );

            /* instantiate an T1 cell */
            t1_candidates.emplace( target_leaves, T1_OUTPUTS( t1_cell.in_phase, false, inv_cbar_to == klut_decomposed.get_constant( false ), 
                                                              false, cbar_to == klut_decomposed.get_constant( false ), 
                                                              sum_to, klut_decomposed.get_constant( false ), klut_decomposed.get_constant( false ), cbar_to, inv_cbar_to ) );
            break;
          }
        }
      }
    }

    for ( auto const& maj3_cut : maj3_cuts )
    {
      auto target_leaves = maj3_cut.first;
      if ( or3_cuts.contains( maj3_cut.first ) )
      {
        if ( !t1_candidates.contains( target_leaves ) )
        {
          /* find a candidate potentially using 2 outputs of an T1 cell  */
          /* ( MAJ3,  OR3 ), check whether the combnation of tt is valid */
          bool maj_phase{ false }, or_phase{ false };
          for ( auto const& t1_cell : standard_T1_cells )
          {
            if ( std::get<0>( maj3_cut.second )._bits[0] == ~t1_cell.carry_truth_table )
            {
              maj_phase = true;
            }
            else if ( std::get<0>( maj3_cut.second )._bits[0] != t1_cell.carry_truth_table )
            {
              continue;
            }

            if ( std::get<0>( or3_cuts[target_leaves] )._bits[0] == ~t1_cell.cbar_truth_table )
            {
              or_phase = true;
            }
            else if ( std::get<0>( or3_cuts[target_leaves] )._bits[0] != t1_cell.cbar_truth_table )
            {
              continue;
            }

            /* managed to find a T1 cell-based implementation      */
            /* figure out the nodes connecting to the output ports */
            klut::node carry_to = klut_decomposed.get_constant( false );
            klut::node inv_carry_to = klut_decomposed.get_constant( false );
            klut::node cbar_to = klut_decomposed.get_constant( false );
            klut::node inv_cbar_to = klut_decomposed.get_constant( false );
            config_t1_ports_connection( klut_decomposed, maj_phase, std::get<1>(          maj3_cut.second ), carry_to, inv_carry_to );
            config_t1_ports_connection( klut_decomposed,  or_phase, std::get<1>(  or3_cuts[target_leaves] ),  cbar_to,  inv_cbar_to );

            /* instantiate an T1 cell */
            t1_candidates.emplace( target_leaves, T1_OUTPUTS( t1_cell.in_phase, inv_carry_to == klut_decomposed.get_constant( false ), inv_cbar_to == klut_decomposed.get_constant( false ), 
                                                              carry_to == klut_decomposed.get_constant( false ), cbar_to == klut_decomposed.get_constant( false ), 
                                                              klut_decomposed.get_constant( false ), carry_to, inv_carry_to, cbar_to, inv_cbar_to ) );
            break;
          }
        }
      }
    }

    fmt::print( "Potentially usage of the T1 cells:\n" );
    for ( auto const& t1_candidate : t1_candidates )
    {
      t1_candidate.second.report();
    }

    /* estimate the gain of implementing parts of the circuits using T1s instead */
    for ( auto it_t1_cands{ t1_candidates.begin() }; it_t1_cands < t1_candidates.end(); ++i1_t1_cands )
    {
      if ( !t1_usage_sanity_check( klut_decomposed, *it_t1_cands ) )
      {
        t1_candidates.erase( it_t1_cands );
      }
    }
    
    /* bind the instantiated T1 cells to the 2-LUT network for csv generation    */
    write_klut_specs_supporting_t1( klut_decomposed, t1_candidates, std::string const& filename )

    continue;
  
    // *** [temporary] GREEDILY ASSIGN A STAGE (sigma) TO EACH ELEMENT ***
    std::unordered_map<klut::signal, glob_phase_t> glob_phase = greedy_assign(klut_decomposed, klut_prim_params, false);
    printUnorderedMap(glob_phase);

    // for (auto n_phases = MIN_N_PHASES; n_phases <= MAX_N_PHASES; ++n_phases)
    for (const auto n_phases : PHASES)
    {
      fmt::print("[i] Mapping with {} phases\n", n_phases);
      // *** IF i = 0, we assign phases with the CP-SAT
      // *** IF i = 1, we assign phases greedily
      // for (auto i = 0; i < 2; ++i)
      for (auto i = 0; i < 1; ++i)
      {
        klut network { klut_decomposed.clone() };
        std::unordered_map<unsigned int, unsigned int> assignment;

        if (i == 0)
        {
          const std::string ilp_cfg_filename = fmt::format("ilp_configs/{}.csv", benchmark);
          fmt::print("\tWriting config {}\n", ilp_cfg_filename);
          write_klut_specs(network, klut_prim_params, glob_phase, ilp_cfg_filename);

          fmt::print("\tCalling OR-Tools\n");
          auto [obj_val, assignment_local, status] = cpsat_macro_opt(ilp_cfg_filename, n_phases);

          if (status == "SUCCESS") // (true) // 
          {
            assignment.insert(std::make_move_iterator(assignment_local.begin()), std::make_move_iterator(assignment_local.end()));
            fmt::print("[i] CP-SAT PHASE ASSIGNMENT: SUCCESS\n");
          }
          else
          {
            fmt::print("[i] CP-SAT PHASE ASSIGNMENT: FAIL\n");
            continue;
          }
        }
        
        // *** IF i = 0, "assignment" has stages assigned by the CP-SAT
        // *** IF i = 1, "assignment" is empty
        greedy_ntk_assign(network, klut_prim_params, n_phases, assignment, true);

        // *** Greedily insert splitters
        splitter_ntk_insertion( network, true);

        // network.foreach_node([&] ( const klut::signal & node ) 
        // {
        //   if ( network.fanout_size( node ) > 1 )
        //   {
        //     assert( network.node_function( node ) == 0x2 );
        //   };
        // });

        fmt::print("[i] FINISHED PHASE ASSIGNMENT\n");

        fmt::print("[i] EXTRACTING PATHS\n");
        std::vector<Path> paths = extract_paths( network, true );
        // auto [DFF_REG, precalc_ndff] = dff_vars(NR, paths, N_PHASES);

        auto total_num_dff = 0u;
        auto file_ctr = 0u;
        auto path_ctr = 0u;
        for (const Path & path : paths)
        {
          fmt::print("\tAnalyzing the path {} out of {}\n", ++path_ctr, paths.size());
          // *** Create binary variables
          auto [DFF_REG, precalc_ndff, required_SA_DFFs] = dff_vars_single_paths(path, network, n_phases);
          total_num_dff += precalc_ndff;
          fmt::print("\t\t\t\t[i]: Precalculated {} DFFs, total #DFF = {}\n", precalc_ndff, total_num_dff);
          
          // *** Generate constraints
          std::vector<Snake> snakes = sectional_snake(path, network, DFF_REG, n_phases, true);

          fmt::print("\tCreated {} snakes\n", snakes.size());
          // *** If there's anything that needs optimization
          if (!snakes.empty())
          {
            std::string cfg_file = fmt::format("ilp_configs/{}_cfgNR_{}.csv", benchmark, file_ctr++);
            write_snakes(snakes, DFF_REG, required_SA_DFFs, cfg_file, n_phases, true);
            auto num_dff = cpsat_ortools(cfg_file);
            // fmt::print("OR Tools optimized to {} DFF\n", num_dff);
            total_num_dff += num_dff;
            fmt::print("\t\t\t\t[i] total CPSAT #DFF = {}\n", total_num_dff);
          }
        }

        // *** Record maximum phase
        uint64_t max_phase = 0u;
        // *** Record number of splitters and total number of DFFs (not only path balancing DFFs)
        uint64_t total_num_spl = 0;
        network.foreach_gate([&](const klut::signal & node)
        {
          NodeData node_data { network.value(node) };
          fmt::print("[Node {}] old max_phase = {}\tnode_data = {}\t", node, max_phase, static_cast<int>(node_data.sigma));
          max_phase = generic_max(max_phase, node_data.sigma);
          fmt::print("new max_phase = {}\n", max_phase);

          auto fo_size = network.fanout_size(node);
          if (fo_size > 1)
          {
            total_num_spl += fo_size - 1;
          }
        });

        network.foreach_po([&](const klut::signal & node)
        {
          NodeData node_data { network.value(node) };
          fmt::print("[PO {}] max_phase = {}, sigma = {}, node #DFF = {}\n", node, max_phase, static_cast<int>(node_data.sigma), ( (max_phase - node_data.sigma) / n_phases ) );
          total_num_dff += (max_phase - node_data.sigma) / n_phases;
          fmt::print("\t\t\t\t[i] total #DFF = {}\n", total_num_dff);
        });

        fmt::print("{} PHASES: #DFF   for {} is {}\n", n_phases, benchmark, total_num_dff);
        int total_area = raw_area + total_num_dff * COSTS_MAP[fDFF] + total_num_spl * COSTS_MAP[fSPL];
        fmt::print("{} PHASES: #AREA  for {} is {}\n", n_phases, benchmark, total_area);
        fmt::print("{} PHASES: #MAX GLOB PHASE for {} is {}\n", n_phases, benchmark, max_phase);

        exp(fmt::format("{}_{}", benchmark, (i==0)?"CPSAT":"GREEDY"), n_phases, total_num_dff, total_area, ( (max_phase - 1) / n_phases + 1 ));
        exp.save();
        exp.table();
      }
    }
  }
  // TODO : now, count #DFFs in extracted paths. Perhaps, the function can be written within the Path object.
  return 0;
}

// Adder 18542 vs 24384


/* 
  7 phases
  |-------------|--------------|--------------|--------|--------------|--------------|------------|---------|
  |   benchmark | #DFF (tight) | #DFF (loose) | Ratio  | area (tight) | area (loose) | Area Ratio |  delay  |
  |-------------|--------------|--------------|--------|--------------|--------------|------------|---------|
  |   int2float |      217     |      217     |  1.00  |     5'136    |     5'136    |     1.00   |    2    |
  |    priority |    3'375     |    3'284     |  1.03  |    45'724    |    45'087    |     1.01   |   12    |
  |       voter |    2'380     |    2'180     |  1.09  |   164'204    |   162'804    |     1.01   |    6    |
  |        c432 |      352     |      342     |  1.03  |     5'186    |     5'116    |     1.01   |    4    |
  |        c880 |      270     |      254     |  1.06  |     6'302    |     6'190    |     1.02   |    3    |
  |       c1355 |       46     |       46     |  1.00  |     4'515    |     4'515    |     1.00   |    2    |
  |       c1908 |      125     |      125     |  1.00  |     3'529    |     3'529    |     1.00   |    2    |
  |       c3540 |      579     |      589     |  0.98  |    16'946    |    17'016    |     1.00   |    3    |
  |    s382.aig |      156     |       89     |  1.75  |     2'880    |     2'411    |     1.19   |    2    |
  |   s5378.aig |      857     |      808     |  1.06  |    22'104    |    21'761    |     1.02   |    2    |
  |  s13207.aig |    2'013     |    1'837     |  1.10  |    43'614    |    42'382    |     1.03   |    3    |
  |-------------|--------------|--------------|--------|--------------|--------------|------------|---------|
*/
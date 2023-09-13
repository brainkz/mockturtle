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

#include <mockturtle/algorithms/multiphase.hpp>

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
template <size_t N, typename T>
using array_map = phmap::flat_hash_map<std::array<uint32_t, N>, T, ArrayHash<N>>;

// Sunmagnetics Technology Library
constexpr std::array<int, 13> COSTS_MAP = {7, 9, 8, 8, 12, 8, 999, 999, 999, 8, 3, 0, 25};

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

  /// @brief Checks whether the truth table matches the carry TT of this T1 cell
  /// @param tt truth table to match
  /// @param phase - set to *true* if the matched to negated output, or *false* if matched to non-negated output
  /// @return *true* if there is a match at all
  bool check_carry( const kitty::dynamic_truth_table & tt, bool & phase ) const
  {
    if ( tt._bits[0] == ~carry_truth_table )
    {
      phase = true;
      return true;
    }
    else if ( tt._bits[0] != carry_truth_table )
    {
      return false;
    }
    phase = false;
    return true;
  }
  /// @brief Checks whether the truth table matches the cbar TT of this T1 cell
  /// @param tt truth table to match
  /// @param phase - set to *true* if the matched to negated output, or *false* if matched to non-negated output
  /// @return *true* if there is a match at all
  bool check_cbar( const kitty::dynamic_truth_table & tt, bool & phase ) const
  {
    if ( tt._bits[0] == ~cbar_truth_table )
    {
      phase = true;
      return true;
    }
    else if ( tt._bits[0] != cbar_truth_table )
    {
      return false;
    }
    phase = false;
    return true;
  }
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
  // uint8_t   sum_truth_table{};
  // uint8_t carry_truth_table{};
  // uint8_t  cbar_truth_table{};
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
    // fmt::print("\t\t      Sum ( Node {0} ): {1:08b} ( 0x{1:02x} )\n",       sum_to,    sum_truth_table );
    // fmt::print("\t\t    Carry ( Node {0} ): {1:08b} ( 0x{1:02x} )\n",     carry_to,  carry_truth_table );
    // fmt::print("\t\tInv carry ( Node {0} ): {1:08b} ( 0x{1:02x} )\n", inv_carry_to, ~carry_truth_table );
    // fmt::print("\t\t     Cbar ( Node {0} ): {1:08b} ( 0x{1:02x} )\n",      cbar_to,   cbar_truth_table );
    // fmt::print("\t\t Inv cbar ( Node {0} ): {1:08b} ( 0x{1:02x} )\n",  inv_cbar_to,  ~cbar_truth_table );
    fmt::print("\t\t      Sum : {}\n", (            has_sum  ? "Node " + std::to_string(       sum_to ) : "N/A" ) );
    fmt::print("\t\t    Carry : {}\n", (          has_carry  ? "Node " + std::to_string(     carry_to ) : "N/A" ) );
    fmt::print("\t\tInv carry : {}\n", ( has_carry_inverted  ? "Node " + std::to_string( inv_carry_to ) : "N/A" ) );
    fmt::print("\t\t     Cbar : {}\n", (           has_cbar  ? "Node " + std::to_string(      cbar_to ) : "N/A" ) );
    fmt::print("\t\t Inv cbar : {}\n", (  has_cbar_inverted  ? "Node " + std::to_string(  inv_cbar_to ) : "N/A" ) );

  }
};

template <typename CutType>
array_map<3, std::tuple<kitty::dynamic_truth_table, uint32_t>> match_cuts( const std::vector<TT3> & template_tts, const klut & ntk, const CutType & cuts)
{
  array_map<3, std::tuple<kitty::dynamic_truth_table, uint32_t>> matching_cuts;

  ntk.foreach_node( [&]( const klut::signal & node ) 
  { 
    const auto idx = ntk.node_to_index( node );
    const auto & node_cut_set = cuts.cuts( idx );
    // fmt::print("[{}] Extracted {} cuts for node {} (idx={})\n", benchmark, node_cut_set.size(), node, idx);
    // if (node_cut_set.size() == 0) { return; }

    /* check if the function of the cut is NPN-equivalent to any of the template_tts */
    for (const auto & cut_entry : node_cut_set)
    {
      const auto tt = cuts.truth_table( *cut_entry );
      if (std::find_if(template_tts.begin(), template_tts.end(), [&](const TT3 & tt){return tt._bits.front() == tt._bits;}) != template_tts.end())
      {
        std::array<uint32_t,3> leaves;
        std::copy(cut_entry->begin(), cut_entry->end(), leaves.begin());
        matching_cuts.emplace( leaves, std::make_tuple(tt, idx) );
      }
    }
  } );
  return matching_cuts;
}

array_map<3, T1_OUTPUTS> find_t1_candidates( 
  const klut & ntk,
  const array_map<3, std::tuple<kitty::dynamic_truth_table, uint32_t>>& xor3_cuts,
  const array_map<3, std::tuple<kitty::dynamic_truth_table, uint32_t>>& maj3_cuts, 
  const array_map<3, std::tuple<kitty::dynamic_truth_table, uint32_t>>& or3_cuts
) 
{
  array_map<3, T1_OUTPUTS> t1_candidates;

  for ( auto const& [target_leaves, xor3_tt_and_node] : xor3_cuts )
  {
    if ( !( maj3_cuts.contains( target_leaves ) &&  // if there's a match with MAJ3
             or3_cuts.contains( target_leaves ) ) ) // and a match with OR3
    {
      continue;
    }
    const auto & [xor3_tt, xor_index] = xor3_tt_and_node;
    const auto & [maj3_tt, maj_index] = maj3_cuts.at(target_leaves);
    const auto & [ or3_tt,  or_index] = or3_cuts.at(target_leaves); 
    /* find a candidate potentially using all 3 outputs of an T1 cell */
    /* check whether the combination of truth table is valid          */
    for ( auto const& t1_cell : standard_T1_cells )
    {
      bool maj_phase{ false }, or_phase{ false };
      if (
        xor3_tt._bits[0] != t1_cell.sum_truth_table ||
        !t1_cell.check_carry(maj3_tt, maj_phase) || 
        !t1_cell.check_cbar(or3_tt, or_phase)
      )
      {
        continue;
      }

      /* managed to find a T1 cell-based implementation      */
      /* figure out the nodes connecting to the output ports */
      klut::node sum_to       = ntk.index_to_node( xor_index );

      klut::node carry_to     = ntk.get_constant( false );
      klut::node inv_carry_to = ntk.get_constant( false );
      config_t1_ports_connection( ntk, maj_phase, maj_index, carry_to, inv_carry_to );

      klut::node cbar_to      = ntk.get_constant( false );
      klut::node inv_cbar_to  = ntk.get_constant( false );
      config_t1_ports_connection( ntk,  or_phase, or_index, cbar_to,  inv_cbar_to );

      /* instantiate an T1 cell */
      t1_candidates.emplace( target_leaves, T1_OUTPUTS( 
        t1_cell.in_phase, 
        inv_carry_to != ntk.get_constant( false ), 
        inv_cbar_to  != ntk.get_constant( false ), 
        sum_to       != ntk.get_constant( false ), 
        carry_to     != ntk.get_constant( false ), 
        cbar_to      == ntk.get_constant( false ), 
        sum_to, carry_to, inv_carry_to, cbar_to, inv_cbar_to ) );
      break;
      /* TODO: think about if it is possible that a design can be implemented under */
      /* different input phases of a T1 cell; if yes, it is better to enumerate all */
      /* the potential input phases and choose the one requires the minimum cost    */
    }
  }

  for ( auto const& [target_leaves, xor3_tt_and_node] : xor3_cuts )
  {
    if ( !maj3_cuts.contains( target_leaves ) ||   // if there's no XOR-MAJ match
          t1_candidates.contains( target_leaves ) ) // or a 3-way match already found
    {
      continue;
    }

    auto [xor3_tt, xor_index] = xor3_tt_and_node;
    auto [maj3_tt, maj_index] = maj3_cuts.at(target_leaves);
    /* find a candidate potentially using 2 outputs of an T1 cell  */
    /* ( XOR3, MAJ3 ), check whether the combnation of tt is valid */
    bool maj_phase{ false };
    for ( auto const& t1_cell : standard_T1_cells )
    {
      if ( xor3_tt._bits[0] != t1_cell.sum_truth_table )
      {
        continue;
      }

      bool maj_phase{ false };
      if (!t1_cell.check_carry(maj3_tt, maj_phase))
      {
        continue;
      }

      /* managed to find a T1 cell-based implementation      */
      /* figure out the nodes connecting to the output ports */
      klut::node sum_to       = ntk.index_to_node( xor_index );
      klut::node carry_to     = ntk.get_constant( false );
      klut::node inv_carry_to = ntk.get_constant( false );
      config_t1_ports_connection( ntk, maj_phase, maj_index, carry_to, inv_carry_to );

      /* instantiate an T1 cell */
      t1_candidates.emplace( target_leaves, T1_OUTPUTS( 
        t1_cell.in_phase, 
        inv_carry_to != ntk.get_constant( false ), 
        false, 
        true, 
        carry_to != ntk.get_constant( false ), 
        false, 
        sum_to, 
        carry_to, 
        inv_carry_to, 
        ntk.get_constant( false ), 
        ntk.get_constant( false ) ) );
      break;
    }
  }

  for ( auto const& [target_leaves, xor3_tt_and_node] : xor3_cuts )
  {
    if ( !or3_cuts.contains( target_leaves ) ||  // if there's no match with OR3
        t1_candidates.contains( target_leaves )) // or a 3-way match already found
    {
      continue;
    }
    auto [xor3_tt, xor_index] = xor3_tt_and_node;
    auto [or3_tt, or_index] = or3_cuts.at(target_leaves); 
    /* find a candidate potentially using 2 outputs of an T1 cell  */
    /* ( XOR3,  OR3 ), check whether the combnation of tt is valid */
    bool or_phase{ false };
    for ( auto const& t1_cell : standard_T1_cells )
    {
      if ( xor3_tt._bits[0] != t1_cell.sum_truth_table )
      {
        continue;
      }

      bool or_phase{ false };
      if (!t1_cell.check_cbar(or3_tt, or_phase))
      {
        continue;
      }

      /* managed to find a T1 cell-based implementation      */
      /* figure out the nodes connecting to the output ports */
      klut::node sum_to      = ntk.index_to_node( xor_index );
      klut::node cbar_to     = ntk.get_constant( false );
      klut::node inv_cbar_to = ntk.get_constant( false );
      config_t1_ports_connection( ntk, or_phase, or_index, cbar_to, inv_cbar_to );

      /* instantiate an T1 cell */
      t1_candidates.emplace( target_leaves, T1_OUTPUTS( 
        t1_cell.in_phase, 
        false, 
        inv_cbar_to != ntk.get_constant( false ), 
        true, 
        false, 
        cbar_to != ntk.get_constant( false ), 
        sum_to, 
        ntk.get_constant( false ), 
        ntk.get_constant( false ), 
        cbar_to, 
        inv_cbar_to ) );
      break;
    }
  }

  for ( auto const& [target_leaves, maj3_tt_and_node] : maj3_cuts )
  {
    if ( !or3_cuts.contains( target_leaves ) ||  // if there's no match with OR3
        t1_candidates.contains( target_leaves )) // or a 3-way match already found
    {
      continue;
    }
    auto [maj3_tt, maj_index] = maj3_tt_and_node;
    auto [or3_tt, or_index] = or3_cuts.at(target_leaves); 
    /* find a candidate potentially using 2 outputs of an T1 cell  */
    /* ( MAJ3,  OR3 ), check whether the combnation of tt is valid */
    for ( auto const& t1_cell : standard_T1_cells )
    {
      bool maj_phase{ false }, or_phase{ false };
      if ( maj3_tt._bits[0] == ~t1_cell.carry_truth_table )
      {
        maj_phase = true;
      }
      else if ( maj3_tt._bits[0] != t1_cell.carry_truth_table )
      {
        continue;
      }

      bool or_phase{ false };
      if (!t1_cell.check_cbar(or3_tt, or_phase))
      {
        continue;
      }

      /* managed to find a T1 cell-based implementation      */
      /* figure out the nodes connecting to the output ports */
      klut::node carry_to     = ntk.get_constant( false );
      klut::node inv_carry_to = ntk.get_constant( false );
      klut::node cbar_to      = ntk.get_constant( false );
      klut::node inv_cbar_to  = ntk.get_constant( false );
      config_t1_ports_connection( ntk, maj_phase, maj_index, carry_to, inv_carry_to );
      config_t1_ports_connection( ntk,  or_phase, or_index,  cbar_to,  inv_cbar_to );

      /* instantiate an T1 cell */
      t1_candidates.emplace( target_leaves, T1_OUTPUTS( 
        t1_cell.in_phase, 
        inv_carry_to != ntk.get_constant( false ), 
        inv_cbar_to != ntk.get_constant( false ), 
        false, 
        carry_to != ntk.get_constant( false ), 
        cbar_to != ntk.get_constant( false ), 
        ntk.get_constant( false ), 
        carry_to, 
        inv_carry_to, 
        cbar_to, 
        inv_cbar_to ) );
      break;
    }
  }

  return t1_candidates;
}

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
    if (it != path.targets.end() && ( ( fanout_sigma < 0 ) || ( ( fanout_sigma >= 0 ) && ( dff.sigma >= static_cast<uint32_t>( fanout_sigma ) ) ) ) )
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

void config_t1_ports_connection( klut const& ntk, const bool output_phase, const uint32_t target_node, klut::node& port, klut::node& inv_port )
{
  if ( !output_phase )
  {
    /*      connect the node to the output port      */
    port = ntk.index_to_node( target_node );
    /* if (1) the target node is an inverter, and    */
    /* (2) its fanin node has more than one fanout,  */
    /* the fanin node of the target node shall be    */
    /* connected to the output port of inverted port */
    if ( ntk.func_lit( target_node ) == 3 )
    {
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
    if ( ntk.func_lit( target_node ) == 3 )
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

  uint32_t gain = get_node_cost( ntk.func_lit( n ) );
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

bool t1_usage_sanity_check( klut& ntk, std::pair<const std::array<uint32_t, 3>, T1_OUTPUTS>& t1_candidate, int64_t& updated_area )
{
  int32_t gain{ 0 };
  const std::array<uint32_t, 3> leaves = std::get<0>( t1_candidate );
  T1_OUTPUTS& t1_outputs = std::get<1>( t1_candidate );
  std::array<uint32_t, 3> roots;
  roots[0] = ntk.node_to_index( t1_outputs.sum_to );
  roots[1] = std::max( ntk.node_to_index( t1_outputs.carry_to ), ntk.node_to_index( t1_outputs.inv_carry_to ) );
  roots[2] = std::max( ntk.node_to_index( t1_outputs.cbar_to  ), ntk.node_to_index( t1_outputs.inv_cbar_to  ) );

  for ( const uint32_t root : roots )
  {
    if ( root != ntk.node_to_index( ntk.get_constant( false ) ) )
    {
      gain += static_cast<int32_t>( deref_node( ntk, leaves, root ) );
    }
  }
  gain -= static_cast<int32_t>( __builtin_popcount( t1_outputs.in_phase ) * ( COSTS_MAP[fNOT] - COSTS_MAP[fDFF] ) );

  if ( gain -= static_cast<int32_t>( COSTS_MAP[fFA] ); gain < 0 )
  {
    for ( const uint32_t root : roots )
    {
      reref_node( ntk, leaves, root );
    }

    return false;
  }
  /* update gate type for the committed T1 cells */
  if ( t1_outputs.has_sum ) { ntk.set_value( t1_outputs.sum_to, NodeData( static_cast<NodeData>( ntk.value( t1_outputs.sum_to ) ).sigma, T1_GATE ).value ); }
  if ( t1_outputs.has_carry ) { ntk.set_value( t1_outputs.carry_to, NodeData( static_cast<NodeData>( ntk.value( t1_outputs.carry_to ) ).sigma, T1_GATE ).value ); }
  if ( t1_outputs.has_carry_inverted ) { ntk.set_value( t1_outputs.inv_carry_to, NodeData( static_cast<NodeData>( ntk.value( t1_outputs.inv_carry_to ) ).sigma, T1_GATE ).value ); }
  if ( t1_outputs.has_cbar ) { ntk.set_value( t1_outputs.cbar_to, NodeData( static_cast<NodeData>( ntk.value( t1_outputs.cbar_to ) ).sigma, T1_GATE ).value ); }
  if ( t1_outputs.has_cbar_inverted ) { ntk.set_value( t1_outputs.inv_cbar_to, NodeData( static_cast<NodeData>( ntk.value( t1_outputs.inv_cbar_to ) ).sigma, T1_GATE ).value ); }

  updated_area -= gain;

  return true;
}

void write_klut_specs_supporting_t1( klut const& ntk, array_map<3, T1_OUTPUTS> const& t1s, std::string const& filename )
{
  std::ofstream spec_file( filename );

  ntk.foreach_gate( [&]( auto const& n ) {
    if ( ntk.is_dangling( n ) )
    {
      /* a dangling node due to the usage of T1 cells */
      return true;
    }

    if ( static_cast<NodeData>( ntk.value( n ) ).type == FA_GATE )
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

    spec_file << fmt::format( "{0},{1},{2}\n", ntk.node_to_index( n ), static_cast<NodeData>( ntk.value( n ) ).type, fmt::join( n_fanins, "|" ) );
    return true;
  } );

  /* write information of the committed T1 cells into  */
  /* the csv file                                      */
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

bool customized_T1_classifier( T1_OUTPUTS const& t1 )
{
  /* pay attention to T1s that use exclusively the    */
  /* output ports of SUM and (inverted) CBAR          */
  // if ( t1.has_sum && ( !t1.has_carry && !t1.has_carry_inverted ) && ( t1.has_cbar || t1.has_cbar_inverted ) )
  // {
  //   return true;
  // }

  /* pay attention to T1s that use more than 3 output */
  /* ports                                            */
  uint8_t num_output_ports_used{ 0u };
  if ( t1.has_sum )
  {
    ++num_output_ports_used;
  }
  if ( t1.has_carry )
  {
    ++num_output_ports_used;
  }
  if ( t1.has_carry_inverted )
  {
    ++num_output_ports_used;
  }
  if ( t1.has_cbar )
  {
    ++num_output_ports_used;
  }
  if ( t1.has_cbar_inverted )
  {
    ++num_output_ports_used;
  }
  if ( num_output_ports_used > 3 )
  {
    return true;
  }
  return false;
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
  if ( argc == 1 )
  {
    PHASES.push_back( 7 );
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
  auto benchmarks1 = epfl_benchmarks( experiments::arithmetic );
  //auto benchmarks1 = epfl_benchmarks( experiments::adder );
  //auto benchmarks2 = iscas_benchmarks( experiments::c432 | experiments::c880 | experiments::c1908 | experiments::c1355 | experiments::c3540 );
  //benchmarks1.insert(benchmarks1.end(), benchmarks2.begin(), benchmarks2.end());

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
  //const std::vector<std::string> ISCAS89_BENCHMARKS {"s382.aig", "s5378.aig", "s13207.aig"};
  //benchmarks1.insert(benchmarks1.end(), ISCAS89_BENCHMARKS.begin(), ISCAS89_BENCHMARKS.end());
  // std::reverse(benchmarks1.begin(), benchmarks1.end());

  // *** LIST ALL CONSIDERED BENCHMARKS ***
  fmt::print("Benchmarks:\n\t{}\n", fmt::join(benchmarks1, "\n\t"));
  
  // *** READ COMPOUND GATE LIBRARY ***
  phmap::flat_hash_map<ULL, Node> GNM_global;
  bool status = LoadFromFile(GNM_global, NODEMAP_BINARY_PREFIX);
  assert(status);

  phmap::flat_hash_map<std::string, LibEntry> entries = read_LibEntry_map(LibEntry_file);

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
    auto _result = decompose_to_klut(res_wo_pb, GNM_global, entries, COSTS_MAP);
    auto klut_decomposed = std::get<0>(_result);
    auto raw_area = std::get<1>(_result);
    fmt::print("Decomposition complete\n");

    // *** ENUMERATE 3-CUTS ***
    cut_enumeration_params ce_params; 
    ce_params.cut_size = 3u;
    const auto cuts = mockturtle::cut_enumeration<klut, true>( klut_decomposed, ce_params );

    /* print enumerated cuts */
    // klut_decomposed.foreach_node( [&]( auto node ) {
    //   auto idx = klut_decomposed.node_to_index( node );
    //   auto & node_cuts = cuts.cuts( idx );
    //   std::cout << node_cuts << "\n";
    // } );

    // *** FIND THOSE CUTS MATCHING THE XOR3/MAJ3/OR3 FUNCTIONS ***
    const auto xor3_cuts = match_cuts( xor3_tts, klut_decomposed, cuts);
    const auto maj3_cuts = match_cuts( maj3_tts, klut_decomposed, cuts);
    const auto  or3_cuts = match_cuts(  or3_tts, klut_decomposed, cuts);

    /* TODO: adopt the assumption that, it is a good deal if an implementation can make use of more than 2 out of the 3 outputs of T1,  */
    /* then we would need four rounds of matching: (1) 3 leaves using all 3 outputs; (2) ...using XOR and MAJ; (3) ...using MAJ and OR; */
    /* (4) ...using XOR and OR. For each 3 leaves found, check validity by deciding which of the 8 T1s to choose                        */
    array_map<3, T1_OUTPUTS> t1_candidates = find_t1_candidates(klut_decomposed, xor3_cuts, maj3_cuts, or3_cuts);

    /* estimate the gain of implementing parts of the circuits using T1s instead */
    auto updated_area{ raw_area };
    // uint32_t num_t1_use_more_than_3{ 0u };
    // uint32_t num_t1_cells{ 0u };

    /* Rewrote this using iterator output */
    for ( auto it_t1_cands{ t1_candidates.begin() }; it_t1_cands != t1_candidates.end(); )
    {
      if ( !t1_usage_sanity_check( klut_decomposed, *it_t1_cands, updated_area ) )
      {
        // update iterator after erasing
        it_t1_cands = t1_candidates.erase( it_t1_cands );
      }
      else
      {
        ++it_t1_cands;
      }
    }

    /* bind the instantiated T1 cells to the 2-LUT network for csv generation    */
    std::string filename = benchmark + "_before_phase_assignment.csv";
    write_klut_specs_supporting_t1( klut_decomposed, t1_candidates, filename );

    // fmt::print( "Usage of the T1 cells:\n" );
    // uint32_t num_t1_cells{ 0u };
    // for ( auto const& t1_candidate : t1_candidates )
    // {
    //   fmt::print( "[{}]\n", ++num_t1_cells );
    //   fmt::print( "\t\tInput : Node {}, Node {}, Node {}\n", t1_candidate.first[0], t1_candidate.first[1], t1_candidate.first[2] );
    //   t1_candidate.second.report();
    // }
    // fmt::print( "\t# T1s : {}, # T1s with more than 3 ouput ports used : {}\n", num_t1_cells, num_t1_use_more_than_3 );

    fmt::print( "Area before : {}, \tArea after : {}, \tRed. : {:>5.2f}%\n", raw_area, updated_area, ( ( static_cast<float>( raw_area ) - static_cast<float>( updated_area ) ) / static_cast<float>( raw_area ) * 100 ) );
    

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
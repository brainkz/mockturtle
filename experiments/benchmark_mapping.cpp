
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

#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/multiphase.hpp>
// #include <mockturtle/algorithms/compound_gate_mapping.hpp>

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

#include <chrono>



template <size_t N, typename T>
using array_map = phmap::flat_hash_map<std::array<klut::node, N>, T, ArrayHash<N>>;

// // Sunmagnetics Technology Library
// constexpr std::array<int,12> COSTS_MAP = {7, 9, 8, 8, 12, 8, 999, 999, 999, 8, 3, 0};
// Sunmagnetics Technology Library
constexpr std::array<int,12> COSTS_MAP = COSTS_SUNMAGNETICS;

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
  phmap::flat_hash_map<std::string, int> & nDFF_global, 
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

  uint32_t append( const uint64_t dff_hash, DFF_registry& DFF_REG )
  {
    DFF_var& dff = DFF_REG.at( dff_hash );
    std::vector<uint64_t>& head_section = sections.back();
    uint64_t& head_hash = head_section.back();
    DFF_var& head_dff = DFF_REG.at( head_hash );

    if ( dff.sigma == head_dff.sigma ) // add to the same section
    {
      head_section.push_back( dff_hash );
    }
    else
    {
      // assert( head_dff.phase - dff.phase == 1 );
      sections.push_back( { dff_hash } );
    }
    return sections.size();
  }

  /* for adding helper variables, which */
  /* are not stored in 'DFF_REG'        */
  void append( DFF_var const& dff )
  {
    sections.emplace_back( dff_hash( dff ) );
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
    if (it != path.targets.end() && ( ( fanout_sigma < 0 ) || ( ( fanout_sigma >= 0 ) && ( dff.sigma >= static_cast<uint32_t>( fanout_sigma ) ) ) ))
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
  fmt::print(iss.str());

  if (line.find("Objective value: ") != std::string::npos) 
  {
      objective_value = std::stoi(line.substr(17));
  } 
  else 
  {
      // Handle missing or incorrect format for objective value
      std::cerr << "Error: Objective value not found or invalid format." << std::endl;
      std::string command = fmt::format("{} {}", LAUNCH_CMD, cfg_name);
      system(command.c_str());
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
          std::string command = fmt::format("{} {}", LAUNCH_CMD, cfg_name);
          system(command.c_str());
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


int cpsat_ortools_union(const std::string & cfg_name, const uint8_t n_phases) 
{
  std::string command = fmt::format("{} {} {} {}", PYTHON_EXECUTABLE, PYTHON_DFF_PLACEMENT_UNION, cfg_name, n_phases);
  fmt::print("Executing command:\n{}\n", command);
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

enum GateTypes : uint8_t
{
  CONST0_FUNC,
  CONST1_FUNC,
  PI_FUNC,
  SPL_FUNC,
  MRG_FUNC,
  DFF_FUNC,
  NOT_FUNC,
  XOR_FUNC,
  AND_FUNC,
  OR_FUNC
};

namespace mockturtle {
template<class Ntk>
void write_bench_with_types( Ntk const& ntk, std::ostream& os)
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_get_constant_v<Ntk>, "Ntk does not implement the get_constant method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_num_pos_v<Ntk>, "Ntk does not implement the num_pos method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_node_function_v<Ntk>, "Ntk does not implement the node_function method" );

  // os << fmt::format( "n{}, gnd\n", ntk.node_to_index( ntk.get_node( ntk.get_constant( false ) ) ) );
  os << fmt::format( "{},{}\n", ntk.node_to_index( ntk.get_node( ntk.get_constant( false ) ) ) , CONST0_FUNC);
  if ( ntk.get_node( ntk.get_constant( false ) ) != ntk.get_node( ntk.get_constant( true ) ) )
  {
    // os << fmt::format( "n{}, vdd\n", ntk.node_to_index( ntk.get_node( ntk.get_constant( true ) ) ) );
    os << fmt::format( "{},{}\n", ntk.node_to_index( ntk.get_node( ntk.get_constant( true ) ) ) , CONST1_FUNC);
  }

  ntk.foreach_pi( [&]( auto const& n ) {
    // os << fmt::format( "INPUT(n{})\n", ntk.node_to_index( n ) );
    os << fmt::format( "{},{}\n", ntk.node_to_index( n ) , PI_FUNC );
  } );

  // for ( auto i = 0u; i < ntk.num_pos(); ++i )
  // {
  //   // os << fmt::format( "OUTPUT(po{})\n", i );
  //   os << fmt::format( "o,o{},PI_FUNCTION\n", i );
  // }


  ntk.foreach_node( [&]( auto const& n ) {
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return; /* continue */

    auto func = ntk.node_function( n );
    std::string children;
    auto first = true;
    ntk.foreach_fanin( n, [&]( auto const& c, auto i ) {
      if ( ntk.is_complemented( c ) )
      {
        kitty::flip_inplace( func, i );
      }
      if ( first )
      {
        first = false;
      }
      else
      {
        // children += ", ";
        children += "|";
      }

      children += fmt::format( "n{}", ntk.node_to_index( ntk.get_node( c ) ) );
    } );
    
    const NodeData & data = ntk.value(n);
    const uint8_t gate_type = static_cast<uint8_t>(data.type);
    const std::string func_hex = kitty::to_hex( func );

    std::string func_str;
    uint8_t func_code;
    if ( gate_type == AA_GATE )
    {
      // it can only be a merger, no other AA gate is supported
      assert( func_hex == "e");
      func_str = "MRG";
      func_code = MRG_FUNC;
    }
    else if ( gate_type == AS_GATE )
    {
      if ( func_hex == "6" && ntk.fanin_size(n) == 2)
      {
        func_str = "XOR";
        func_code = XOR_FUNC;
      }
      else if ( func_hex == "1" && ntk.fanin_size(n) == 1)
      {
        func_str = "NOT";
        func_code = NOT_FUNC;
      }
      else if ( func_hex == "2" && ntk.fanin_size(n) == 1)
      {
        fmt::print("Detected a DFF {}. This should not happen in a decomposed KLUT", n);
        func_str = "DFF";
        func_code = DFF_FUNC;
      }
      else 
      {
        throw std::runtime_error(fmt::format("Unsupported AS gate with function 0x{} and fanin {}", func_hex, ntk.fanin_size(n)));
      }
    }
    else if ( gate_type == SA_GATE )
    {
      if ( func_hex == "8" && ntk.fanin_size(n) == 2)
      {
        func_str = "AND";
        func_code = AND_FUNC;
      }
      else if ( func_hex == "e" && ntk.fanin_size(n) == 2)
      {
        func_str = "OR";
        func_code = OR_FUNC;
      }
      else 
      {
        throw std::runtime_error(fmt::format("Unsupported SA gate with function 0x{} and fanin {}", func_hex, ntk.fanin_size(n)));
      }
    }
    // os << fmt::format( "n{} = {} ({})\n",
    //                   ntk.node_to_index( n ),
    //                   func_str, children);
    os << fmt::format( "{},{},{}\n",
                      ntk.node_to_index( n ),
                      func_code, children);

    // os << fmt::format( "n{} = LUT 0x{} ({})\n",
    //                    ntk.node_to_index( n ),
    //                    kitty::to_hex( func ), children);
  } );

  /* outputs */
  // ntk.foreach_po( [&]( auto const& s, auto i ) {
  //   if ( ntk.is_constant( ntk.get_node( s ) ) )
  //   {
  //     os << fmt::format( "po{} = {}\n",
  //                        i,
  //                        ( ntk.constant_value( ntk.get_node( s ) ) ^ ntk.is_complemented( s ) ) ? "vdd" : "gnd" );
  //   }
  //   else
  //   {
  //     // os << fmt::format( "po{} = LUT 0x{} (n{})\n",
  //     //                    i,
  //     //                    ntk.is_complemented( s ) ? 1 : 2,
  //     //                    ntk.node_to_index( ntk.get_node( s ) ) );
  //     os << fmt::format( "po{} = DFF (n{})\n",
  //                        i, ntk.node_to_index( ntk.get_node( s ) ) );
  //   }
  // } );

  os << std::flush;
}
}

int main(int argc, char* argv[])  //
{
  using namespace experiments;
  using namespace mockturtle;

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
    const std::vector<std::tuple<std::string, std::string, std::string>> benchmarks { {
      // Arithmetic for T1 cells
      {DEFAULT_FOLDER, "adder",       "aig"},
      {DEFAULT_FOLDER, "sin",         "aig"},
      {DEFAULT_FOLDER, "voter",       "aig"},
      {DEFAULT_FOLDER, "square",      "aig"},
      {DEFAULT_FOLDER, "multiplier",  "aig"},
      {DEFAULT_FOLDER, "log2",        "aig"},
      {DEFAULT_FOLDER, "c7552",       "aig"},
      {DEFAULT_FOLDER, "c6288",       "aig"},

      // ASPDAC 2024
      {DEFAULT_FOLDER, "int2float", "aig"},
      {DEFAULT_FOLDER, "priority",  "aig"},
      {DEFAULT_FOLDER, "voter",     "aig"},
      {DEFAULT_FOLDER, "c432",      "aig"},
      {DEFAULT_FOLDER, "c880",      "aig"},
      {DEFAULT_FOLDER, "c1908",     "aig"},
      {DEFAULT_FOLDER, "c3540",     "aig"},
      {DEFAULT_FOLDER, "c1355",     "aig"},
      {ISCAS89_FOLDER, "s13207",    "aig"},
      {ISCAS89_FOLDER, "s5378",     "aig"},

      // Arxiv 2401.06411
      {ISCAS89_FOLDER, "s382",  "aig"},
      {ISCAS89_FOLDER, "s27",   "aig"},
      {ISCAS89_FOLDER, "s298",  "aig"},
      {ISCAS89_FOLDER, "s382",  "aig"},
      {ISCAS89_FOLDER, "s526",  "aig"},
      {DEFAULT_FOLDER, "c3540", "aig"},
      {DEFAULT_FOLDER, "c5315", "aig"},
      {DEFAULT_FOLDER, "c7552", "aig"},
      {DEFAULT_FOLDER, "c6288", "aig"},
      // {"!!!unknown", "Counter64", "???"},
      // {"!!!unknown, "Shoup", "???"},
      
      // Original ISVLSI
      {OPENCORES_FOLDER, "simple_spi-gates",    "blif"},
      {OPENCORES_FOLDER, "des_area-gates",      "blif"},
      {OPENCORES_FOLDER, "pci_bridge32-gates",  "blif"},
      {OPENCORES_FOLDER, "spi-gates",           "blif"},
      {OPENCORES_FOLDER, "mem_ctrl-gates",      "blif"},

      // VLSI-SoC 2023
      {DEFAULT_FOLDER, "sin",       "aig"},
      {DEFAULT_FOLDER, "cavlc",     "aig"},
      {DEFAULT_FOLDER, "dec",       "aig"},
      {DEFAULT_FOLDER, "int2float", "aig"},
      {DEFAULT_FOLDER, "priority",  "aig"},
      {DEFAULT_FOLDER, "c499",      "aig"},
      {DEFAULT_FOLDER, "c880",      "aig"},
      {DEFAULT_FOLDER, "c1908",     "aig"},
      {DEFAULT_FOLDER, "c3540",     "aig"},
      {DEFAULT_FOLDER, "c5315",     "aig"},
      {DEFAULT_FOLDER, "c7552",     "aig"},
    } };


    // // *** LIST ALL CONSIDERED BENCHMARKS ***
    // fmt::print("Benchmarks:\n\t{}\n", fmt::join(benchmarks, "\n\t"));

    // *** READ COMPOUND GATE LIBRARY ***
    phmap::flat_hash_map<ULL, Node> GNM_global;
    bool load_status = LoadFromFile(GNM_global, NODEMAP_BINARY_PREFIX);
    assert(load_status);

    phmap::flat_hash_map<std::string, LibEntry> entries = read_LibEntry_map(LibEntry_file);

  #pragma endregion benchmark_parsing

  phmap::flat_hash_set<std::string> completed {};

  // *** START PROCESSING BECNHMARKS ***
  for ( auto const& [folder, benchmark, type] : benchmarks )
  {

    if (completed.contains(benchmark))
    {
      fmt::print( "[i] Benchmark {} has already been processed. Skipping...\n", benchmark );
        continue;
    }

    fmt::print( "[i] processing {}\n", benchmark );

    #pragma region load network
    // *** LOAD NETWORK INTO MIG ***
    mig ntk_original;
    if (type == "blif")
    {
      fmt::print("USING THE BLIF READER\n");

      std::string abc_command = fmt::format("\"{}\" -c \"read_blif {}{}.blif\" -c strash -c \"write_aiger temp.aig\" ", ABC_EXECUTABLE, folder, benchmark);
      std::system(abc_command.c_str());

      if ( lorina::read_aiger( "temp.aig", aiger_reader( ntk_original ) ) != lorina::return_code::success )
      {
        fmt::print("Failed to read {}\n", benchmark);
        continue;
      }
    }
    // else if (type == "bench") 
    // {
    //   fmt::print("USING THE BENCH READER\n");
    //   const std::string path = fmt::format("{}{}", folder, benchmark);
    //   if ( lorina::read_bench( path, bench_reader( ntk_original ) ) != lorina::return_code::success )
    //   {
    //     fmt::print("Failed to read {}\n", path);
    //     continue;
    //   }
    // }
    else if ( type == "aig" )
    {
      fmt::print( "USING THE AIGER READER\n" );
      const std::string path = fmt::format("{}{}.aig", folder, benchmark);
      if ( lorina::read_aiger( path, aiger_reader( ntk_original ) ) != lorina::return_code::success )
      {
        fmt::print("Failed to read {}\n", path);
        continue;
      }
    }
    else
    {
      fmt::print("Unknown benchmark file type\n");
      fmt::print("Folder {}\n", folder);
      fmt::print("Benchmark {}\n", benchmark);
      fmt::print("Type {}\n", type);
      continue;
    }
    #pragma endregion

    #pragma region mapping with compound gates 
    // *** MAP, NO NEED FOR RETIMING/PATH BALANCING ***
    fmt::print("Started mapping {}\n", benchmark);
    auto [res_wo_pb, st_wo_pb] = map_wo_pb(ntk_original, tech_lib, false); //benchmark, true, nDFF_global, total_ndff_w_pb, total_area_w_pb, cec_w_pb 
    fmt::print("Finished mapping {}\n", benchmark);
    #pragma endregion

    #pragma region decomposition of the mapped network into a klut
    // *** DECOMPOSE COMPOUND GATES INTO PRIMITIVES, REMOVE DFFS, REPLACE OR GATES WITH CB WHERE POSSIBLE ***
    auto _result = decompose_to_klut(res_wo_pb, GNM_global, entries, COSTS_MAP);
    auto klut_decomposed = std::get<0>(_result);
    auto raw_area = std::get<1>(_result);
    fmt::print("Decomposition complete\n");
    #pragma endregion

    const auto paths = extract_paths(klut_decomposed, true);
    for (const auto & path : paths)
    {
      const auto threads = path.path_threads(klut_decomposed);
      
    }
    
    // print in a digestible format

    // std::string output_filename { fmt::format(".\\csv\\{}_compound.csv", benchmark) };
    // std::ofstream file(output_filename);
    // if (file.is_open()) {
    //     write_bench_with_types(klut_decomposed, file);
    // } else {
    //     std::cerr << "Error writing the output file " << output_filename << std::endl;
    // }
    completed.emplace(benchmark);
  }
  return 0;
}

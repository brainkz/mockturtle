
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

    fmt::print("Compatible TTs:\n");
    for (auto i = 0u; i < xor3_tts.size(); ++i)
    {
      fmt::print("\t[{}]:\n", i);
      fmt::print("\t\tXOR3: {0:08b}={0:02x}={0:d}\n", xor3_tts[i]._bits);
      fmt::print("\t\tMAJ3: {0:08b}={0:02x}={0:d}\n", maj3_tts[i]._bits);
      fmt::print("\t\t OR3: {0:08b}={0:02x}={0:d}\n",  or3_tts[i]._bits);
    }
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
      // experiments::adder | experiments::div  | 
      auto benchmarks1 = epfl_benchmarks( experiments::sin  | experiments::multiplier );
      // auto benchmarks1 = epfl_benchmarks( experiments::int2float | experiments::priority | experiments::voter);
      // auto benchmarks2 = iscas_benchmarks( experiments::c432 | experiments::c880 | experiments::c1908 | experiments::c1355 | experiments::c3540 );
      // benchmarks1.insert(benchmarks1.end(), benchmarks2.begin(), benchmarks2.end());

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
      phmap::flat_hash_map<ULL, Node> GNM_global;
      bool status = LoadFromFile(GNM_global, NODEMAP_BINARY_PREFIX);
      assert(status);

      phmap::flat_hash_map<std::string, LibEntry> entries = read_LibEntry_map(LibEntry_file);

    #pragma endregion benchmark_parsing

    // *** START PROCESSING BECNHMARKS ***
    for ( auto const& benchmark : benchmarks1 )
    {
      fmt::print( "[i] processing {}\n", benchmark );

      #pragma region load network
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

      // TODO: To be replaced with Mingfei's code
      #pragma region cut enumeration to find TTs that could be shared with a T1 cell
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
          auto qq = *cut_entry;
          const auto tt = cuts.truth_table( *cut_entry );
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
          // std::ostringstream oss;
          // oss << *cut_entry;
          // fmt::print("[{}] \tCut {} with tt 0b{}=0x{}\n", benchmark, oss.str(), kitty::to_binary( tt ), kitty::to_hex( tt ));
        }
        // flat_hash_set<std::array<uint32_t, 3>> arraySet;
        // auto tt = cuts.truth_table( node_cut_set );
        // std::cout << node_cut_set << " -> " <<  << "\n";
        // cuts.truth_table( node_cut_set );
      } );

      for (const auto & [leaves, xor3_tuple] : xor3_cuts) 
      {
        const auto xor3_tt = std::get<0>(xor3_tuple);
        // const auto xor3_root = std::get<1>(xor3_tuple);
        if (maj3_cuts.contains(leaves)) 
        {
          auto maj3_tuple = maj3_cuts[leaves];
          auto maj3_tt = std::get<0>(maj3_tuple);
          auto maj3_root = std::get<1>(maj3_tuple);
          // cut_intersection.push_back(xor3_nodes);
          // const auto & node_cut_set = cuts.cuts( maj3_root );

          fmt::print("[{0}]\tCut [{1}]\n\t XOR3 tt 0b{2:08b}=0x{2:02x}\n\t MAJ3 tt 0b{3:08b}=0x{3:02x}\n", benchmark, fmt::join(leaves, ","), xor3_tt._bits[0], maj3_tt._bits[0]);
          // for (const auto & cut_entry : node_cut_set)
          // {
          //   const auto tt = cuts.truth_table( *cut_entry );
          //   std::ostringstream oss;
          //   oss << *cut_entry;
          // }
        }
      }

      // return 0;

      // for (const auto & cut : xor3_cuts)
      // {
      //   if ( std::find( maj3_cuts.begin(), maj3_cuts.end(), cut ) != maj3_cuts.end() )
      //   {
      //     auto tt = cuts.truth_table( node_cut_set );
      //     std::cout << node_cut_set << " -> " <<  << "\n";
      //     cuts.truth_table( node_cut_set );
      //   }
      // }

      #pragma endregion

      // TODO: can be removed, in principle
      // *** [temporary] GREEDILY ASSIGN A STAGE (sigma) TO EACH ELEMENT *** 
      // std::unordered_map<klut::signal, glob_phase_t> glob_phase = greedy_assign(klut_decomposed, klut_prim_params, false); // printUnorderedMap(glob_phase);

      // for (auto n_phases = MIN_N_PHASES; n_phases <= MAX_N_PHASES; ++n_phases)
      for (const auto n_phases : PHASES)
      {
        fmt::print("[i] Mapping with {} phases\n", n_phases);
        // *** IF i = 0, we assign phases with the CP-SAT
        // *** IF i = 1, we assign phases greedily (currently unused)
        // for (auto i = 0; i < 2; ++i)
        for (auto i = 0; i < 1; ++i)
        {
          klut network { klut_decomposed.clone() };
          phmap::flat_hash_map<unsigned int, unsigned int> assignment;

          // *** IF i = 0, "assignment" has stages assigned by the CP-SAT
          // *** IF i = 1, "assignment" is empty
          if (i == 0)
          {
            const std::string ilp_cfg_filename = fmt::format("ilp_configs/{}.csv", benchmark);
            fmt::print("\tWriting config {}\n", ilp_cfg_filename);
            write_klut_specs(network, ilp_cfg_filename);

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
          greedy_ntk_assign(network, n_phases, assignment, true );

          // *** Greedily insert splitters
          splitter_ntk_insertion( network, true);

          // network.foreach_node([&] ( const klut::signal & node ) {if ( network.fanout_size( node ) > 1 ){assert( network.node_function( node ) == 0x2 );};});

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

#if false
  int main() 
  {
    /*
    // *** READ COMPOUND GATE LIBRARY ***
    const std::vector<std::vector<UI>> sets_of_levels { { {0,0,0,0}, {0,0,0,1}, {0,0,0,2}, {0,0,1,1}, {0,0,1,2}, {0,1,1,1}, {0,1,1,2}, {0,1,2,2}, {0,1,2,3} } }; //  {0,1,1,3},

    auto start_time_csv = std::chrono::high_resolution_clock::now();

    phmap::flat_hash_map<ULL, Node> GNM_global = read_global_gnm( sets_of_levels, NODEMAP_PREFIX );

    auto end_time_csv = std::chrono::high_resolution_clock::now();

    // Calculate the duration in milliseconds
    auto ms_csv = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_csv - start_time_csv);
    fmt::print("CSV Runtime: {}ms\n", ms_csv.count());

    const auto cost_it = std::max_element(GNM_global.begin(), GNM_global.end(), [&](const auto & pair1, const auto & pair2) { return pair1.second.cost < pair2.second.cost; });
    const auto & [hash_cost, node_cost] = *cost_it;
    fmt::print("The most expensive node: {} {}\n{}\n", hash_cost, node_cost.to_str(), node_cost.to_stack(GNM_global));

    const auto depth_it = std::max_element(GNM_global.begin(), GNM_global.end(), [&](const auto & pair1, const auto & pair2) 
    {
      return pair1.second.depth < pair2.second.depth; 
    });
    const auto & [hash_depth, node_depth] = *depth_it;
    fmt::print("The deepest node: {} {}\n{}\n", hash_depth, node_depth.to_str(), node_depth.to_stack(GNM_global));
    */
    const std::string DATABASE_PREFIX { "../GNM/GNM_global" };

    // // Dump the map to a binary file
    // SaveToFile(GNM_global, DATABASE_PREFIX);
    // fmt::print("Saved to {}\n", DATABASE_PREFIX);
    auto start_time_dat = std::chrono::high_resolution_clock::now();

    // Load the map from the binary file
    phmap::flat_hash_map<uint64_t, Node> loadedMap;
    bool status = LoadFromFile(loadedMap, DATABASE_PREFIX);
    if (!status)
    {
      fmt::print("READING FAILED");
      return 1;
    };

    auto end_time_dat = std::chrono::high_resolution_clock::now();

    auto ms_dat = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_dat - start_time_dat);

    fmt::print("DAT Runtime: {}ms\n", ms_dat.count());
    /*

    for (const auto& [hash, old_node] : GNM_global) 
    {
      if (loadedMap.find(hash) == loadedMap.end())
      {
        fmt::print("Entry at {} is not found\n", hash);
        continue;
      };
      const Node & new_node = loadedMap.at(hash);
      if (old_node != new_node)
      {
        fmt::print("Entries at {} are not equivalent\n", hash);
        fmt::print("\tOld: {}\n", old_node.to_str());
        fmt::print("\tNew: {}\n", new_node.to_str());
        continue;
      };
    }

    */
   
    auto start_time_par = std::chrono::high_resolution_clock::now();

    // Load the map from the binary file
    phmap::flat_hash_map<uint64_t, Node> loadedMapPar = ParallelLoadFromFiles(DATABASE_PREFIX);

    auto end_time_par = std::chrono::high_resolution_clock::now();

    auto ms_par = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_par - start_time_par);

    fmt::print("PAR Runtime: {}ms\n", ms_par.count());

    for (const auto& [hash, old_node] : loadedMap) 
    {
      if (loadedMapPar.find(hash) == loadedMapPar.end())
      {
        fmt::print("Entry at {} is not found\n", hash);
        continue;
      };
      const Node & new_node = loadedMapPar.at(hash);
      if (old_node != new_node)
      {
        fmt::print("Entries at {} are not equivalent\n", hash);
        fmt::print("\tOld: {}\n", old_node.to_str());
        fmt::print("\tNew: {}\n", new_node.to_str());
        continue;
      };
    }

    

    return 0;
  }

#endif
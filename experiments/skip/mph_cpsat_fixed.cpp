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
#include <thread>

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
#include <mockturtle/views/mph_view.hpp>
#include <mockturtle/algorithms/functional_reduction.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>

// mockturtle/algorithms/mig_algebraic_rewriting.hpp

#include <mockturtle/io/genlib_utils.hpp>

// #include <mockturtle/utils/GNM_global.hpp> // GNM global is stored here

#include <mockturtle/utils/misc.hpp>

#include <experiments.hpp>

constexpr uint8_t NUM_PHASES { 4u };

typedef mockturtle::mph_view<mockturtle::klut_network, NUM_PHASES> mph_klut;
typedef uint64_t node_t;

constexpr std::array<int, 12> COSTS_MAP = COSTS_SUNMAGNETICS;

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

std::tuple<int, phmap::flat_hash_map<unsigned int, unsigned int>, std::string>  cpsat_macro_opt(const std::string & cfg_name, uint8_t n_phases) 
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
  phmap::flat_hash_map<unsigned int, unsigned int> key_value_pairs;

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

// Function to read unordered_map from CSV file
phmap::flat_hash_map<std::string, int> readCSV(const std::string& filename) 
{
  std::ifstream infile(filename);             // Open the input file stream
  phmap::flat_hash_map<std::string, int> map;   // Create the unordered_map
  
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

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, double, double, double, double, double, double, double, double, double> exp( "mapper", "benchmark", "N_PHASES", "#DFF", "area", "delay", "mapping", "phase_assgn", "spl_insertion", "DFF_insertion", "total");

  fmt::print( "[i] processing technology library\n" );

  // library to map to technology
  std::vector<gate> gates;
  std::ifstream inputFile( DATABASE_PATH );
  if ( lorina::read_genlib( inputFile, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  mockturtle::tech_library_params tps; // tps.verbose = true;
  tech_library<NUM_VARS, mockturtle::classification_type::p_configurations> tech_lib( gates, tps );

  #pragma region benchmark_parsing
    // *** BENCHMARKS OF INTEREST ***
    // auto benchmarks1 = epfl_benchmarks( experiments::int2float | experiments::priority | experiments::voter);
    // auto benchmarks1 = epfl_benchmarks( experiments::adder );
    // auto benchmarks2 = iscas_benchmarks( experiments::c432 | experiments::c880 | experiments::c1908 | experiments::c1355 | experiments::c3540 );

    auto benchmarks1 = all_benchmarks( 
      experiments::int2float | 
      experiments::priority |
      experiments::voter  |
      experiments::c432 |
      experiments::c880 |
      experiments::c1908  |
      experiments::c3540  |
      experiments::c1355 
    );
    // benchmarks1.insert(benchmarks1.end(), benchmarks2.begin(), benchmarks2.end());

    // *** READ COMPOUND GATE LIBRARY ***
    phmap::flat_hash_map<ULL, Node> GNM_global;
    bool status = LoadFromFile(GNM_global, NODEMAP_BINARY_PREFIX);
    assert(status);
    
    phmap::flat_hash_map<std::string, LibEntry> entries = read_LibEntry_map(LibEntry_file);

  #pragma endregion

  // *** START PROCESSING BECNHMARKS ***
  for ( auto const& benchmark : benchmarks1 )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    // *** LOAD NETWORK INTO MIG ***
    mockturtle::mig_network ntk_original;
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
      std::string path = fmt::format("{}/{}", ISCAS89_FOLDER, benchmark);
      fmt::print("READING {}\n", path);

      // std::filesystem::path currentDir = std::filesystem::current_path();
      // // Convert the path to a string and print it
      // std::string currentDirStr = currentDir.string();
      // std::cout << "Current Working Directory: " << currentDirStr << std::endl;

      if ( lorina::read_aiger( path, aiger_reader( ntk_original ) ) != lorina::return_code::success )
      {
        fmt::print("Failed to read {}\n", benchmark);
        continue;
      }
    }
    else // regular benchmark
    {
      fmt::print("USING THE AIGER READER\n");
      if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( ntk_original ) ) != lorina::return_code::success )
      {
        fmt::print("Failed to read {}\n", benchmark);
        continue;
      }
      // convert_klut_to_graph<mig>(ntk_original, temp_klut);
    }

    std::chrono::high_resolution_clock::time_point mapping_start = std::chrono::high_resolution_clock::now();

    // *** MAP, NO NEED FOR RETIMING/PATH BALANCING AT THIS TIME ***
    fmt::print("Started mapping {}\n", benchmark);
    auto [mapped_ntk, mapper_stats] = map_wo_pb(ntk_original, tech_lib, false); 
    fmt::print("Finished mapping {}\n", benchmark);

    // *** DECOMPOSE COMPOUND GATES INTO PRIMITIVES, REMOVE DFFS, REPLACE OR GATES WITH CB WHERE POSSIBLE ***
    auto _result = decompose_to_klut<NUM_PHASES>(mapped_ntk, GNM_global, entries, COSTS_MAP);
    auto klut_decomposed = std::get<0>(_result);
    auto raw_area = std::get<1>(_result);
    fmt::print("Decomposition complete\n");

    // Stop the timer
    std::chrono::high_resolution_clock::time_point mapping_end = std::chrono::high_resolution_clock::now();

    // Calculate elapsed time
    std::chrono::duration<double> mapping_seconds = mapping_end - mapping_start;
  
    fmt::print("[i] Mapping with {} phases\n", NUM_PHASES);
    // *** IF i = 0, we assign phases with the CP-SAT
    // *** IF i = 1, we assign phases greedily
    // for (auto i = 0; i < 2; ++i)
    for (auto i = 0; i < 1; ++i)
    {
      mph_klut network { klut_decomposed.clone() };
      phmap::flat_hash_map<unsigned int, unsigned int> assignment;

      std::chrono::high_resolution_clock::time_point pa_start = std::chrono::high_resolution_clock::now();

      if (i == 0)
      {
        const std::string ilp_cfg_filename = fmt::format("ilp_configs/{}.csv", benchmark);
        fmt::print("\tWriting config {}\n", ilp_cfg_filename);
        write_klut_specs(network, ilp_cfg_filename);

        fmt::print("\tCalling OR-Tools\n");
        auto [obj_val, assignment_local, status] = cpsat_macro_opt(ilp_cfg_filename, NUM_PHASES);

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
      
      assign_stages(network, NUM_PHASES, assignment, false );

      // Stop the timer
      std::chrono::high_resolution_clock::time_point pa_end = std::chrono::high_resolution_clock::now();

      // Calculate elapsed time
      std::chrono::duration<double> pa_seconds = pa_end - pa_start;


      // Start the timer
      std::chrono::high_resolution_clock::time_point spl_start = std::chrono::high_resolution_clock::now();

      // *** Greedily insert splitters
      splitter_ntk_insertion( network, false );

      // Stop the timer
      std::chrono::high_resolution_clock::time_point spl_end = std::chrono::high_resolution_clock::now();

      // Calculate elapsed time
      std::chrono::duration<double> spl_seconds = spl_end - spl_start;

      fmt::print("[i] FINISHED PHASE ASSIGNMENT\n");

      // Start the timer
      std::chrono::high_resolution_clock::time_point dff_start = std::chrono::high_resolution_clock::now();


      fmt::print("[i] EXTRACTING PATHS\n");
      std::vector<Path<NUM_PHASES>> paths = extract_paths( network, false );
      // auto [DFF_REG, precalc_ndff] = dff_vars(NR, paths, N_PHASES);

      std::vector<int> local_num_dff( paths.size() );
      std::vector<std::thread> threads;
      threads.reserve(paths.size());

      int idx = 0;
      for (const Path<NUM_PHASES> & path : paths)
      {
        threads.emplace_back(process_path<NUM_PHASES>, std::ref(local_num_dff), idx, std::ref(path), std::ref(network), benchmark);
        idx++;
      }
      // Wait for all threads to finish
      for (auto& thread : threads) {
          thread.join();
      }

      fmt::print("{}\n", fmt::join(local_num_dff, ","));
      auto total_num_dff = std::accumulate(local_num_dff.begin(), local_num_dff.end(), 0);


      // Stop the timer
      std::chrono::high_resolution_clock::time_point dff_end = std::chrono::high_resolution_clock::now();

      // Calculate elapsed time
      std::chrono::duration<double> dff_seconds = dff_end - dff_start;

      // *** Record maximum phase
      uint64_t max_stage = 0u;
      // *** Record number of splitters and total number of DFFs (not only path balancing DFFs)
      uint64_t total_num_spl = 0;
      network.foreach_gate([&](const klut::signal & node)
      {
        // fmt::print("[Node {}] old max_stage = {}\tnode_data = {}\t", node, max_stage, static_cast<int>(node_data.stage));
        max_stage = generic_max(max_stage, network.get_stage(node));
        // fmt::print("new max_stage = {}\n", max_stage);

        auto fo_size = network.fanout_size(node);
        if (fo_size > 1)
        {
          total_num_spl += fo_size - 1;
        }
      });

      network.foreach_po([&](const klut::signal & node)
      {
        // fmt::print("[PO {}] max_stage = {}, sigma = {}, node #DFF = {}\n", node, max_stage, static_cast<int>(node_data.stage), ( (max_stage - node_data.stage) / n_phases ) );
        total_num_dff += (max_stage - network.get_stage(node)) / NUM_PHASES;
        // fmt::print("\t\t\t\t[i] total #DFF = {}\n", total_num_dff);
      });

      // fmt::print("{} PHASES: #DFF   for {} is {}\n", n_phases, benchmark, total_num_dff);
      int total_area = raw_area + total_num_dff * COSTS_MAP[fDFF] + total_num_spl * COSTS_MAP[fSPL];
      // fmt::print("{} PHASES: #AREA  for {} is {}\n", n_phases, benchmark, total_area);
      // fmt::print("{} PHASES: #MAX GLOB PHASE for {} is {}\n", n_phases, benchmark, max_stage);

      exp(fmt::format("{}_{}", benchmark, (i==0)?"CPSAT":"GREEDY"), NUM_PHASES, total_num_dff, total_area, ( (max_stage - 1) / NUM_PHASES + 1 ),
      mapping_seconds.count(), pa_seconds.count(), spl_seconds.count(), dff_seconds.count(),
      mapping_seconds.count() + pa_seconds.count() + spl_seconds.count() + dff_seconds.count());
      exp.save();
      exp.table();
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
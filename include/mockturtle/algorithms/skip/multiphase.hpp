
#pragma once

#include <mockturtle/io/genlib_utils.hpp>
#include <mockturtle/algorithms/nodes.hpp>
#include <mockturtle/utils/misc.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/rsfq_view.hpp>
#include <mockturtle/views/fanout_view.hpp>
#include <mockturtle/views/mph_view.hpp>

typedef uint32_t stage_t;

typedef mockturtle::klut_network klut;
typedef uint64_t node_t;


const std::vector<std::string> GATE_TYPE { "PI", "AA", "AS", "SA", "T1" }; //, "PO"

template <uint8_t NUM_PHASES>
struct Path
{
  std::set<klut::signal> sources;   // AS/SA gates
  std::set<klut::signal> internals; // AA    gates
  std::set<klut::signal> targets;   // AS/SA gates
  Path(const std::set<klut::signal> & _sources, const std::set<klut::signal>& _internals, const std::set<klut::signal>& _targets)
    : sources(_sources), internals(_internals), targets(_targets) {}
  Path() : sources({}), internals({}), targets({}) {}
  
  void absorb(Path<NUM_PHASES> & other)
  {
    sources.insert(other.sources.begin(), other.sources.end());
    internals.insert(other.internals.begin(), other.internals.end());
    targets.insert(other.targets.begin(), other.targets.end());
  }

  void print() const
  {
    fmt::print(format());
  }

  std::string format() const
  {
    return fmt::format("Path from [{}]\n\tvia [{}]\n\tto [{}]\n", fmt::join(sources, ","), fmt::join(internals, ","), fmt::join(targets, ","));
  }

  void print_bfs(const mockturtle::mph_view<mockturtle::klut_network, NUM_PHASES> & ntk) const
  {
    print();
    std::deque<klut::signal> queue { targets.begin(), targets.end() };
    while (!queue.empty())
    {
      const klut::signal & node = queue.front();
      queue.pop_front();

      const std::vector<klut::signal> fanins { preds(node, ntk) };

      fmt::print("\t[{} {} {}]\tσ={}\tfanins : {{ {} }} \n", GATE_TYPE[ntk.get_type(node)], _kind(node), node, ntk.get_stage(node), fmt::join(fanins, ", "));

      for (const klut::signal & fanin : fanins)
      {
        queue.push_back(fanin);
      }
    }
    return;
  }

  std::vector<klut::signal> preds(const klut::signal & sig, const klut & ntk) const
  {
    if ( sources.count(sig) != 0)
    {
      return {};
    }

    std::vector<klut::signal> predecessors;

    ntk.foreach_valid_fanin(sig, [&](const klut::signal & parent)
    {
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

  std::string _kind( const klut::signal & node ) const
  {
    if ( targets.count(node) > 0 )  
    { 
      return "Target";  
    }
    else if ( internals.count(node) > 0 )  
    { 
      return "Internal";  
    }
    else if ( sources.count(node) > 0 )  
    { 
      return "Source";  
    }
    else 
    {
      throw;
    }
  }

  /// @brief Find paths of signals in Path object from target signals to their sources.
  /// This function performs a depth-first search starting from target signals and follows
  /// their predecessors until no more predecessors are found, creating threads of signals.
  /// 
  /// @param ntk The logic network.
  /// @return A vector of vectors, where each inner vector represents a thread of signals
  ///         from a target signal to its sources.
  std::vector<std::vector<klut::signal>> path_threads(const klut & ntk, bool verbose = false) const
  {
    // Initialize two vectors to store incomplete and completed threads
    std::vector<std::vector<klut::signal>> incomplete_threads;
    std::vector<std::vector<klut::signal>> done_threads;

    // Start a thread for each target signal
    for (const klut::signal & tgt : targets)
    {
      incomplete_threads.push_back({ tgt });
    }

    DEBUG_PRINT("[path_threads] CREATED {} incomplete_threads\n", incomplete_threads.size());
    for (auto thread: incomplete_threads)
    {
      DEBUG_PRINT("[path_threads]\tincomplete thread: [{}]\n", fmt::join(thread, ","));
    }

    // Continue until there are no incomplete threads left
    while (!incomplete_threads.empty())
    {
      // Take the last thread from the incomplete_threads vector
      std::vector<klut::signal> thread = incomplete_threads.back();
      incomplete_threads.pop_back();

      DEBUG_PRINT("[path_threads] processing thread: [{}]\n", fmt::join(thread, ","));

      // Get the gate signal at the end of the current thread
      const klut::signal & gate = thread.back();

      // Find predecessors of the gate signal in the network
      std::vector<klut::signal> predecessors { preds(gate, ntk) };

      DEBUG_PRINT("[path_threads]\tpredecessors: [{}]\n", fmt::join(predecessors, ","));

      // If there are no predecessors, this thread is complete
      if (predecessors.empty())
      {
        done_threads.push_back( thread );
        for (auto thread: done_threads)
        {
          DEBUG_PRINT("[path_threads] thread complete: [{}]\n", fmt::join(thread, ","));
        }
      }
      else
      {
        // If there are predecessors, extend the thread with each predecessor
        for ( const klut::signal & predecessor : predecessors )
        {
          incomplete_threads.push_back( thread );
          incomplete_threads.back().push_back(predecessor);
          DEBUG_PRINT("[path_threads]\tnew thread: [{}]\n", fmt::join(incomplete_threads.back(), ","));
        }
      }
    }

    // Return the completed threads
    return done_threads;
  } 
};


/// @brief Structure representing the potential DFF location. Uniquely defined by fanin, fanout and stage
struct DFF_var 
{
  klut::signal fanin;
  klut::signal fanout;
  uint32_t stage;
  std::unordered_set<uint64_t> parent_hashes;

  DFF_var(klut::signal _fanin, klut::signal _fanout, uint32_t _stage, std::unordered_set<uint64_t> _parent_hashes = {})
      : fanin(_fanin), fanout(_fanout), stage(_stage), parent_hashes(_parent_hashes) {}

  DFF_var(const DFF_var& other)
      : fanin(other.fanin), fanout(other.fanout), stage(other.stage), parent_hashes(other.parent_hashes) {}

  DFF_var( klut::signal _index, uint32_t _stage, std::unordered_set<uint64_t> _parent_hashes = {} )
      : fanin( 0u ), fanout( _index ), stage( _stage ), parent_hashes( _parent_hashes ) {}

  DFF_var( uint64_t dff_hash, std::unordered_set<uint64_t> _parent_hashes = {} )
      : fanin( (uint64_t)dff_hash >> 40 ), fanout( (uint64_t)( dff_hash << 24 ) >> 40 ), stage( (uint64_t)( dff_hash & 0xFFFF ) ), parent_hashes( _parent_hashes ) {}

  std::string str() const
  {
    if ( fanin == 0 )
    {
      return fmt::format( "gate_{}_{}", fanout, stage );
    }

    return fmt::format("var_{}_{}_{}", fanin, fanout, stage);
  }
};

uint64_t dff_hash(klut::signal _fanin, klut::signal _fanout, uint32_t _stage)
{
  return ( (uint64_t)_fanin << 40 ) | ( (uint64_t)_fanout << 16 ) | _stage;
}

uint64_t dff_hash( DFF_var const& dff )
{
  return ( (uint64_t)dff.fanin << 40 ) | ( (uint64_t)dff.fanout << 16 ) | dff.stage;
}

/// @brief Enhanced map of the DFF variables for easy tracking of DFFs
struct DFF_registry
{
  phmap::flat_hash_map<uint64_t, DFF_var> variables;

  DFF_var & at(node_t _fanin, node_t _fanout, stage_t _stage)
  {
    return variables.at( dff_hash(_fanin, _fanout, _stage) );
  } 
  DFF_var & at(uint64_t _hash)
  {
    return variables.at( _hash );
  } 
  uint64_t add(node_t _fanin, node_t _fanout, stage_t _phase, std::unordered_set<uint64_t> _parent_hashes = {})
  {
    uint64_t _hash = dff_hash(_phase, _fanout, _fanin);
    DFF_var temp { _fanin, _fanout, _phase, _parent_hashes };
    variables.emplace(_hash, temp);
    return _hash;
  } 

  std::string str(uint64_t hash, bool negated = false) const
  {
    const DFF_var & dff = variables.at(hash);
    if (negated)
    {
      return fmt::format( "var_{}_{}_{}.Not()", dff.fanin, dff.fanout, dff.stage );
    }
    return fmt::format( "var_{}_{}_{}", dff.fanin, dff.fanout, dff.stage );
  }
};



template <uint8_t NUM_PHASES>
struct Snake
{
  std::deque<std::vector<uint64_t>> sections;

  Snake(): sections({}) {}
  Snake( const uint64_t head ): sections({ { head } }) {}
  Snake( const std::deque<std::vector<uint64_t>> _sections ): sections(_sections) {}
  Snake( const Snake & _other ): sections(_other.sections) {}

  bool append(const uint64_t dff_hash, DFF_registry &DFF_REG)
  {
    DFF_var & dff = DFF_REG.at( dff_hash );

    std::vector<uint64_t> & head_section = sections.back();
    uint64_t & head_hash = head_section.back();
    DFF_var & head_dff = DFF_REG.at( head_hash );
    if (dff.stage == head_dff.stage) // add to the same section
    {
      head_section.push_back( dff_hash );
      return false;
    }
    else
    {
      // assert( head_dff.phase - dff.phase == 1 );
      sections.push_back( { dff_hash } );
      if (sections.size() > NUM_PHASES)
      {
        sections.pop_front();
      }
      return true;
    }
  }
};

template <uint8_t NUM_PHASES>
void write_snakes(const std::vector<Snake<NUM_PHASES>> & snakes, DFF_registry & DFF_REG, const std::vector<uint64_t> & required_SA_DFFs, const std::string cfg_name, bool verbose = false)
{
  std::ofstream spec_file (cfg_name);

  for (const auto & snake : snakes)
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
    if (vars_bucket.size() == NUM_PHASES)
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
/// @param NUM_PHASES 
/// @param verbose 
/// @return 
template <uint8_t NUM_PHASES>
std::vector<Snake<NUM_PHASES>> sectional_snake(const Path<NUM_PHASES> & path, const mph_klut & ntk,  DFF_registry & DFF_REG, const bool verbose = false)
{
  std::vector<Snake<NUM_PHASES>> out_snakes; 
  std::vector<Snake<NUM_PHASES>> stack;
  
  DEBUG_PRINT("[i]: Starting extraction of worms \n");
  // get all DFFs 
  for (const auto & [hash, dff]: DFF_REG.variables)
  {
    const auto [fo_stage, fo_type] = ntk.get_stage_type(dff.fanout);
    auto fanout_stage = fo_stage - ( fo_type == AS_GATE );
    auto it = std::find(path.targets.begin(), path.targets.end(), dff.fanout);
    if (it != path.targets.end() && ( ( fanout_stage < 0 ) || ( ( fanout_stage >= 0 ) && ( dff.stage >= static_cast<uint32_t>( fanout_stage ) ) ) ))
    {
      stack.emplace_back( hash );
    }
  }
  
  while (!stack.empty())
  {
    DEBUG_PRINT("[i] Stack size is {} \n", stack.size());
    auto snake = stack.back();
    stack.pop_back();

    DEBUG_PRINT("\t[i] The snake has {} sections\n", snake.sections.size());
    uint64_t hash = snake.sections.back().back();
    DFF_var & dff = DFF_REG.at( hash );

    // fmt::print("\tCurrent worm size {}, between phases {} and {} \n", worm.size(), DFF_REG.str(worm.front()), DFF_REG.str(worm.back()));


    DEBUG_PRINT("\t\t[i] The DFF {} has {} parents\n", DFF_REG.at( hash ).str(),  dff.parent_hashes.size() );

    bool returned_current_snake = false;
    for (const uint64_t parent_hash : dff.parent_hashes)
    {
      auto snake_copy = snake; 
      DEBUG_PRINT("\t\t[i] Advancing towards fanin {}\n", DFF_REG.at( parent_hash ).str() );
      bool status = snake_copy.append(parent_hash, DFF_REG);
      DEBUG_PRINT((status) ? "\t\t\tAdded new section!\n" :"\t\t\tExtended existing section!\n"  );
      DEBUG_PRINT("\t\t\tThe new length is {}\n", snake_copy.sections.size() );
      
      stack.push_back( snake_copy );
      if (status && !returned_current_snake && snake_copy.sections.size() == NUM_PHASES)
      {
        DEBUG_PRINT("\t\tAdding the snake to the output\n");
        out_snakes.push_back(snake);
        returned_current_snake = true;
      }
    }
  }
  return out_snakes;
}


template <uint8_t NUM_PHASES>
class multiphase_balancing_impl
{
public:
  multiphase_balancing_impl( 
    mockturtle::binding_view<mockturtle::klut_network> src, 
    phmap::flat_hash_map<ULL, Node> nodemap,
    phmap::flat_hash_map<std::string, LibEntry> entries,
    const std::array<int, 12> COSTS_MAP 
  )


public:
  /// @brief Assigns stages to nodes based on stage assignment. If the assignment is not found, assigns a stage greedily
  /// @param ntk - multiphase view of the network
  /// @param phase_assignment - mapping from node index to assigned stage
  /// @param verbose 
  void assign_stages(mockturtle::mph_view<mockturtle::klut_network, NUM_PHASES> & ntk, const phmap::flat_hash_map<unsigned int, unsigned int> & phase_assignment, const bool verbose = false)
  {
    mockturtle::topo_view<klut> ntk_topo ( ntk );

    ntk_topo.foreach_node([&] ( const klut::signal & node ) 
    {
      if ( ntk_topo.is_constant( node ) )
      {
        return;
      }

      if ( ntk_topo.is_pi( node ) )
      {
        auto ct = phase_assignment.count(node);
        uint32_t stage = (ct != 0) ? phase_assignment.at( node ) : 0;
        ntk.set_stage_type(node, stage, AS_GATE);

        DEBUG_PRINT("PI {} placed at ɸ=0 [S=0, φ=0]\n", node);
        return;
      }

      // if (verbose) fmt::print("{} GATE {}:\n", GATE_TYPE.at(node_params.type), node);

      uint8_t node_type = ntk.get_type( node );
      DEBUG_PRINT("{} GATE {}:\n", GATE_TYPE.at(node_type), node);

      auto ct = phase_assignment.count(node);
      if (ct != 0) // there is a precalculated phase assignment
      {
        ntk.set_stage_type(node, phase_assignment.at( node ), node_type);
      }
    });
  }

  /// @brief compares stages of two klut nodes. 
  /// @param a - first klut signal to be compared
  /// @param b - first klut signal to be compared
  /// @param ntk - multiphase view of the network
  /// @return 
  bool phase_ntk_comparison(const klut::signal & a, const klut::signal & b, const mockturtle::mph_view<mockturtle::klut_network, NUM_PHASES> & ntk )
  {
    if (ntk.is_constant(a))
    {
      return true;
    }
    else if (ntk.is_constant(b))
    {
      return false;
    }
    return ntk.get_stage(a) < ntk.get_stage(b);
  }

  /// @brief Explicit insertion of the splitters. Creates a simple chain (only needed for assigning the stages)
  /// @param ntk - multiphase network
  /// @param verbose - debugging flag
  void splitter_ntk_insertion(mockturtle::mph_view<mockturtle::klut_network, NUM_PHASES> & ntk, const bool verbose = false)
  {
    // Lambda function for comparing the phases of two signals in the network.
    auto phase_comp = [&](const klut::signal & a, const klut::signal & b)
    {
      return phase_ntk_comparison(a, b, ntk);
    };

    // Create a view of the network that provides access to fanout information.
    auto ntk_fo = mockturtle::fanout_view<klut>(ntk);

    auto init_size = ntk.size();
    // For each node in the fanout view:
    ntk_fo.foreach_node([&](const klut::signal & node)
    {
      if ( ntk_fo.is_dangling( node ) || ntk_fo.is_constant(node) )
      {
        return;
      }

      // Get the number of fanouts for the current node.
      uint32_t fo_size{ 0u };
      std::vector<klut::signal> fanouts;
      fanouts.reserve(fo_size);
      ntk_fo.foreach_fanout( node, [&]( auto const& fo_node ) {
        if ( !ntk_fo.is_dangling( fo_node ) )
        {
          ++fo_size;
          fanouts.push_back( fo_node );
        }
      } );
      ntk._storage->nodes[node].data[0].h1 = fo_size;

      DEBUG_PRINT("\t[NODE {}] FANOUT SIZE = {}\n", node, fo_size);
      // If the current node is a constant or it has fanout ≤ 1, skip to the next node.
      if ( fo_size <= 1 )
      {
        return;
      }
      fmt::print("Processing node {:>5d} out of {:>5d}\r", node, init_size);

      // Fix the fanout count (bugged fanouts_size()?)
      ntk._storage->nodes[node].data[0].h1 = fanouts.size();

      // Sort the fanouts using the phase comparison function.
      std::sort(fanouts.begin(), fanouts.end(), phase_comp);
      DEBUG_PRINT("\t[NODE {}] SORTED FANOUTS:\n", node);
      if (verbose) 
      {
        printVector(fanouts, 2);
      }

      // Create [fo_size - 1] splitter nodes.
      klut::signal last_spl = node;
      std::vector<klut::signal> splitters;
      splitters.reserve( fanouts.size() - 1 );
      for (auto it = fanouts.begin(); it < fanouts.end() - 1; it++)
      {
        DEBUG_PRINT("\t\t[NODE {}] LAST SPL: {}\n", node, last_spl);

        // Create a new splitter node connected to 'last_spl'.
        DEBUG_PRINT("\t\t[NODE {}] CREATING SPL FOR {}\n", node, *it);
        DEBUG_PRINT("\t\t[NODE {}] LAST_SPL {} FANOUT BEFORE: {}\n", node, last_spl, ntk.fanout_size(last_spl));
        const klut::signal spl = ntk.explicit_buffer(last_spl);
        ntk.set_stage_type(spl, ntk.get_stage(*it), AA_GATE);
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
          if ( ntk.is_dangling( pred_it->data ) )
          {
            continue;
          }

          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PRED {}\n", node, last_spl, pred_it->data);
          if (pred_it->data == node)
          {
            DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] RECORDING OLD PREDS\n", node, last_spl);

            DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS BEFORE\n", node, last_spl);
            pred_it->data = spl;                             // Replace the connection with the splitter.
            DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS AFTER\n", node, last_spl);

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
            
            // splitter has only one fanin, can break the loop
            break;
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
        if ( ntk.is_dangling( pred_it->data ) )
        {
          continue;
        }

        if (pred_it->data == node)
        {
          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] RECORDING OLD PREDS\n", node, last_spl);

          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS BEFORE\n", node, last_spl);
          // for (const auto& entry : ntk._storage->nodes[fanouts.back()].children)  { DEBUG_PRINT("\t\t\t\t{}\n", entry.data); }
          pred_it->data = last_spl;                            // Replace the connection with the last splitter.
          DEBUG_PRINT("\t\t\t[NODE {}, LAST_SPL {}] PREDS AFTER\n", node, last_spl);

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
        }
      }

      // debug printing
      if (verbose)
      {
        auto fo_ntk { mockturtle::fanout_view( ntk ) };
        for (const auto & node : fanouts)
        {
          DEBUG_PRINT("\t\t\t node: {}\n", node);
          
          fo_ntk.foreach_fanin( node, [&](const klut::signal & fi_node)
          {
            if ( fo_ntk.is_dangling( fi_node ) )
            {
              return;
            }

            DEBUG_PRINT("\t\t\t\t fanin : {}\n", fi_node);
          });
          fo_ntk.foreach_fanout( node, [&](const klut::signal & fo_node)
          {
            if ( fo_ntk.is_dangling( fo_node ) )
            {
              return;
            }

            DEBUG_PRINT("\t\t\t\t fanout: {}\n", fo_node);
          });
        }
        for (const auto & node : splitters)
        {
          DEBUG_PRINT("\t\t\t spl : {}\n", node);
          
          fo_ntk.foreach_fanin( node, [&](const klut::signal & fi_node)
          {
            if ( fo_ntk.is_dangling( fi_node ) )
            {
              return;
            }

            DEBUG_PRINT("\t\t\t\t fanin : {}\n", fi_node);
          });
          fo_ntk.foreach_fanout( node, [&](const klut::signal & fo_node)
          {
            if ( fo_ntk.is_dangling( fo_node ) )
            {
              return;
            }

            DEBUG_PRINT("\t\t\t\t fanout: {}\n", fo_node);
          });
        }
      }
      // Ensure that the current node's fan-out count is now 1 (since all other fanouts have been replaced by splitters).
      assert(ntk._storage->nodes[node].data[0].h1 == 1);
    });
  }

  std::vector<klut::signal> get_fanins(const mockturtle::mph_view<mockturtle::klut_network, NUM_PHASES>  & ntk, const klut::signal & fo_node, const phmap::flat_hash_map<klut::signal, std::vector<klut::signal>> & cut_leaves)
  {
    if (ntk.get_type(fo_node) == T1_GATE)
    {
      // get leaves of the 3-cut
      return cut_leaves.at(fo_node);
    }
    else
    {
      std::vector<klut::signal> fanins;
      ntk.foreach_fanin(fo_node, [&](const klut::signal & fi_node)
      {
        fanins.push_back(fi_node);
      });
    }
  }

  /* the version of 'get_fanins' that does not distinguish T1 gates from others */
  std::vector<klut::signal> get_fanins( klut const& ntk, klut::signal const& n )
  {
    std::vector<klut::signal> fanins;
    ntk.foreach_fanin( n, [&]( klut::signal const& ni )
    {
      fanins.push_back( ni );
    } );
  }

  void write_klut_specs( const mockturtle::mph_view<mockturtle::klut_network, NUM_PHASES>  & ntk, const std::string & filename, const bool verbose = false )
  {
    std::ofstream spec_file( filename );

    spec_file << "PI";
    ntk.foreach_pi( [&]( const auto & node ){
      spec_file << "," << node;
    } );
    spec_file << "\n";

    ntk.foreach_gate( [&]( auto const& n ) 
    {
      std::vector<klut::node> n_fanins;
      ntk.foreach_fanin( n, [&n_fanins]( auto const& ni ) 
      {
        n_fanins.push_back( ni );
      } );

      spec_file << fmt::format( "{0},{1},{2}\n", n, ntk.get_type(n), fmt::join( n_fanins, "|" ) );
      // spec_file << fmt::format( "{0},{1},{2}\n", ntk.node_to_index( n ), ntk.get_type(n), fmt::join( n_fanins, "|" ) );
      return true;
    } );
  }

  std::vector<Path<NUM_PHASES>> extract_paths(const mockturtle::mph_view<mockturtle::klut_network, NUM_PHASES> & ntk, bool verbose = false)
  {
    DEBUG_PRINT("\t[i] ENTERED FUNCTION extract_paths\n");
    std::vector<Path<NUM_PHASES>> paths;

    ntk.foreach_node([&](const klut::signal & fo_node)
    {
      DEBUG_PRINT("\t\t[i] PROCESSING NODE {}\n", fo_node);
      if (ntk.is_constant(fo_node) || ntk.is_pi(fo_node))
      {
        DEBUG_PRINT("\t\t\t[NODE {}] the node is CONSTANT/PI\n", fo_node);
        return;
      }
      if ( ntk.get_type(fo_node) == AA_GATE ) 
      {
        DEBUG_PRINT("\t\t\t[NODE {}] the node is AA, skipping\n", fo_node);
        return;
      }

      // at this point, the node should be AS/SA/T1
      DEBUG_PRINT("\t\t[NODE {}] the node is AS/SA/T1, continuing...\n", fo_node);

      ntk.foreach_fanin(fo_node, [&](const klut::signal & fi_node)
      {
        DEBUG_PRINT("\t\t\t[NODE {}] processing fanin {}\n", fo_node, fi_node);
        // Create a path object with only a target
        Path<NUM_PHASES> node_path;
        node_path.targets.emplace( fo_node );

        std::vector<klut::signal> stack { fi_node };
        
        DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}] created stack\n", fo_node, fi_node);
        
        std::set<klut::signal> seen;
        while (!stack.empty())
        {
          DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}] stack contents:\n", fo_node, fi_node);
          if (verbose) 
          {
            printVector(stack, 4);
          }

          const klut::signal & n = stack.back();
          stack.pop_back();

          DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}]: Analyzing node {}\n", fo_node, fi_node, n);

          // A constant does not have any effect on the DFF placement, we can skip it
          if ( ntk.is_constant( n ) )
          {
            continue;
          }        
          const auto n_type = ntk.get_type(n);
          // Found a source of the path, add to sources, do not continue traversal
          if ( ntk.is_pi(n) || n_type == AS_GATE || n_type == SA_GATE || n_type == T1_GATE )
          {
            DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}]: node {} is a source \n", fo_node, fi_node, n);
            node_path.sources.emplace( n );
          }
          // Found AA gate, add to internal nodes, add parents for further traversal
          else if ( n_type == AA_GATE )
          {
            DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}]: node is INTERNAL adding fanins \n", fo_node, fi_node, n);
            node_path.internals.emplace( n );

            ntk.foreach_fanin(n, [&](const klut::signal & sig){
              stack.push_back( sig );
            });
          }
          else
          {
            DEBUG_PRINT("\t\t\t[NODE {}][FANIN {}]: Signal {}: {} is not recognized \n", fo_node, fi_node, n, GATE_TYPE.at( n_type ));
            throw "Unsupported case";
          }
          seen.emplace( n );
        }

        // Identify overlapping paths
        std::vector<size_t> to_merge;
        for (size_t i = 0u; i < paths.size(); ++i)
        {
          Path<NUM_PHASES> & known_paths = paths[i];
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

  void merge_overlapping_paths( std::vector<Path<NUM_PHASES>>& paths, Path<NUM_PHASES>& current_path )
  {
    /* identify overlapping paths that can be merged */
    std::vector<size_t> index_to_merge;
    for ( auto i{ 0u }; i < paths.size(); ++i )
    {
      auto known_path = paths[i];
      if ( haveCommonElements( known_path.sources, current_path.sources ) )
      {
        index_to_merge.push_back( i );
      }
    }

    /* merge identified overlapping paths */
    for ( auto it{ index_to_merge.rbegin() }; it != index_to_merge.rend(); ++it )
    {
      auto idx = *it;
      current_path.absorb( paths[idx] );
      paths.erase( paths.begin() + idx );
    }
  }

  void buildPath(const mockturtle::mph_view<mockturtle::klut_network, NUM_PHASES> & ntk, Path<NUM_PHASES> & node_path, std::vector<klut::signal> stack, const bool verbose = false)
  {
    std::set<klut::signal> seen;
    while (!stack.empty())
    {
      DEBUG_PRINT("\t\t\tStack contents:\n");
      if (verbose) 
      {
        printVector(stack, 4);
      }

      const klut::signal & n = stack.back();
      stack.pop_back();

      DEBUG_PRINT("\t\t\tAnalyzing node {}\n", n);

      // A constant does not have any effect on the DFF placement, we can skip it
      if ( ntk.is_constant( n ) )
      {
        continue;
      }        
      const auto n_type = ntk.get_type(n);
      // Found a source of the path, add to sources, do not continue traversal
      if ( ntk.is_pi(n) || n_type == AS_GATE || n_type == SA_GATE || n_type == T1_GATE )
      {
        DEBUG_PRINT("\t\t\tnode {} is a source \n", n);
        node_path.sources.emplace( n );
      }
      // Found AA gate, add to internal nodes, add parents for further traversal
      else if ( n_type == AA_GATE )
      {
        DEBUG_PRINT("\t\t\tnode is INTERNAL adding fanins \n", n);
        node_path.internals.emplace( n );

        ntk.foreach_fanin(n, [&](const klut::signal & sig){
          stack.push_back( sig );
        });
      }
      else
      {
        DEBUG_PRINT("\t\t\tSignal {}: {} is not recognized \n", n, GATE_TYPE.at( n_type ));
        throw "Unsupported case";
      }
      seen.emplace( n );
    }
  }

  /// @brief Create binary variables for DFF placement in a given path
  /// @param path - a path object to insert DFFs into
  /// @param NR - unordered_map of NtkNode objects 
  /// @param verbose - prints debug messages if set to *true*
  /// @return 
  std::tuple<DFF_registry, uint64_t, std::vector<uint64_t>> dff_vars_single_paths(const Path<NUM_PHASES> & path, const mockturtle::mph_view<mockturtle::klut_network, NUM_PHASES> & ntk, bool verbose = false)
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
      // fixing the C++17 bug with structured bindings
      auto _node_tuple = stack.back();
      const klut::signal fo_node          = std::get<0>(_node_tuple);
      const uint64_t earliest_child_hash  = std::get<1>(_node_tuple);
      stack.pop_back();
      const auto [fo_stage, fo_type] = ntk.get_stage_type(fo_node);

      // AS gates and T1 gate are clocked, so one needs to start one stage earlier
      uint32_t latest_sigma = fo_stage - (fo_type == AS_GATE || fo_type == T1_GATE);
      DEBUG_PRINT("[DFF] Analyzing child: {}({})[{}]\n", GATE_TYPE.at(fo_type), fo_node, (int)fo_stage);

      ntk.foreach_fanin(fo_node, [&](const klut::signal & fi_node)
      {
        const auto [fi_stage, fi_type] = ntk.get_stage_type(fi_node);
        /* if fi_node is an AA gate, it is allowed to put a DFF at the phase that fi_node is assigned to */
        uint32_t earliest_sigma = fi_stage + (fi_type != AA_GATE);

        DEBUG_PRINT("\t[DFF] Analyzing parent: {}({})[{}]\n", GATE_TYPE.at(fi_type), fi_node, (int)fi_stage);

        // check if the chain is straight - #DFF is just floor(delta-phase), no need to create the dff vars
        if (fo_type != AA_GATE && fi_type != AA_GATE)
        {
          // special case when an AS gate feeds directly into SA gate
          if (fo_stage == fi_stage)
          {
            DEBUG_PRINT("\t[DFF] Straight chain: AS{} -> SA{}\n", fi_node, fo_node);
            // do nothing, no additional DFFs needed
            assert(fo_type == SA_GATE && fi_type == AS_GATE && ntk.fanout_size(fi_node) == 1);
          }
          else
          {
            DEBUG_PRINT("\t[DFF] Straight chain: {}[{}] -> {}[{}]\n", GATE_TYPE.at(fi_type), (int)fi_stage, GATE_TYPE.at(fo_type), (int)fo_stage);
            // straight chain, just floor the difference!
            precalc_ndff += (fo_stage - fi_stage - 1)/NUM_PHASES + (fo_type == SA_GATE); //extra DFF before SA gate
          }
          return;
        }

        DEBUG_PRINT("\t[DFF] Non-straight chain: {}[{}] -> {}[{}]\n", GATE_TYPE.at(fi_type), (int)fi_stage, GATE_TYPE.at(fo_type), (int)fo_stage);
        std::vector<uint64_t> out_hashes;
        DEBUG_PRINT("\tAdding new DFFs [reg size = {}]\n", DFF_REG.variables.size());

        for (stage_t stage = earliest_sigma; stage <= latest_sigma; ++stage)
        {
          uint64_t new_hash = DFF_REG.add(fi_node, fo_node, stage);
          out_hashes.push_back(new_hash);
          DEBUG_PRINT("\tAdded new DFFs at phase {} [reg size = {}]\n", stage, DFF_REG.variables.size());
        }
        DEBUG_PRINT("\tConnecting new DFFs\n");
        for (auto i = 1u; i < out_hashes.size(); ++i)
        {
          DFF_var & dff = DFF_REG.at( out_hashes[i] );
          dff.parent_hashes.emplace(out_hashes[i-1]);
        }
        if (fo_type == SA_GATE)
        {
          // ensure that the SA gate is placed at least one stage after the fanin 
          //(to leave space for a DFF)
          assert( !out_hashes.empty() );
          // The last DFF in the chain is required
          required_SA_DFFs.push_back(out_hashes.back());
        }
        // if there are DFFs, the earliest_hash is the first hash in the chain
        // otherwise, it is the earliest_hash of the previous chain
        uint64_t earliest_hash = (out_hashes.empty()) ? earliest_child_hash : out_hashes.front();
        // if the node is internal, connect with the fanout phase
        if (fo_type == AA_GATE && !out_hashes.empty() && earliest_hash != 0 && earliest_child_hash != 0)
        {
          DFF_var & child_dff = DFF_REG.at( earliest_child_hash );
          DEBUG_PRINT("\tPrior node is {}[{}]\n", child_dff.str(), (int)child_dff.stage); 
          // assert(child_dff.fanin == fo_id);
          child_dff.parent_hashes.emplace( out_hashes.back() );
        }
        if (fi_type == AA_GATE)
        {
          stack.emplace_back( fi_node, earliest_hash );
          DEBUG_PRINT("\tEmplacing {}({})[{}], {}\n", GATE_TYPE.at(fi_type), fi_node, (int)fi_stage, (earliest_hash!=0)?DFF_REG.at(earliest_hash).str():"");
        }
      });
    }
    return std::make_tuple(DFF_REG, precalc_ndff, required_SA_DFFs);
  }

  void process_path(std::vector<int> & local_num_dff, const int idx, const Path<NUM_PHASES>& path, const mockturtle::mph_view<mockturtle::klut_network, NUM_PHASES> & network, const std::string benchmark)
  {
    bool verbose = true;
    DEBUG_PRINT("\tAnalyzing the path\n");

    // *** Create binary variables
    auto [DFF_REG, precalc_ndff, required_SA_DFFs] = dff_vars_single_paths(path, network, NUM_PHASES);
    local_num_dff[idx] += precalc_ndff;
    DEBUG_PRINT("\t\t\t\t[i]: Precalculated {} DFFs\n", precalc_ndff);

    // *** Generate constraints
    std::vector<Snake<NUM_PHASES>> snakes = sectional_snake<NUM_PHASES>(path, network, DFF_REG, verbose);
    DEBUG_PRINT("\tCreated {} snakes\n", snakes.size());

    // *** If there's anything that needs optimization
    if (!snakes.empty())
    {
      std::string cfg_file = fmt::format("ilp_configs/{}_cfgNR_mt_{}.csv", benchmark, idx);
      write_snakes<NUM_PHASES>(snakes, DFF_REG, required_SA_DFFs, cfg_file, verbose);
      auto num_dff = cpsat_ortools(cfg_file);
      // DEBUG_PRINT("OR Tools optimized to {} DFF\n", num_dff);
      local_num_dff[idx] += num_dff;
      DEBUG_PRINT("\t\t\t\t[i] total CPSAT #DFF = {}\n", local_num_dff[idx]);
    }
  }
private:
  mockturtle::binding_view<mockturtle::klut_network> src, 
  phmap::flat_hash_map<uint64_t, Node> nodemap,
  phmap::flat_hash_map<std::string, LibEntry> entries,
  const std::array<int, 12> COSTS_MAP 
};


template<uint8_t NUM_PHASES>
mockturtle::mph_view<mockturtle::klut_network, NUM_PHASES> multiphase_mapping( 
  mockturtle::binding_view<mockturtle::klut_network> src, 
  phmap::flat_hash_map<ULL, Node> nodemap,
  phmap::flat_hash_map<std::string, LibEntry> entries,
  const std::array<int, 12> COSTS_MAP, 
  const bool verbose = false
  )
{

}
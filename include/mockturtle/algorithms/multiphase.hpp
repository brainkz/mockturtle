
#pragma once

#include <mockturtle/io/auxiliary_genlib.hpp>
#include <mockturtle/algorithms/nodes.hpp>
#include <mockturtle/utils/misc.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/rsfq_view.hpp>
#include <mockturtle/views/fanout_view.hpp>


typedef uint32_t glob_phase_t;

typedef mockturtle::klut_network klut;
typedef mockturtle::xag_network   xag;
typedef mockturtle::xmg_network   xmg;
typedef mockturtle::mig_network   mig;
typedef mockturtle::aig_network   aig;
typedef uint64_t node_t;

enum GateType : uint8_t {
    PI_GATE = 0u,
    AA_GATE = 1u,
    AS_GATE = 2u,
    SA_GATE = 3u,
    T1_GATE = 4u
};
const std::vector<std::string> GATE_TYPE { "PI", "AA", "AS", "SA", "T1" }; //, "PO"

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

std::tuple<klut, int64_t> decompose_to_klut(mockturtle::binding_view<klut> src, phmap::flat_hash_map<ULL, Node> nodemap, phmap::flat_hash_map<std::string, LibEntry> entries, const std::array<int, 12> COSTS_MAP, bool verbose = false)
{
  phmap::flat_hash_map<klut::signal, klut::signal> src2tgt;
  int64_t area = 0;
  
  klut tgt;
  src.foreach_pi( [&]( const klut::signal & src_pi ) 
  {
    klut::signal tgt_pi = tgt.create_pi();
    tgt.set_value(tgt_pi, NodeData(0, PI_GATE).value);
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
    DEBUG_PRINT("Processing node {0} ({1} out of {2})\r", src_node, ++ctr, num_node_src);
    if (ctr == num_node_src)
    {
      DEBUG_PRINT("\n");
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

    // Get the topological order of internal nodes (i.e., nodes inside the cell)
    std::vector<ULL> topo_order;
    root.topo_sort(nodemap, topo_order, true);

    std::vector<klut::signal> tgt_fanins;
    src.foreach_fanin(src_node, [&](const auto & src_fanin)
    {
      tgt_fanins.push_back(src2tgt.at(src_fanin));
      // fmt::print("\tRecorded fanin of src_node {0}:\n\t\tsrc_fanin: {1}\n\t\ttgt_fanin: {2}\n", src_node, src_fanin, tgt_fanins.back().data);
    } );

    phmap::flat_hash_map<ULL, klut::signal> node2tgt;

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
        //  nothing else to do, no new constraints, no contribution to objective function
      }
      else if (node.last_func == fNOT)
      {
        klut::signal parent_tgt = node2tgt.at(node.parent_hashes.front());
        klut::signal tgt_sig = tgt.create_not( parent_tgt );
        tgt.set_value(tgt_sig, NodeData(0, AS_GATE).value);
        node2tgt.emplace(hash, tgt_sig);
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
    tgt.foreach_fanin(xor_sig, [&] (const klut::signal & fanin)
    {
      // if the OR gate is before XOR, do not replace it with the CB
      auto it = std::find(OR_replacement_candidates.begin(), OR_replacement_candidates.end(), fanin);
      if (it != OR_replacement_candidates.end())
      {
        OR_replacement_candidates.erase(it);
      }
    });
  }

  for (klut::signal const & sig : OR_replacement_candidates)
  {
    tgt.set_value(sig, NodeData(0, AA_GATE).value);
    area += COSTS_MAP[fCB];
    area -= COSTS_MAP[fOR];
    DEBUG_PRINT("REPLACED OR WITH CB = {}\n", area);
  }

  src.foreach_po([&](auto const & src_po)
  {
    tgt.create_po(src2tgt.at(src_po));
  } );
  
  return std::make_tuple(tgt, area);
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

/// @brief Assigns stages to nodes based on stage assignment. If the assignment is not found, assigns a stage greedily
/// @param ntk 
/// @param n_phases 
/// @param phase_assignment 
/// @param verbose 
void greedy_ntk_assign(const klut & ntk, const uint8_t n_phases, const phmap::flat_hash_map<unsigned int, unsigned int> & phase_assignment, const bool verbose = false)
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
      uint32_t sigma = (ct != 0) ? phase_assignment.at( node ) : 0;
      NodeData node_data {sigma, AS_GATE};
      ntk.set_value(node, node_data.value);

      DEBUG_PRINT("PI {} placed at ɸ=0 [S=0, φ=0]\n", node);
      return;
    }

    // if (verbose) fmt::print("{} GATE {}:\n", GATE_TYPE.at(node_params.type), node);

    uint8_t node_type = NodeData( ntk.value( node ) ).type;
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

/// @brief compares stages of two klut nodes. 
/// @param a - first klut signal to be compared
/// @param b - first klut signal to be compared
/// @param ntk - klut network
/// @return 
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


// Function to insert splitter nodes in a KLUT network.
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


/// @brief Structure representing the potential DFF location uniquely defined by fanin, fanout and stage
struct DFF_var 
{
  klut::signal fanin;
  klut::signal fanout;
  uint32_t sigma;
  std::unordered_set<uint64_t> parent_hashes;

  DFF_var(klut::signal _fanin, klut::signal _fanout, uint32_t _sigma, std::unordered_set<uint64_t> _parent_hashes = {})
      : fanin(_fanin), fanout(_fanout), sigma(_sigma), parent_hashes(_parent_hashes) {}

  DFF_var(const DFF_var& other)
      : fanin(other.fanin), fanout(other.fanout), sigma(other.sigma), parent_hashes(other.parent_hashes) {}

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
  phmap::flat_hash_map<uint64_t, DFF_var> variables;

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
    DFF_var temp { _fanin, _fanout, _phase, _parent_hashes };
    variables.emplace(_hash, temp);
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


std::vector<klut::signal> get_fanins(const klut & ntk, const klut::signal & fo_node, const phmap::flat_hash_map<klut::signal, std::vector<klut::signal>> & cut_leaves)
{
  if (NodeData(ntk.value(fo_node)).type == T1_GATE)
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

void write_klut_specs(const klut & ntk, const std::string & filename)
{
  auto ntk_fo = mockturtle::fanout_view<klut, false>( ntk );
  std::ofstream spec_file (filename);
  ntk.foreach_node( [&]( const klut::signal & sig ) 
  {
    std::vector<klut::node> fanouts;
    ntk_fo.foreach_fanout(sig, [&](const klut::node fo) { fanouts.push_back(fo);  });
    std::vector<klut::node> fanins;
    ntk_fo.foreach_fanin (sig, [&](const klut::node fi) {  fanins.push_back(fi);  });

    // TODO: Infer from the ntk, not Primitive object
    const NodeData & node_data {ntk.value( sig )};
    uint8_t prim_type = ( ntk.fanin_size( sig ) == 0 ) ? PI_GATE : node_data.type;
    spec_file << fmt::format("{0},{1},{2},{3}\n", sig, prim_type, fmt::join(fanins, "|"), fmt::join(fanouts, "|"));
  });
}

void write_klut_specs_T1(const klut & ntk, 
const phmap::flat_hash_map<klut::signal, std::vector<klut::signal>> & cut_leaves, 
const phmap::flat_hash_map<klut::signal, klut::signal> & representatives, 
const std::string & filename)
{
  auto ntk_fo = mockturtle::fanout_view<klut, false>( ntk );
  std::ofstream spec_file (filename);
  ntk.foreach_node( [&]( const klut::signal & sig ) 
  {
    if (NodeData(ntk.value(sig)).type == T1_GATE)
    {
      if (representatives.at(sig) != sig)
      {
        return;
      }
      // TODO: also need to record the bound nodes
      // need to ensure the sigmas given to T1 outputs are equal
    }
    std::vector<klut::node> fanouts;
    ntk_fo.foreach_fanout(sig, [&](const klut::node fo) { fanouts.push_back(fo);  });
    std::vector<klut::node> fanins = get_fanins(ntk, sig, cut_leaves);
    // ntk_fo.foreach_fanin (sig, [&](const klut::node fi) {  fanins.push_back(fi);  });

    const NodeData & node_data {ntk.value( sig )};
    uint8_t prim_type = ( ntk.fanin_size( sig ) == 0 ) ? PI_GATE : node_data.type;
    spec_file << fmt::format("{0},{1},{2}\n", sig, prim_type, fmt::join(fanins, "|"));
  });
}


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
    if ( fo_node_data.type == AA_GATE ) 
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
        if ( ntk.is_pi(n) || n_data.type == AS_GATE || n_data.type == SA_GATE || n_data.type == T1_GATE )
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

void buildPath(const klut & ntk, Path & node_path, std::vector<klut::signal> stack, const phmap::flat_hash_map<klut::signal, std::vector<klut::signal>> & cut_leaves, const bool verbose = false)
{
  std::set<klut::signal> seen;
  while (!stack.empty())
  {
    DEBUG_PRINT("\t\t\tStack contents:\n");
    printVector(stack, 4);

    const klut::signal & n = stack.back();
    stack.pop_back();

    DEBUG_PRINT("\t\t\tAnalyzing node {}\n", n);

    // A constant does not have any effect on the DFF placement, we can skip it
    if ( ntk.is_constant( n ) )
    {
      continue;
    }        
    const NodeData n_data { ntk.value(n) };
    // Found a source of the path, add to sources, do not continue traversal
    if ( ntk.is_pi(n) || n_data.type == AS_GATE || n_data.type == SA_GATE || n_data.type == T1_GATE )
    {
      DEBUG_PRINT("\t\t\tnode {} is a source \n", n);
      node_path.sources.emplace( n );
    }
    // Found AA gate, add to internal nodes, add parents for further traversal
    else if ( n_data.type == AA_GATE )
    {
      DEBUG_PRINT("\t\t\tnode is INTERNAL adding fanins \n", n);
      node_path.internals.emplace( n );

      ntk.foreach_fanin(n, [&](const klut::signal & sig)
      {
        stack.push_back( sig );
      });
    }
    else
    {
      DEBUG_PRINT("\t\t\tSignal {}: {} is not recognized \n", n, GATE_TYPE.at( n_data.type ));
      throw "Unsupported case";
    }
    seen.emplace( n );
  }
}

/// @brief same as extract_paths but with the support of T1 cells. Uses helper functions
/// @param ntk 
/// @param verbose 
/// @return 
std::vector<Path> extract_paths_t1(const klut & ntk, const phmap::flat_hash_map<klut::signal, std::vector<klut::signal>> & cut_leaves, const phmap::flat_hash_map<klut::signal, klut::signal> & representatives, bool verbose = false)
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
    if ( fo_node_data.type == AA_GATE ) 
    {
      DEBUG_PRINT("\t\t\t[NODE {}] the node is AA, skipping\n", fo_node);
      return;
    }

    // at this point, the node should be AS/SA/T1
    DEBUG_PRINT("\t\t[NODE {}] the node is AS/SA/T1, continuing...\n", fo_node);

    // If the gate is T1:
    //   - all leaves of the cut belong to the same path due to the transverse constraint (unlike AS/SA)
    //   - only a single node among T1 outputs needs to be considered. Other outputs would represent the same path
    if (fo_node_data.type == T1_GATE)
    {
      // Check if the node is the main representative of the T1 cell
      if (representatives.at(fo_node) != fo_node)
      {
        // if not, skip, the T1 cell was (or will be) considered by another node
        return;
      }

      Path node_path;
      node_path.targets.emplace( fo_node );
      std::vector<klut::signal> stack { get_fanins(ntk, fo_node, cut_leaves) };
      buildPath(ntk, node_path, stack, cut_leaves, true);
      paths.push_back(node_path);
    }
    else if (fo_node_data.type == AS_GATE || fo_node_data.type == SA_GATE)
    {
      Path node_path;
      node_path.targets.emplace( fo_node );
      ntk.foreach_fanin( fo_node, [&] (const klut::signal & fi_node) 
      {
        std::vector<klut::signal> stack { fi_node };
        buildPath(ntk, node_path, stack, cut_leaves, true);
        paths.push_back(node_path);
      });
    }
    else
    {
      DEBUG_PRINT("\t\t\tSignal {}: {} is not recognized \n", fo_node, GATE_TYPE.at( fo_node_data.type ));
      throw "Unsupported case";
    }
  });

  // TODO: NEED FUNCTION TO MERGE THE PATHS
  // // Identify overlapping paths
  // std::vector<size_t> to_merge;
  // for (size_t i = 0u; i < paths.size(); ++i)
  // {
  //   Path & known_paths = paths[i];
  //   // merge if there are sources in common
  //   if( haveCommonElements( known_paths.sources, node_path.sources) )
  //   {
  //     to_merge.push_back(i);
  //   }
  // }

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
    // fixing the C++17 bug with structured bindings
    auto _node_tuple = stack.back();
    const klut::signal fo_node          = std::get<0>(_node_tuple);
    const uint64_t earliest_child_hash  = std::get<1>(_node_tuple);
    stack.pop_back();
    NodeData fo_data { ntk.value( fo_node ) };

    // AS gates and T1 gate are clocked, so one needs to start one stage earlier
    uint32_t latest_sigma = fo_data.sigma - (fo_data.type == AS_GATE);
    DEBUG_PRINT("[DFF] Analyzing child: {}({})[{}]\n", GATE_TYPE.at(fo_data.type), fo_node, (int)fo_data.sigma);


    ntk.foreach_fanin(fo_node, [&](const klut::signal & fi_node)
    {
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



/// @brief Create binary variables for DFF placement in a given path (with support of T1 cells)
/// @param path - a path object to insert DFFs into
/// @param NR - unordered_map of NtkNode objects 
/// @param n_phases - # of phases
/// @param verbose - prints debug messages if set to *true*
/// @return 
std::tuple<DFF_registry, uint64_t, std::vector<uint64_t>> dff_vars_single_paths_T1(const Path & path, const klut & ntk, const uint8_t n_phases, bool verbose = false)
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
    // fixing the C++17 bug with structured bindings
    auto _node_tuple = stack.back();
    const klut::signal fo_node          = std::get<0>(_node_tuple);
    const uint64_t earliest_child_hash  = std::get<1>(_node_tuple);
    stack.pop_back();
    NodeData fo_data { ntk.value( fo_node ) };

    // AS gates and T1 gate are clocked, so one needs to start one stage earlier
    uint32_t latest_sigma = fo_data.sigma - (fo_data.type == AS_GATE || fo_data.type == T1_GATE);
    DEBUG_PRINT("[DFF] Analyzing child: {}({})[{}]\n", GATE_TYPE.at(fo_data.type), fo_node, (int)fo_data.sigma);

    ntk.foreach_fanin(fo_node, [&](const klut::signal & fi_node)
    {
      NodeData fi_data { ntk.value( fi_node ) };
      uint32_t earliest_sigma = fi_data.sigma + (fi_data.type != AA_GATE);

      DEBUG_PRINT("\t[DFF] Analyzing parent: {}({})[{}]\n", GATE_TYPE.at(fi_data.type), fi_node, (int)fi_data.sigma);

      // check if the chain is straight - #DFF is just floor(delta-phase), no need to create the DFF vars
      if (  (fo_data.type == AS_GATE && fo_data.type == SA_GATE) // the gate should be AS or SA
            && fi_data.type != AA_GATE )                         // the fanin to the gate should be AS/SA/T1
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
      
      // DFF locations from the earliest to the latest
      std::vector<uint64_t> dff_hashes;
      DEBUG_PRINT("\tAdding new DFFs [reg size = {}]\n", DFF_REG.variables.size());

      for (glob_phase_t sigma = earliest_sigma; sigma <= latest_sigma; ++sigma)
      {
        uint64_t new_hash = DFF_REG.add(fi_node, fo_node, sigma);
        dff_hashes.push_back(new_hash);
        DEBUG_PRINT("\tAdded new DFFs at phase {} [reg size = {}]\n", sigma, DFF_REG.variables.size());
      }
      DEBUG_PRINT("\tConnecting new DFFs\n");
      for (auto i = 1u; i < dff_hashes.size(); ++i)
      {
        DFF_var & dff = DFF_REG.at( dff_hashes[i] );
        dff.parent_hashes.emplace(dff_hashes[i-1]);
      }

      // TODO : make sure that the inverted fanins of the T1 cell are placed at least one stage away from the fanin
      // ensure that the SA gate is placed at least one stage after the fanin 
      if (fo_data.type == SA_GATE)
      {
        //(there has to be at least one DFF)
        assert( !dff_hashes.empty() );
        // The last DFF in the chain is required
        required_SA_DFFs.push_back(dff_hashes.back());
      }
      // if there are DFFs, the earliest_hash is the first hash in the chain
      // otherwise, it is the earliest_hash of the previous chain
      uint64_t earliest_hash = (dff_hashes.empty()) ? earliest_child_hash : dff_hashes.front();
      // if the node is internal, connect with the fanout phase
      if (fo_data.type == AA_GATE && !dff_hashes.empty() && earliest_hash != 0 && earliest_child_hash != 0)
      {
        DFF_var & child_dff = DFF_REG.at( earliest_child_hash );
        DEBUG_PRINT("\tPrior node is {}[{}]\n", child_dff.str(), (int)child_dff.sigma); 
        // assert(child_dff.fanin == fo_id);
        child_dff.parent_hashes.emplace( dff_hashes.back() );
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
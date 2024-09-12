from collections import namedtuple
from typing import Dict, FrozenSet, List, Tuple, Union

import networkx as nx
from ortools.sat.python import cp_model as cp

edgeVars = namedtuple("edgeVars", ["En", "First", "Last", "CQ", "Setup", "Hold"])


AS = ("PI", "DFF", "NOT", "XOR")
SA = ("AND", "OR")
AA = ("MRG", "SPL")

MAX_SIGMA: int = 1000
NUM_PHASES: int = 4
T_CLK: int = 100000
T_PHASE = T_CLK // NUM_PHASES

# time unit is 1e-14 s
GATES = {
    "PI": {"Setup": 440, "Hold": 390, "CQ": 790},
    "XOR": {"Setup": 690, "Hold": 550, "CQ": 720},
    "NOT": {"Setup": 440, "Hold": 690, "CQ": 930},
    "DFF": {"Setup": 440, "Hold": 390, "CQ": 790},
    "AND": {"Delay": 570},
    "OR": {"Delay": 570},
    "MRG": {"Delay": 570},
    "SPL": {"Delay": 660},
}

# Helper function to perform DFS and collect valid paths


def dfs_collect_paths(G, node, path, result):
    attr = G.nodes[node]["attr"]
    print(f"Extending path {path}")
    print(f"\tNode: {node}")
    print(f"\tType: {attr}")

    # Check if this node is an end node and not the start node
    if node != path[0] and attr in AS:
        final_path = tuple(path) + (node,)
        result[frozenset(final_path)] = [final_path]
        # result.append(final_path)
        return

    # Only traverse further if the current node is AA
    if attr not in AS:
        for neighbor in G.successors(node):
            dfs_collect_paths(G, neighbor, path + [node], result)


def get_paths(G):
    # Find all nodes with attribute (1) or (2) as start/end points
    start_end_nodes = [n for n, attr in G.nodes(data="attr") if attr in AS]
    print(f"AS nodes: {start_end_nodes}")

    # Store all valid paths
    all_paths: Dict[FrozenSet, List[Tuple]] = {}

    # Perform DFS from each start node
    for start_node in start_end_nodes:
        print(f"Analyzing start node {start_node}")
        for neighbor in G.successors(start_node):
            dfs_collect_paths(G, neighbor, [start_node], all_paths)
            print(all_paths)

    final_paths: Dict[FrozenSet, List[Tuple]] = {}

    # Merge intersecting paths here
    while all_paths:
        nodes, paths = all_paths.popitem()
        print(f"Popped: nodes = {nodes}, paths = {paths}")

        for other_nodes, other_paths in all_paths.items():
            print(f"Checking intersection with other_nodes = {other_nodes}")

            # Check for intersection
            if not nodes.isdisjoint(other_nodes):
                print(f"Intersection found with {other_nodes}, merging")
                new_nodes = frozenset(nodes.union(other_nodes))
                all_paths[new_nodes] = paths + other_paths
                del all_paths[other_nodes]
                print(f"Merged nodes = {new_nodes}, merged paths = {all_paths[new_nodes]}")
                break
        else:
            # The path does not intersect with anything
            print(f"No intersection found, adding to final paths: {nodes} -> {paths}")
            final_paths[nodes] = paths

    return final_paths


def add_FLEn(model: cp.CpModel, u: int, v: int, type_u: str, type_v: str):
    First = model.NewIntVar(0, MAX_SIGMA, f"First_Phase_{u}_{v}")
    Last = model.NewIntVar(0, MAX_SIGMA, f"Last_Phase_{u}_{v}")

    if (type_u in AA) or (type_v in AA) or (type_u in AS and type_v in SA):
        model.Add(First <= Last)
    else:
        model.Add(First < Last)

    if (type_u in AA) and (type_v in AA):
        En = model.NewBoolVar(f"Has_DFF_{u}_{v}")
    else:
        En = model.NewConstant(1)
    return First, Last, En


def add_nDFF(model: cp.CpModel, u: int, v: int, En: cp.IntVar, type_u: str, type_v: str, hold_safe: bool):
    NUM_PHASES_EFF = NUM_PHASES - hold_safe

    Factor = (type_u in AA) + (type_v in (AA + SA))
    numerator = model.NewIntVar(0, MAX_SIGMA, f'NUM_{u}_{v}')
    model.Add(numerator == (Last - First - 1 + Factor * NUM_PHASES_EFF))

    nDFF = model.NewIntVar(0, MAX_SIGMA, f"nDFF_{u}_{v}")
    # The edge is completely unclocked
    if (type_u in AA) and (type_v in AA):
        # Only enabled if En is true, i.e., if there is a DFF along the edge
        model.AddDivisionEquality(nDFF, numerator, NUM_PHASES_EFF).OnlyEnforceIf(En)
        # zero otherwise
        model.Add(nDFF == 0).OnlyEnforceIf(En.Not())
    else:
        # No conditionals needed
        model.AddDivisionEquality(nDFF, numerator, NUM_PHASES_EFF)
    return nDFF


def add_CQ(model: cp.CpModel, u: int, v: int, nDFF: cp.IntVar, type_u: str, pred_types: Tuple[str] = tuple()):
    # if ever needed, the CQ time is DFF CQ time
    if type_u in AA + ("DFF",):
        CQ = model.NewConstant(GATES["DFF"]["CQ"])
    elif type_u in AS:
        viable_CQ = [GATES[type_u]["CQ"], GATES["DFF"]["CQ"]]
        CQ = model.NewIntVarFromDomain(cp.Domain.FromValues(viable_CQ), f'CQ_{u}_{v}')
        model.Add(CQ == GATES[type_u]["CQ"]).OnlyEnforceIf(nDFF == 0)
        model.Add(CQ == GATES["DFF"]["CQ"]).OnlyEnforceIf(nDFF > 0)
    elif type_u in SA:
        if not pred_types:
            CQ = model.NewConstant(GATES["DFF"]["CQ"])
        else:
            pred_CQ = max(GATES[pred_type]["CQ"] for pred_type in pred_types)
            viable_CQ = [pred_CQ + GATES[type_u]["Delay"], GATES["DFF"]["CQ"]]
            CQ = model.NewIntVarFromDomain(cp.Domain.FromValues(viable_CQ), f'CQ_{u}_{v}')
            model.Add(CQ == GATES[type_u]["CQ"]).OnlyEnforceIf(nDFF == 0)
            model.Add(CQ == GATES["DFF"]["CQ"]).OnlyEnforceIf(nDFF > 0)
    return CQ


def add_SetupHold(model: cp.CpModel, u: int, v: int, nDFF: cp.IntVar, type_u: str, type_v: str):
    if type_u in AA:
        # if ever needed, the setup time is DFF setup time
        if type_v in AA + SA + ("DFF",):
            Setup = model.NewConstant(GATES["DFF"]["Setup"])
            Hold = model.NewConstant(GATES["DFF"]["Hold"])
        elif type_v in AS:
            viable_Setup = [GATES[type_v]["Setup"], GATES["DFF"]["Setup"]]
            Setup = model.NewIntVarFromDomain(cp.Domain.FromValues(viable_Setup), f'Setup_{u}_{v}')
            model.Add(Setup == GATES[type_v]["Setup"]).OnlyEnforceIf(nDFF == 0)
            model.Add(Setup == GATES["DFF"]["Setup"]).OnlyEnforceIf(nDFF > 0)

            viable_Hold = [GATES[type_v]["Hold"], GATES["DFF"]["Hold"]]
            Hold = model.NewIntVarFromDomain(cp.Domain.FromValues(viable_Hold), f'Hold_{u}_{v}')
            model.Add(Hold == GATES[type_v]["Hold"]).OnlyEnforceIf(nDFF == 0)
            model.Add(Hold == GATES["DFF"]["Hold"]).OnlyEnforceIf(nDFF > 0)

        return Setup, Hold
    else:
        # Setup time is not needed for edges starting with AS gates
        # Return some dummy values
        return model.NewConstant(0), model.NewConstant(0)


hold_safe_short_paths = True

if __name__ == "__main__":

    # Model.
    model = cp.CpModel()

    # Initialize a directed graph
    G = nx.DiGraph()

    # Parse the CSV file and build the graph
    with open("ilp_config.csv", "r") as f:
        for line in f:
            gate_id_str, func, fanins_str, attr_str = line.split(",")
            gate_id = int(gate_id_str)
            fanins = tuple(int(q) for q in fanins_str.split("|"))
            attr = int(attr_str)

            # Add the gate node to the graph, with its attribute
            G.add_node(gate_id, attr=attr, func=func)

            # Add edges from fanins to the current gate
            for fanin in fanins:
                G.add_edge(fanin, gate_id)

    independent_paths = get_paths(G).items()
    cost_fun = []

    edge_vars: Dict[Tuple[int, int], List[Union[cp.IntVar]]] = {}

    # For each edge in a path, create appropriate variables and constraints
    for nodes, threads in independent_paths.items():
        edge_vars: Dict[Tuple[int, int], Dict[str, cp.IntVar]] = {}
        for thread in threads:
            path_threads = zip(thread[:-1], thread[1:])
            for e in path_threads:
                if e in edge_vars:
                    continue

                u, v = e
                type_u = G.nodes[u]["attr"]
                type_v = G.nodes[v]["attr"]
                # Variables associated with this edge

                First, Last, En = add_FLEn(model, type_u, type_v)
                nDFF = add_nDFF(model, u, v, En, type_u, type_v, hold_safe=hold_safe_short_paths)

                CQ = add_CQ(model, u, v, nDFF, type_u)
                Setup, Hold = add_SetupHold(model, u, v, nDFF, type_u, type_v)

                cost_fun.append(nDFF)

                edge_vars[e] = edgeVars(En, First, Last, CQ, Setup, Hold)

            for i, e1 in enumerate(path_threads[:-1]):
                En_i, First_i, Last_i, CQ_i, Setup_i, Hold_i = edge_vars[e1]
                for j, e2 in enumerate(path_threads[i + 1:]):
                    En_j, First_j, Last_j, CQ_j, Setup_j, Hold_j = edge_vars[e2]
                    all_En = [edge_vars[e].En.Not() for e in path_threads[i + 1: j]] + [En_i, En_j]

                    En_ij = model.NewBoolVar(f"En_{e1}_{e2}")
                    model.AddBoolAnd(all_En).OnlyEnforceIf(En_ij)

                    delayBetween = 0
                    for u, v in path_threads[i: j]:
                        type_v = G.nodes[v]["attr"]
                        delayBetween += GATES[type_v]["Delay"]

                    # Diff = model.NewIntVar(0, MAX_SIGMA, f"Diff_{e1}_{e2}")
                    model.Add(First_j > Last_i).OnlyEnforceIf(En_ij)

                    prod = model.NewIntVar(0, NUM_PHASES * T_PHASE, f"prod_{e1}_{e2}")
                    model.AddMultiplicationEquality(prod, [T_PHASE, First_j - Last_i])

                    model.Add(CQ_i + delayBetween + Setup_j <= prod).OnlyEnforceIf(En_ij)

        model.Minimize(sum(cost_fun))

        # Solves and prints out the solution.
        solver = cp.CpSolver()
        solver.parameters.max_time_in_seconds = 600.0
        print(f'Starting Macro ILP')
        status = solver.Solve(model)
        print(status)
        print(f'Solve status: {solver.StatusName(status)}')
        if (solver.StatusName(status) in ("OPTIMAL", "FEASIBLE")):
            print(f'Objective value: {solver.ObjectiveValue()}')

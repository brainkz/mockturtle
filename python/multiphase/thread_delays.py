import sys
from collections import namedtuple
from pprint import pprint
from typing import Dict, FrozenSet, List, Tuple, Union  # noqa

edgeVars = namedtuple("edgeVars", ["En", "First", "Last", "nDFF", "hasDFF", "CQ", "Setup", "Hold"])

import networkx as nx
from ortools.sat.python import cp_model as cp

AS = ("PI", "DFF", "NOT", "XOR")
SA = ("AND", "OR")
AA = ("MRG", "SPL")

MAX_SIGMA: int = 1000
NUM_PHASES: int = 1

# time unit is 1e-13 s
UNIT = 1e-13

T_PHASE = 189
T_CLK: int = T_PHASE * NUM_PHASES

GATES = {
    "PI": {"Setup": 44, "Hold": 39, "CQ": 79},
    "XOR": {"Setup": 69, "Hold": 55, "CQ": 72},
    "NOT": {"Setup": 44, "Hold": 69, "CQ": 93},
    "DFF": {"Setup": 44, "Hold": 39, "CQ": 79},
    "AND": {"Delay": 57},
    "OR": {"Delay": 57},
    "MRG": {"Delay": 57},
    "SPL": {"Delay": 66},
}

# Helper function to perform DFS and collect valid paths


def dfs_collect_paths(G, node, path, result):
    _type = G.nodes[node]["type"]
    print(f"Extending path {path}")
    print(f"\tNode: {node}")
    print(f"\tType: {_type}")

    # Check if this node is an end node and not the start node
    if node != path[0] and _type in (AS + SA):
        final_path = tuple(path) + (node,)
        result[frozenset(final_path)] = [final_path]
        # result.append(final_path)
        return

    # Only traverse further if the current node is AA
    if _type in AA:
        for neighbor in G.successors(node):
            dfs_collect_paths(G, neighbor, path + [node], result)


def get_paths(G):
    # Find all AS/SA nodes as start/end points
    start_end_nodes = [n for n, _type in G.nodes(data="type") if _type in (AS + SA)]
    print(f"AS/SA nodes: {start_end_nodes}")

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


def add_FLEn(model: cp.CpModel, u: int, v: int, type_u: str, type_v: str, e_prev: Tuple[Union[None, Tuple[int, int]]]):
    First = model.NewIntVar(0, MAX_SIGMA, f"First_Phase_{u}_{v}")
    Last = model.NewIntVar(0, MAX_SIGMA, f"Last_Phase_{u}_{v}")

    if (type_u in AA) or (type_v in AA) or (type_u in AS and type_v in SA):
        model.Add(First <= Last)
    else:
        model.Add(First < Last)

    if (type_u in AA) and (type_v in AA):
        En = model.NewBoolVar(f"isClocked_{u}_{v}")
        model.Add(First == edge_vars[e_prev].Last).OnlyEnforceIf(En.Not())
        model.Add(Last == edge_vars[e_prev].Last).OnlyEnforceIf(En.Not())
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
        # # Only enabled if En is true, i.e., if there is a DFF along the edge
        # # zero otherwise
        # model.AddDivisionEquality(nDFF, numerator, NUM_PHASES_EFF).OnlyEnforceIf(En)
        # model.Add(nDFF == 0).OnlyEnforceIf(En.Not())

        # REFORMULATING DIVISION due to error with AddDivisionEquality
        model.Add(nDFF * NUM_PHASES_EFF <= numerator).OnlyEnforceIf(En)
        model.Add(numerator < (nDFF + 1) * NUM_PHASES_EFF).OnlyEnforceIf(En)
        # # zero otherwise
        model.Add(nDFF == 0).OnlyEnforceIf(En.Not())
    else:
        # No conditionals needed
        # model.AddDivisionEquality(nDFF, numerator, NUM_PHASES_EFF)

        # REFORMULATING DIVISION due to error with AddDivisionEquality
        model.Add(nDFF * NUM_PHASES_EFF <= numerator)
        model.Add(numerator < (nDFF + 1) * NUM_PHASES_EFF)

    hasDFF = model.NewBoolVar(f"hasDFF_{u}_{v}")
    model.Add(nDFF == 0).OnlyEnforceIf(hasDFF.Not())
    model.Add(nDFF > 0).OnlyEnforceIf(hasDFF)
    return nDFF, hasDFF


def add_CQ(model: cp.CpModel, u: int, v: int, hasDFF: cp.IntVar, type_u: str, pred_types: Tuple[str] = tuple()):
    # if ever needed, the CQ time is DFF CQ time
    if type_u in AA + ("DFF",):
        CQ = model.NewConstant(GATES["DFF"]["CQ"])  #
        # viable_CQ = [GATES[type_u]["CQ"], GATES["DFF"]["CQ"]]
        # CQ = model.NewIntVarFromDomain(cp.Domain.FromValues(viable_CQ), f'CQ_{u}_{v}')
    elif type_u in AS:
        viable_CQ = [GATES[type_u]["CQ"], GATES["DFF"]["CQ"]]
        CQ = model.NewIntVarFromDomain(cp.Domain.FromValues(viable_CQ), f'CQ_{u}_{v}')
        # CQ = model.NewIntVar(*sorted(viable_CQ), f'CQ_{u}_{v}')
        model.Add(CQ == GATES[type_u]["CQ"]).OnlyEnforceIf(hasDFF.Not())
        model.Add(CQ == GATES["DFF"]["CQ"]).OnlyEnforceIf(hasDFF)
    elif type_u in SA:
        if not pred_types:
            CQ = model.NewConstant(GATES["DFF"]["CQ"])
        else:
            pred_CQ = max(GATES[pred_type]["CQ"] for pred_type in pred_types)
            viable_CQ = [pred_CQ + GATES[type_u]["Delay"], GATES["DFF"]["CQ"]]
            CQ = model.NewIntVarFromDomain(cp.Domain.FromValues(viable_CQ), f'CQ_{u}_{v}')
            model.Add(CQ == GATES[type_u]["CQ"]).OnlyEnforceIf(hasDFF.Not())
            model.Add(CQ == GATES["DFF"]["CQ"]).OnlyEnforceIf(hasDFF)
    return CQ


def add_SetupHold(model: cp.CpModel, u: int, v: int, hasDFF: cp.IntVar, type_u: str, type_v: str):
    if type_u in AA:
        # if ever needed, the setup time is DFF setup time
        if type_v in AA + SA + ("DFF",):
            Setup = model.NewConstant(GATES["DFF"]["Setup"])
            Hold = model.NewConstant(GATES["DFF"]["Hold"])
        elif type_v in AS:
            viable_Setup = [GATES[type_v]["Setup"], GATES["DFF"]["Setup"]]
            Setup = model.NewIntVarFromDomain(cp.Domain.FromValues(viable_Setup), f'Setup_{u}_{v}')
            model.Add(Setup == GATES[type_v]["Setup"]).OnlyEnforceIf(hasDFF.Not())
            model.Add(Setup == GATES["DFF"]["Setup"]).OnlyEnforceIf(hasDFF)

            viable_Hold = [GATES[type_v]["Hold"], GATES["DFF"]["Hold"]]
            Hold = model.NewIntVarFromDomain(cp.Domain.FromValues(viable_Hold), f'Hold_{u}_{v}')
            model.Add(Hold == GATES[type_v]["Hold"]).OnlyEnforceIf(hasDFF.Not())
            model.Add(Hold == GATES["DFF"]["Hold"]).OnlyEnforceIf(hasDFF)

        return Setup, Hold
    else:
        # Setup time is not needed for edges starting with AS gates
        # Return some dummy values
        return model.NewConstant(0), model.NewConstant(0)


def solve_model(model: cp.CpModel):
    solver = cp.CpSolver()
    solver.parameters.cp_model_presolve = False  # Turn off presolve to preserve all constraints
    solver.parameters.log_search_progress = True  # Log solver steps
    solver.parameters.max_time_in_seconds = 600.0
    print(f'Starting Macro ILP')
    status = solver.Solve(model)
    print(status)
    assert solver.StatusName() in ("OPTIMAL", "FEASIBLE")
    return solver, status


hold_safe_short_paths = (NUM_PHASES > 1)

if __name__ == "__main__":

    # Initialize a directed graph
    G = nx.DiGraph()

    G.add_node("A", type="DFF", letter="A")  # 2
    G.add_node("B", type="AND", letter="B")  # 3
    G.add_node("C", type="NOT", letter="C")  # 4

    G.add_node("a", type="SPL", letter="a")  # 5
    G.add_node("b", type="MRG", letter="b")  # 6
    G.add_node("c", type="SPL", letter="c")  # 7
    G.add_node("d", type="MRG", letter="d")  # 8
    G.add_node("e", type="SPL", letter="e")  # 9

    G.add_node("W", type="AND", letter="W")  # 10
    G.add_node("X", type="DFF", letter="X")  # 11
    G.add_node("Y", type="XOR", letter="Y")  # 12
    G.add_node("Z", type="OR", letter="Z")  # 13

    G.add_edge("A", "d")
    G.add_edge("d", "Z")
    G.add_edge("B", "a")
    G.add_edge("a", "W")
    G.add_edge("a", "b")
    G.add_edge("b", "c")
    G.add_edge("c", "d")
    G.add_edge("c", "e")
    G.add_edge("e", "X")
    G.add_edge("e", "Y")
    G.add_edge("C", "b")

    # # Parse the CSV file and build the graph
    # with open("ilp_config.csv", "r") as f:
    #     for line in f:
    #         gate_id_str, func, fanins_str, attr_str = line.split(",")
    #         gate_id = int(gate_id_str)
    #         fanins = tuple(int(q) for q in fanins_str.split("|"))
    #         attr = int(attr_str)

    #         # Add the gate node to the graph, with its attribute
    #         G.add_node(gate_id, type=attr, func=func)

    #         # Add edges from fanins to the current gate
    #         for fanin in fanins:
    #             G.add_edge(fanin, gate_id)

    independent_paths = get_paths(G)
    cost_fun = []

    # For each edge in a path, create appropriate variables and constraints
    for nodes, threads in independent_paths.items():

        model = cp.CpModel()

        edge_vars: Dict[Tuple[int, int], edgeVars] = {}
        inter_vars: Dict[Tuple[int, int], cp.IntVar] = {}
        for thread in threads:

            thread_edges = tuple(zip(thread[:-1], thread[1:]))
            prev_edges = (None, ) + thread_edges[:-1]

            for e, e_prev in zip(thread_edges, prev_edges):
                if e in edge_vars:
                    continue
                u, v = e
                type_u = G.nodes[u]["type"]
                type_v = G.nodes[v]["type"]
                # Variables associated with this edge

                First, Last, En = add_FLEn(model, u, v, type_u, type_v, e_prev)
                nDFF, hasDFF = add_nDFF(model, u, v, En, type_u, type_v, hold_safe=hold_safe_short_paths)

                CQ = add_CQ(model, u, v, hasDFF, type_u)
                Setup, Hold = add_SetupHold(model, u, v, hasDFF, type_u, type_v)

                cost_fun.append(nDFF)

                edge_vars[e] = edgeVars(En, First, Last, nDFF, hasDFF, CQ, Setup, Hold)

            # print('Starting outer check')
            # solve_model(model)
            # print('Passed outer check')

            for i, ei in enumerate(thread_edges[:-1]):
                En_i, First_i, Last_i, nDFF_i, hasDFF_i, CQ_i, Setup_i, Hold_i = edge_vars[ei]
                for j, ej in enumerate(thread_edges[i + 1:], i + 1):
                    En_j, First_j, Last_j, nDFF_j, hasDFF_j, CQ_j, Setup_j, Hold_j = edge_vars[ej]
                    all_En = [edge_vars[e].En.Not() for e in thread_edges[i + 1: j]] + [En_i, En_j]

                    all_En_Inv = [edge_vars[e].En for e in thread_edges[i + 1: j]] + [En_i.Not(), En_j.Not()]

                    print(f"Analyzing {''.join(ei)} and {''.join(ej)}")
                    print(f"\tedges in between: {[''.join(e) for e in thread_edges[i + 1: j]]}")

                    En_ij = model.NewBoolVar(f"En_{ei}_{ej}")
                    model.AddBoolAnd(all_En).OnlyEnforceIf(En_ij)
                    model.AddBoolOr(all_En_Inv).OnlyEnforceIf(En_ij.Not())
                    inter_vars[ei, ej] = En_ij

                    print(f"En_ij: {En_ij}, all_En: {all_En}, all_En_Inv: {all_En_Inv}")

                    delayBetween = 0
                    for u, v in thread_edges[i: j]:
                        type_v = G.nodes[v]["type"]
                        delayBetween += GATES[type_v]["Delay"]

                    Diff = model.NewIntVar(0, MAX_SIGMA, f"Diff_{ei}_{ej}")
                    model.Add(Diff == First_j - Last_i)
                    model.Add(Diff > 0).OnlyEnforceIf(En_ij)
                    model.Add(Diff <= NUM_PHASES).OnlyEnforceIf(En_ij)
                    prod = model.NewIntVar(0, MAX_SIGMA * T_PHASE, f"prod_{ei}_{ej}")
                    model.AddMultiplicationEquality(prod, [T_PHASE, Diff])
                    model.Add(CQ_i + delayBetween + Setup_j <= prod).OnlyEnforceIf(En_ij)

                    print(f"T_PHASE: {T_PHASE}, CQ_i: {CQ_i}, delayBetween: {delayBetween}, Setup_j: {Setup_j}, total: {CQ_i + delayBetween + Setup_j}")

                    # print("*" * 10 + 'STARTING INTER CHECK')
                    # solve_model(model)
                    # print("*" * 10 + 'PASSED INTER CHECK')
            # sys.exit(0)

        model.Minimize(sum(cost_fun))

        # # Solves and prints out the solution.
        # solver = cp.CpSolver()
        # solver.parameters.cp_model_presolve = False  # Turn off presolve to preserve all constraints
        # solver.parameters.log_search_progress = True  # Log solver steps
        # solver.parameters.max_time_in_seconds = 600.0
        # print(f'Starting Macro ILP')
        # status = solver.Solve(model)
        # print(status)
        solver, status = solve_model(model)
        print(f'Solve status: {solver.StatusName(status)}')
        if (solver.StatusName(status) in ("OPTIMAL", "FEASIBLE")):
            print(f'Objective value: {solver.ObjectiveValue()}')
            for (u, v), var in edge_vars.items():
                type_u = G.nodes[u]["type"]
                type_v = G.nodes[v]["type"]
                print("*" * 10 + f" {u} ({type_u}) -> {v} ({type_v})")
                print(solver.Values((var.En, var.First, var.Last, var.nDFF)))
            for ((ui, vi), (uj, vj)), var in inter_vars.items():
                type_ui = G.nodes[ui]["type"]
                type_vi = G.nodes[vi]["type"]
                type_uj = G.nodes[uj]["type"]
                type_vj = G.nodes[vj]["type"]
                print(f"e_i: {ui} ({type_ui}) -> {vi} ({type_vi})")
                print(f"e_j: {uj} ({type_uj}) -> {vj} ({type_vj})")
                print(solver.Value(var))

import os
import sys
import traceback
from collections import defaultdict
from itertools import combinations
from pprint import pprint
from typing import Any, Dict, FrozenSet, List, Set, Tuple, Union

import networkx as nx
from ortools.sat.python import cp_model

AS = ("PI", "DFF", "NOT", "XOR")
SA = ("AND", "OR")
AA = ("MRG", "SPL")

# Helper function to perform DFS and collect valid paths


def dfs_collect_paths(G, node, path, result):
    attr = G.nodes[node]['attr']
    print(f"Extending path {path}")
    print(f"\tNode: {node}")
    print(f"\tType: {attr}")

    # Check if this node is an end node (attr 1 or 2) and not the start node
    if node != path[0] and attr in (AS + SA):
        final_path = tuple(path) + (node,)
        result[frozenset(final_path)] = [final_path]
        # result.append(final_path)
        return

    # Only traverse further if the current node is AA
    if attr in AA:
        for neighbor in G.successors(node):
            dfs_collect_paths(G, neighbor, path + [node], result)


def get_paths(G):
    # Find all nodes with attribute (1) or (2) as start/end points
    start_end_nodes = [n for n, attr in G.nodes(data='attr') if attr in (AS + SA)]
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


MAX_SIGMA: int = 1000
NUM_PHASES: int = 4

GATE_TYPES = {
    "PI": {
        "sync_in": False,
        "sync_out": True,
        "C-Q": 0,
        "Setup": 0,
    },
    "DFF": {
        "sync_in": False,
        "sync_out": True,
        "C-Q": 9.1e-12,
        "Setup": 4.5e-12,
    },
    "NOT": {
        "sync_in": False,
        "sync_out": True,
        "C-Q": 10.0e-12,
        "Setup": 4.8e-12,
    },
    "XOR": {
        "sync_in": False,
        "sync_out": True,
        "C-Q": 7.9e-12,
        "Setup": 7.3e-12,
    },
    "AND": {
        "sync_in": True,
        # "sync_in_tol" : 5e-12, # IDEA: consider tolerance to SFQ pulse mismatch
        "sync_out": False,
        "Delay": 6.2e-12
    },
    "OR": {
        "sync_in": True,
        # "sync_in_tol" : 5e-12, # IDEA: consider tolerance to SFQ pulse mismatch
        "sync_out": False,
        "Delay": 6.3e-12
    },
    "MRG": {
        "sync_in": False,
        # "sync_in_tol" : math("inf"), # IDEA: consider tolerance to SFQ pulse mismatch
        "sync_out": False,
        "Delay": 6.3e-12
    },
    "SPL": {
        "sync_in": False,
        # "sync_in_tol" : math("inf"), # IDEA: consider tolerance to SFQ pulse mismatch
        "sync_out": False,
        "Delay": 8.0e-12
    },
}

if __name__ == "__main__":

    # Model.
    model = cp_model.CpModel()

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

    edge_vars: Dict[Tuple[int, int], List[Union[cp_model.IntVar]]] = {}
    # For each edge in a path, create a variable
    for nodes, paths in independent_paths.items():
        edges = {e for path in paths for e in zip(path[:-1], path[1:])}
        for e in edges:
            u, v = e
            type_u = G.nodes[u]['attr']
            type_v = G.nodes[v]['attr']
            # Variables associated with this edge

            # En = model.NewBoolVar(f"Has_Clocked_Gate_{u}_{v}") # moving to cases
            First = model.NewIntVar(0, MAX_SIGMA, f"First_Phase_{u}_{v}")
            Last = model.NewIntVar(0, MAX_SIGMA, f"Last_Phase_{u}_{v}")

            if type_u in AA or type_v in AA:
                model.Add(First <= Last)
            else:
                # if both are clocked, cannot have the same sigma
                model.Add(First < Last)

            if type_u not in AA or type_v not in AA:
                # if edge has clocked gate, En==1
                En = model.NewConstant(1)
            else:
                En = model.NewBoolVar(f"Has_DFF_{u}_{v}")

            # hasDFF = model.NewBoolVar(f"Has_DFF_{u}_{v}")

            # Constraints associated with this edge
            if type_u in (AS + SA) and type_v in AA:
                numerator = model.NewIntVar(0, MAX_SIGMA, f'NUM_{u}_{v}')
                model.Add(numerator == (Last - First + NUM_PHASES - 1))

                # No conditionals are needed here
                nDFF = model.NewIntVar(0, MAX_SIGMA, f"nDFF_{u}_{v}")
                model.AddDivisionEquality(nDFF, numerator, NUM_PHASES)

                Setup = model.NewConstant(0)  # Not applicable here

                if type_u == "DFF":
                    # The gate is DFF, CQ delay is always the same
                    C_Q = model.NewConstant(GATE_TYPES["DFF"]["C-Q"])
                elif type_u in AS:
                    # delays can differ
                    possible_CQ = [
                        GATE_TYPES[type_u]["C-Q"],
                        GATE_TYPES["DFF"]["C-Q"],
                    ]
                    C_Q = model.NewIntVarFromDomain(cp_model.Domain.FromValues(possible_CQ), f'TCQ_{u}_{v}')
                    model.Add(C_Q == GATE_TYPES[type_u]["C-Q"]).OnlyEnforceIf(nDFF == 0)
                    model.Add(C_Q == GATE_TYPES["DFF"]["C-Q"]).OnlyEnforceIf(nDFF > 0)

                elif type_u in SA:
                    # Need to determine the C-Q delay of the predecessor
                    pred_CQ = 0
                    for p in G.predecessors(u):
                        pred_type = G.nodes[p]["attr"]
                        pred_CQ = max(pred_CQ, GATE_TYPES[type_u]["C-Q"])

                    # delays can differ
                    possible_CQ = [
                        pred_CQ + GATE_TYPES[type_u]["Delay"],  # no DFF, use pred + AND/OR delay
                        GATE_TYPES["DFF"]["C-Q"],  # has DFF, use just DFF
                    ]
                    C_Q = model.NewIntVarFromDomain(cp_model.Domain.FromValues(possible_CQ), f'TCQ_{u}_{v}')
                    model.Add(C_Q == pred_CQ + GATE_TYPES[type_u]["Delay"]).OnlyEnforceIf(nDFF == 0)
                    model.Add(C_Q == GATE_TYPES["DFF"]["C-Q"]).OnlyEnforceIf(nDFF > 0)

            elif type_u in AA and type_v in AA:
                numerator = model.NewIntVar(0, MAX_SIGMA, f'NUM_{u}_{v}')
                model.Add(numerator == (Last - First + 2 * NUM_PHASES - 1))

                nDFF = model.NewIntVar(0, MAX_SIGMA, f"nDFF_{u}_{v}")

                # Only enabled if En is true, i.e., there are any DFFs along the path
                model.AddDivisionEquality(nDFF, numerator, NUM_PHASES).OnlyEnforceIf(En)
                # zero otherwise
                model.Add(nDFF == 0).OnlyEnforceIf(En.Not())

                C_Q = model.NewIntVarFromDomain(cp_model.Domain.FromValues([0, GATE_TYPES["DFF"]["C-Q"]]), f'TCQ_{u}_{v}')
                model.Add(C_Q == GATE_TYPES["DFF"]["C-Q"]).OnlyEnforceIf(En)
                model.Add(C_Q == 0).OnlyEnforceIf(En.Not())

                Setup = model.NewIntVarFromDomain(cp_model.Domain.FromValues([0, GATE_TYPES["DFF"]["Setup"]]), f'Setup_{u}_{v}')
                model.Add(Setup == GATE_TYPES["DFF"]["Setup"]).OnlyEnforceIf(En)
                model.Add(Setup == 0).OnlyEnforceIf(En.Not())

            elif type_u in AA and type_v in AS:
                numerator = model.NewIntVar(0, MAX_SIGMA, f'NUM_{u}_{v}')
                model.Add(numerator == (Last - First + NUM_PHASES - 1))

                # No conditionals are needed here
                nDFF = model.NewIntVar(0, MAX_SIGMA, f"nDFF_{u}_{v}")
                model.AddDivisionEquality(nDFF, numerator, NUM_PHASES)

                C_Q = model.NewConstant(0)

                if type_v == "DFF":
                    Setup = model.NewConstant(GATE_TYPES["DFF"]["Setup"])
                elif type_v in AS:
                    possible_Setup = [
                        GATE_TYPES[type_v]["Setup"],
                        GATE_TYPES["DFF"]["Setup"],
                    ]
                    Setup = model.NewIntVarFromDomain(cp_model.Domain.FromValues(possible_Setup), f'Setup_{u}_{v}')
                    model.Add(Setup == GATE_TYPES["DFF"]["Setup"]).OnlyEnforceIf(nDFF > 0)
                    model.Add(Setup == GATE_TYPES[type_v]["Setup"]).OnlyEnforceIf(nDFF == 0)

            elif type_u in AA and type_v in SA:
                numerator = model.NewIntVar(0, MAX_SIGMA, f'NUM_{u}_{v}')
                model.Add(numerator == (Last - First + 2 * NUM_PHASES - 1))

                # No conditionals are needed here
                nDFF = model.NewIntVar(0, MAX_SIGMA, f"nDFF_{u}_{v}")
                model.AddDivisionEquality(nDFF, numerator, NUM_PHASES)

                C_Q = model.NewConstant(0)
                Setup = model.NewConstant(GATE_TYPES["DFF"]["Setup"])

            # The rest are isolated cases, should we include them?
            elif type_u in (AS, SA) and type_v in AS:
                numerator = model.NewIntVar(0, MAX_SIGMA, f'NUM_{u}_{v}')
                model.Add(numerator == (Last - First - 1))

                # No conditionals are needed here
                nDFF = model.NewIntVar(0, MAX_SIGMA, f"nDFF_{u}_{v}")
                model.AddDivisionEquality(nDFF, numerator, NUM_PHASES)

                C_Q = model.NewConstant(0)
                Setup = model.NewConstant(0)

            elif type_u in (AS, SA) and type_v in SA:
                numerator = model.NewIntVar(0, MAX_SIGMA, f'NUM_{u}_{v}')
                model.Add(numerator == (Last - First + NUM_PHASES - 1))

                # No conditionals are needed here
                nDFF = model.NewIntVar(0, MAX_SIGMA, f"nDFF_{u}_{v}")
                model.AddDivisionEquality(nDFF, numerator, NUM_PHASES)

                C_Q = model.NewConstant(0)
                Setup = model.NewConstant(0)

            cost_fun.append(nDFF)

            edge_vars[e] = {
                "En": En,
                "First": First,
                "Last": Last,
                "C-Q": C_Q,
                "Setup": Setup,
            }
    # import builtins

    # LOG_FILE = "ilp_log.log"

    # def print(*args, **kwargs):
    #     with open("ilp_log.log", "w") as f:
    #         f.write()

if False:

    EPS = 1e-6

    # TYPES = {"0": "PI",
    #          "1": "AA",
    #          "2": "AS",
    #          "3": "SA",
    #          "4": "FA"}

    # TYPES = {"0": "AA", "1": "AS", "2": "SA"}

    class Primitive:
        __slots__ = ['sig', 'type', 'fanins', 'in_neg']

        def __init__(self, _sig: int, _type: str, _fanins: list, _in_neg: list = []) -> None:
            self.sig = _sig
            self.type = _type if _fanins else 'PI'
            self.fanins = _fanins
            self.in_neg = _in_neg

        def __repr__(self) -> str:
            return f"Primitive(sig={self.sig}, type='{self.type}', fanins={self.fanins}, in_neg={self.in_neg}"

    def parse_specs(filename: str) -> dict[int, Primitive]:
        items = {}
        print(f"Reading {filename}")
        with open(filename, 'r') as f:
            line_iter = iter(f)

            # first line is for PIs
            pi_line = next(line_iter)
            for _sig_str in pi_line.split(',')[1:]:
                _sig = int(_sig_str)
                items[_sig] = Primitive(_sig, "PI", [])

            for line in line_iter:
                line = line.strip()

                _all_sig_str, _type_str, _fanin_str = line.split(',')
                print(f"Processing line: {line}")
                _fanins = [int(_fanin_sig) for _fanin_sig in _fanin_str.split('|') if _fanin_sig]
                _sig = int(_all_sig_str)
                items[_sig] = Primitive(_sig, _type_str, _fanins)
                print(f"\tCreating: {items[_sig]}")
        return items

    if __name__ == "__main__":

        # Model.
        model = cp_model.CpModel()

        NPH = int(sys.argv[1])
        cfg_path = sys.argv[2]

        base_name, _ = os.path.splitext(cfg_path)
        log_filename = base_name + ".log"

        # Save the original sys.stdout
        original_stdout = sys.stdout

        os.makedirs(os.path.dirname(log_filename), exist_ok=True)
        log_file = open(log_filename, 'w')
        # Redirect sys.stdout to the log file
        sys.stdout = log_file

        try:
            print(f"Parsing specs file {cfg_path}")
            all_signals = parse_specs(cfg_path)
            # print(all_signals)

            print(f'Finished creating [all_signals]')
            pprint(all_signals)

            # max_phase = max(g.phase for g in all_signals.values() if g.phase is not None)

            # sigma_bounds = (0, max_phase + NPH)
            sigma_bounds = (0, 1000)

            print(f"Setting up the model {cfg_path}")
            fanout = defaultdict(int)
            Sigma = {}
            fanin_constr = {}
            ctr = 0
            for g in all_signals.values():
                if g.type == 'PI':
                    print(f"Creating PI SIGMA_{g.sig}")
                    Sigma[g.sig] = model.NewIntVar(0, NPH - 1, f"SIGMA_{g.sig}")
                else:
                    print(f"Creating {g.type} gate SIGMA_{g.sig}")
                    Sigma[g.sig] = model.NewIntVar(*sigma_bounds, f"SIGMA_{g.sig}")
                    for p_sig in g.fanins:
                        fanout[p_sig] += 1

            print(f'Finished creating [Sigma]')
            # pprint(Sigma)

            fanin_diff = {}
            Delta = {}
            # IMPORTANT: PIs are assumed to be at clock stage 0 but with an arbitrary phase
            expr = []
            already_processed = []
            for g in all_signals.values():
                # print(f"The type of the gate {g.sig} is {g.type}")
                if g.type == 'PI':
                    continue
                elif g.type == 'AA':
                    # any PI is the earliest fanin, with zero phase and zero clock stage
                    # IMPORTANT: only 2-input CB is supported
                    assert (len(g.fanins) == 2)
                    a_sig, b_sig = g.fanins
                    a = all_signals[a_sig]
                    b = all_signals[b_sig]

                    max_diff = model.NewIntVar(*sigma_bounds, f'max_{a.sig}_{b.sig}')
                    model.AddMaxEquality(max_diff, [Sigma[g.sig] - Sigma[a.sig],
                                                    Sigma[g.sig] - Sigma[b.sig]])

                    div_val = model.NewIntVar(*sigma_bounds, f'div_{a.sig}_{b.sig}')
                    model.AddDivisionEquality(div_val, max_diff, NPH)
                    expr.append(div_val)

                    # print(f"AA constraint: {Sigma[a.sig].Name()} <= {Sigma[g.sig].Name()}")
                    # print(f"AA constraint: {Sigma[b.sig].Name()} <= {Sigma[g.sig].Name()}")

                    model.Add(Sigma[a.sig] <= Sigma[g.sig])
                    model.Add(Sigma[b.sig] <= Sigma[g.sig])

                # Regular gate - use the standard definition
                elif g.type == 'AS':
                    for i, p_sig in enumerate(g.fanins):
                        # print(f"AS constraint: {Sigma[p_sig].Name()} < {Sigma[g.sig].Name()}")
                        model.Add(Sigma[p_sig] < Sigma[g.sig])

                        # expr.append( Sigma[g.sig] - Sigma[p_sig] )

                        delta = model.NewIntVar(*sigma_bounds, f'delta_{p_sig}_{g.sig}')
                        model.Add(delta == (Sigma[g.sig] - Sigma[p_sig]))

                        div_val = model.NewIntVar(*sigma_bounds, f'div_{p_sig}_{g.sig}')
                        model.AddDivisionEquality(div_val, delta, NPH)
                        expr.append(div_val)

                elif g.type == 'SA':

                    for i, p_sig in enumerate(g.fanins):
                        p = all_signals[p_sig]
                        if (p.type == 'AS') and (fanout[p.sig] == 1):
                            model.Add(Sigma[p_sig] <= Sigma[g.sig])
                            # print(f"SA constraint: {Sigma[p_sig].Name()} <= {Sigma[g.sig].Name()} ({p.type}) ({fanout[p.sig]})")
                        else:
                            model.Add(Sigma[p_sig] < Sigma[g.sig])
                            # print(f"SA constraint: {Sigma[p_sig].Name()} <  {Sigma[g.sig].Name()} ({p.type}) ({fanout[p.sig]})")

                        # expr.append( Sigma[g.sig] - Sigma[p_sig] )

                        delta = model.NewIntVar(*sigma_bounds, f'delta_sa_{p_sig}_{g.sig}')
                        model.Add(delta == (Sigma[g.sig] - Sigma[p_sig] + NPH - 1))

                        div_val = model.NewIntVar(*sigma_bounds, f'div_{p_sig}_{g.sig}')
                        model.AddDivisionEquality(div_val, delta, NPH)
                        expr.append(div_val)

                elif g.type == 'FA':
                    # Make sure the T1 cell has not yet been processed
                    if g.sig in already_processed:
                        continue

                    # find siblings
                    cell_idx = -1
                    for i, sibling_list in enumerate(all_t1_cells):
                        if g.sig in sibling_list.values():
                            cell_idx = i
                            break
                    else:
                        raise ValueError(f"T1 cell output\n {g} \n\tis not represented in all_t1_cells")

                    assert (len(g.fanins) == 3)
                    # sig_idx = [p_sig for p_sig in g.fanins]
                    # TODO : check fanins
                    # TODO :    if the fanin is AA, add DFF/NOT to cost and [++SIGMA_EFF] (not accurate but too complex otherwise)
                    # TODO :    if the fanin is AS/SA and negated, [++SIGMA_EFF]
                    incr_Sigma = [0, 0, 0]
                    for i, (p_sig, p_neg) in enumerate(zip(g.fanins, g.in_neg)):
                        p = all_signals[p_sig]
                        if (p.type == 'AA') or (p.type in ('AS', 'SA') and p_neg):
                            incr_Sigma[i] += 1

                        # TODO : create vars [ABS_MIN=min(A,B,C)], [MED=median(A,B,C)-ABS_MIN], [MAX=max(A,B,C)-ABS_MIN]
                    ABS_MIN = model.NewIntVar(*sigma_bounds, f'absmin_{g.sig}')
                    model.AddMinEquality(ABS_MIN, [(Sigma[p_sig] + incr_Sigma[i]) for i, p_sig in enumerate(g.fanins)])

                    ABS_MAX = model.NewIntVar(*sigma_bounds, f'absmax_{g.sig}')
                    model.AddMaxEquality(ABS_MAX, [(Sigma[p_sig] + incr_Sigma[i]) for i, p_sig in enumerate(g.fanins)])

                    ABS_MED = model.NewIntVar(*sigma_bounds, f'absmed_{g.sig}')
                    model.Add(ABS_MED == (Sigma[g.fanins[0]] + incr_Sigma[0]) + (Sigma[g.fanins[1]] + incr_Sigma[1]) + (Sigma[g.fanins[2]] + incr_Sigma[2]) - ABS_MIN - ABS_MAX)

                    # IMPORTANT: this is the stage of the XOR cell only!!!
                    # IMPORTANT: TODO: ADD the stages of other outputs
                    T1_XOR_SIGMA = model.NewIntVar(*sigma_bounds, f'sigma_xor_{g.sig}')
                    model.Add(T1_XOR_SIGMA > ABS_MIN + 2)
                    model.Add(T1_XOR_SIGMA > ABS_MED + 1)
                    model.Add(T1_XOR_SIGMA > ABS_MAX)

                    MED_MIN_diff_TMP = model.NewIntVar(*sigma_bounds, f'mod_med_min_diff_tmp_{g.sig}')
                    model.Add(MED_MIN_diff_TMP == ABS_MED - ABS_MIN)
                    MED_MIN_diff = model.NewIntVar(*sigma_bounds, f'mod_med_min_diff_{g.sig}')
                    model.AddModuloEquality(MED_MIN_diff, MED_MIN_diff_TMP, NPH)

                    MAX_MED_diff_TMP = model.NewIntVar(*sigma_bounds, f'mod_max_med_diff_tmp_{g.sig}')
                    model.Add(MAX_MED_diff_TMP == ABS_MAX - ABS_MED)
                    MAX_MED_diff = model.NewIntVar(*sigma_bounds, f'mod_max_med_diff_{g.sig}')
                    model.AddModuloEquality(MAX_MED_diff, MAX_MED_diff_TMP, NPH)

                    MAX_MIN_diff_TMP = model.NewIntVar(*sigma_bounds, f'mod_max_min_diff_tmp_{g.sig}')
                    model.Add(MAX_MIN_diff_TMP == ABS_MAX - ABS_MIN)
                    MAX_MIN_diff = model.NewIntVar(*sigma_bounds, f'mod_max_min_diff_{g.sig}')
                    model.AddModuloEquality(MAX_MIN_diff, MAX_MIN_diff_TMP, NPH)

                    FLOOR_DIFF_MIN_TMP = model.NewIntVar(-1000, 1000, f'floor_diff_min_tmp_{g.sig}')
                    model.Add(FLOOR_DIFF_MIN_TMP == T1_XOR_SIGMA - ABS_MIN - 1)
                    FLOOR_DIFF_MIN = model.NewIntVar(*sigma_bounds, f'floor_diff_min_{g.sig}')
                    model.AddDivisionEquality(FLOOR_DIFF_MIN, FLOOR_DIFF_MIN_TMP, NPH)

                    FLOOR_DIFF_MED_TMP = model.NewIntVar(*sigma_bounds, f'floor_diff_med_tmp_{g.sig}')
                    model.Add(FLOOR_DIFF_MED_TMP == T1_XOR_SIGMA - ABS_MED - 1)
                    FLOOR_DIFF_MED = model.NewIntVar(*sigma_bounds, f'floor_diff_med_{g.sig}')
                    model.AddDivisionEquality(FLOOR_DIFF_MED, FLOOR_DIFF_MED_TMP, NPH)

                    FLOOR_DIFF_MAX_TMP = model.NewIntVar(*sigma_bounds, f'floor_diff_max_tmp_{g.sig}')
                    model.Add(FLOOR_DIFF_MAX_TMP == T1_XOR_SIGMA - ABS_MAX - 1)
                    FLOOR_DIFF_MAX = model.NewIntVar(*sigma_bounds, f'floor_diff_max_{g.sig}')
                    model.AddDivisionEquality(FLOOR_DIFF_MAX, FLOOR_DIFF_MAX_TMP, NPH)

                    MIN_FLAG = model.NewBoolVar(f'min_flag_{g.sig}')
                    model.Add(MED_MIN_diff == 0).OnlyEnforceIf(MIN_FLAG)
                    model.Add(FLOOR_DIFF_MIN == 0).OnlyEnforceIf(MIN_FLAG)

                    MED_FLAG = model.NewBoolVar(f'med_flag_{g.sig}')
                    model.Add(MAX_MED_diff == 0).OnlyEnforceIf(MED_FLAG)
                    model.Add(FLOOR_DIFF_MED == 0).OnlyEnforceIf(MED_FLAG)

                    # place
                    for g_sig in sibling_list.values():
                        model.Add(Sigma[g_sig] == T1_XOR_SIGMA)

                    # Cost function
                    expr.append(FLOOR_DIFF_MIN)
                    expr.append(FLOOR_DIFF_MED)
                    expr.append(FLOOR_DIFF_MAX)
                    expr.append(MIN_FLAG)
                    expr.append(MED_FLAG)

                    already_processed.extend(sibling_list.values())
                else:
                    raise ValueError(f"{g}\nUnsupported gate type")

            model.Minimize(sum(expr))

            # Solves and prints out the solution.
            solver = cp_model.CpSolver()
            solver.parameters.max_time_in_seconds = 600.0
            print(f'Starting Macro ILP')
            status = solver.Solve(model)
            print(status)
            print(f'Solve status: {solver.StatusName(status)}')
            if (solver.StatusName(status) in ("OPTIMAL", "FEASIBLE")):
                print(f'Objective value: {solver.ObjectiveValue()}')
                for sig, sigma in Sigma.items():
                    g = all_signals[sig]
                    print(f'{sigma.Name()[6:]} :{solver.Value(sigma)}')

            # Restore sys.stdout to the original value
            sys.stdout = original_stdout
            # Close the log file

            print(f'Solve status: {solver.StatusName(status)}')
            if (solver.StatusName(status) in ("OPTIMAL", "FEASIBLE")):
                print(f'Objective value: {solver.ObjectiveValue()}')
                for sig, sigma in Sigma.items():
                    g = all_signals[sig]
                    print(f'{sigma.Name()[6:]} :{solver.Value(sigma)}')

        except Exception as e:
            print(f'An error occurred: {str(e)}')
            traceback.print_exc(file=log_file)
            os.system(f"code {log_filename}")
        finally:
            log_file.close()

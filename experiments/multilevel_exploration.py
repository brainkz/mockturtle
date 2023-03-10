import networkx as nx
from itertools import combinations
from pprint import pprint
from collections import Counter, defaultdict
from pickle import dump, load
import builtins
import sys
import logging
from heapq import nlargest
from sortedcontainers import SortedList
from joblib import Parallel, delayed
import asyncio 

LEVELS = sorted((0, 0, 0, 0))
LOGFILE = 'results_' + '_'.join(str(lvl) for lvl in LEVELS) + '.log'

for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
logging.basicConfig(filename = LOGFILE, format='%(asctime)s : %(message)s', level=logging.DEBUG)
# , handlers=[ logging.FileHandler(LOGFILE), logging.StreamHandler(sys.stdout) ] 

ONES = 0xFFFF

COSTS = {   'MERGE' : 8,   'MERGE3': 11,  'XOR'   : 7,  'NOT'   : 9,  'DFF'   : 7,  'SPL'   : 7, }

MAXLVL = 5
    
# Func, cost, depth, xorable flag
PI = [(0xFF00, 0, LEVELS[0], True), (0xF0F0, 0, LEVELS[1], True), (0xCCCC, 0, LEVELS[2], True), (0xAAAA, 0, LEVELS[3], True), (0xFFFF, 0, 0, True), (0x0000, 0, 0, True)]

c_any  = {}
d_any  = {}
c_xor  = {}
d_xor  = {}
for func, *_ in PI:
    c_any[func] = {'cost' : 0, 'depth' : 0, 'xorable' : True, 'edges' : set()}
    d_any[func] = {'cost' : 0, 'depth' : 0, 'xorable' : True, 'edges' : set()}
    c_xor[func] = {'cost' : 0, 'depth' : 0, 'xorable' : True, 'edges' : set()}
    d_xor[func] = {'cost' : 0, 'depth' : 0, 'xorable' : True, 'edges' : set()}

nspl_dupe_cost_cache = {}

G = nx.DiGraph()
G.add_nodes_from(PI, nspl = 0, edges = set())

def count_spl(n1, n2):
    '''
    Assumes the "edges" attribute is the set of edges leading to each of n1 and n2
    '''
    if   (n1, n2) in nspl_dupe_cost_cache:
        return nspl_dupe_cost_cache[n1, n2]
    elif (n2, n1) in nspl_dupe_cost_cache:
        return nspl_dupe_cost_cache[n2, n1]
    # TODO: add nodes efficiently
    all_edges = G.nodes[n1]['edges'].union(G.nodes[n2]['edges'])
    count_spl = defaultdict(int)
    seen_sa = defaultdict(int)
    for u,v in all_edges:
        count_spl[u] += 1
        if G.edges[u,v]['type'] in ('AND', 'OR'):
            seen_sa[u] += 1
    nspl = sum((c-1 for c in count_spl.values()))
    overhead = 0
    seen_sa = defaultdict(int)
    for sa_node, count in seen_sa.items():
        if count > 1:
            preds = list(G.pred[sa_node])
            if len(preds) == 1: #NOT or DFF
                gate_type = G.edges[preds[0], sa_node]['type']
                overhead += COSTS[gate_type] * (count - 1)
            elif len(preds) == 2: # XOR
                assert(gate_type == 'XOR')
                overhead += COSTS[gate_type] * (count - 1)
                nspl += count - 1 # need to do twice more splittings for duplicating an XOR gate
    nspl_dupe_cost_cache[n1, n2] = nspl, overhead
    return nspl, overhead

'''
# Unused for now
def count_spl(n1, n2):
    if   (n1, n2) in nspl_dupe_cost_cache:
        return nspl_dupe_cost_cache[n1, n2]
    elif (n2, n1) in nspl_dupe_cost_cache:
        return nspl_dupe_cost_cache[n2, n1]
    nspl = 0
    seen = set()
    stack = [n1, n2]
    seen_sa = Counter()
    while stack:
        n = stack.pop()
        for p in G.pred[n]:
            if p not in seen:
                seen.add(p)
                stack.append(p)
                if G.edges[p,n]['type'] in ('AND', 'OR'):
                    seen_sa[p] += 1
            else:
                nspl += 1
    overhead = 0
    for sa_node, count in seen_sa.items():
        if count > 1:
            preds = list(G.pred[sa_node])
            if len(preds) == 1: #NOT or DFF
                gate_type = G.edges[preds[0], sa_node]['type']
                overhead += COSTS[gate_type] * (count - 1)
            elif len(preds) == 2: # XOR
                assert(gate_type == 'XOR')
                overhead += COSTS[gate_type] * (count - 1)
                nspl += count - 1 # need to do twice more splittings for duplicating an XOR gate
    nspl_dupe_cost_cache[n1, n2] = nspl, overhead
    return nspl, overhead
'''

def add_node(parents, func, edge_type, cost, lvl, xorable, nspl, force_add = False):
    node = (func, cost, lvl, xorable)
    edges = set()
    for p in parents:
        edges.update(G.nodes[p]['edges'])
        edges.add((p, node))
    attr = {'cost' : cost, 'depth' : lvl, 'xorable' : xorable, 'nspl' : nspl, 'edges' : edges}
    if force_add:
        G.add_node(node, **attr)
        for n in parents:
            G.add_edge(n, node, type = edge_type)
    if func not in c_any: # no such node exists
        G.add_node(node, **attr)
        for n in parents:
            G.add_edge(n, node, type = edge_type)
        c_any[func] = attr
        d_any[func] = attr
        if xorable:
            c_xor[func] = attr
            d_xor[func] = attr
    else:
        add_flag = False
        if c_any[func]['cost'] > cost or (c_any[func]['cost'] >= cost and c_any[func]['depth'] > lvl): # the new node has best cost 
            add_flag = True
            c_any[func] = attr
        if d_any[func]['depth'] > lvl or (d_any[func]['depth'] >= lvl and d_any[func]['cost'] > cost): # the new node has best depth 
            add_flag = True
            d_any[func] = attr
        if xorable:
            if func not in c_xor: # no xorable node existed before
                add_flag = True    
                c_xor[func] = attr
                d_xor[func] = attr
            else:
                if c_xor[func]['cost'] > cost or (c_xor[func]['cost'] >= cost and c_xor[func]['depth'] > lvl): # the new node has best cost among xorables
                    add_flag = True    
                    c_xor[func] = attr
                if d_xor[func]['depth'] > lvl or (d_xor[func]['depth'] >= lvl and d_xor[func]['cost'] > cost): # the new node has best depth among xorables
                    add_flag = True
                    d_xor[func] = attr
        if add_flag:
            G.add_node(node, **attr)
            for n in parents:
                G.add_edge(n, node, type = edge_type)
                
async def req(url):
    async with httpx.AsyncClient() as client:
        r = await client.get(url)

    return r.status_code

async def main():
    tasks = []

    for url in URLS:
        tasks.append(asyncio.create_task(req(url)))

    results = await asyncio.gather(*tasks)

    print(dict(zip(URLS, results)))


if __name__ == "__main__":
    asyncio.run(main())

def _one_iteration(n1, i, all_nodes, max_lvl):
    # normal loop
    out_dict = _process_row(n1, i, all_nodes, max_lvl)
    # yield the result
    return out_dict

def _process_row(n1, i, all_nodes, max_lvl):
    out_dict = defaultdict(SortedList)
    f1, c1, d1, x1 = n1
    for n2 in all_nodes[i+1:]:
        (f2, c2, d2, x2) = n2
        func = (f1 | f2)
        if func != f1 and func != f2: #the function should make sense
            # the function is unknown or a known implementation has larger delay
            # func is xorable AND (no known xorable implementation exists or known xorable has worse depth)
            xorable = (f1 ^ f2) == func
            preliminary_cost = c1 + c2
            if xorable and (                                # func is xorable AND
                (func not in d_xor) or                      # (no known xorable implementation exists OR
                (d_xor[func]['depth'] > max_lvl or          #  (the known xorable implementation has worse depth OR
                    (d_xor[func]['depth'] == max_lvl and       #   it has the same depth AND
                    preliminary_cost <= d_xor[func]['cost'])  #   its cost is not better)
                    )                                          # )
                ):
                out_dict[func].add((preliminary_cost, xorable, n1, n2))
            elif not xorable and (                          # func is not xorable AND
                (func not in d_any) or                      # (no known implementation exists OR
                (d_any[func]['depth'] > max_lvl or          #  (the known implementation has worse depth OR
                    (d_any[func]['depth'] == max_lvl and       #   it has the same depth AND
                    preliminary_cost <= d_any[func]['cost'])  #   its cost is not better)
                    )                                          # )
                ):
                out_dict[func].add((preliminary_cost, xorable, n1, n2))
    print(f"Finished gathering candidate node pairs from {i}", end = '\r')
    return out_dict
                
def merge_tree(min_lvl):
    max_lvl = min_lvl + 1
    merged = set()
    while True:
        old_len = len(G)
        all_nodes = [(f1, c1, d1, x1) for (f1, c1, d1, x1) in G if min_lvl <= d1 <= max_lvl and f1 != 0 and f1 != 0xFFFF]
        # total_pairs = len(all_nodes) * (len(all_nodes) - 1) // 2
        # costs_any = [d_any[func]['cost'] for func in range(0xFFFF+1)]
        # costs_xor = [d_xor[func]['cost'] for func in range(0xFFFF+1)]
        # depths_any = [d_any[func]['depth'] for func in range(0xFFFF+1)]
        # depths_xor = [d_xor[func]['depth'] for func in range(0xFFFF+1)]
        
        print(f'Gathering the candidate node pairs from {len(all_nodes)} nodes')
        all_parents = Parallel(n_jobs=40)(delayed(_process_row)(n1, i, all_nodes, max_lvl) for i, n1 in enumerate(all_nodes))
        
        # coros = [_one_iteration(n1, i, all_nodes, max_lvl) for i, n1 in enumerate(all_nodes)]
        # loop = asyncio.get_event_loop()
        # all_parents = loop.run_until_complete(asyncio.gather(*coros))
        
        parents = defaultdict(SortedList)
        for i, out_dict in enumerate(all_parents):
            print(f'\t Stitching candidate node pairs #{i+1:,} out of {len(all_parents):,}                              ', end = '\r') #)
            for func, lnodes in out_dict.items():
                parents[func].update(lnodes)
                
        
        '''
        for i, n1 in enumerate(all_nodes):
            print(f'\t Filtering node pairs #{i+1:,} out of {len(all_nodes):,}                              ', end = '\r') #)
            f1, c1, d1, x1 = n1
            for n2 in all_nodes[i+1:]:
                (f2, c2, d2, x2) = n2
                func = (f1 | f2)
                if func != f1 and func != f2: #the function should make sense
                    # the function is unknown or a known implementation has larger delay
                    # func is xorable AND (no known xorable implementation exists or known xorable has worse depth)
                    xorable = (f1 ^ f2) == func
                    preliminary_cost = c1 + c2
                    if xorable and (                                # func is xorable AND
                        (func not in d_xor) or                      # (no known xorable implementation exists OR
                        (d_xor[func]['depth'] > max_lvl or          #  (the known xorable implementation has worse depth OR
                         (d_xor[func]['depth'] == max_lvl and       #   it has the same depth AND
                          preliminary_cost <= d_xor[func]['cost'])  #   its cost is not better)
                         )                                          # )
                        ):
                        parents[func].add((preliminary_cost, xorable, n1, n2))
                    elif not xorable and (                          # func is not xorable AND
                        (func not in d_any) or                      # (no known implementation exists OR
                        (d_any[func]['depth'] > max_lvl or          #  (the known implementation has worse depth OR
                         (d_any[func]['depth'] == max_lvl and       #   it has the same depth AND
                          preliminary_cost <= d_any[func]['cost'])  #   its cost is not better)
                         )                                          # )
                        ):
                        parents[func].add((preliminary_cost, xorable, n1, n2))
            '''
            # TODO : add dfs tree to attr. Then, use induced_subgraph to count splitters
            # TODO : sort valid_nodes
        print('\n')
        # print({k:nlargest(3,v) for k,v in parents.items()})
        for i, (func, candidates) in enumerate(parents.items()):
            # using depth-based dicts to make sure funcs with larger cost but smaller depth still get reviewed
            best_cost_any = d_any[func]['cost'] if func in d_any else float('inf')
            best_cost_xor = d_xor[func]['cost'] if func in d_xor else float('inf')
            print(f'\t Finding best node pairs #{i+1:,} out of {len(parents):,}                              ', end = '\r') #)
            # TODO: Consider only creating a list of new_nodes and adding only the best instead of add_node right away
            for preliminary_cost, xorable, n1, n2 in candidates:
                # Not xorable but has better preliminary cost
                if not xorable and (preliminary_cost <= best_cost_any):
                    nspl, dupe_gate_cost = count_spl(n1, n2)
                    new_nspl = nspl - G.nodes[n1]['nspl'] - G.nodes[n2]['nspl']
                    cost = preliminary_cost  + COSTS['SPL'] * new_nspl + dupe_gate_cost
                    add_node((n1, n2), func, 'CB', cost, max_lvl, xorable, nspl)
                    best_cost_any = min(best_cost_any, cost)
                # Xorable and has better preliminary cost among xorables
                elif xorable and (preliminary_cost <= best_cost_xor):
                    nspl, dupe_gate_cost = count_spl(n1, n2)
                    new_nspl = nspl - G.nodes[n1]['nspl'] - G.nodes[n2]['nspl']
                    cost = preliminary_cost  + COSTS['SPL'] * new_nspl + dupe_gate_cost
                    add_node((n1, n2), func, 'CB', cost, max_lvl, xorable, nspl)
                    best_cost_xor = min(best_cost_any, cost)
                # best cost for xorables and non-xorables already found
                else:
                    break
        print(f'\n')
        if old_len == len(G):
            break
        logging.info('\n')
        
        
        # for i, (n1, n2) in enumerate(combinations(all_nodes, 2)):
        #     (f1, c1, d1, x1), (f2, c2, d2, x2) = n1, n2
        #     func = f1 | f2
        #     func_same = (func == f1 or func == f2)
        #     if func_same:
        #         continue
            
        #     func_seen = func in d_xor
        #     func_worse_depth = func_seen and max_lvl > d_xor[func]['depth']
            
        #     preliminary_cost = c1 + c2 + COSTS['MERGE']
        #     func_worse_cost  = func_seen and preliminary_cost > c_xor[func]['cost']
            
        #     if func_worse_depth and func_worse_cost:
        #         continue
            
        #     # logging.info(f'\t Merging node pair #{i+1:05d} out of {total_pairs}') #, end = '\r'
        #     print(f'\t Merging node pair #{i+1:,} out of {total_pairs:,}                              ', end = '\r') #)
        #     nspl1 = G.nodes[n1]['nspl']
        #     nspl2 = G.nodes[n2]['nspl']
        #     nspl, dupe_gate_cost = count_spl(n1, n2)
        #     new_nspl = nspl - nspl1 - nspl2
        #     xorable = (func == (f1 ^ f2))
        #     cost = preliminary_cost  + COSTS['SPL'] * new_nspl + dupe_gate_cost
        #     add_node((n1, n2), func, 'CB', cost, max_lvl, xorable, nspl)
        
def as_gates(min_src_lvl, max_src_lvl, tgt_lvl):
    all_nodes = [(n1, nspl) for n1, nspl in G.nodes(data = 'nspl') if min_src_lvl <= n1[2] <= max_src_lvl and n1[0] != 0 and n1[0] != 0xFFFF]
    
    for i, (n1, nspl1) in enumerate(all_nodes):
        f1, c1, d1, x1 = n1
        dff_cost = c1 + COSTS['DFF']
        not_func = ONES ^ f1
        not_cost = c1 + COSTS['NOT']
        add_node([n1], f1      , 'DFF', dff_cost, tgt_lvl, True, nspl1, force_add = True) # add DFF in any case to make sure the TT is available for next stage
        add_node([n1], not_func, 'NOT', not_cost, tgt_lvl, True, nspl1) #, force_add = True (No need to add inverter if better function exists)
        
        # logging.info(f'\t Analyzing node #{i+1:04d} out of {len(all_nodes)}') #, end = '\r'
        print(f'\t Analyzing node #{i+1:,} out of {len(all_nodes):,}                              ', end = '\r') #
        if not x1:
            continue
        
        for n2, nspl2 in all_nodes[(i+1):]:
            f2, c2, d2, x2 = n2
            if not x2: # both should be xorable
                continue
            
            func = f1 ^ f2
            func_seen        = func in d_xor
            func_worse_depth = func_seen and tgt_lvl > d_xor[func]['depth']
            
            preliminary_cost = c1 + c2 + COSTS['XOR']
            func_worse_cost  = func_seen and preliminary_cost > c_xor[func]['cost']
            
            if func_worse_depth and func_worse_cost:
                continue
            
            nspl, dupe_gate_cost = count_spl(n1, n2)
            new_nspl = nspl - nspl1 - nspl2
            cost = preliminary_cost + COSTS['SPL'] * new_nspl + dupe_gate_cost
            add_node([n1, n2], func, 'XOR', cost, tgt_lvl, True, nspl) #, force_add = True
        print('\n')
    logging.info('\n')
    
    # for n1, nspl in all_nodes:
    #     dff_func, c1, d1, x1 = n1
    #     dff_cost = c1 + COSTS['DFF']
    #     not_func = ONES ^ dff_func
    #     not_cost = c1 + COSTS['NOT']
    #     add_node([n1], dff_func, 'DFF', dff_cost, tgt_lvl, True, nspl, force_add = True)
    #     add_node([n1], not_func, 'NOT', not_cost, tgt_lvl, True, nspl, force_add = True)

    # # XOR gates
    # for (n1, nspl1), (n2, nspl2) in combinations(all_nodes, 2):
    #     f1, c1, d1, x1 = n1
    #     f2, c2, d2, x2 = n2
    #     if not (x1 and x2): # both should be xorable
    #         continue
    #     func = f1^f2
    #     nspl, dupe_gate_cost = count_spl(n1, n2)
    #     new_nspl = nspl - nspl1 - nspl2
    #     cost = c1 + c2 + COSTS['XOR'] + COSTS['SPL'] * new_nspl + dupe_gate_cost
    #     add_node([n1, n2], func, 'XOR', cost, tgt_lvl, True, nspl, force_add = True)
        
def sa_gates(src_lvl, tgt_lvl):
    all_nodes = [(n1, nspl) for n1, nspl in G.nodes(data = 'nspl') if src_lvl == n1[2] and n1[0] != 0 and n1[0] != 0xFFFF]
    total_pairs = len(all_nodes) * (len(all_nodes) - 1) // 2
    for i, ((n1, nspl1), (n2, nspl2)) in enumerate(combinations(all_nodes, 2)):
        (f1, c1, d1, x1), (f2, c2, d2, x2) = n1, n2
        and_func = f1 & f2
        or_func  = f1 | f2
        # breakpoint()
        and_func_same         = (and_func == f1 or and_func == f2)
        or_func_same          = ( or_func == f1 or  or_func == f2)
        if and_func_same and or_func_same:
            continue
            
        and_func_seen         = and_func in d_xor
        and_func_worse_depth  = and_func_seen and tgt_lvl > d_xor[and_func]['depth']
        or_func_seen          =  or_func in d_xor
        or_func_worse_depth   =  or_func_seen and tgt_lvl > d_xor[ or_func]['depth']
        
        preliminary_cost = c1 + c2 + COSTS['MERGE']
        and_func_worse_cost  = and_func_seen and preliminary_cost > c_xor[and_func]['cost']
        or_func_worse_cost   =  or_func_seen and preliminary_cost > c_xor[ or_func]['cost']
        
        if and_func_worse_depth and or_func_worse_depth and and_func_worse_cost and or_func_worse_cost:
            continue
        
        # logging.info(f'\t Processing AND/OR for node pair #{i+1:05d} out of {total_pairs}') #, end = '\r'
        print(f'\t Processing AND/OR for node pair #{i+1:,} out of {total_pairs:,}                              ', end = '\r')
        
        nspl, dupe_gate_cost = count_spl(n1, n2)
        new_nspl = nspl - nspl1 - nspl2
        cost = preliminary_cost + COSTS['SPL'] * new_nspl + dupe_gate_cost
        if not (and_func_worse_depth and and_func_worse_cost):
            add_node([n1, n2], and_func, 'AND', cost, tgt_lvl, True, nspl)
        if not ( or_func_worse_depth and  or_func_worse_cost):
            add_node([n1, n2],  or_func,  'OR', cost, tgt_lvl, True, nspl)

if True: 
    lvl = 0
    # for lvl in range(3, 16, 3):
    while True:
        logging.info(f'Processing level {lvl // 3}')
        # First merger tree
        merge_tree(lvl)
        logging.info(f'\t After merge tree :')
        logging.info(f'\t\t{len(G)} nodes')
        logging.info(f'\t\t{len(c_any)} truth tables')
        logging.info(f'\t\t{len(c_xor)} xorable truth tables')
        if len(c_xor) == (ONES + 1):
            break
        # DFF/NOT gates
        as_gates(lvl, lvl+1, lvl+2)
        logging.info(f'\t After AS gates :')
        logging.info(f'\t\t{len(G)} nodes')
        logging.info(f'\t\t{len(c_any)} truth tables')
        logging.info(f'\t\t{len(c_xor)} xorable truth tables')
        if len(c_xor) == (ONES + 1):
            break
        # AND/OR gates
        sa_gates(lvl+2, lvl+3)
        logging.info(f'\t After SA gates :')
        logging.info(f'\t\t{len(G)} nodes')
        logging.info(f'\t\t{len(c_any)} truth tables')
        logging.info(f'\t\t{len(c_xor)} xorable truth tables')
        if len(c_xor) == (ONES + 1):
            break
        lvl += 3
else:
    merge_tree(0)
    logging.info(f'\t After merge tree :')
    logging.info(f'\t\t{len(G)} nodes')
    logging.info(f'\t\t{len(c_any)} truth tables')
    logging.info(f'\t\t{len(c_xor)} xorable truth tables')

    min_src_lvl = 0
    max_src_lvl = 1
    tgt_lvl = 2
    all_nodes = [(n1, nspl) for n1, nspl in G.nodes(data = 'nspl') if min_src_lvl <= n1[2] <= max_src_lvl]
    for n1, nspl in all_nodes:
        dff_func, c1, d1, x1 = n1
        dff_cost = c1 + COSTS['DFF']
        not_func = ONES ^ dff_func
        not_cost = c1 + COSTS['NOT']
        add_node([n1], dff_func, 'DFF', dff_cost, tgt_lvl, True, nspl, force_add = True)
        add_node([n1], not_func, 'NOT', not_cost, tgt_lvl, True, nspl, force_add = True)

    # XOR gates
    for (n1, nspl1), (n2, nspl2) in combinations(all_nodes, 2):
        f1, c1, d1, x1 = n1
        f2, c2, d2, x2 = n2
        if not (x1 and x2): # both should be xorable
            continue
        func = f1^f2
        nspl, dupe_gate_cost = count_spl(n1, n2)
        new_nspl = nspl - nspl1 - nspl2
        cost = c1 + c2 + COSTS['XOR'] + COSTS['SPL'] * new_nspl + dupe_gate_cost
        add_node([n1, n2], func, 'XOR', cost, tgt_lvl, True, nspl, force_add = True)
        
    logging.info(f'\t After AS gates :')
    logging.info(f'\t\t{len(G)} nodes')
    logging.info(f'\t\t{len(c_any)} truth tables')
    logging.info(f'\t\t{len(c_xor)} xorable truth tables')

    # AND/OR gates
    src_lvl = 2
    tgt_lvl = 3
    all_nodes = [(n1, nspl) for n1, nspl in G.nodes(data = 'nspl') if src_lvl == n1[2]]

    for (n1, nspl1), (n2, nspl2) in combinations(all_nodes, 2):
        f1, c1, d1, x1 = n1
        f2, c2, d2, x2 = n2
        nspl, dupe_gate_cost = count_spl(n1, n2)
        new_nspl = nspl - nspl1 - nspl2
        cost = c1 + c2 + COSTS['MERGE'] + COSTS['SPL'] * new_nspl + dupe_gate_cost
        
        and_func = f1 & f2
        or_func = f1 | f2
        add_node([n1, n2], and_func, 'AND', cost, tgt_lvl, True, nspl)
        add_node([n1, n2],  or_func,  'OR', cost, tgt_lvl, True, nspl)
        
    logging.info(f'\t After SA gates :')
    logging.info(f'\t\t{len(G)} nodes')
    logging.info(f'\t\t{len(c_any)} truth tables')
    logging.info(f'\t\t{len(c_xor)} xorable truth tables')
        

        
    # logging.info(list(G.nodes(data = True)))
    
        
    


        



from ortools.sat.python import cp_model
from collections import defaultdict
import sys
import os
from pprint import pprint
# import builtins

# LOG_FILE = "ilp_log.log"

# def print(*args, **kwargs):
#     with open("ilp_log.log", "w") as f:
#         f.write()

EPS = 1e-6

TYPES = { "0": "PI", 
          "1": "AA", 
          "2": "AS", 
          "3": "SA", 
          "4": "FA"}

# TYPES = {"0": "AA", "1": "AS", "2": "SA"}

class T1_cell:
  __slots__ = ['xor3', 'maj3', 'min3', 'or3', 'nor3', 'inputs3']
  def __init__(self, _outputs: dict[str, int], _inputs: list[tuple[int, bool]] ) -> None:
    self.outputs = _outputs
    self.inputs = _inputs
    

class Primitive:
  __slots__ = ['sig', 'type', 'fanins', 'siblings']
  def __init__(self, _sig: int, _type : str, _fanins : list, _siblings : dict = {}, _in_neg : list) -> None:
    self.sig      = _sig
    self.type     = _type if _fanins else 'PI'
    self.fanins   = _fanins
    self.siblings = _siblings
    self.in_neg = _in_neg

  def __repr__(self) -> str:
    if self.siblings:
      return f"Primitive(sig={self.sig}, type='{self.type}', fanins={self.fanins}, siblings={self.siblings}"
    else:
      return f"Primitive(sig={self.sig}, type='{self.type}', fanins={self.fanins}"
      
def parse_specs(filename: str) -> tuple[dict[int, Primitive], list[dict[str, int]]]:
  items = {}
  t1_cells = []
  with open(filename, 'r') as f:
    for line in f:
      line = line.strip()
      
      _sig_str, _type_str, _fanin_str = line.split(',')
      _type = TYPES[_type_str]
      
      if "|" in _sig_str:
        assert(_type_str == "4")
        _fanins = []
        _in_neg = []
        for _fanin_sig in _fanin_str.split('|'):
          if _fanin_sig.startswith('~'):
            _fanins.append(int(_fanin_sig[1:]))
            _in_neg.append(True)
          else:
            _fanins.append(int(_fanin_sig))
            _in_neg.append(False)
        _siblings = {}
        for _key, _sig in zip( ['XOR','MAJ','MIN','OR','NOR'],  _sig_str.split('|')):
          if _sig != '0':
            _siblings[_key] = (int(_sig), True)
            items[g.sig] = Primitive(_sig, _type, _fanins, _in_neg)
        assert('XOR' in _siblings)
        t1_cells.append(_siblings)
      else:
        assert(_type_str != "4")
        _fanins  = [int(_fanin_sig) for _fanin_sig in _fanin_str.split('|') if _fanin_sig]
        _sig  = int(_sig_str)
        _type = TYPES[_type_str]
        items[g.sig] = Primitive(_sig, TYPES[_type_str], _fanins)
  return items, t1_cells

if __name__ == "__main__":
  
  # Model.
  model = cp_model.CpModel()
  
  NPH = int(sys.argv[1])
  cfg_path = sys.argv[2]
  
  
  # Split the filename and its extension
  base_name, _ = os.path.splitext(cfg_path)
  # Create the new filename with the ".log" extension
  log_filename = base_name + ".log"
  # Save the original sys.stdout
  original_stdout = sys.stdout
  # Open the log file in write mode
  log_file = open(log_filename, 'w')
  # Redirect sys.stdout to the log file
  sys.stdout = log_file

  print(f"Parsing specs file {cfg_path}")
  all_signals, all_t1_cells = parse_specs(cfg_path)
  # print(all_signals)
  
  print(f'Finished creating [all_signals]')
  pprint(all_signals)
  
  max_phase = max(g.phase for g in all_signals.values() if g.phase is not None)
  
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
  pprint(Sigma)
  
  fanin_diff = {}
  Delta = {}
  # IMPORTANT: PIs are assumed to be at clock stage 0 but with an arbitrary phase
  expr = []
  already_processed = []
  for g in all_signals.values():
    # print(f"The type of the gate {g.sig} is {g.type}")
    if g.type == 'AA':
      # any PI is the earliest fanin, with zero phase and zero clock stage
      # IMPORTANT: only 2-input CB is supported
      assert(len(g.fanins) == 2)
      a_sig, b_sig = g.fanins
      a = all_signals[a_sig]
      b = all_signals[b_sig]
      
      max_diff = model.NewIntVar(*sigma_bounds, f'max_{a.sig}_{b.sig}')
      model.AddMaxEquality(max_diff, [ Sigma[g.sig] - Sigma[a.sig], 
                      Sigma[g.sig] - Sigma[b.sig] ])
      
      div_val = model.NewIntVar(*sigma_bounds, f'div_{a.sig}_{b.sig}')
      model.AddDivisionEquality(div_val, max_diff, NPH)
      expr.append( div_val )
      
      # print(f"AA constraint: {Sigma[a.sig].Name()} <= {Sigma[g.sig].Name()}")
      # print(f"AA constraint: {Sigma[b.sig].Name()} <= {Sigma[g.sig].Name()}")
      
      model.Add( Sigma[a.sig] <= Sigma[g.sig] )
      model.Add( Sigma[b.sig] <= Sigma[g.sig] )
      
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
        expr.append( div_val )
        
    elif g.type == 'SA':
      
      for i, p_sig in enumerate(g.fanins):
        p = all_signals[p_sig]
        if (p.type == 'AS') and (fanout[p.sig] == 1):
          model.Add(Sigma[p_sig] <= Sigma[g.sig])
          # print(f"SA constraint: {Sigma[p_sig].Name()} <= {Sigma[g.sig].Name()} ({p.type}) ({fanout[p.sig]})")
        else: 
          model.Add(Sigma[p_sig] <  Sigma[g.sig])
          # print(f"SA constraint: {Sigma[p_sig].Name()} <  {Sigma[g.sig].Name()} ({p.type}) ({fanout[p.sig]})")
          
        # expr.append( Sigma[g.sig] - Sigma[p_sig] )
          
        delta = model.NewIntVar(*sigma_bounds, f'delta_sa_{p_sig}_{g.sig}')
        model.Add(delta == (Sigma[g.sig] - Sigma[p_sig] + NPH - 1))
        
        div_val = model.NewIntVar(*sigma_bounds, f'div_{p_sig}_{g.sig}')
        model.AddDivisionEquality(div_val, delta, NPH)
        expr.append( div_val )
        
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
        
        
      assert(len(g.fanins) == 3)
      # sig_idx = [p_sig for p_sig in g.fanins]
        # TODO : check fanins
        # TODO :    if the fanin is AA, add DFF/NOT to cost and [++SIGMA_EFF] (not accurate but too complex otherwise)
        # TODO :    if the fanin is AS/SA and negated, [++SIGMA_EFF]
      incr_Sigma = [0, 0, 0]
      for i, (p_sig, p_neg) in enumerate(zip(g.fanins, g.in_neg)):
        p = all_signals[p_sig]
        if (p.type == 'AA') or (p.type in ('AS', 'SA') and p.p_neg):
          incr_Sigma[i] += 1
      
        # TODO : create vars [ABS_MIN=min(A,B,C)], [MED=median(A,B,C)-ABS_MIN], [MAX=max(A,B,C)-ABS_MIN]
      ABS_MIN = model.NewIntVar(*sigma_bounds, f'absmin_{g.sig}')
      model.AddMinEquality(ABS_MIN, [(Sigma[p_sig]+incr_Sigma[i]) for i, p_sig in enumerate(g.fanins)])
      
      ABS_MAX = model.NewIntVar(*sigma_bounds, f'absmax_{g.sig}')
      model.AddMaxEquality(ABS_MAX, [(Sigma[p_sig]+incr_Sigma[i]) for i, p_sig in enumerate(g.fanins)])
      
      ABS_MED = model.NewIntVar(*sigma_bounds, f'absmed_{g.sig}')
      model.Add(ABS_MED == (Sigma[g.fanins[0]]+incr_Sigma[0]) + (Sigma[g.fanins[1]]+incr_Sigma[1]) + (Sigma[g.fanins[2]]+incr_Sigma[2]) - ABS_MIN - ABS_MAX)
      
      # IMPORTANT: this is the stage of the XOR cell only!!!
      # IMPORTANT: TODO: ADD the stages of other outputs
      T1_XOR_SIGMA        = model.NewIntVar(*sigma_bounds, f'sigma_xor_{g.sig}')
      # T1_CARRY_SIGMA      = model.NewIntVar(*sigma_bounds, f'sigma_carry_{g.sig}')
      # T1_CBAR_SIGMA       = model.NewIntVar(*sigma_bounds, f'sigma_cbar_{g.sig}')
      # T1_CARRY_NOT_SIGMA  = model.NewIntVar(*sigma_bounds, f'sigma_carry_not_{g.sig}')
      # T1_CBAR_NOT_SIGMA   = model.NewIntVar(*sigma_bounds, f'sigma_cbar_not_{g.sig}')
      model.Add(T1_XOR_SIGMA > ABS_MIN + 2)
      model.Add(T1_XOR_SIGMA > ABS_MED + 1)
      model.Add(T1_XOR_SIGMA > ABS_MAX    )
      
      MED_MIN_diff = model.NewIntVar(*sigma_bounds, f'mod_med_min_diff_{g.sig}')
      model.AddModuloEquality(MED_MIN_diff, ABS_MED - ABS_MIN, NPH)
      
      MAX_MED_diff = model.NewIntVar(*sigma_bounds, f'mod_max_med_diff_{g.sig}')
      model.AddModuloEquality(MAX_MED_diff, ABS_MAX - ABS_MED, NPH)
      
      MAX_MIN_diff = model.NewIntVar(*sigma_bounds, f'mod_max_min_diff_{g.sig}')
      model.AddModuloEquality(MAX_MIN_diff, ABS_MAX - ABS_MIN, NPH)
      
      FLOOR_DIFF_MIN = model.NewIntVar(*sigma_bounds, f'floor_diff_min_{g.sig}')
      model.AddDivisionEquality(FLOOR_DIFF_MIN, T1_XOR_SIGMA-ABS_MIN-1, NPH)

      FLOOR_DIFF_MED = model.NewIntVar(*sigma_bounds, f'floor_diff_med_{g.sig}')
      model.AddDivisionEquality(FLOOR_DIFF_MED, T1_XOR_SIGMA-ABS_MED-1, NPH)

      FLOOR_DIFF_MAX = model.NewIntVar(*sigma_bounds, f'floor_diff_max_{g.sig}')
      model.AddDivisionEquality(FLOOR_DIFF_MAX, T1_XOR_SIGMA-ABS_MAX-1, NPH)
      
      MIN_FLAG = model.NewBoolVar(f'min_flag_{g.sig}')
      model.AddBoolAnd(MED_MIN_diff == 0, FLOOR_DIFF_MIN == 0).OnlyEnforceIf(MIN_FLAG)
      
      MED_FLAG = model.NewBoolVar(f'med_flag_{g.sig}')
      model.AddBoolAnd(MAX_MED_diff == 0, FLOOR_DIFF_MED == 0).OnlyEnforceIf(MED_FLAG)
      
      # place 
      for g_sig in sibling_list.values():
        model.Add(Sigma[g_sig] == T1_XOR_SIGMA)
      
      # Cost function
      expr.append( FLOOR_DIFF_MIN )
      expr.append( FLOOR_DIFF_MED )
      expr.append( FLOOR_DIFF_MAX )
      expr.append( MIN_FLAG )
      expr.append( MED_FLAG )
      
      already_processed.extend(sibling_list.values())
      
  model.Minimize(sum(expr))

  # Solves and prints out the solution.
  solver = cp_model.CpSolver()
  print(f'Starting Macro ILP')
  solver.parameters.max_time_in_seconds = 300.0
  status = solver.Solve(model)
  print(f'Solve status: {solver.StatusName(status)}')
  if (solver.StatusName(status) in ("OPTIMAL", "FEASIBLE")):
    print(f'Objective value: {solver.ObjectiveValue()}')
    for sig, sigma in Sigma.items():
      g = all_signals[sig]
      print(f'{sigma.Name()[6:]} :{solver.Value(sigma)}')
      
  # Restore sys.stdout to the original value
  sys.stdout = original_stdout
  # Close the log file
  log_file.close()
  
  print(f'Solve status: {solver.StatusName(status)}')
  if (solver.StatusName(status) in ("OPTIMAL", "FEASIBLE")):
    print(f'Objective value: {solver.ObjectiveValue()}')
    for sig, sigma in Sigma.items():
      g = all_signals[sig]
      print(f'{sigma.Name()[6:]} :{solver.Value(sigma)}')

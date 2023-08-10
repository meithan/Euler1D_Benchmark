# Utility to automaticlly run all benchmarks
from subprocess import run

# ----------------------------------------------------

def run_command(cmd, echo=False, check=True):
  if echo: print(cmd)
  result = run(cmd, capture_output=True, check=True, shell=True, text=True)
  return result

def run_benchmark(setup_cmds, test_cmd, num_runs):
  for cmd in setup_cmds:
    run_command(cmd)
  times = []
  for i in range(num_runs):
    result = run_command(test_cmd)
    elapsed = float(result.stdout.strip())
    print("{:.3f} ".format(elapsed), end="", flush=True)
    times.append(elapsed)
  print()
  return times

# ----------------------------------------------------

# Define (single) benchmarks here
# (name, [commands to run before the tests eg compilation], command to run test)
# The only output of the test must be the execution time
benchmarks = [
  ("C/C++", ["rm -f Euler1D_C", "g++ -O3 Euler1D.cpp -o Euler1D_C"], "./Euler1D_C"),
  ("Fortran", ["rm -f Euler1D_F", "gfortran -O3 Euler1D.f90 -o Euler1D_F"], "./Euler1D_F"),
  ("Java", ["rm -f Euler1D.class", "javac Euler1D.java"], "java Euler1D"),
  ("Python3.8", [], "python3.8 Euler1D.py"),
  ("Python3.8+Numpy", [], "python3.8 Euler1D_numpy.py"),
  ("Python3.11", [], "python3.11 Euler1D.py"),
  ("Python3.11+Numpy", [], "python3.11 Euler1D_numpy.py"),
  ("Julia", [], "julia Euler1D.jl"),
  # Julia+LoopVec run manually in Julia REPL
  ("Rust", ["rustc -C opt-level=3 -o Euler1D_Rust Euler1D.rs"], "./Euler1D_Rust"),
]

num_runs = 10

for name, setup_cmds, test_cmd in benchmarks:
  print("> {}".format(name))
  times = run_benchmark(setup_cmds, test_cmd, num_runs)
  # print(",".join([name]+["{:.3f}".format(x) for x in times]))

# ----------------------------------------------------

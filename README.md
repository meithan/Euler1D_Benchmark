# Euler1D Benchmark
A comparison of various programming languages solving a 1D hydrodynamics problem.

The programs implement a simple finite-difference solver (the [Lax-Friedrichs method](https://en.wikipedia.org/wiki/Lax%E2%80%93Friedrichs_method)) for the 1D [Euler equations](https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics)) (inviscid, compressible hydrodynamics) in various languages to compare execution speed and syntax.

The benchmark is to solve the [Sod Shock Tube](https://en.wikipedia.org/wiki/Sod_shock_tube), with NX = 5000 points and a CFL parameter of 0.9. The test's normal disk output is used to verify correctness, but is skipped for the benchmarks.

Execution time is measured internally by each program by comparing the platform's wall clock or CPU clock at the start and end of the main program, and printed to the terminal, in seconds, as the sole program output. For `Julia+LoopVec` the code was run manually within the Julia REPL after the running `using LoopVectorization` and `include("Euler1D_opt.jl)` (so as to compile/optimize as much as possible before the tests are run), but the program still measures its own run time internally; the results of doing this are in line with those obtained with the `@benchmark` macro from [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl). Ten runs are carried out with each problem implementation to get an idea of the inherent variability (which is generally found to be small).

Current implementations are:

- C/C++
- Fortran 90 
- Java
- Python, using native nested lists and `for` loops
- Python with Numpy arrays and vectorized operations
- Julia, native
- Julia, leveraging the [LoopVectorizations.jl](https://github.com/JuliaCI/BenchmarkTools.jl) library
- Rust

The tests were executed on a personal computer with an Intel Core i7-9700F processor running Gentoo Linux with kernel 6.1.31. Compiler/interpreter versions were gcc 12.3.1 for GCC (C/C++ and Fortran 90), OpenJDK 17.0.6 for Java, CPython 3.8.17 and 3.11.4 for Python (with Numpy 1.24.4), 1.8.5 for Julia and 1.69.1 for Rust.

Optimization flags were used where available (e.g. `O3` for gcc and `opt-level=3` for Rust); see the `run_benchmarks.py` script for details. The [LoopVectorizations.jl](https://github.com/JuliaCI/BenchmarkTools.jl) library was used for the `Julia+LoopVec` benchmark to vectorize/optimize the main loops (by simply prepending the `@turbo` macro); thanks to Luis Arcos ([LAlbertoA](https://github.com/LAlbertoA)) for the tip! The two Python benchmarks were executed with both Python 3.8 and 3.11; it turns out that 3.11 is substantially faster!

It was my first time writing Julia and Rust code so those implementations might be a bit rough. If you know how to further optimize any of them to make the comparison more fair, please let me know!

## Results

The results of the benchmark are presented in the following plots. The bar heights are the average execution time averaged over 10 runs for each case, normalized to the execution time of the fasest benchmark implementation (so far, Fortran). In the first plot the native Python implementations are shown on a separate scale.

#### **Linear** Y scale
![Lin scale](https://github.com/meithan/Euler1D_Benchmark/blob/main/benchmark_lin.png)


#### **Logarithmic** Y scale
![Log scale](https://github.com/meithan/Euler1D_Benchmark/blob/main/benchmark_log.png)


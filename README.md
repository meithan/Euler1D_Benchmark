# Euler1D Benchmark
A comparison of various programming languages solving a 1D hydrodynamics problem.

The programs implement a simple finite-difference solver (the [Lax-Friedrichs method](https://en.wikipedia.org/wiki/Lax%E2%80%93Friedrichs_method)) for the 1D [Euler equations](https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics)) (inviscid, compressible hydrodynamics) in various languages to compare execution speed and syntax.

The benchmark is to solve the [Sod Shock Tube](https://en.wikipedia.org/wiki/Sod_shock_tube), with NX = 5000 points and a CFL parameter of 0.9. Ten runs, without disk output, are carried out with each version and the execution time measured with the `time` command.

Current implementations are:

- C
- Fortran (90)
- Python, using native nested lists and `for` loops
- Python with Numpy arrays and vectorized operations
- Julia
- Rust

## Results

The results of the benchmark are presented in the following plots:

![Lin scale](https://github.com/meithan/Euler1D_Benchmark/blob/main/benchmark_lin.png)
(**Linear** Y scale)


![Log scale](https://github.com/meithan/Euler1D_Benchmark/blob/main/benchmark_log.png)
(**Logarithmic** Y scale)

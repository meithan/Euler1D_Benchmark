# Euler1D Benchmark
A comparison of various programming languages solving a 1D hydrodynamics problem.

The programs implement a simple finite-difference solver (the [Lax-Friedrichs method](https://en.wikipedia.org/wiki/Lax%E2%80%93Friedrichs_method)) for the 1D [Euler equations](https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics)) (inviscid, compressible hydrodynamics) in various languages to compare execution speed and syntax.

The benchmark is to solve the [Sod Shock Tube](https://en.wikipedia.org/wiki/Sod_shock_tube), with NX = 5000 points and a CFL parameter of 0.9. Ten runs, without disk output, are carried out with each version and the execution time measured with the `time` command.

Current implementations are:

- C/C++ (gcc 12.2.1)
- Fortran 90 (gcc 12.2.1)
- Python, using native nested lists and `for` loops (CPython 3.8 and 3.11)
- Python with Numpy arrays and vectorized operations (CPython 3.8 and 3.11)
- Julia (1.8.5)
- Rust (1.66.1)

Python 3.8 and 3.11 were also tried for both Python implementations; 3.11 is significantly faster! it's my first time writing Julia and Rust code, so if you know how to further optimize any of them, please let me know!

## Results

The results of the benchmark are presented in the following plots. The bars show the execution time averaged over 10 runs (on a computer with an Intel Core i7-9700F running Gentoo Linux with kernel 6.1.28)

#### **Linear** Y scale
![Lin scale](https://github.com/meithan/Euler1D_Benchmark/blob/main/benchmark_lin.png)


#### **Logarithmic** Y scale
![Log scale](https://github.com/meithan/Euler1D_Benchmark/blob/main/benchmark_log.png)


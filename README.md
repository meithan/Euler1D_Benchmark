# Euler1D Benchmark
A comparison of various programming languages solving a 1D hydrodynamics problem.

The programs implement a simple finite-difference solver (the [Lax-Friedrichs method](https://en.wikipedia.org/wiki/Lax%E2%80%93Friedrichs_method)) for the 1D [Euler equations](https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics)) (inviscid, compressible hydrodynamics) in various languages to compare execution speed and syntax.

The benchmark is to solve the [Sod Shock Tube](https://en.wikipedia.org/wiki/Sod_shock_tube), with NX = 5000 points and a CFL parameter of 0.9. Ten runs, without disk output, are carried out with each version and the execution time measured with the `time` command.

Current implementations are:

- C/C++
- Fortran 90 
- Python, using native nested lists and `for` loops
- Python with Numpy arrays and vectorized operations
- Julia
- Rust

The tests were executed on a personal computer with an Intel Core i7-9700F processor running Gentoo Linux with kernel 6.1.28. Compiler/interpreter versions were gcc 12.2.1 for C/C++ and Fortran 90, CPython 3.8 and 3.11 for Python, 1.8.5 for Julia and 1.66.1 for Rust. Optimization flags were used where available (`O3` for gcc and `opt-level=3` for Rust).

Python 3.8 and 3.11 were tried for both Python implementations (3.11 is significantly faster!). It was my first time writing Julia and Rust code so those implementations might be a bit rough. If you know how to further optimize any of them to make the comparison more fair, please let me know!

## Results

The results of the benchmark are presented in the following plots. The bar heights are the average execution times over 10 runs in each case, while the numbers above them are execution time ratios vs the C/C++ implementation.

#### **Linear** Y scale
![Lin scale](https://github.com/meithan/Euler1D_Benchmark/blob/main/benchmark_lin.png)


#### **Logarithmic** Y scale
![Log scale](https://github.com/meithan/Euler1D_Benchmark/blob/main/benchmark_log.png)


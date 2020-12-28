```@meta
CurrentModule = Optini
```

# Optini

Optini is a prototype numerical optimization package for pedagogical and research purposes. 
It primarily references Nocedal & Wright's [*Numerical Optimization*](https://www.springer.com/gp/book/9780387303031)
and Kochenderfer & Wheeler's [*Algorithms for Optimization*](https://mitpress.mit.edu/books/algorithms-optimization).
The main purpose of Optini is to provide an extensible but generic interface for non-linear 
optimization algorithms and simplify the prototype process for new methods. As such, it is 
a useful scratchpad for experimental algorithms, and should not be considered production-ready.
For state-of-the-art Julia libraries in this space, we recommend [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl)
and [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl) instead. 

## Installation

Optini is only compatible with Julia 1.5 and above. If you would like to experiment or 
contribute to the package, you can install it from your Julia REPL:
```
]add https://github.com/lhnguyen-vn/Optini.jl
```
# Univariate Optimization

Optini provides a collection of popular bracketing and root-finding algorithms including 
[`Bisection`](@ref), [`Fibonacci`](@ref), [`GoldenSection`](@ref), and [`QuadraticFit`](@ref).
The [`optimize`](@ref) call in this case requires the lower and upper bounds of the search 
interval, for which Optini provides a convenient [`bracket`](@ref) function to find an 
initial bracket around the minimum. 

## References

```@docs
optimize(::Function, ::Real, ::Real)
bracket
Bisection
Fibonacci
GoldenSection
QuadraticFit
```
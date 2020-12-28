# Tutorials

Optini is designed to be easy to use, with a minimal API that allows you to try out 
different algorithms. Here you can find a few examples on how to get started with the 
package.

## Rosenbrock Function

In this tutorial, we use the Rosenbrock function from $\mathbb{R}^2$ to $\mathbb{R}$ as the 
objective to minimize. Rosenbrock has only one global minimum, but it is located in a narrow
valley that resembles the shape of a banana. Finding this trench is quite trivial, but it 
can be rather difficult to converge to the minimum as the bottom of the valley is flat. We 
first define the function as follows:

```jldoctest rosenbrock
julia> rosenbrock(x; a=1, b=100) = (a - x[1])^2 + b*(x[2] - x[1]^2)^2
rosenbrock (generic function with 1 method)
```

We then select an appropriate algorithm. At the moment, Optini provides support for generic
line search and trust region methods, with convenient aliases for popular algorithms. In 
this instance, we will use the famous gradient descent method:

```jldoctest rosenbrock
julia> using Optini

julia> alg = GradientDescent()
Line Search Algorithm:
  • Direction: Steepest()
  • Initial step length guess: StaticInitial{Float64}(0.001)
  • Step length optimization method: BacktrackingLineSearch{Float64,Float64}(0.0001, 0.5)
```

Optini informs us that the default `GradientDescent` is a `LineSearch` method using the 
`Steepest` direction, a static initial step length `0.001`, and `BacktrackingLineSearch` to 
search for the appropriate step. All of these parameters are modular, and you can refer to 
the documentation on [Line Search Algorithms](@ref) to select the appropriate methods or 
customize your own.

Once we have defined the algorithm, the other necessary ingredient is an initial point to 
kickstart our process. Here, we choose the origin:

```jldoctest rosenbrock
julia> x = zeros(2)
2-element Array{Float64,1}:
 0.0
 0.0
```

Finally, we simply call the `optimize` function:

```jldoctest rosenbrock
julia> optimize(rosenbrock, x; alg)
Solution Summary:
  • Converged: true
  • Total Iterations: 65897
  • Minimizer: [0.9999999999988823, 0.9999999999977602]
  • Minimum: 1.2511394230708787e-24
```

Considering the analytical minimum is `0` attained at `[1, 1]`, the solution is by and large
accurate. The `optimize` interface also includes a few solver options, such as 
absolute tolerance for the gradient norm or the maximum number of iterations. For instance, 
if we want to save the trace of the solver:

```jldoctest rosenbrock
julia> sol = optimize(rosenbrock, x; alg, savetrace=true);

julia> sol.metadata[:trace][1:20_000:end]
4-element Array{Optini.FirstOrderState{Array{Float64,1},Float64},1}:
 Optini.FirstOrderState{Array{Float64,1},Float64}([0.0, 0.0], 1.0, [-2.0, 0.0])
 Optini.FirstOrderState{Array{Float64,1},Float64}([0.9998977748147321, 0.9997951509901917], 1.0466723905409629e-8, [-4.083139381724218e-5, -8.181785220440219e-5])
 Optini.FirstOrderState{Array{Float64,1},Float64}([0.9999999653270628, 0.9999999305153797], 1.2041376497141562e-15, [-1.384702634028865e-8, -2.774942498007249e-8])
 Optini.FirstOrderState{Array{Float64,1},Float64}([0.9999999999882364, 0.9999999999764259], 1.3860259910853339e-22, [-4.742206627604874e-12, -9.392486788328824e-12])
```
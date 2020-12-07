@testset "Univariate Bracketing" begin
    f(x) = 2x
    @test bracket(f) === nothing

    g(x) = -(x + sin(x))*exp(-x^2)
    br = bracket(g)
    @test br !== nothing
    @test br[1] < 0.67956 < br[2]
end

@testset "Univariate Optimization" begin
    f(x) = x ≤ 3 ? (x-2.0)^2 : 2*log(x-2.0) + 1
    for alg in [Fibonacci(), GoldenSection(), Bisection(), QuadraticFit()]
        sol = optimize(f, 0f0, 6f0; alg, maxiter=40)
        # Test convergence
        @test sol.converged
        # Test result
        @test sol.minimizer ≈ 2
        # Test type stability
        @test typeof(sol.minimizer) == Float32
    end
    # Test with integer input
    sol = optimize(f, 0, 6; maxiter=40)
    @test sol.converged
    @test sol.minimizer ≈ 2
    @test typeof(sol.minimizer) == Float64
end
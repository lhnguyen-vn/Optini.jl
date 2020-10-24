@testset "Univariate Bracketing" begin
    f(x) = 2x
    @test bracket(f) === nothing

    g(x) = -(x + sin(x))*exp(-x^2)
    br = bracket(g)
    @test br !== nothing
    @test br[1] < 0.67956 < br[2]
end

@testset "Univariate Optimization" begin
    f(x) = x ≤ 3 ? (x-2)^2 : 2*log(x-2) + 1
    for alg in [Fibonacci(), GoldenSection(), Bisection(), QuadraticFit()]
        sol = optimize(f, 0, 6; alg)
        @test sol.converged
        @test sol.minimizer ≈ 2
    end
end
let
    f(x) = (x[1] - 1)^2 + 2(x[2] - 1)^2
    x = [0.0, 0.0]
    state = Optini.FirstOrderState(x, f(x), Zygote.gradient(f, x)[1])
    p = -state.∇f

    @testset "Line Search Initial Step Length" begin
        inits = [
            StaticInitial(1f0), 
            PreviousDecreaseInitial(1f0), 
            QuadraticInitial(1f0)
        ]

        for init in inits
            α = init(state, p) 
            @test typeof(α) == Float32 # type stability
            @test α == 1f0 # first iteration
            Optini.update!(init, state, p, 0.1f0)
        end

        next_x = [0.2, 0.4]
        next_state = Optini.FirstOrderState(next_x, f(next_x), Zygote.gradient(f, next_x)[1])
        next_p = -next_state.∇f
        steps = [1f0, 0.24038461f0, 0.39817306f0] # expected initial step lengths

        for (init, step) in zip(inits, steps)
            α = init(next_state, next_p)
            @test typeof(α) == Float32 # type stability
            @test α == step # check initial step
        end
    end

    @testset "Line Search Algorithms" begin
        init = StaticInitial(1f0) # initial step length
        linesearches = [
            StaticLineSearch(; init),
            ExactLineSearch(; init),
            BacktrackingLineSearch(; init),
            InterpolationLineSearch(; init),
            StrongWolfeLineSearch(; init)
        ]

        steps = [1f0, 0.2777777f0, 0.5f0, 0.2777778f0, 0.5f0]

        for (ls, step) in zip(linesearches, steps)
            α = ls(f, state, p)
            @test typeof(α) == Float32 # type stability
            @test α == step # check step length
        end
    end
end
let
    f(x) = x⋅x
    g(x) = 2x
    x = [1.0, -1.0]
    state = Optini.FirstOrderState(x, f(x), g(x))
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

        next_x = [0.8, -0.8]
        next_state = Optini.FirstOrderState(next_x, f(next_x), g(next_x))
        next_p = -next_state.∇f
        steps = [1f0, 0.15625f0, 0.2840625f0] # expected initial step lengths

        for (init, step) in zip(inits, steps)
            α = init(next_state, next_p)
            @test typeof(α) == Float32 # type stability
            @test α == step # check initial step
        end
    end

    @testset "Line Search Methods" begin
        linesearches = [
            StaticLineSearch(),
            ExactLineSearch(reltol=1e-7),
            BacktrackingLineSearch(),
            InterpolationLineSearch(),
            StrongWolfeLineSearch()
        ]

        for (iter, ls) in enumerate(linesearches)
            α = ls(f, state, p, 1f0)
            @test typeof(α) == Float32 # type stability
            if iter == 1
                @test α == 1f0 # check step length
            else
                @test α == 0.5f0 # check step length
            end
        end
    end
end
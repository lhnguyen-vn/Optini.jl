"""
    StaticLineSearch

`StaticLineSearch` directly uses the initial guess as the new step length.
"""
struct StaticLineSearch end

(sls::StaticLineSearch)(f, state, p, α₀) = α₀
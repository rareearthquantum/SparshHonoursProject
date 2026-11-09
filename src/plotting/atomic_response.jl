"""
    plot_0sd_polarisation(result)

Return a plot of `abs2.(result.P)` against time for a local atomic run.
"""
function plot_0sd_polarisation(result)
    polarisation_intensity = abs2.(result.P)

    t_vec = result.t_vec
    #title = (title==:default) ? "" : title
    fig = plot(; xlabel="\$t\$", ylabel="\$|\\mathcal{P}|^2\$",
        title="Squared magnitude of the polarisation density against time",
        size=PLOT_PANEL_SIZE)

    plot!(t_vec, polarisation_intensity; label=false)

    return fig
end

"""
    plot_0sd_pol_complex_parts(result)

Return a plot of the ensemble-averaged Bloch components
`u = 2real(P)` and `v = -2imag(P)` against time for a local run.
"""
function plot_0sd_pol_complex_parts(result)
    Pu = 2*real(result.P)
    Pv = -2*imag(result.P)

    t_vec = result.t_vec

    fig = plot(; xlabel="\$t\$", ylabel="Bloch component",
        title="Ensemble averaged Bloch components against time",
        size=PLOT_PANEL_SIZE)

    plot!(t_vec, Pu; label="averaged u")
    plot!(t_vec, Pv; label="averaged v")

    return fig
end

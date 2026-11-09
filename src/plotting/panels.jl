"""
    plots_for_jevon_2sd_echo(result)

Return a two-panel figure with transverse-summed temporal intensity
slices and echo efficiency for a two-pulse 2D result.
"""
function plots_for_jevon_2sd_echo(result)
    heatmap = plot_2d_tz_ysum_heatmap(result; operation=abs2, title="Pulse intensity")
    superimposed_zlines_against_t = plot_2d_superimposed_lines_against_t(result; title="Pulse intensity slices in \$z\$", zoom=true, nslices=10, ylabel="\$I/I_{max}\$")
    echo_efficiency = plot_2sd_echo_efficiency(result)

    return _compose_plots(superimposed_zlines_against_t, echo_efficiency;
        layout=(1, 2))
end

"""
    plots_for_jevon_2sd_SIT(result)

Return a six-panel 2D propagation figure with intensity heatmaps,
temporal slices, normalised energy, and transverse slices of time-summed
field and intensity. The fixed layouts use up to 30 transverse slices.
"""
function plots_for_jevon_2sd_SIT(result)
    heatmap = plot_2d_tz_ysum_heatmap(result; operation=abs2, title="Pulse intensity")
    #area = plot_2d_z_ysum_line(result; operation=real, title="Total pulse area")
    pulse_profile_heatmap = plot_2d_zy_pulseprofile_heatmap(result; operation=abs2, title="Pulse intensity profile time slice")
    energy = plot_2d_z_ysum_line(result; operation=abs2,
        title="Total pulse energy", ylabel="\$E/E(0)\$", unit_scale=true)
    superimposed = plot_2d_superimposed_z_ysum_line(result; operation=abs2, nslices=20, title="Pulse intensity slices in \$z\$", zoom=true, ylabel="\$I/I_{max}\$")
    superimposed_area_yslices = plot_2d_z_yslices_line(result; operation=real, nslices=30, title="Pulse area slices in \$y\$ against \$z\$", ylabel="\$A/A_{max}\$")
    superimposed_intensity_yslices = plot_2d_z_yslices_line(result; operation=abs2, nslices=30, title="Pulse intensity slices in \$y\$ against \$z\$", ylabel="\$I/I_{max}\$")

    return _compose_plots(
        heatmap, pulse_profile_heatmap, superimposed_area_yslices,
        superimposed, energy, superimposed_intensity_yslices,
        layout=(2, 3))
end

"""
    plots_for_jevon_1sd_echo(result)

Return a three-panel 1D echo figure with an intensity heatmap, temporal
profiles, and echo efficiency. Requires `Ny == 1` and two pulses.
"""
function plots_for_jevon_1sd_echo(result)
    heatmap = plot_1d_tz_heatmap(result; operation=abs2, title="Pulse intensity")
    superimposed_zlines_against_t = plot_1d_superimposed_lines_against_t(result; title="Pulse intensity slices in \$z\$", zoom=true, nslices=10, ylabel="\$I/I_{max}\$")
    echo_efficiency = plot_1sd_echo_efficiency(result)

    return _compose_plots(heatmap, superimposed_zlines_against_t, echo_efficiency;
        layout=(1, 3))
end

"""
    plots_for_jevon_1sd_SIT(result)

Return a four-panel 1D propagation figure with an intensity heatmap,
temporal profiles, normalised energy, and pulse area. Requires `Ny == 1`
and nonzero input energy and area.
"""
function plots_for_jevon_1sd_SIT(result)
    heatmap = plot_1d_tz_heatmap(result; operation=abs2, title="Pulse intensity")
    superimposed_zlines_against_t = plot_1d_superimposed_lines_against_t(result; title="Pulse intensity slices in \$z\$", zoom=true, nslices=10, ylabel="\$I/I_{max}\$")
    energy = plot_2d_z_ysum_line(result; operation=abs2,
        title="Total pulse energy", ylabel="\$E/E(0)\$", unit_scale=true)
    area = plot_2d_z_ysum_line(result; operation=real,
        title="Pulse area", ylabel="\$A/A(0)\$", unit_scale=true)
    #superimposed_tlines_against_z = plot_1d_superimposed_lines_against_z(result; operation=abs2, nslices=10, title="Pulse intensity slices in \$t\$")

    return _compose_plots(
        heatmap, energy,
        superimposed_zlines_against_t, area,
        layout=(2, 2))
end

"""
    plots_for_jevon_0sd(result)

Return a two-panel local-response figure showing polarisation intensity
and ensemble-averaged Bloch components.
"""
function plots_for_jevon_0sd(result)
    pol_abs2 = plot_0sd_polarisation(result)
    pol_real_n_imag = plot_0sd_pol_complex_parts(result)


    return _compose_plots(pol_abs2, pol_real_n_imag; layout=(2, 1))
end

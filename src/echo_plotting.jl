# Plotting and animation helpers.

using Plots

function heatmap_matrix(A)
    # A is stored as time × z. Plots wants z × time for heatmap(time, z, data).
    # Use transpose, not adjoint ('), so complex quantities are not conjugated.
    return transpose(A)
end

function plot_propagation(result; plot_size=(1400, 1000), z_index=nothing)
    time_vec = result.time_vec
    z_vec = result.z_vec
    Omega = result.Omega
    P = result.P

    if z_index === nothing
        z_index = max(1, length(z_vec))
    end

    heatmap_P_abs2 = heatmap(
        time_vec, z_vec, heatmap_matrix(abs2.(P));
        title="abs2 of polarisation", xlabel="t", ylabel="z", c=:viridis,
    )
    heatmap_P_real = heatmap(
        time_vec, z_vec, heatmap_matrix(real.(P));
        title="real part of polarisation", xlabel="t", ylabel="z", c=:viridis,
    )
    heatmap_P_imag = heatmap(
        time_vec, z_vec, heatmap_matrix(imag.(P));
        title="imaginary part of polarisation", xlabel="t", ylabel="z", c=:viridis,
    )

    heatmap_Omega_abs2 = heatmap(
        time_vec, z_vec, heatmap_matrix(abs2.(Omega));
        title="abs2 of Omega", xlabel="t", ylabel="z", c=:viridis,
    )
    heatmap_Omega_real = heatmap(
        time_vec, z_vec, heatmap_matrix(real.(Omega));
        title="real part of Omega", xlabel="t", ylabel="z", c=:viridis,
    )
    heatmap_Omega_imag = heatmap(
        time_vec, z_vec, heatmap_matrix(imag.(Omega));
        title="imaginary part of Omega", xlabel="t", ylabel="z", c=:viridis,
    )

    plot_P_abs2 = plot(time_vec, abs2.(P[:, z_index]); title="abs2(P) at z=$(round(z_vec[z_index], digits=2))")
    plot_P_real = plot(time_vec, real.(P[:, z_index]); title="real(P) at z=$(round(z_vec[z_index], digits=2))")
    plot_P_imag = plot(time_vec, imag.(P[:, z_index]); title="imag(P) at z=$(round(z_vec[z_index], digits=2))")

    return plot(
        heatmap_P_abs2, heatmap_Omega_abs2, plot_P_abs2,
        heatmap_P_real, heatmap_Omega_real, plot_P_real,
        heatmap_P_imag, heatmap_Omega_imag, plot_P_imag;
        size=plot_size,
        layout=(3, 3),
    )
end

function animate_field_and_polarisation(result; filename="field_and_polarisation.gif", fps=30)
    time_vec = result.time_vec
    z_vec = result.z_vec
    Omega = result.Omega
    P = result.P

    omega_ylim = extrema(abs2.(Omega))
    p_ylim = extrema(abs2.(P))

    anim = @animate for i in eachindex(time_vec)
        p1 = plot(
            z_vec,
            abs2.(Omega[i, :]);
            ylabel="abs2(Ω)",
            ylim=omega_ylim,
            legend=false,
            title="t = $(round(time_vec[i], digits=3))",
        )

        p2 = plot(
            z_vec,
            abs2.(P[i, :]);
            xlabel="z",
            ylabel="abs2(P)",
            ylim=p_ylim,
            legend=false,
        )

        plot(p1, p2; layout=(2, 1), size=(700, 600))
    end

    return gif(anim, filename; fps=fps)
end

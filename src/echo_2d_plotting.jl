using Plots



function animate_2d_echo(result; filename="echo_2d.gif", fps=30)
    field_intensity = abs2.(result.Omega)
    z_vec = result.z_vec
    y_vec = result.y_vec
    t_vec = result.time_vec

    anim = @animate for i in eachindex(t_vec)
        field_plot = heatmap(y_vec, z_vec, field_intensity[i, :, :];
            ylabel="|Ω|²",
            clims=extrema(field_intensity),
            legend=false,
            title="t = $(round(t_vec[i]; sigdigits=4))",
            c=:viridis
        )
    end

    return gif(anim, filename; fps)
end
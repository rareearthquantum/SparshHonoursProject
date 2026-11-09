using JLD2
using Dates

"""
    save_result(filename, result, elapsed_seconds)

Write `result` to the JLD2 file `filename`, including configuration,
elapsed seconds, coordinate vectors, `Omega`, and `P`. The parent
directory must already exist; an existing file at this path is overwritten.
"""
function save_result(filename, result, elapsed_seconds)
    jldsave(
        filename;
        config=Dict(key => getfield(result.cfg, key) for key in fieldnames(typeof(result.cfg))),
        elapsed_seconds=elapsed_seconds,
        detunings=collect(result.detunings),
        t_vec=collect(result.t_vec),
        z_vec=collect(result.z_vec),
        y_vec=collect(result.y_vec),
        Omega=result.Omega,
        P=result.P,
    )
end

"""
    save_data(result, elapsed, subdir_name; parameter_info="placeholder", timestamp="placeholder")

Create `data/subdir_name` under the project root and save the simulation
with [`save_result`](@ref). `elapsed` is in seconds. Placeholder metadata
is replaced by a parameter label and current timestamp. Print the output
directory after saving.
"""
function save_data(result, elapsed, subdir_name; parameter_info="placeholder", timestamp="placeholder")
    (parameter_info == "placeholder") && (parameter_info=make_param_info(result.cfg);)
    (timestamp == "placeholder") && (timestamp=Dates.format(now(), dateformat"yyyymmdd-HHMMSS-sss");)

    data_output_dir = joinpath(dirname(dirname(@__DIR__)), "data", subdir_name)
    mkpath(data_output_dir)
    data_path = joinpath(data_output_dir, parameter_info * "_$(timestamp)_" * ".jld2")
    save_result(data_path, result, elapsed)

    path_elems = split(data_path, "/")
    path_from_projroot = path_elems[end-2] * path_elems[end-1]

    println("Saved jld2 data to directory .../PROJECT_ROOT/" * path_from_projroot)
end;

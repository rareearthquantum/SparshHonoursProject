import Pkg
Pkg.activate(dirname(@__DIR__))

const PROJECT_ROOT = dirname(@__DIR__)

# Keep GR headless when this script is run on a compute node.
get!(ENV, "GKSwstype", "100")

using Dates
using JLD2
using Printf
using TOML

include("../src/input_pulse_methods.jl")
include("../src/echo_propagation.jl")
include("../src/echo_parameter_sweep.jl")

# Parameters not overridden in the profile retain these values.
const BASE_CONFIG = EchoConfig()
const DEFAULT_SWEEP_PROFILE = joinpath(PROJECT_ROOT, "scripts", "parameters", "echo_parameter_sweep_default.toml")

function main()
    profile_path = isempty(ARGS) ? DEFAULT_SWEEP_PROFILE : abspath(first(ARGS))
    profile = load_sweep_profile(profile_path)
    profile_base_config = make_config(BASE_CONFIG, profile.base_config)
    combinations = sweep_combinations(profile.sweep_parameters)
    save_data = lowercase(get(ENV, "ECHO_SWEEP_SAVE_DATA", "false")) in ("1", "true", "yes")

    output_root = if length(ARGS) < 2
        joinpath(PROJECT_ROOT, "data", "echo_sweeps")
    else
        abspath(ARGS[2])
    end
    timestamp = Dates.format(now(), dateformat"yyyymmdd-HHMMSS-sss")
    job_suffix = haskey(ENV, "SLURM_JOB_ID") ? "_job-$(ENV["SLURM_JOB_ID"])" : ""
    output_dir = joinpath(output_root, "$(profile.profile_name)_$(timestamp)$(job_suffix)")
    mkpath(output_dir)

    println("Sweep profile: $(profile.profile_path)")
    println("Running $(length(combinations)) echo propagation simulation(s)")
    println("Output directory: $output_dir")
    println("Saving JLD2 data: $save_data")
    if !save_data
        println("JLD2 saving is off by default; set ECHO_SWEEP_SAVE_DATA=true to enable it.")
    end

    for (index, parameters) in enumerate(combinations)
        cfg = make_config(profile_base_config, parameters)
        label = run_label(index, parameters)
        println("\n[$index/$(length(combinations))] $label")

        result = nothing
        elapsed = @elapsed result = run_propagation(cfg)
        println(@sprintf("Propagation completed in %.3f seconds", elapsed))

        figure = plot_propagation(result)
        plot_path = joinpath(output_dir, label * ".png")
        savefig(figure, plot_path)
        println("Saved $plot_path")

        if save_data
            data_path = joinpath(output_dir, label * ".jld2")
            save_result(data_path, result, parameters, elapsed)
            println("Saved $data_path")
        end

        # Make memory from one simulation available before starting the next.
        figure = nothing
        result = nothing
        GC.gc()
    end

    println("\nSweep complete: $output_dir")
    return output_dir
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

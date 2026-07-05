function config_values(cfg::EchoConfig)
    names = fieldnames(EchoConfig)
    return NamedTuple{names}(Tuple(getfield(cfg, name) for name in names))
end

function named_tuple_from_dict(dict)
    pairs_tuple = Tuple(Symbol(key) => value for (key, value) in pairs(dict))
    return (; pairs_tuple...)
end

function normalize_sweep_parameters(raw::NamedTuple)
    normalized = map(values(raw)) do value
        value isa AbstractVector && return collect(value)
        return [value]
    end
    return NamedTuple{keys(raw)}(Tuple(normalized))
end

function make_config(base::EchoConfig, overrides::NamedTuple)
    unknown = setdiff(keys(overrides), fieldnames(EchoConfig))
    isempty(unknown) || error("Unknown EchoConfig field(s): $(join(unknown, ", "))")
    return EchoConfig(; merge(config_values(base), overrides)...)
end

function load_sweep_profile(filename::AbstractString)
    profile = TOML.parsefile(filename)
    haskey(profile, "sweep_parameters") || error("Sweep profile must define a [sweep_parameters] table: $filename")

    base_config = haskey(profile, "base_config") ? named_tuple_from_dict(profile["base_config"]) : NamedTuple()
    sweep_parameters = normalize_sweep_parameters(named_tuple_from_dict(profile["sweep_parameters"]))
    profile_name = splitext(basename(filename))[1]

    return (
        base_config=base_config,
        sweep_parameters=sweep_parameters,
        profile_name=profile_name,
        profile_path=abspath(filename),
    )
end

function sweep_combinations(parameters::NamedTuple)
    names = keys(parameters)
    isempty(names) && error("SWEEP_PARAMETERS must contain at least one parameter")
    all(!isempty, values(parameters)) || error("Every swept parameter needs at least one value")

    return [NamedTuple{names}(Tuple(values)) for values in Iterators.product(parameters...)]
end

function run_label(index::Integer, parameters::NamedTuple)
    sorted_parameters = sort(collect(pairs(parameters)); by=pair -> string(first(pair)))
    parameter_name = join(("$(key)=$(value)" for (key, value) in sorted_parameters), "_")
    return @sprintf("run-%04d_%s", index, parameter_name)
end

function save_result(filename, result, swept_parameters, elapsed_seconds)
    # rotate_grid is deliberately omitted: it is large and can be reconstructed
    # from detunings and time_vec.
    jldsave(
        filename;
        config=config_values(result.cfg),
        swept_parameters=swept_parameters,
        elapsed_seconds=elapsed_seconds,
        detunings=collect(result.detunings),
        time_vec=collect(result.time_vec),
        z_vec=collect(result.z_vec),
        Omega=result.Omega,
        P=result.P,
    )
end

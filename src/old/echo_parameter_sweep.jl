using Printf

function config_values(cfg::EchoConfig)
    names = fieldnames(EchoConfig)
    return NamedTuple{names}(Tuple(getfield(cfg, name) for name in names))
end

function make_config(base::EchoConfig, overrides::NamedTuple)
    haskey(overrides, :Zf) &&
        (overrides = merge(overrides, (Zf=meters(overrides.Zf),)))
    haskey(overrides, :alpha) &&
        (overrides = merge(overrides, (alpha=coupling_si(overrides.alpha),)))
    return EchoConfig(; merge(config_values(base), overrides)...)
end

function sweep_combinations(parameters::NamedTuple)
    return [NamedTuple{keys(parameters)}(values) for values in Iterators.product(parameters...)]
end

function run_label(index::Integer, parameters::NamedTuple)
    sorted_parameters = sort(collect(pairs(parameters)); by=pair -> string(first(pair)))
    parameter_name = join(("$(key)=$(value)" for (key, value) in sorted_parameters), "_")
    return @sprintf("run-%04d_%s", index, parameter_name)
end

function save_result(filename, result, swept_parameters, elapsed_seconds)
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

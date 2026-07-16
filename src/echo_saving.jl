using JLD2

function save_result(filename, result, elapsed_seconds)
    jldsave(
        filename;
        config=Dict(key => getfield(result.cfg, key) for key in fieldnames(typeof(result.cfg))),
        elapsed_seconds=elapsed_seconds,
        detunings=collect(result.detunings),
        time_vec=collect(result.time_vec),
        z_vec=collect(result.z_vec),
        Omega=result.Omega,
        P=result.P,
    )
end

function save_sweep_result(filename, result, swept_parameters, elapsed_seconds) #unfinished
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

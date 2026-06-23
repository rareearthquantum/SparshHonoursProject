using TOML

include(srcdir("constants.jl"))

function required_value(table::AbstractDict, key::String, location::String)
    haskey(table, key) || throw(ArgumentError("missing '$location.$key' in simulation config"))
    return table[key]
end

function numeric_value(table::AbstractDict, key::String, location::String)
    value = required_value(table, key, location)
    value isa Real || throw(ArgumentError("'$location.$key' must be numeric"))
    return Float64(value)
end

function integer_value(table::AbstractDict, key::String, location::String)
    value = required_value(table, key, location)
    value isa Integer || throw(ArgumentError("'$location.$key' must be an integer"))
    return Int(value)
end

function range_value(table::AbstractDict, key::String, location::String)
    value = required_value(table, key, location)
    (value isa AbstractVector && length(value) == 2 && all(x -> x isa Real, value)) ||
        throw(ArgumentError("'$location.$key' must be a two-element numeric array"))
    range = (Float64(value[1]), Float64(value[2]))
    range[1] < range[2] || throw(ArgumentError("'$location.$key' must be increasing"))
    return range
end

function validate_grid(grid::AbstractDict)
    N = (
        integer_value(grid, "Nz", "grid"),
        integer_value(grid, "Nd", "grid"),
        integer_value(grid, "Nt", "grid"),
        integer_value(grid, "Ny", "grid")
    )
    all(n -> n >= 2, N) || throw(ArgumentError("all grid sizes must be at least 2"))
    return N
end

function dimensionless_pulses(pulses)
    return [
        (
            (
                center=numeric_value(pulse, "time_center", "pulses[$i]"),
                width=numeric_value(pulse, "time_width", "pulses[$i]"),
                area=numeric_value(pulse, "time_area", "pulses[$i]")
            ),
            (
                center=numeric_value(pulse, "y_center", "pulses[$i]"),
                width=numeric_value(pulse, "y_width", "pulses[$i]"),
                area=numeric_value(pulse, "y_area", "pulses[$i]")
            )
        ) for (i, pulse) in enumerate(pulses)
    ]
end

function physical_pulses(pulses, scales)
    return [
        (
            (
                center=numeric_value(pulse, "time_center_s", "pulses[$i]") / scales.t,
                width=numeric_value(pulse, "time_width_s", "pulses[$i]") / scales.t,
                # Pulse area is integral(Omega dt), so it is dimensionless.
                area=numeric_value(pulse, "time_area_rad", "pulses[$i]")
            ),
            (
                center=numeric_value(pulse, "y_center_m", "pulses[$i]") / scales.y,
                width=numeric_value(pulse, "y_width_m", "pulses[$i]") / scales.y,
                # The spatial Gaussian is unitless when its area is divided by Y0.
                area=numeric_value(pulse, "y_area_m", "pulses[$i]") / scales.y
            )
        ) for (i, pulse) in enumerate(pulses)
    ]
end

function validate_pulses(pulse_params)
    isempty(pulse_params) && throw(ArgumentError("the config must contain at least one [[pulses]] table"))
    for (i, pulse) in enumerate(pulse_params)
        pulse[1].width > 0 || throw(ArgumentError("pulses[$i].time_width must be positive"))
        pulse[2].width > 0 || throw(ArgumentError("pulses[$i].y_width must be positive"))
    end
    return pulse_params
end

function load_simulation_config(path::AbstractString)
    data = TOML.parsefile(path)
    mode = Symbol(required_value(data, "mode", "root"))
    name = String(get(data, "name", splitext(basename(path))[1]))
    grid = required_value(data, "grid", "root")
    domain = required_value(data, "domain", "root")
    model = required_value(data, "model", "root")
    pulses = required_value(data, "pulses", "root")
    N = validate_grid(grid)

    if mode == :dimensionless
        scales = (; z=1.0, t=1.0, y=1.0)
        ranges = (
            range_value(domain, "z", "domain"),
            range_value(domain, "t", "domain"),
            range_value(domain, "y", "domain")
        )
        alpha = numeric_value(model, "alpha", "model")
        beta = numeric_value(model, "beta", "model")
        pulse_params = validate_pulses(dimensionless_pulses(pulses))
        physical = nothing
        axis_labels = (; z="z", t="t", y="y")
    elseif mode == :physical_si
        scale_table = required_value(data, "scales", "root")
        scales = (
            z=numeric_value(scale_table, "z_m", "scales"),
            t=numeric_value(scale_table, "time_s", "scales"),
            y=numeric_value(scale_table, "y_m", "scales")
        )
        all(x -> x > 0, scales) || throw(ArgumentError("all physical scales must be positive"))

        ranges_si = (
            range_value(domain, "z_m", "domain"),
            range_value(domain, "t_s", "domain"),
            range_value(domain, "y_m", "domain")
        )
        ranges = (
            (ranges_si[1][1] / scales.z, ranges_si[1][2] / scales.z),
            (ranges_si[2][1] / scales.t, ranges_si[2][2] / scales.t),
            (ranges_si[3][1] / scales.y, ranges_si[3][2] / scales.y)
        )

        carrier_frequency_hz = numeric_value(model, "carrier_frequency_hz", "model")
        coupling_hz_per_m = numeric_value(model, "coupling_hz_per_m", "model")
        coefficients = dimensionless_coefficients(;
            carrier_frequency_hz,
            coupling_hz_per_m,
            z_scale_m=scales.z,
            time_scale_s=scales.t,
            y_scale_m=scales.y
        )
        alpha, beta = coefficients.alpha, coefficients.beta
        pulse_params = validate_pulses(physical_pulses(pulses, scales))
        physical = (;
            carrier_frequency_hz,
            coupling_hz_per_m,
            beta_m=coefficients.beta_m,
            ranges_si
        )
        axis_labels = (; z="z (m)", t="t (s)", y="y (m)")
    else
        throw(ArgumentError("unsupported config mode '$mode'; use 'dimensionless' or 'physical_si'"))
    end

    return (; name, mode, N, ranges, pulse_params, alpha, beta, scales, physical, axis_labels)
end

function output_grids(config, grids)
    z_grid, t_grid, y_grid = grids
    config.mode == :physical_si || return grids
    return (
        z_grid .* config.scales.z,
        t_grid .* config.scales.t,
        y_grid .* config.scales.y
    )
end

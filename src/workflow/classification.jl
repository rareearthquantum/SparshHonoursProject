const PULSE_AREA_RTOL = 1e-8
const PULSE_AREA_ATOL = 1e-10

"""
    _area_is(area, target)

Return whether a pulse area matches `target` using the shared
`PULSE_AREA_RTOL` and `PULSE_AREA_ATOL` tolerances.
"""
_area_is(area, target) = isapprox(area, target; rtol=PULSE_AREA_RTOL, atol=PULSE_AREA_ATOL)

"""
    classify_simulation(cfg) -> Symbol

Infer the physical test represented by `cfg`. The classification uses the pulse
sequence, coupling (`alpha`), diffraction (`beta`), and simulated spatial axes.
It throws an `ArgumentError` when the parameters do not describe one of the
supported tests.
"""
function classify_simulation(cfg)
    npulses = length(cfg.pulses)
    npulses in (1, 2) ||
        throw(ArgumentError("expected one or two pulses, got $npulses"))

    cfg.alpha >= 0 || throw(ArgumentError("alpha must be non-negative"))
    cfg.beta >= 0 || throw(ArgumentError("beta must be non-negative"))

    areas = [pulse[1].area for pulse in cfg.pulses]
    all(>(0), areas) || throw(ArgumentError("pulse areas must be positive; got $(areas)"))

    if npulses == 2
        is_echo = _area_is(areas[1], pi / 2) && _area_is(areas[2], pi)
        is_echo || throw(ArgumentError(
            "a two-pulse echo test requires pulse areas pi/2 and pi; got $(areas)",
        ))

        if cfg.Nz > 1 && iszero(cfg.alpha)
            throw(ArgumentError(
                "a propagating echo test requires alpha > 0 so the echo field can be emitted",
            ))
        end

        return :echo
    end

    # With no propagation axis, only the local atomic response is being tested.
    cfg.Nz == 1 && return :zero_dimensional

    has_medium_interaction = !iszero(cfg.alpha)
    has_diffraction = cfg.Ny > 1 && !iszero(cfg.beta)

    if !has_medium_interaction
        has_diffraction || throw(ArgumentError(
            "the configuration has neither medium interaction nor transverse diffraction",
        ))
        return :diffraction
    end

    medium_test = areas[1] > pi ? :sit : :attenuation

    if has_diffraction
        return medium_test == :sit ? :sit_with_diffraction : :attenuation_with_diffraction
    end

    return medium_test
end

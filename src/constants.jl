const tera, giga, mega, kilo, centi, milli, micro, nano, pico =
    1e12, 1e9, 1e6, 1e3, 1e-2, 1e-3, 1e-6, 1e-9, 1e-12

const SPEED_OF_LIGHT_M_PER_S = 299_792_458.0
const VACUUM_PERMEABILITY_H_PER_M = 1.256_637_062_12e-6

"Paraxial diffraction coefficient 1/(2k), in metres."
function diffraction_coefficient_m(carrier_frequency_hz::Real)
    omega = 2pi * carrier_frequency_hz
    k = omega / SPEED_OF_LIGHT_M_PER_S
    return 1 / (2k)
end

"Maxwell prefactor multiplying a physical polarization, not a coherence sum."
function maxwell_polarization_prefactor(carrier_frequency_hz::Real)
    omega = 2pi * carrier_frequency_hz
    return VACUUM_PERMEABILITY_H_PER_M * omega * SPEED_OF_LIGHT_M_PER_S / 2
end

"""
Convert physical Maxwell--Bloch coefficients to the dimensionless coefficients
used by the numerical kernel.

`coupling_hz_per_m` is the effective coefficient multiplying the normalized
dimensionless coherence/polarization. It must already include atom density and
dipole moment. It is deliberately distinct
from `maxwell_polarization_prefactor`, which multiplies an SI polarization.
"""
function dimensionless_coefficients(;
    carrier_frequency_hz::Real,
    coupling_hz_per_m::Real,
    z_scale_m::Real,
    time_scale_s::Real,
    y_scale_m::Real
)
    alpha = coupling_hz_per_m * z_scale_m * time_scale_s
    beta_m = diffraction_coefficient_m(carrier_frequency_hz)
    beta = beta_m * z_scale_m / y_scale_m^2
    return (; alpha, beta, beta_m)
end

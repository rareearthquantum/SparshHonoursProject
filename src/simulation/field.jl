"""
    field_2d!(rabi_frequency_derivative, rabi_frequency_spectrum, z, p)

Fill the spectral field derivative with `im * alpha * polarisation_spectrum
* rotfactor`, using `p = (alpha, polarisation_spectrum, rotfactor)`.
The field spectrum and `z` are accepted for the stepping interface.
Mutates `rabi_frequency_derivative` and returns `nothing`.
"""
function field_2d!(rabi_frequency_derivative, rabi_frequency_spectrum, z, p)
    alpha, polarisation_spectrum, rotfactor = p

    @. rabi_frequency_derivative = im * alpha * polarisation_spectrum * rotfactor

    return nothing
end

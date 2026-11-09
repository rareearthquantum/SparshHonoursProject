struct PulseParams
    center::Float64
    width::Float64
    area::Float64
end

"""
Construct Gaussian or top-hat pulse parameters, converting keyword values
to `Float64`. For Gaussians, `width` is the standard deviation; for top hats,
it is the full support width. Validation is performed by `EchoConfig`.
"""
PulseParams(; center, width, area) =
    PulseParams(Float64(center), Float64(width), Float64(area))

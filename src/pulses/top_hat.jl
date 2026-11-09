"""
Return an area-normalised rectangular pulse as `Float64`. The value is
`area / width` when `abs(u - center) <= width / 2`, including the endpoints,
and zero otherwise. `width` must be positive.
"""
top_hat_pulse(u, center, width, area)::Float64 = ( abs(u-center) <= width/2 ) ? area/width : 0

"""
Evaluate an area-normalised top hat at `u` using the centre, positive full
width, and area stored in `pulseparams`. Return `Float64`; the support
includes both endpoints.
"""
top_hat_pulse(u,  pulseparams::PulseParams)::Float64 = top_hat_pulse(u, pulseparams.center, pulseparams.width, pulseparams.area)

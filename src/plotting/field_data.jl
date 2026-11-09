"""
    _normalise_by_max!(values)

Divide every entry of `values` by its maximum in place and return the
same array. Requires a nonzero maximum for finite normalised values.
"""
function _normalise_by_max!(values)
    values ./= maximum(values)
    return values
end

"""
    _field_2d_summed_over_y(result, operation)

Apply `operation` elementwise to `result.Omega[time, z, y]` and sum
over the transverse axis. Return a time-by-z matrix; no grid-spacing
factor or normalisation is applied.
"""
_field_2d_summed_over_y(result, operation) =
    dropdims(sum(operation.(result.Omega); dims=3); dims=3)

"""
    _field_1d(result, operation)

Remove the singleton transverse axis of `result.Omega` and apply
`operation` elementwise, returning a time-by-z matrix. Requires `Ny == 1`.
"""
_field_1d(result, operation) = operation.(dropdims(result.Omega; dims=3))

"""
    zoomingin(vector; threshold=0.001)

Return inclusive row bounds around entries at least `threshold` in any
column of the matrix `vector`, padded by one row within array bounds.
Requires at least one above-threshold entry for meaningful ordered bounds.
"""
function zoomingin(vector; threshold=0.001)

    first_above = lastindex(vector[:, 1])
    last_above = firstindex(vector[:, 1])

    for i in 1:length(vector[1, :])
        viewvec = @view vector[:, i]
        temp = findfirst(>=(threshold), viewvec)
        if temp !== nothing && temp < first_above
            first_above = temp
        end

        temp = findlast(>=(threshold), viewvec)
        if temp !== nothing && temp > last_above
            last_above = temp
        end
    end

    startindex = max(first_above - 1, firstindex(vector[:, 1]))
    endindex = min(last_above + 1, lastindex(vector[:, 1]))

    return (startindex, endindex)
end

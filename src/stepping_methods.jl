function euler_step!(v,f,h,i)
	v[i+1] = v[i] + h*f(i)
end

function adams_bashforth_2step!(
	v::AbstractArray,
	f::Function,
	h::Float64,
	i::Int64
)
	v[i+1] = v[i] + (h/2)*(3*f(i) - f(i-1))
end

function adams_bashforth_3step!(
	v::AbstractArray,
	f::Function,
	h::Float64,
	i::Int64
)
	v[i+1] = v[i] + (h/12)*(23f(i) - 16f(i-1) + 5f(i-2))
end

function adams_bashforth_4step!(
	v::AbstractArray,
	f::Function,
	h::Float64,
	i::Int64
)
	v[i+1] = v[i] + (h/24)*(55f(i) - 59f(i-1) + 37f(i-2) - 9f(i-3))
end

function adams_moulton_2step!(
	v::AbstractArray,
	f::Function,
	fp,
	h::Float64,
	i::Int64
)
	v[i+1] = v[i] + (h/12)*(5*fp + 8*f(i) - f(i-1))
end

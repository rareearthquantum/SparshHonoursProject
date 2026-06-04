function pulse(
	t::Float64,
	t0::Float64,
	width::Float64,
	A::Float64
)
	A * exp(-((t-t0)/width)^2)
end

function Ein_t(t::Float64)
	p1c, p1w, p1A = 5.0, 1.0, 1.0
	p1 = pulse(t,p1c,p1w,p1A)
	#p2 = pulse(t,18.0,1.0,1.0)
	return p1
end

function Ein_y(y::Float64)
	p_c, p_w, p_A = 0.0, 1.0, 1.0
	p = pulse(y,p_c,p_w,p_A)
	return p
end


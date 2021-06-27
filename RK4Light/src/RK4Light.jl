module RK4Light

	export RK4Problem, rk4
	
	struct RK4Problem
		f
		u0
		tspan
		p
	end
	
	function rk4(prob::RK4Problem; savenum = 10^5)
		u = zeros(savenum)
		t = zeros(savenum)
		u[1] = prob.u0
		t[1] = prob.tspan[1]
		
		tstep = (prob.tspan[2] - prob.tspan[1])/savenum
		
		for i in 2:savenum
			t[i] = t[i-1] + tstep
			
			k1 = prob.f(u[i-1], t[i-1], prob.p)
			k2 = prob.f(u[i-1] + tstep * k1/2, t[i-1] + tstep/2, prob.p)
			k3 = prob.f(u[i-1] + tstep * k2/2, t[i-1] + tstep/2, prob.p)
			k4 = prob.f(u[i-1] + tstep * k3, t[i], prob.p)
			
			u[i] = u[i-1] + 1/6 * tstep * (k1 + 2k2 + 2k3 + k4)
			
			@show k1, k2, k3, k4
		end
		
		(t, u)
	end
	
end # module

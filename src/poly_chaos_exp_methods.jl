# Uncertainty Propagation Solvers Methods

"""
Polynomial Chaos Expansion Methods
"""

####################################
#######  Basis and Sampling ########
####################################

function generate_polynomials!(prop::PCEPropagator3)
    # unpack options
    opts = prop.opts
    uncert_type = opts.uncert_type  # length(uncert_type == n)
    @assert length(uncert_type) == prop.n
    d = opts.d
    # build corresponding multi-ortho-poly basis (can be mixed)
    mop_list = Vector{AbstractOrthoPoly}([])                 # store uncertainty types
    for i=1:length(uncert_type)
        type = uncert_type[i]
        if type == :gaussian      # polynomials are Hermite (probablist's Hermite)
            op = GaussOrthoPoly(d)
        elseif type == :uniform   # polynomials are Legendre (check Uniform01OrthoPoly)
            op = LegendreOrthoPoly(d)
        elseif type == :gamma     # polynomials are Laguerre
            op = GammaOrthoPoly(d)
        elseif type == :beta      # polynomials are Jacobi
            op = Beta01OrthoPoly(d)
        else
            prinln("Unknown or unsupported type of Uncertainty")
        end
        push!(mop_list, op)
    end
    mop = MultiOrthoPoly(mop_list, d)       # build multivariate-ortho-poly-basis
    prop.mop = mop                     # add multivariate-ortho basis in struct
    prop.P = mop.dim                   # Get dimension of the basis
    return nothing
end

function draw_samples!(prop::PCEPropagator3)
    # Draw corresponding samples for propagation
    opts = prop.opts
    uncert_type = opts.uncert_type     # length(uncert_type == n)
    n = prop.n                         # get system dimension
    M = opts.M                         # number of samples to draw
    x0 = prop.x0                       # get initial point
    ξ = zeros(n ,M)                    # initialize sample matrix
    for i=1:n
        type = uncert_type[i]
        if type == :gaussian      # polynomials are Hermite (probablist's Hermite)
            D = Normal(x0[i], sqrt(Q0[i, i]))  # independence of the coordinates assumed here
        elseif type == :uniform   # polynomials are Legendre (check Uniform01OrthoPoly)
            D = Uniform(x0[i]-0.05, x0[i]+0.05)  # modify bounds I guess
        elseif type == :gamma     # polynomials are Laguerre
            D = Gamma(α, θ)         # modify parameters
        elseif type == :beta      # polynomials are Jacobi
            D = Beta(α, β)          # modify parameters
        else
            prinln("Unknown or unsupported type of Uncertainty")
        end
        ξ_line = rand(D, M)
        #ξ_col = randn(M)
        ξ[i, :] = ξ_line
    end
    prop.ξ = ξ                      # add to the propagator struct for future use
    return nothing
end


function propagate_samples!(prop::PCEPropagator3)
    ξ = prop.ξ                    # samples drawn
    opts = prop.opts    # unpack options
    t0 = opts.t0        # get start date
    tf = opts.tf        # get end date
    dt = opts.dt        # get time step propagation
    n = prop.n          # get system dimension (will be replaced with dist dimensions n+p later)
    M = opts.M          # get number of samples to draw
    model = prop.model  # dynamic function to use for propagation (rk4 for now)
    f = prop.f          # dynamic function to use for propagation (rk4 for now)
    t_list = t0:dt:tf    # time span
	prop.nbr_steps = length(t_list)
    #samples = zeros(n, length(t_list), M)  # initialize samples, pre-allocated
    samples = []
    for i=1:1:M
        t_sim, Z = rk4(f, ξ[:, i], 0.0, dt, [t0; tf])  # propagate using RK4 (see using TrajOpt later)
        push!(samples, Z)                              # each matrix contains the history of one sample
    end
    prop.ξ_propagated  = samples
    return nothing
end


function rk4(f, y_0, p, dt, t_span)
    T = t_span[1]:dt:t_span[end]
    y = zeros(length(T), length(y_0))
    if length(y_0) == 1
        y[1, :] = [y_0]
    else
        y[1, :] = y_0
    end
    for i=1:1:length(T)-1
        t = T[i]
        y_star = y[i, :]
        k1 = f(t, y_star, p)
        y1 = y_star+k1*dt/2 #intermediate evaluation value
        k2 = f(t+dt/2, y1, p)
        y2 = y_star+k2*dt/2
        k3 = f(t+dt/2, y2, p)
        y3 = y_star+k3*dt
        k4 = f(t+dt, y3, p)
        m = (k1+2*k2+2*k3+k4)/6 #slope average
        y[i+1, :] = y_star + m*dt
    end
    return T, y'
end


# need to work on rollout using TrajOptimization here
function rollout!(model::AbstractModel, ) where {T}
	n = solver.n
	m = solver.m
	N = solver.N

	X_traj = solver.X_traj
	U_traj = solver.U_traj
	Δt_traj = solver.Δt_traj
	K_traj = solver.K_traj
	d_traj = solver.d_traj

	X = [SVector{n,T}(zeros(n)) for k=1:N]
	U = [SVector{m,T}(zeros(m)) for k=1:N-1]
	X[1] = x0

	δx = SVector{n,T}(zeros(n))
	δu = SVector{m,T}(zeros(m))

    for k = 1:solver.N-1
        δx = TO.state_diff(solver.model, X[k], X_traj[k])
		δu = 0.0*d_traj[k] ################################################ PROBLEM HERE
		mul!(δu, K_traj[k], δx, 1.0, 1.0)
        U[k] = U_traj[k] + δu
		zk = TO.KnotPoint(X[k], U[k], Δt_traj[k])
		X[k+1] = TO.discrete_dynamics(TO.RK4, solver.model, zk)
    end
	return X
end


####################################
######Coefficients Computation######
####################################

function compute_Phi_matrix!(prop::PCEPropagator3)
    # Uses Samples to compute Φ Matrix once
    opts = prop.opts
    M = opts.M
    mop = prop.mop
    Q0 = prop.Q
    x0 = prop.x0
    P = prop.P
    n = prop.n
    ξ = prop.ξ      # get initial sampled points
    Phi = zeros(M, P)
    for i=1:M
        for j=1:1:P
            vec = [(ξ[k, i]-x0[k])/sqrt(Q0[k, k]) for k=1:1:n]
            res = PolyChaos.evaluate(mop.ind[j, :], vec, mop)[1]
            Phi[i, j] = res
        end
    end
    prop.Φ = Phi   # update structure
    prop.pinv_Φ = pinv(Phi)
    return nothing
end

function compute_coeff_PCE_step(prop::PCEPropagator3, t::T) where {T}
    #function computes PCE coefficients for time step t
    #t is the time step we want to look at t belongs to
    opts = prop.opts
    ξ_propagated = prop.ξ_propagated
    P = prop.P
    n = prop.n
    M = opts.M
    pinv_Φ = prop.pinv_Φ
    C = zeros(P,n) #contains coefficients for one time step
    for i=1:1:n
        vec = [ξ_propagated[j][i, t] for j=1:M]
        C[:, i] = (pinv_Φ*vec)'
    end
    return C
end

function compute_coeff_PCE_full!(prop::PCEPropagator3)
    opts = prop.opts
    n = prop.n
    P = prop.P
    t0 = opts.t0
    tf = opts.tf
    dt = opts.dt
    ξ_propagated = prop.ξ_propagated
    nbr_steps = prop.nbr_steps
    #C = zeros(P, n, nbr_steps)
    C = []
    for j = 1:1:nbr_steps
        c = compute_coeff_PCE_step(prop, j)
        push!(C, c)                    #each matrix in C is associated with a time step
        @show(j)
    end
    prop.C = C
    return nothing
end

function compute_moments_PCE!(prop::PCEPropagator3)
    # Compute the first two moments from PCE coefficients
    C = prop.C                                      # get coefficients
    n = prop.n                                      # get system dimension
    P = prop.P                                      # get mop dimensions
    nbr_steps = prop.nbr_steps                      # get nbr time steps
    m = zeros(n, nbr_steps)                         # initialize mean
    var = zeros(n, nbr_steps)                       # initialize variance
    for j=1:1:nbr_steps
        m[:, j] = C[j][1, :]                            # extract mean
        for i = 1:1:n
            var[i, j] = sum(C[j][k, i]^2 for k=2:1:P)  # compute variance
        end
    end
    prop.mean = m                                      # update propagator
    prop.var = var                                  # update propagator
    return nothing
end

function PCE_solve!(prop::PCEPropagator3)
    # Solve Routine for Polynomial Chaos Expansion methods
    generate_polynomials!(prop)         # get multivariate orthogonal basis, stored in prop.mop
    draw_samples!(prop)                 # draw samples, stored in prop.ξ
    propagate_samples!(prop)            # propagate samples through model, stored in prop.ξ_propagated
    compute_Phi_matrix!(prop)           # compute Φ Matrix, stored in prop.Φ, used for coefficients computation
    compute_coeff_PCE_full!(prop)       # compute PCE coefficients at each time step
    compute_moments_PCE!(prop)          # get the probabilistic moments from PCE coefficients
end

############################
##### Visualization ########
############################

#reconstruct density distributions (detail for x position)
function kernel_density_estimation_plot(prop::PCEPropagator3; Nr::Int=1000, dim::Int=1, t::Int=1) where {T}
	# function plotting the density function using PCE results
	# Nr = number of samples for reconstruction
	# t = time step at which we compute the density function for dimension dim
	P = prop.P                      # get mop dimension
	x0 = prop.x0                    # get initial state
	Q0 = prop.Q                     # get initial uncertainty on the state
	n = prop.n                      # get system dimension
	mop = prop.mop                  # get multi-ortho-poly basis
	C = prop.C                      # get the computed coefficients during PCE prop
	Ξ = zeros(length(x0), Nr)       # generate samples (necesaary for reconstruction of density)
	D = MvNormal(x0, sqrt(Q0))
	rand!(D, Ξ)
	ũ_x = zeros(Nr, 1)
	for p=1:P
		ũ_x = ũ_x + [C[t][p, dim]*PolyChaos.evaluate(mop.ind[p, :], Ξ[:, i], mop)[1] for i=1:Nr]
	end
	pdf = kde(ũ_x[:])
	plt = plot()
	plot!(pdf.x, pdf.density)
	display(plt)
	return nothing
end


function sig_pce(var_pce)
    n, t = size(var_pce)
    V = zeros(n, n, t)
    for j = 1:1:t
        V[:, :, j] = Diagonal(var_pce[:, j])
    end
    VV = sig(V)
    return VV
end




"""
Linear Covariance Methods
"""

function lincov_prop(x0, Q0, Z, t0, tf, dt, model_lincov)
    #Z is nominal trajectory
    #dt should be same on nominal and desired step here
    #model_lincov is the model modified to integrate forwardDiff
    #E does not make much sense here so far.
    T = t0:dt:tf
    n = length(x0)
    E = zeros(n, length(T))
    QQ = zeros(n, n, length(T))
    x1 = x0
    E[:, 1] = x0
    QQ[:, :, 1] = Q0
    Q1 = Q0
    for i=1:1:length(T)-1
        t = T[i]
        p = [45*pi/180]
        f(u) = model_lincov(t, u, p)
        A = ForwardDiff.jacobian(f, Z[:, i])
        Φ = exp(A*dt)
        x = Φ*x1
        Q = Φ*Q1*Φ'
        E[:, i+1] = x
        QQ[:, :, i+1] = Q
        x1 = x
        Q1 = Q
    end
    return T, QQ
end



"""
Set-based Unscented Propagation Methods
"""

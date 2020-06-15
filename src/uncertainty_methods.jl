# Uncertainty Propagation Solvers Methods


"""
Polynomial Chaos Expansion Methods
"""

####################################
#######  Basis and Sampling ########
####################################

function generate_polynomials(prop::PCEPropagator)
    # unpack options
    opts = prop.opts
    uncert_type = opts.uncert_type  # length(uncert_type == n)
    @assert length(uncert_type) == prop.n
    d = opts.d
    # build corresponding multi-ortho-poly basis (can be mixed)
    mop_list = []                 # store uncertainty types
    for i=1:length(uncert_type)
        type = uncert_type[i]
        if type == :gaussian      # polynomials are Hermite (probablist's Hermite)
            op = GaussianOrthoPoly(d)
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
    mop = MultiOrthoPoly(mop, d)       # build multivariate-ortho-poly-basis
    return mop
end

function draw_samples(prop::PCEPropagator)
    # Draw corresponding samples for propagation
    opts = prop.opts
    uncert_type = opts.uncert_type     # length(uncert_type == n)
    n = prop.n                         # get system dimension
    M = opts.M                         # number of samples to draw
    ξ = zeros(n ,M)                    # initialize sample matrix
    for i=1:n
        type = uncert_type[i]
        if type == :gaussian      # polynomials are Hermite (probablist's Hermite)
            D = Normal(0.0, 1.0)
        elseif type == :uniform   # polynomials are Legendre (check Uniform01OrthoPoly)
            D = Uniform(-1.0, 1.0)  # modify bounds I guess
        elseif type == :gamma     # polynomials are Laguerre
            D = Gamma(α, θ)         # modify parameters
        elseif type == :beta      # polynomials are Jacobi
            D = Beta(α, β)          # modify parameters
        else
            prinln("Unknown or unsupported type of Uncertainty")
        end
        ξ_col = rand(D, M)
        ξ[i, :] = ξ_col
    end
    return ξ                    # add to the propagator struct for future use
end


function propagate_samples(prop::PCEPropagator)
    ξ = prop.ξ                    # samples drawn
    opts = prop.opts    # unpack options
    t0 = opts.t0        # get start date
    tf = opts.tf        # get end date
    dt = opts.dt        # get time step propagation
    n = prop.n          # get system dimension (will be replaced with dist dimensions n+p later)
    M = opts.M          # get number of samples to draw
    t_list = t0:dt:tf    # time span
    samples = zeros(n, length(t_list), M)  # initialize samples, pre-allocated
    for i=1:1:M
        t_sim, Z = rk4(model, ξ[:, i], dt, [t0; tf])  # propagate using RK4 (see using TrajOpt later)
        samples[:, :, i] = Z
    end
    return samples      # should be added to PCEPropagator struct
end


####################################
######Coefficients Computation######
####################################

function compute_Phi_matrix(Ms, P, mop, samples, x0, Q0, n, ξ)
    # Uses Samples to compute Φ Matrix once
    opts = prop.opts
    M = opts.M
    mop = prop.mop
    samples = prop.ξ_prop
    Q0 = prop.Q
    P = prop.P
    Phi = zeros(M, P)
    for i=1:M
        for j=1:1:P
            vec = [(ξ[k, i]-x0[k])/sqrt(Q0[k, k]) for k=1:1:n]
            res = PolyChaos.evaluate(mop.ind[j, :], vec, mop)[1]
            Phi[i, j] = res
        end
    end
    prop.Φ = Phi   # update structure
    return nothing
end

function compute_coeff_PCE_step(samples, t, A, n, P) #t is the time step we want to look at t belongs to
    #function computes PCE coefficients for time step t
    C = zeros(P,n) #contains coefficients for one time step
    for i=1:1:n
        C[:, i] = (A*samples[i, t, :])'
    end
    return C
end

function compute_coeff_PCE_full(samples, Phi, t0, tf, dt, P)
    n, t, Ms = size(samples)
    C = zeros(P, n, t)
    A = pinv(Phi)
    T = t0:dt:tf
    for j = 1:1:length(T)
        c = compute_coeff_PCE_step(samples, j, A, n, P)
        C[:, :, j] = c
        @show(j)
    end
    return C
end

function mean_var_PCE(CC, T, n, P)
    t = length(T)
    m = zeros(n, t)
    var = zeros(n, t)
    for j=1:1:t
        m[:, j] = CC[1, :, j]
        for i = 1:1:n
            var[i, j] = sum(CC[k, i, j]^2 for k=2:1:P)
        end
    end
    return m, var
end

function simulation_PCE(t0, tf, dt, Ms, type, D, d, x0, Q0, model, p)
    #returns the coefficients of the PC expansion for all time steps specified
    #Ms is number of samples (non-intrusive PCE method here)
    #D is the distribution type we want (linked with "type")
    n = length(x0)
    nrec = 100
    mop = generate_poly(nrec, d, n, type)
    DD = MvNormal(x0, Q0/9)
    ξ = rand!(DD, zeros(n, Ms))
    samples = generate_samples_PCE(t0, tf, dt, ξ, Ms, model, n, p) #counting
    P = mop.dim #number of poly in the family
    Phi = compute_Phi_matrix(Ms, P, mop, samples, x0, Q0, n, ξ)
    CC = compute_coeff_PCE_full(samples, Phi, t0, tf, dt, P)
    T = t0:dt:tf
    m, var = mean_var_PCE(CC, T, n, P)
    return T, m, var
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

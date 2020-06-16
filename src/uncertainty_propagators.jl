# Uncertainty Propagation Solvers Propagation

abstract type AbstractUncertaintyPropagator end  # abstract type for dispatch

"""
Polynomial Chaos Expansion Structures
"""

# Funnel structures for multipliers version of funnel computation

@with_kw mutable struct PCEOptions{T}

	"Type of Solve."
	solve_type::Symbol=:samples  #:galerkin(analytical method)

	"Numnber of Samples."
	M::Int=1000                       #Number of sampled trajectories

	"Degree of Polynomial Expansion."
	d::Int=3                          #Max degree of terms on expansion

	"Type of Uncertainty Considered, Distributions for each param."
	uncert_type::Vector{Symbol}            #vector of types, for each parameters

	"Experiment Starting Date."
	t0::T=0.0                              #secondes

	"Experiment Ending Date."
	tf::T=25.0                             #secondes

	"Propagation Dynamics Time Step"       # Will be changed in the future
	dt::T=1e-2

end


mutable struct PCEPropagator3{T} <: AbstractUncertaintyPropagator
    # Polynomial Chaos Expansion Solver
	n::Int                                # State size
	m::Int                                # Control size
	p::Int                                # Disturbance size
	#N::Int                               # Number of time steps
	model::TO.AbstractModel               # Dynamics model of the system (TrajOpt struct)
	f::Any                                # Regular Dynamics Function (not adapted to TO)
	#model_dist::TO.AbstractModel         # Dynamics model of the system with distrubance

	nbr_steps::Int                        # Nbr of steps (integration)

	x0::Vector{T}                         # Initial system state
	Q::Matrix{T}                  # Initial Uncertainty Region State (uncertainty around initial state)
	W::Matrix{T}                  # Disturbances Sets (other uncertainties, not related to the state)

	# Variables generated during Solve, useful throughout (need to be initialized to 0)
	ξ::Matrix{T}                  # Samples (not needed if have follower)
	ξ_propagated::Vector{Matrix{T}}       # Propagated samples
	mop::MultiOrthoPoly                   # Multivariate Orthogonal Basis
	Φ::Matrix{T}                  # Φ Matrix cached, reused several times
	pinv_Φ::Matrix{T}             # pinv_Φ Matrix cached in here too

	# PCE Coefficients to be stored
	C::Vector{Matrix{T}}                  # Contains computed PCE coefficients at all time steps

	# Get Moments from coefficients
	mean::Matrix{T}                          # Lists of means
	var::Matrix{T}                        # Lists of variances

	P::Int                                # Number of polynomial in mop (depends on n and d)
	opts::PCEOptions{T}           	      # Propagator Options
end

# PCEPropagator constructor requiring minimial number of inputs
function PCEPropagator3(model::TO.AbstractModel, f::Any, x0::Vector{T},
	 					Q::AbstractMatrix{T}, W::AbstractMatrix{T},
						opts::PCEOptions{T}) where {T}
	n, m = size(model)
	p, p = size(W)

	nbr_steps = 0
	ξ = zeros(n, n)
	ξ_propagated = [zeros(n, n)]
	op1 = GaussOrthoPoly(opts.d)
	mop = MultiOrthoPoly([op1, op1], 2)
	Φ = zeros(n, n)
	pinv_Φ = zeros(n ,n)
	C = [zeros(n, n)]
	mean = zeros(n, n)
	var = zeros(n, n)
	P = 0
	propagator = PCEPropagator3{T}(n, m, p, model, f, nbr_steps, x0, Q, W,
								  ξ,
								  ξ_propagated,
								  mop,
								  Φ,
								  pinv_Φ,
								  C,
								  mean,
								  var,
								  P,
								  opts)
	return propagator
end




#= FOR REFERENCE USE
# Solver constructor that requires minimal number of inputs.
function MultiplierFunnelSolver(
	data::Funnel73Data,
	opts::MultiplierFunnelOptions{T}
	) where {T}

	n,m = size(data.model)
	p = data.p
	N = data.N

	num_var = opts.exp_opts.disturbance ? n + p : n
	@polyvar X[1:n]                                 # needed
	@polyvar XW[1:num_var]                          # needed later
	# basis = get_basis(X,opts.d,opts.basis_type)   # not needed in multiplier for now
	# basisW = get_basis(XW,opts.d,opts.basis_type) # not needed in multiplier for now
	# b = length(basis)                             # not needed in multiplier for now
	# bW = length(basisW)                           # not needed in multiplier for now

	#ρ = zeros(T,N)                                  # initial ρ values
	ρ = opts.ρ_ini                                   #chosen initial value for ρ, important for multiplier version
	σ = [0.0*X[1] for i=1:N]                  # initial multipliers

	# Function needed right now - needed for Taylor Expansion, will interface differently
	function f(x, u)
	    return [u[1]*cos(x[3]); u[1]*sin(x[3]); u[2]]
	end
	# Need to deal with arrays/Sarrays really here

	solver = MultiplierFunnelSolver{n,m,T}(
		data.n,
		data.m,
		data.p,
		data.N,
		data.model,
		data.model_dist,
		data.X_traj,
		data.U_traj,
		data.Δt_traj,
		data.A_traj,
		data.B_traj,
		data.K_traj,
		data.d_traj,
		data.P_traj,
		opts.M1,
		opts.W,
		X,
		XW,
		f,
		ρ,
		σ,
		opts
		)
	return solver
end

=#






"""
Linear Covariance Analysis Structures
"""









"""
Set-based Uscented Propagation Structures
"""

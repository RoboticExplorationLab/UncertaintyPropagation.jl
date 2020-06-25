# Set-based Unscented Transform Method
# See paper "SAMPLE-BASED ROBUST UNCERTAINTY PROPAGATION FOR ENTRY VEHICLES"
# AAS 2020, Derollez and Manchester

"""
Set-based Unscented Transform Propagation Structures
RUT : Robust Unscented Transform Propagator
"""

abstract type AbstractUncertaintyPropagator end

@with_kw mutable struct RUTOptions{T}

	"Points Extraction Scheme."
	extract_scheme::Symbol=:minimal_boundary

	"Number of Point to Extract."      # directly related to previous item
	nbr_points::Int=:4

	"Minimum Volume Ellipsoid Solver."
	ell_solver::Symbol=:Mosek         #:DRN, :SCS

	"Presence of Disturbances."       # boolean
	disturbances::Bool=false          # if false, only consider state uncertainty

	"Integration Scheme when propagating points."
	int_scheme::typeof(TO.RK3)=RK3            # can also use RK4 and so on

	"Experiment Starting Date."
	t0::T=0.0                              #secondes

	"Experiment Ending Date."
	tf::T=10.0                             #secondes

	"Number of Knot Points."
	N::Int=100
end


mutable struct RUTPropagator{T} <: AbstractUncertaintyPropagator
    # Polynomial Chaos Expansion Solver
	n::Int                                # State size
	m::Int                                # Control size
	p::Int                                # Disturbance size
	N::Int                                # Number of Knot Points
	M::Int                                # Number of points to extract (each time step)
	d::Int                              # Ellipsoid space dension
	model::TO.AbstractModel               # Dynamics model of the system (TrajOpt struct)
	#f::Any                                # Regular Dynamics Function (not adapted to TO)
	#model_dist::TO.AbstractModel         # Dynamics model of the system with distrubance

	x0::Vector{T}                         # Initial system state
	Q::Matrix{T}                  # Initial Uncertainty Region State (uncertainty around initial state)
	W::Matrix{T}                  # Disturbances Sets (other uncertainties, not related to the state)
	dt::T                         # Integration Time Step

	# Containers
	A::Matrix{T}                          # Current A matrix
	b::Vector{T}                          # Current b vector
	A_list::Vector{Matrix{T}}             # List of all A matrices
	b_list::Vector{Vector{T}}             # List of all b vectors
	points::Matrix{T}                     # current extracted points (each column is a point)
	points_list::Vector{Matrix{T}}        # List of all extracted points at all time

	propagated_points::Matrix{T}          # Current points after propagation
	propagated_points_list::Vector{Matrix{T}} # All propagated points

	# Get Moments from coefficients
	#mean::Matrix{T}                          # Lists of means
	#var::Matrix{T}                           # Lists of variances

	opts::RUTOptions{T}           	         # RUT Propagator Options
end

# RUTPropagator constructor requiring minimial number of inputs
function RUTPropagator(model::TO.AbstractModel, x0::Vector{T},
	 					Q::AbstractMatrix{T}, W::AbstractMatrix{T},
						opts::RUTOptions{T}) where {T}

	# get system dimension
	n, m = size(model)
	p, p = size(W)

	# get propagation information
	M = opts.nbr_points                     # get nbr of points (extraction step)

	# get dimension of ellipsoid space (depends on disturbances)
	disturbances = opts.disturbances
	disturbances ? d = n : d = n + p

	# get the integration step based on number of knot points
	t0 = opts.t0                            # start time
	tf = opts.tf                            # end time
	N = opts.N                              # get nbr of knot points
	dt = (tf-t0)/(N-1)                      # get integration time step

	# Initialize containers

	A0 = inv(sqrt(Q))                        # transform ellipsoid representation
    b0 = -A0*x0                              # transform ellipsoid representation

	# Need to do a trick when d!n (here works only if n==d)
	
	A = A0
	b = b0
	A_list = [zeros(d, d) for i=1:N]
	b_list = [zeros(d) for i=1:N]
	points = zeros(d, M)
	points_list = [zeros(d, M) for j=1:N]

	propagated_points = zeros(d, M)
	propagated_points_list = [zeros(d, M) for j=1:N]

	propagator = RUTPropagator{T}(n,
								  m,
								  p,
								  N,
								  M,
								  d,
								  model,
								  x0,
								  Q,
								  W,
								  dt,
								  A,
								  b,
								  A_list,
								  b_list,
								  points,
								  points_list,
								  propagated_points,
								  propagated_points_list,
								  opts)

	return propagator
end

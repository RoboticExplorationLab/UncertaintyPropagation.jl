# Experiment on Duffing Oscillator

include("includes.jl")
include("uncertainty_propagators.jl")
include("uncertainty_methods.jl")

# Define model

function duffing(t,u,p)
    # p = α, β, δ, γ, ω
    du = zeros(2)
    α = -1.0 #p[1]
    β = 1.0 #p[2]
    δ = 0.2 #p[3]
    γ = 0.1 #p[4]
    ω = 1.0 #p[5]
    du[1] = u[2]
    du[2] =  γ*cos(t)-δ*u[2]-α*u[1]-β*u[1]^3 #-0.01*(u[1]-1.0)-0.01*(u[2]-0.0) #+F(t)[1]
    return du
end

struct Duffing{T} <: AbstractModel
    α::T
    β::T
    δ::T
    #γ::T
    #ω::T
end

Base.size(::Duffing) = 2,1

function TO.dynamics(model::Duffing, x, u)
    α = model.α
    β = model.β
    δ = model.δ
    #γ = model.γ
    #ω = model.ω
    ẋ[1] = x[2]
    ẋ[2] = u[1] - δ*x[2]- α*x[1] - β*x[1]^3
    return ẋ
end

# Model and discretization
α = -1.0
β = 1.0
δ = 0.2

# Model and Discretization
model = Duffing(α, β, δ)
n,m = size(model)
f = duffing
x0 = [0.1; 0.1]            # initial state
U0 = [0.01; 0.01]
Q0 = Diagonal(U0.^2)       # initial uncertainty on state (matrix)
S = [1.0; 1.0]             # scaling on coordinates (not used for now, see later)
W = Diagonal(ones(2))      # uncertainty matrix for other parameters

solve_type = :samples      # :samples or :galerkin (coming soon)
M = 2000                   # number of samples
d = 2                      # PCE polynomial degree expansion
uncert_type = [:uniform, :uniform]   # uncertainty type on state coordinates
tf = 10.0
dt = 1e-2
pce_options = PCEOptions(solve_type=solve_type,
						 M=M,
						 d=d,
						 dt=dt,
						 uncert_type=uncert_type,
						 tf=tf)


# Run Polynomial Chaos Expansion propagator on model
pce_propagator = PCEPropagator3(model, f, x0, Q0, W, pce_options)
PCE_solve!(pce_propagator)

mean = pce_propagator.mean
var = pce_propagator.var

Plots.plot(mean[1, :], mean[2, :])
Plots.plot!(mean[2, :])

# Undisturbed trajectory for comparison with the mean
t_sim, Z = rk4(f, x0, 0.0, dt, [0.0,tf])
plot!(Z[1, :], Z[2, :])

# Kernel Density Estimation of the probability density function in dimension "dim"
# method using samples to reconstruct pdf
kernel_density_estimation_plot(pce_propagator; Nr=3000, dim=1, t=600)

# generate_polynomials!(pce_propagator)         # get multivariate orthogonal basis, stored in prop.mop
# pce_propagator.P
# pce_propagator.mop
# show(pce_propagator.mop)
#
# draw_samples!(pce_propagator)                 # draw samples, stored in prop.ξ
# pce_propagator.ξ
# propagate_samples!(pce_propagator)            # propagate samples through model, stored in prop.ξ_propagated
# pce_propagator.ξ_propagated
# compute_Phi_matrix!(pce_propagator)           # compute Φ Matrix, stored in prop.Φ, used for coefficients computation
# pce_propagator.Φ
# pce_propagator.pinv_Φ
# compute_coeff_PCE_full!(pce_propagator)       # compute PCE coefficients at each time step
#
# compute_moments_PCE!(pce_propagator)
# compute_coeff_PCE_step(pce_propagator, 1)

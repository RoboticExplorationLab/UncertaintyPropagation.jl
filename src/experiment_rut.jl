# Test Experiment on RUT method

include("imports.jl")             # get useful packages
include("dynamics.jl")            # get Duffing dynamics
include("rut_struct.jl")          # get rut structures
include("rut_methods.jl")         # get rut methods

# Define Model
α = -1.0
β = 1.0
δ = 0.2

# Model and Discretization
model = Duffing(α, β, δ)
n,m = size(model)
x0 = [0.1; 0.1]            # initial state
U0 = [0.01; 0.01]
Q0 = Diagonal(U0.^2)       # initial uncertainty on state (matrix)
S = [1.0; 1.0]             # scaling on coordinates (not used for now, see later)
W = Diagonal(ones(0))      # uncertainty matrix for other parameters

tf = 5.0
N = 51
opts = RUTOptions(extract_scheme=:minimal_boundary,
                  nbr_points=4,
                  ell_solver=:mosek,
                  disturbances=false,
                  int_scheme=RK4,
                  t0=0.0,
                  tf=tf,
                  N=N)

rut_propagator = RUTPropagator(model, x0, Q0, W, opts)
rut_propagate!(rut_propagator)


plot_center_trajectory(rut_propagator)
animate_center_trajectory(rut_propagator)
animate_ellipses(rut_propagator)
rut_propagator.A_list

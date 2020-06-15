# Experiment on Duffing Oscillator

include("includes.jl")

T = Float64

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

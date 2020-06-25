# Experiment Structures containing information about an experiment

mutable struct ExperimentData{T}
	N::Int	                  # Number of time steps
	t0::T                     # Starting time
	tf::T                     # Final time
	x0::AbstractVector{T}     # Initial state of the system
	#xf::AbstractVector{T}     # Goal state of the system
end

function ExperimentData(;
	N::Int=21,
	t0::T = 0.0
	tf::T=3.0,
	x0::AbstractVector{T}=SVector{3,T}([1., 4., pi/2]),
	) where T
	return ExperimentData{T}(N,t0,tf,x0)
end

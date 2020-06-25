#TODO document


#PCE
#Remove PCEPorpagator.f and interface with TrajOpt and rollout
#See if TrajOpt can handle non-autonomous systems and see how to deal with it
#get nbr steps into the propagator (reused many times)
#check other density estimation methods for reconstruction after compute coefficients
#Try to use the already generated ξ to do the reconstruction job
#Adding other visualization tools (projected ellipsoids plots)

#Next Step is Set-based Unscented Propagation

# Deal with Control input if controlled trajectory (see one_step_rollout function)
# in the Robust Unscented Transform (RUT method)

# Fix the problem with the first step of rut (initial condition not saved right now)
# Check ellipses on the propagation on RUT method

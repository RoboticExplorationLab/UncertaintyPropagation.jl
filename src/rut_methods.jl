# Robust (Set-based) Unscented Transform Porpagation Methods

# Point Extraction Schemes on the Ellipsoid (different schemes = information)

function ellipse2points!(prop::RUTPropagator)
	# Extractin point scheme
    # unpack data
    A = prop.A                        # get matrix representing ellipsoid
    b = prop.b                        # get vector related to ellipsoid's center
	d = prop.d                        # get ellipsoid dimension
    points2 = zeros(d, 2*d+1)
    M = -inv(A)*b                     # get center of the ellipsoid
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    for i =1:d
        points2[:, 2*i-1] = M + W[:, i]
        points2[:, 2*i] = M - W[:, i]
        #@show(points2[:, 2*i])
    end
    points2[:, 2*d+1] = M
	prop.points = points2
    return nothing
end

function ellipse2points_s!(prop::RUTPropagator)
    A = prop.A                        # get matrix representing ellipsoid
    b = prop.b                        # get vector related to ellipsoid's center
    d = prop.d                        # get ellipsoid dimension
    points2 = zeros(d, 2*d+1)
    M = -inv(A)*b                     # get center of the ellipsoid
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    for i =1:d
        σ = W[:, i]/(norm[:, i])
        points2[:, 2*i-1] = M + σ
        points2[:, 2*i] = M - σ
    end
	prop.points = points2
	return nothing
end


function ellipse2points6!(prop::RUTPropagator)
    A = prop.A                        # get matrix representing ellipsoid
    b = prop.b                        # get vector related to ellipsoid's center
    d = prop.d                        # get ellipsoid dimension
    points2 = zeros(d, 6*d+1)
    M = -inv(A)*b                     # get ellipsoid center (coordinates)
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    for i =1:d
        points2[:, 6*i-1] = M + W[:, i]
        points2[:, 6*i-2] = M - W[:, i]
        points2[:, 6*i-3] = M + 0.5*W[:, i]
        points2[:, 6*i-4] = M - 0.5*W[:, i]
        points2[:, 6*i-5] = M + 0.8*W[:, i]
        points2[:, 6*i] = M - 0.8*W[:, i]
    end
    points2[:, 6*d+1] = M
	prop.points = points2
    return nothing
end

function extract_points!(prop::RUTPropagator)
	# general function for extracting points
	opts = prop.opts
	extract_scheme = opts.extract_scheme
	if extract_scheme == :minimal_boundary
		ellipse2points!(prop)
	elseif extract_scheme == :axes_6
		nothing
	elseif extract_scheme == :allover
		nothing
	else
		println("Unknown points extraction scheme entered")
	end
	return nothing
end


function min_volume_covering_ellipsoid!(prop::RUTPropagator, k::Int)
    # Function fitting an ellipsoid to a set of points
    # Minimum Volume Covering Ellipsoid to be obtained
	# k is the time step at which it happens

    # unpack data from propagator
    opts = prop.opts                         # get options
    d = prop.d                               # get ellipsoid dimension
    nbr_points = prop.M                      # get nbr of propagated points
    ell_solver = opts.ell_solver             # get type of ellipsoid solver
	X = prop.propagated_points               # get propagated points

    # set the convex optimization solver
    s = MathOptInterface.LogDetConeTriangle(n)
    #model = Model(with_optimizer(Mosek.Optimizer, QUIET= true, MSK_DPAR_INTPNT_CO_TOL_DFEAS=10^(-20), MSK_DPAR_INTPNT_CO_TOL_PFEAS=10^(-20), MSK_DPAR_INTPNT_CO_TOL_MU_RED = 10^(-20)))
    model = Model(Mosek.Optimizer)
    @variable(model, A[1:d, 1:d], PSD)
    @variable(model, b[1:d])
    @variable(model , t)
    @objective(model, Max, t)
    @constraint(model, con[i = 1:nbr_points], [1.0; A*X[:, i]+b] in SecondOrderCone())
    V = [A[i, j] for j in 1:1:d for i in 1:1:j] #vectorize form of the matrix
    @constraint(model, [t;1.0;V] in s)

    # Add another type of constraint
    #ϵ = 1e-1
    #I = Diagonal(ones(n))
    #B = A-ϵ*I
    #@SDconstraint(model, -B ⪰ 0)
    #@show(B)
    #W = [-B[i, j] for j in 1:1:n for i in 1:1:n]
    #@constraint(model, [0.0;1.0;W] in s)
    #@show(con)
    #MathOptInterface.TimeLimitSec() = 0.5

    # Solve the convex optimization problem
    JuMP.optimize!(model)
    status = JuMP.termination_status(model)                 # get problem status
    A = JuMP.value.(A)                                      # get A value
    b = JuMP.value.(b)                                      # get b value
    prop.A = A                                              # update current A
    prop.b = b                                              # update current b
    return nothing
end


function propagate_samples!(prop::RUTPropagator, k::Int)
	# k is the time step
	M = prop.M          # get number of sample points to propagate
	for i = 1:M
		one_step_rollout!(prop,k,i)
	end
	return nothing
end

function one_step_rollout!(prop::RUTPropagator, k::Int, i::Int)
	# Propagate sample point i forward in time by one time step, at time k
	# get options
	opts = prop.opts
	# get dimensions
	n = prop.n
	m = prop.m
	p = prop.p
	N = prop.N
	dt = prop.dt
	model = prop.model                    # get dynamical system's model
	int_scheme = opts.int_scheme          # get integrator scheme (RK3, RK4...)
	#disturbances = opts.disturbances
	points = prop.points                  # get current points
	# get specific point of interest
	point = points[:, i]
	# propagate the point forward in time by one time step
	z_point = TO.KnotPoint(point, [0.0], dt)
	propagated_point = Array(TO.discrete_dynamics(int_scheme, model, z_point))  # NEED TO DEAL WITH CONTROL
	prop.propagated_points[:, i] = propagated_point               # each column of the matrix is a point
	return nothing
end


function rut_step!(prop::RUTPropagator, k::Int)
	# Entire Step of Robust Unscented Transform Method at time step k
	A = prop.A                                   # get current A
	b = prop.b                                   # get current b
	extract_points!(prop)                        # extract point on ellipsoids according to extract_scheme
	propagate_samples!(prop, k)                  # propagate samples forward one step
	min_volume_covering_ellipsoid!(prop, k)      # run covering ellipsoid code
	# Add values to the prop containers
	prop.A_list[k] = prop.A
	prop.b_list[k] = prop.b
	prop.points_list[k] = prop.points
	prop.propagated_points_list[k] = prop.propagated_points
	return nothing
end


function rut_propagate!(prop::RUTPropagator)
	N = prop.N
	for j=1:N-1
		rut_step!(prop::RUTPropagator, j)
	end
	return nothing
end

# Visualization tools (on Ellipsoids). Will be migrated later and generalized to
# be applicable to all the different UncertaintyPropagator for comparison

# Need to select the right dimensions to plot against each other and then the right
# way to project ellipsoids in the right planes

function plot_center_trajectory(prop::RUTPropagator)
	A_list = prop.A_list
	b_list = prop.b_list
	N = prop.N
	centers_list = [(-inv(A_list[i])*b_list[i]) for i=1:N]
	plt = plot()
	Plots.plot!([centers_list[i][1] for i=1:N],
					[center_list[i][2] for i=1:N])
	display(plt)
	return nothing
end

function animate_center_trajectory(prop::RUTPropagator)
	A_list = prop.A_list
	b_list = prop.b_list
	N = prop.N
	centers_list = [(-inv(A_list[i])*b_list[i]) for i=1:N]
	plt = plot()
	anim = @animate for k=1:1:N
		Plots.scatter!([centers_list[k][1]],
					[center_list[k][2]])
		display(plt)
	end
	return nothing
end


function animate_ellipses(prop::RUTPropagator)
	N = prop.N
	A_list = prop.A_list                  # get ellipsoids' shapes
	b_list = prop.b_list                  # get ~ellipsoids' center info
	plt = plot()
	anim = @animate for j=1:1:N
		angles = 0.0:0.01:2*π
		B = zeros(2, length(angles))
		for i=1:length(angles)
			B[:, i] = [cos(angles[i]) - b_list[j][1], sin(angles[i]) - b_list[j][2]]
			ellipse  = A_list[j][1:2, 1:2] \ B
			Plots.plot!(ellipse[1, :], ellipse[2, :], legend = false, color = :red, linewidth = 3.0)
			display(plt)
		end
	end
	return anim
end


# anim = @animate for j=1:1:200
#     angles = 0.0:0.01:2*pi
#     B = zeros(2, length(angles))
#     for i = 1:1:length(angles)
#         B[:, i] = [cos(angles[i]) - blist[1, j], sin(angles[i]) - blist[2, j]]
#         #B[:, i] = [cos(angles[i]), sin(angles[i])]
#     end
#     ellipse  = Alist[1:2, 1:2, j] \ B
#     #Plots.plot(sol, vars=(1,2), legend = false)
#     Plots.scatter(centerlist[1, 1:5:end], centerlist[2, 1:5:end], linewidth = 2.0, color = :blue, label = false)
#     scatter!([centerlist[1, j]],[centerlist[2, j]] )
#     Plots.plot!(ellipse[1, :], ellipse[2, :], legend = false, color = :red, linewidth = 3.0)
#     xlims!(0.0, 0.7)
#     ylims!(0.0, 0.7)
#     #scatter!(XX[1, :, j], XX[2, :, j])
#     #plot_traj_center(centerlist)
#     #scatter!(E[1, :, j], E[2, :, j])
#     xlabel!("position")
#     ylabel!("velocity")
#     title!("Ellipse propagation step=$j")
# end
# gif(anim, "duff_prop_2.gif", fps = 20)
#
# #ELLIPSE FITTING AT DIFFERENT STEPS
# anim = @animate for j=1:1:200
#     angles = 0.0:0.01:2*pi
#     B = zeros(2, length(angles))
#     for i = 1:1:length(angles)
#         B[:, i] = [cos(angles[i]), sin(angles[i])]
#     end
#     ellipse  = Alist[1:2, 1:2, j+1] \ B
#     if j ==1
#         plot(ellipse[1, :], ellipse[2, :], legend = false, linewidth = 3.0, linestyle = :dash)
#     else
#         plot!(ellipse[1, :], ellipse[2, :], legend = false, linewidth = 3.0, linestyle = :dash)
#     end
#     #scatter!(XX[1, :, j+1], XX[2, :, j+1])
#     title!("Ellipse fitting t=$((j-1)*0.01)")
#     xlabel!("position")
#     ylabel!("velocity")
# end
# gif(anim, "ellipse_fit.gif", fps = 10)




# OLD STYLE (WORKING)
# function rk4(f, y_0, p, dt, t_span)
#     T = t_span[1]:dt:t_span[end]
#     y = zeros(length(T), length(y_0))
#     if length(y_0) == 1
#         y[1, :] = [y_0]
#     else
#         y[1, :] = y_0
#     end
#     for i=1:1:length(T)-1
#         t = T[i]
#         y_star = y[i, :]
#         k1 = f(t, y_star, p)
#         y1 = y_star+k1*dt/2 #intermediate evaluation value
#         k2 = f(t+dt/2, y1, p)
#         y2 = y_star+k2*dt/2
#         k3 = f(t+dt/2, y2, p)
#         y3 = y_star+k3*dt
#         k4 = f(t+dt, y3, p)
#         m = (k1+2*k2+2*k3+k4)/6 #slope average
#         y[i+1, :] = y_star + m*dt
#     end
#     return T, y'
# end
#
# function prop_points_rk(X, t, dt_rk, p, model, dt_e)
#     m = length(X[1, :])
#     Xnew = zeros(size(X))
#     for i=1:1:m
#         t_sim, Z = rk4(model, X[:, i], p, dt_rk, [t, t+dt_e])#integration2(dyna_coeffoff_inplace!, X[:, i], dt)
#         #rk4(dyna_coeffoff, X[:, i], u, 0.001, [0.0, dt])
#         #@show(i)
#         Xnew[:, i] = Z[:, end]
#     end
#     return Xnew
# end

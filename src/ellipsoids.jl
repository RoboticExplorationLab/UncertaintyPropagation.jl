mutable struct Ellipsoid{T}
    Q::AbstractMatrix{T} # Symmetric matrix defining the ellipsoid
    c::AbstractVector{T} # Center of the ellipsoid
    n::Int               # Dimension of the ellipsoid
end

function Ellipsoid(Q::AbstractMatrix{T}, c::AbstractVector{T}=zeros(T,size(Q)[1])) where T
    n = size(Q)[1]
    return Ellipsoid{T}(Q,c,n)
end

function plot_ellipsoid_2D(ell::Ellipsoid;
    dims::Vector{Int}=[1,2], L::Int=200, linewidth::T=3.0, color::Symbol=:blue) where T
    plt = plot()
    plot_ellipsoid_2D!(ell, dims=dims, L=L, linewidth=linewidth, color=color)
    display(plt)
    return nothing
end

function plot_ellipsoid_2D(ells::Vector{E};
    dims::Vector{Int}=[1,2], L::Int=200, linewidth::T=3.0, color::Symbol=:blue) where {E,T}
    plt = plot()
    for ell in ells
        plot_ellipsoid_2D!(ell, dims=dims, L=L, linewidth=linewidth, color=color)
    end
    display(plt)
    return nothing
end

function plot_ellipsoid_2D!(ell::Ellipsoid;
    dims::Vector{Int}=[1,2], L::Int=200, linewidth::T=3.0, color::Symbol=:blue) where T
    # Plot the 2D ellipse corresponding to the equation 1/2*x'*Q*x = 1 centered on point c.
    # If Q and c are higher dimensional than 2 we project on the correct dimensions dims.
    # L is the number of points on the ellipse
    Q = ell.Q
    c = ell.c
    Θ = [2*pi*i/L for i=0:L]
    X = [[cos(θ), sin(θ)] for θ in Θ]
	idims = setdiff(1:ell.n, dims)
	J = Q[ dims, dims]
	L = Q[idims, dims]
	K = Q[idims,idims] + 1e-3*I
    @show K
	Q_shur = J - L'*(K\L)
	# https://math.stackexchange.com/questions/2431159/how-to-obtain-the-equation-of-the-projection-shadow-of-an-ellipsoid-into-2d-plan
    X_ell = [x ./ sqrt(0.5*x'Q_shur*x) for x in X]
    plot!(
        [x[1].+c[dims[1]] for x in X_ell],
        [x[2].+c[dims[2]] for x in X_ell],
        title="Ellipsoid projected on dimensions "*string(dims),
		linewidth=linewidth,
        xlabel="dim "*string(dims[1]),
        ylabel="dim "*string(dims[2]),
		color=color,
        aspect_ratio=1,)
    scatter!(
        [c[dims[1]]], [c[dims[2]]],
        legend=false,
		color=:black
        )
    return nothing
end

function plot_samples_2D!(samples::Vector{Z};
    dims::Vector{Int}=[1,2]) where Z
    # Plot the 2D samples projected on dimensions dims.
    scatter!(
        [s[dims[1]] for s in samples],
		[s[dims[2]] for s in samples],
		color=:red,
		markersize=3.,
        legend=false,
		xlabel="dim "*string(dims[1]),
		ylabel="dim "*string(dims[2]),
		aspect_ratio=1,
        )
    return nothing
end

function plot_monte_carlo_samples_2D!(samples::Vector{Vector{T}}, xk::AbstractVector{T};
    dims::Vector{Int}=[1,2]) where T
    # Plot the 2D samples projected on dimensions dims.
	scatter!(
		[s[dims[1]]+xk[dims[1]] for s in samples],
		[s[dims[2]]+xk[dims[2]] for s in samples],
		color=:red,
		markersize=3.,
		legend=false,
		xlabel="dim "*string(dims[1]),
		ylabel="dim "*string(dims[2]),
		aspect_ratio=1,
		)
	return nothing
end



# # Test struct
# R_ell = rand(5,5)
# Q_ell = R_ell'*R_ell
# Ellipsoid(Q_ell,rand(5),5)
# Ellipsoid(Q_ell,rand(5))
# Ellipsoid(Q_ell)

# # Test plot
# R_ell = randn(5,5)
# Q_ell = R_ell'*R_ell
# c_ell = [1., 2, 3, 4, 5]
# ell = Ellipsoid(Q_ell,c_ell)
# ell2 = Ellipsoid(Q_ell)
# ell3 = Ellipsoid(Diagonal(ones(5))./50)
# ells = [ell, ell2, ell3]
#
# dims_ell = [1,3]
# L_ell = 200
# plot_ellipsoid_2D(ell, dims=dims_ell, L=L_ell)
# plot_ellipsoid_2D(ells, dims=[1,2], L=L_ell)

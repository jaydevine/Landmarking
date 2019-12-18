#-------------------------------------------------------------------------------------------------------------------------------------------------
# Julia script for landmark neural net.
#-------------------------------------------------------------------------------------------------------------------------------------------------
# Import packages.
using Pkg
using CSV
using Flux
using Flux: @epochs
using Flux.Tracker
using LinearAlgebra
using StatsBase
using DataFrames
using Plots
#using CuArrays # If you want to use the GPU instead of CPU. A bit tedious to set up, but things will run much faster.
#using CUDAnative# If you want to use the GPU.
#-------------------------------------------------------------------------------------------------------------------------------------------------
# Some important notes:

# 1) You need to perform a single GPA on the registration derived landmarks (e.g., using geomorph in R),
# then export the Procrustes shape coordinates for import into Julia.

# 2) If you want to retain size information, extract the centroid size vector from the original GPA.

# 3) Because you've trained (or will be training) your network on a certain set of configurations, future work
# will require you to project the registration derived landmarks into this training shape space. Our script
# Training_Projection.R will show you how to do this.

# 4) Let's assume you have all of your data in /path/to/landmarks/.
#-------------------------------------------------------------------------------------------------------------------------------------------------
# Read in datasets.
# Small deformation (ANIMAL) dataset:
small_single = CSV.read("C:/path/to/landmarks/Small_Single_and_Manual_Train.csv")
# Large deformation (SyN) dataset:
# large_single = CSV.read("C:/path/to/landmarks/Large_Single_and_Manual_Train.csv")

# Read in our test set indices for subsetting
split_index = CSV.read("C:/path/to/landmarks/Test_Set_Index.csv")

# x is a matrix of the automated examples; y is a matrix of the manual examples. Here, we have 218 homologous automated
# and manual configurations for training and testing.
x = Matrix(small_single[1:218,2:end])'
y = Matrix(small_single[219:end,2:end])'

# Define training data. We're using N=171 for training.
x_train = reshape(collect(skipmissing(x[:,collect(1:218)[setdiff(1:end, collect(skipmissing(split_index[:,1])))]])), size(x)[1],171)
y_train = reshape(collect(skipmissing(y[:,collect(1:218)[setdiff(1:end, collect(skipmissing(split_index[:,1])))]])), size(y)[1],171)

# Define testing data. We're using N=47 for testing.
x_test = reshape(collect(skipmissing(x[:,collect(skipmissing(split_index[:,1]))])), size(x)[1], 47)
y_test = reshape(collect(skipmissing(y[:,collect(skipmissing(split_index[:,1]))])), size(y)[1], 47)

# If you want to plot the loss, alter this however you wish.
# col_vector = fill("red", 218, 1)
# col_vector[collect(skipmissing(split_index[:,1]))] .= "black"
#-------------------------------------------------------------------------------------------------------------------------------------------------
# Run Bookstein's (1989) TPS interpolation internally. See https://github.com/anj1/ThinPlateSplines.jl.git for more detail.

struct ThinPlateSpline
	λ  # Stiffness.
	x1 # control points
	Y  # Homogeneous control point coordinates
	Φ  # TPS kernel
	d  # Affine component
	c  # Non-affine component
end

is_zero(r::AbstractFloat) = abs(r)<eps(r)
is_zero(r) = false

tps_basis(r::T) where T<:Any  = is_zero(r) ? zero(T) : r*r*log(r)

my_norm(a) = sqrt(sum(a.^2))

# x: matrix of size KxD
tps_kernel(x) = [tps_basis(my_norm(x[i,:] - x[j,:])) for i=1:size(x,1),j=1:size(x,1)]

# find solution of tps transformation
# compute_affine: compute affine component
# (required for some operations but takes additional time.)
function tps_solve(x,y,λ,compute_affine=true)
	# D: number of dimensions
	# K: number of data points
	K,D = size(x)

	# homogeneous coordinates
	X=cat(dims=2,ones(K,1),x)
	Y=cat(dims=2,ones(K,1),y)

	# compute TPS kernel
	Φ = tps_kernel(x)

	# full QR decomposition
	Q,r = qr(X)
	q1 = Q[:,1:(D+1)]
	q2 = Q[:,(D+2):end]

	# warping coefficients
	c = q2*inv(UniformScaling(λ) + q2'*Φ*q2)*q2'*Y

	# affine component
	d = compute_affine ?  r\(q1'*(Y - Φ*c)) : []
	return ThinPlateSpline(λ,x,Y,Φ,d,c)
end

# Thin-plate spline bending energy at minimum
tps_energy(tps::ThinPlateSpline) = tps.λ*tr(tps.c*tps.Y')

#-------------------------------------------------------------------------------------------------------------------------------------------------
# Solve the automated-manual TPS interpolation, extract the bending energy, and store it in a vector.
function tps_loss(x, y)
bend = zeros(218)
  for i = 1:size(x)[2]
      tar = convert(Array{Float64,2}, reshape(collect(x)[:,i], 204, 1))
      tar2 = Flux.chunk(tar, 68)
      tar3 = convert(Array{Float64,2}, permutedims(reshape(hcat(tar2...), (length(tar2[1]), length(tar2)))))
      ref = convert(Array{Float64,2}, reshape(collect(y)[:,i], 204, 1))
      ref2 = Flux.chunk(ref, 68)
      ref3 = convert(Array{Float64,2}, permutedims(reshape(hcat(ref2...), (length(ref2[1]), length(ref2)))))
      # We tested different stiffness coefficients -- 0.001 worked best.
      tps = tps_solve(tar3, ref3, 0.001)
      metric = tps_energy(tps)
      bend[i] = metric
      #println(bend[i])
  end
return mean(bend)
end

# Let's define our architecture. We're going to use 204 units (i.e., 68 landmarks in 3 dimensions = 204 inputs) with
# relative linear units as activation. Make sure they are not in the final layer.
model = Chain(
    Dense(204, 204, relu),
    Dense(204, 204, relu),
    Dense(204, 204, relu),
    Dense(204, 204),
    identity)
penalty() = norm(params(model))
opt = ADAM()

# For plotting:
cb = function ()
  @show(loss(x,y))
end

# Jointly define our training data.
data = [(x_train,y_train)]

# For the first part of the loss, we take the sqrt of Flux's mse, as it improves performance over mse alone.
# For the second part of the loss, we use our previously defined tps_loss function.
loss(x, y) =  sqrt(Flux.mse(model(x), y)) + tps_loss(model(x),y)

# Train the model over 10,000 epochs using a) our loss, b) our dense architecture with RelU, c) the training data,
# and d) the Adam optimizer
@epochs 10000 Flux.train!(loss, Flux.params(model), data, opt, cb = cb)

# If you want to plot things, edit this code below. It's a bit rough...
# l = @layout [a ; b]
# p1 = plot([collect(1:47)], vec(sqrt.(mean(collect(model(x_test) - y_test).^2, dims = 1))), label=["Optimized", "Unoptimized"])
# plot!([collect(1:47)], vec(sqrt.(mean(collect((x_test) - y_test).^2, dims = 1))), label=["Optimized", "Unoptimized"])
# p2 = Plots.scatter([collect(1:47)], vec(sqrt.(mean(collect(model(x_test) - y_test).^2, dims = 1))), color = col_vector, label = ["Train", "Test"])
# plot(p1,p2, layout = l)
# save figure: png(plot(p1,p2, layout = l), "/path/to/landmarks/training.png")
# Plots.scatter(vec(collect(model(x_train))), vec(y_train), label=["Test", "Training"])
# Plots.scatter!(vec(collect(model(x_test))), vec(y_test), label=["Test", "Training"])

# Evaluate the model on your test data.
CSV.write("C:/path/to/landmarks/Test_Landmarks_Optimized.csv", DataFrame(collect(transpose(model(x_test)))))

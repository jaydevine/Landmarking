#-------------------------------------------------------------------------------------------------------------------------------------------------
# Julia script for landmark neural net.
#-------------------------------------------------------------------------------------------------------------------------------------------------
# Import packages.
using Pkg
using CSV
using Flux # Install Flux v0.9.0 via Pkg.add(Pkg.PackageSpec(;name="Flux", version="0.9.0")).
using Flux: @epochs
using LinearAlgebra
using StatsBase
using DataFrames
using Plots
using BSON
using BSON: @save
using BSON: @load
# using CuArrays # If you want to use the GPU instead of CPU. You'll also want CUDAnative. Note, however, that the script would need
# to be adapted to use CUDA arrays.

#-------------------------------------------------------------------------------------------------------------------------------------------------
# The only notes you need to read, and the only variables you need to alter, are between these lines.
# 1. NOTES:

# A) GPA and project your training/testing landmark data into the same space. Follow https://github.com/jaydevine/Landmarking/blob/master/R/GPA_and_Projection.R
# for instructions on how to do this.

# B) If you want to retain size information, extract the centroid size vector from the original GPA. In general, registration-based
# size information (e.g., volumes, linear distances, centroid size) is more reliable than the shape information.

# C) Let's assume you have all of your data in C:/path/to/wd/.

# 2. VARIABLES:

# A) Read in datasets. First n rows as automated data and last n rows as manual data.
Train = CSV.read("C:/path/to/wd/SyN_MAN_Train_Orp_Revised.csv")
Test = CSV.read("C:/path/to/wd/SyN_Test_Orp_Revised.csv")

# B) Split dataset up, convert to matrix, and transpose. In this example, x is a matrix of N=170 automated training examples and y is a matrix of
# N=170 manual examples. x2 is a matrix of N=47 automated testing examples. All landmark configurations are homologous.
x = Matrix(Train[1:170,2:end])'
y = Matrix(Train[171:end,2:end])'
x2 = Matrix(Test[1:end,2:end])'

# C) Define training and testing data
x_Train = reshape(collect(x), size(x)[1], size(x)[2])
y_Train = reshape(collect(y), size(y)[1], size(y)[2])
x_Test = reshape(collect(x2), size(x2)[1], size(x2)[2])

# D) Define whether your data are in two or three (K) dimensions.
K = 3
#----------------------------------------------------------------------------------------------------------------------------
# Determine how many landmarks (P) there are.
P = trunc(Int, size(x_Train)[1]/K)

# Define loss functions.
# First, we define an RMSE loss, as it improves performance over MSE alone.
function RMSE(x,y)
 	return sqrt(sum((x .- y).^2)) * 1 // size(x,1)
end

# Second, we define a TPS loss. See https://github.com/anj1/ThinPlateSplines.jl.git for more detail.
Is_Zero(r::AbstractFloat) = abs(r)<eps(r)
Is_Zero(r) = false
TPS_Basis(r::T) where T<:Any  = Is_Zero(r) ? zero(T) : r*r*log(r)
My_Norm(a) = sqrt(sum(a.^2))

# Define TPS_Kernel for x.
TPS_Kernel(x) = [TPS_Basis(My_Norm(x[i,:] - x[j,:])) for i=1:size(x,1),j=1:size(x,1)]

# Solve the TPS interpolation and return an energy:
function TPS_Solve(x,y,λ)
	Bend = zeros(size(x,2))
    # Number of shape dimensions is P landmarks * K dimensions
	for i = 1:size(x,2)
        X = x[:,i]
        Y = y[:,i]
        X = convert(Array{Float64,2}, reshape(collect(X), P, K))
        Y = convert(Array{Float64,2}, reshape(collect(Y), P, K))

		# Create homogeneous coordinates.
		X_Hom = cat(dims=2,ones(P,1),X)
		Y_Hom = cat(dims=2,ones(P,1),Y)

		# Compute TPS kernel.
		Φ = TPS_Kernel(X)

		# QR decomposition.
		Q,R = qr(X_Hom)
		Q1 = Q[:,1:(K+1)]
		Q2 = Q[:,(K+2):end]

		# Calculate warping coefficients.
		C = Q2*inv(UniformScaling(λ) + Q2'*Φ*Q2)*Q2'*Y_Hom

		# Compute bending energy at a minimum and store in vector.
		Energy = λ*tr(C*Y_Hom')
		Bend[i] = Energy
	end
    # Calculate mean bending energy.
	return mean(Bend)
end

# Define losses together: 0.001 worked best.
function Loss(x,y)
	return RMSE(Model(x), y) + TPS_Solve(Model(x),y,0.001)
end

#-------------------------------------------------------------------------------------------------------------------------------------------------
# Define the architecture. We're going to use P*K units in each layer with relative linear units as activation.
# Make sure they are not in the final layer.
Model = Chain(
    Dense(P*K, P*K, relu),
    Dense(P*K, P*K, relu),
    Dense(P*K, P*K, relu),
    Dense(P*K, P*K),
    identity)

# Track parameters used to calculate gradient of loss function and use ADAM as optimizer.
Ps = Flux.params(Model)
Opt = ADAM()

# Jointly define training data.
Data = [(x_Train,y_Train)]

# Train the model over 10,000 epochs using our shape loss, densely connected architecture with RelU, the training data,
# and the Adam optimizer
@epochs 10000 Flux.train!(Loss, Ps, Data, Opt)

# Save the model and model weights with BSON package to avoid retraining your network.
Weights = Tracker.data.(params(Model));
@save "Model.bson" Model # look at saving here: https://pkg.julialang.org/docs/Flux/QdkVy/0.9.0/saving/
@save "Weights.bson" Weights

# If you want to load the parameters back into the model to avoid retraining, define the model above and skip the training step. Then,
@load "Weights.bson" Weights
Flux.loadparams!(Model, Weights)

# Evaluate the model on test data.
#CSV.write("C:/path/to/wd/Test_Predictions.csv", DataFrame(collect(transpose(Model(x_Test)))))

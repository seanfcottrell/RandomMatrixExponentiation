using LinearAlgebra

using Manifolds

using Plots; pyplot()

using Distributions


#Visualization of density of exponentiated matrices sampled from Colored Gaussian Symmetric Ensemble


d = 2
N = 10000

 

e = range(0, stop = 15, length = 200)

 

ℙd = SymmetricPositiveDefinite(d) #Define the SPD manifold

 

Imat = Matrix(I, d, d) # this will be a matrix of bools

Imat_float = float.(Imat) # this will convert it to floats

Sigma_1 = Imat_float # point on SPD manifold = Identity Matrix

 

E = zeros(d, d, N)

S_1 = zeros(d, d, N)

Eig_1 = zeros(d, 1, N)

 

 

for i = 1:N

   

     A = randn(d,d) # generate gaussian random matrix

     B = transpose(A) # transpose A

     E[:, :, i] = (A+B)/2 # E is GOE matrix

   

end

 

for j = 1:N

   

    S_1[:, :, j] = exp(ℙd, Sigma_1, E[:, :, j])

    Eig_1[:, :, j] = eigvals(S_1[:, :, j])

   

end

 

E1 = Eig_1[:]

 

#Now add the analytic eigenvalue density for d = 1

 

c_H = (2*pi)^(-1/2) #normalizing constant

# for d = 1, the product of gammas = 1

 

g(e) = c_H * exp(-((log(e)^2)/2 + log(e)))

# for d =1 there is no spacing between eigenvalues

 

Eigenvalues_S_1_Distribution = Plots.histogram(E1,

    bins = 400,

    normalize=:pdf, 

    xlims = (0,15),

    xlabel = "Value",

    ylabel = "Density",

    title = "Eigenvalue Density",

    label = "Sample Data"

)
 

plot!(g, xlims = (0,15), label = "Analytic Density Function")



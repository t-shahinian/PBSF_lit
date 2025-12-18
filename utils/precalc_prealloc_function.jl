## MAYBE LATER: RENAME INTO SOMETHING BETTER BECAUSE THIS HOLDS ALL THE UTILS FOR THE SIMULATION CLEANUP


#let precalc_alloc_function equal a variable called pre_calc_alloc_values so i can access it
#let initialization_function equal a variable called initialization_values so i can access it 

### meeting: 1: address the ASKs
### i didn't update the output, will do that as i clean up the rest of the code and see what is needed

### function in julia which does pre-calculation and pre-allocation

### make p, q, n a function of x and y so we don't need to input it (make this funciton ahve as few inputs as possible)
### add a check to make sure dimensions of x and y match (dimension mismatching check)
## use spzeros instead of sparse when large matrix (stick with hwt dr zhang has, but keep this in mind)
## make a clean version of the jupyter and keep the old one so i can compare
## dimension checking, if we expect positive, but they put negative, put a warning

# create new jupyter file with my updated code to replace the old one (just need to pull the original and start rewriting in there/replacing. keep the old one
# for reference)

## i added p q and n beacuse they are specified in the .jld file. unsure if we want to get rid of that ASK
function precalc_alloc_function(coords,X,Y,ϕ_sam,K; μβ=nothing,μΛ= nothing,VΛ=nothing,Vβ=nothing,aΣ=2,bΣ=nothing,m=10, N_sam = 20000)    # Assign defaults inside the function after knowing p, q, K
    # initializing n, p, q and dimension checking (only checked for n, see if there's a way to check p and q???)
    n = size(X, 1)
    p = size(X, 2)
    q = size(Y, 2)

    if size(Y, 1) != n
        error("Dimension mismatch: X is $(size(X,1))×$(size(X,2)) but Y is $(size(Y,1))×$(size(Y,2)). Number of rows must match.")
    end

    # you can now safely use n, p, q
    println("n = $n, p = $p, q = $q")
    
    # setting initial values if they aren't given
    if μβ === nothing
        μβ = fill(0.0, p, q)
    end

    if bΣ === nothing
        bΣ = fill(1, q)
    end

    if μΛ === nothing
        μΛ = fill(0.0, K, q)
    end

## Precalculation I ##
    # 2.A: Construct maximin ordering
    # Convert Julia matrix to R
    @rput coords
    # Use order_maxmin_exact from GPvecchia
    R"""
    ordered_indices <- order_maxmin_exact(coords)
    """
    # Get results back to Julia
    @rget ordered_indices;

    # reorder data #
    Y_ord = Y[ordered_indices, :];
    X_ord = X[ordered_indices, :];
    coords_ord = coords[ordered_indices, :];
    # F_ord = F0[ordered_indices, :]; # not needed here, used in model comparison step?? maybe delete (commented out because not needed here)

    # 2.B: build nearest neighbor, might take some time 
    # ? m = 15
    # ? K = 2 # No. factors
    NN = BuildNN(coords_ord, m);

## Precalculation II ##
    # 2.C: Calculate Cholesky decompositions V_lambda and V_beta (?? i think it's just initializing)
    # priors #
    # ? K = 2;
    # ? μβ = fill(0.0, p, q); 

    # if we aren't given v lambda and v beta, we prespecify the inverses (since that is what we need). otherwise we must calculate the breakdown (? ASK // how to calculate if it's specified?)
        # need to add the code to calculate if that's true
    if VΛ === nothing
        inv_VΛ = zeros(Float64, K, K)
        inv_LΛ = sparse(1:K, 1:K, fill(0.1, K));
    else ## ASK check this part. don't know exactly how this is to be used
        LΛ = cholesky(VΛ).L # gives lower cholesky decomp
        inv_LΛ = inv(LΛ)
        inv_VΛ = inv(VΛ)
    end

    if Vβ === nothing
        inv_Vr = zeros(Float64, p, p) # is this beta? Vr?? ASK
        inv_Lr = sparse(1:p, 1:p, fill(0.1, p));
    else # ASK HERE 
        # User provided Vβ → compute Cholesky and inverse Cholesky
        Lr = cholesky(Vβ).L      # lower-triangular Cholesky factor
        inv_Lr = inv(Lr)         # inverse Cholesky factor
        inv_Vr = inv(Vβ)         # full precision matrix (optional, if needed)
    end

    # ? inv_Vr = zeros(Float64, p, p); #inv_Vr[1, 1] = 1; ??? from here, determine wha tthe original v lambda and v beta are, so we can put those as initial values in the function definition
    # ? μΛ = fill(0.0, K, q);
    # ? inv_VΛ = zeros(Float64, K, K);
    # ? aΣ = 2;
    # ? bΣ = fill(1, q); # later will want to add something to estimate sigma//initial starting point, here we jsut assumed something 
    # ? inv_Lr = sparse(1:p, 1:p, fill(0.1, p));
    # ? inv_LΛ = sparse(1:K, 1:K, fill(0.1, K));
    # ? ϕ_sam = [4.0, 6.0]; #[6.0, 9.0] #practical range from 3000 to 150

    # 2.D: A and D #
    nIndx = length(NN.nnIndx)
    A = [Array{Float64}(undef, nIndx) for i in 1:K]
    D = [Array{Float64}(undef, n) for i in 1:K]
    I_A = [spzeros(n, n) for i in 1:K] # 1-A
    A_new = [Array{Float64}(undef, nIndx) for i in 1:K]
    D_new = [Array{Float64}(undef, n) for i in 1:K]
    I_A_new = [spzeros(n, n) for i in 1:K]
        
    nnIndx_col = vcat(NN.nnIndx, 1:n) # Index of columns (for getting the original nxn matrix back; we stored it as spzeros so it takes less space)
    nnIndx_row = zeros(Int64, 0) # Index of rows initialization
    for i in 2:m
        nnIndx_row = vcat(nnIndx_row, fill(i, i-1))
    end
    nnIndx_row = vcat(nnIndx_row, repeat((m + 1):n, inner = m), 1:n)
    dim_invD = n*q;

## Pre-allocation ##
# Pre-allocation for MCMC samples
    # ? N_sam = 20000;

    # Pre-allocation for F updates
    nsam = (n * q) + (K * n);
    Xtilde = spzeros(Float64, (q+K)*n, n*K);
    Ytilde = Array{Float64}(undef, nsam);
    F_sam = Array{Float64, 2}(undef, n, K);
    Fqr = qr(zeros(Float64, n, K)) 
    inv_sqrt_Σ_diag = Vector{Float64}(undef, q)
    invD_ele = Vector{Float64}(undef, n*q)
    invD = spzeros(Float64, n*q, n*q)
    invΣhalf = spzeros(Float64, q, q)
    F_m = zeros(Float64, 1, K);
    v = Array{Float64}(undef, nsam); # storage samples from standard normal

    # Pre-allocation for γ, Σ updates
    Ystar = vcat(Y_ord, inv_Lr * μβ, fill(0.0, K, q)); # (NTotal+p+K) by q matrix
    Xstar = vcat(hcat(X_ord, F_sam), hcat(inv_Lr, spzeros(p, K)), hcat(spzeros(K, p), sparse(1:K, 1:K, 1.0)));
    μγstar = vcat(μβ, μΛ); #invVγstar = fill(0.0, p + K, p + K);
    invVγstar = cholesky(sparse(1.0I, p+K, p+K)); # doesn't fine in-place update for this (changed this to 1.0I)
    u = Array{Float64}(undef, (p + K) * q);  # Pre-allocate space for random samples;
    Y_Xm = spzeros(n + p + K, q); # store the residual
    bstar = fill(0.0, q); astar = aΣ + 0.5 * (n);

# specifying what to return which will be used in future steps/functions
    return (Fqr=Fqr,F_m=F_m,F_sam=F_sam,NN.nnIndx=NN.nnIndx,NN.nnDist=NN.nnDist,NN.nnIndxLU=NN.nnIndxLU,μγstar=μγstar,invVγstar=invVγstar,u=u,astar=astar,bstar=bstar,Y_Xm=Y_Xm,bΣ=bΣ,Xstar=Xstar,Ystar=Ystar,v=v,D=D,I_A=I_A,nnIndx_row=nnIndx_row,nnIndx_col=nnIndx_col,Xtilde=Xtilde,Ytilde=Ytilde,A=A,D=D,Y_ord=Y_ord,X_ord=X_ord,coords_ord=coords_ord,NN=NN, K=K, n=n,p=p,q=q, inv_sqrt_Σ_diag=inv_sqrt_Σ_diag,invD_ele=invD_ele,invD=invD,ϕ_sam=ϕ_sam)
end

## add more to the above return as i find necessary


function initialization(pre_calc_alloc_values)
    ## Initalization (some are optional) ##
    β0 = (pre_calc_alloc_values.X_ord'pre_calc_alloc_values.X_ord)\(pre_calc_alloc_values.X_ord'pre_calc_alloc_values.Y_ord);
    Residuals = pre_calc_alloc_values.Y_ord - pre_calc_alloc_values.X_ord*β0;
    #γ_sam = vcat((X'X)\(X'Y), fill(0.0, K, q));
    reordered_result = reorder_svd_by_spatial_range(Residuals, pre_calc_alloc_values.coords_ord, pre_calc_alloc_values.K, 1000, 11) # ?? do we want to change 1000 and 11? figure out
    γ_sam = vcat((pre_calc_alloc_values.X_ord'pre_calc_alloc_values.X_ord)\(pre_calc_alloc_values.X_ord'pre_calc_alloc_values.Y_ord), reordered_result);
    Σ_sam = [var(Residuals[:, j]) for j in 1:size(Residuals, 2)];

    # force the initial column of F to start from the noisy factor 
    #γ_sam = vcat((X'X)\(X'Y), Λ[[2, 1], :]);
    #Σ_sam = [0.5, 1, 0.4, 2, 0.3, 2.5, 3.5, 0.45, 1.5, 0.5]; #fill(1.0, q);
    return (Σ_sam=Σ_sam, γ_sam=γ_sam)
end




function results(result_path::AbstractString)
    # Ensure the directory exists
    if !isdir(result_path)
        mkpath(result_path)
        println("Created directory: $result_path")
    end

    # Build full file paths
    gamma_file = joinpath(result_path, "γ_sam.csv")
    sigma_file = joinpath(result_path, "Σ_sam.csv")
    F_file = joinpath(result_path, "F_sam.csv")

    # Save initial zero data
    writedlm(gamma_file, fill(0.0, 1, pre_calc_alloc_values.q), ", ")
    writedlm(sigma_file, 0.0, ", ")
    writedlm(F_file, fill(0.0, 1, pre_calc_alloc_values.K), ", ")

    println("Files saved to: $result_path")
end


# function for sampling f, function for sampling gamma and sigma, and then a function that wraps all of this chunk together
function sample_f (pre_calc_alloc_values,initialization_values,l)
        # Build the matrix D_Sigma_o^{1/2} #
    # Compute inverse square root in-place
    @. pre_calc_alloc_values.inv_sqrt_Σ_diag = 1 / sqrt(initialization_values.Σ_sam)
    # Efficient broadcasting for invD_ele
    pre_calc_alloc_values.invD_ele .= repeat(pre_calc_alloc_values.inv_sqrt_Σ_diag, inner=n)
    # Update sparse matrices
    pre_calc_alloc_values.invD = sparse(1:(pre_calc_alloc_values.n*pre_calc_alloc_values.q), 1:(pre_calc_alloc_values.n*pre_calc_alloc_values.q), pre_calc_alloc_values.invD_ele)
    invΣhalf = sparse(1:pre_calc_alloc_values.q, 1:pre_calc_alloc_values.q, 1 ./ sqrt.(initialization_values.Σ_sam)) # ??? is this using the same inverse sigma as in preallocation? it doesn't seem like it
    
    if l == 1
        for k in 1:pre_calc_alloc_values.K
            getAD(pre_calc_alloc_values.coords_ord, pre_calc_alloc_values.NN.nnIndx, pre_calc_alloc_values.NN.nnDist, pre_calc_alloc_values.NN.nnIndxLU, pre_calc_alloc_values.ϕ_sam[k, l], 0.5, pre_calc_alloc_values.A[k], pre_calc_alloc_values.D[k]);
            pre_calc_alloc_values.I_A[k] = sparse(pre_calc_alloc_values.nnIndx_row, pre_calc_alloc_values.nnIndx_col, vcat(-pre_calc_alloc_values.A[k], ones(pre_calc_alloc_values.n)));
        end
    end
    pre_calc_alloc_values.Ytilde .= vcat(pre_calc_alloc_values.invD * vec(pre_calc_alloc_values.Y_ord - pre_calc_alloc_values.X_ord * initialization_values.γ_sam[1:p, :]), zeros(pre_calc_alloc_values.K * pre_calc_alloc_values.n));
    pre_calc_alloc_values.Xtilde .= vcat(kron(invΣhalf * sparse(transpose(initialization_values.γ_sam[(pre_calc_alloc_values.p + 1):(pre_calc_alloc_values.p + pre_calc_alloc_values.K), :])), 
                            sparse(1:pre_calc_alloc_values.n, 1:pre_calc_alloc_values.n, ones(pre_calc_alloc_values.n))),
             blockdiag([Diagonal(1 ./ sqrt.(pre_calc_alloc_values.D[k])) * pre_calc_alloc_values.I_A[k] for k in 1:pre_calc_alloc_values.K]...));
       
    # use LSMR to generate sample of F # 
    randn!(pre_calc_alloc_values.v)  # Fills v with standard normal samples
    elapsed_time = @elapsed begin
    pre_calc_alloc_values.F_sam .= reshape(lsmr(pre_calc_alloc_values.Xtilde, collect(pre_calc_alloc_values.Ytilde) + pre_calc_alloc_values.v), :, pre_calc_alloc_values.K);
    end
    
    # Print the elapsed time
    #println("Iteration $l: Time taken for lsmr step = $elapsed_time seconds")
    
    pre_calc_alloc_values.F_m .= mean(pre_calc_alloc_values.F_sam, dims = 1);
    ProgressMeter.F_sam .-= pre_calc_alloc_values.F_m;

    # Perform thin QR on the tall slice of F_sam in-place
    pre_calc_alloc_values.Fqr = qr!(pre_calc_alloc_values.F_sam)  # Note the ! for in-place modification

    # Assign scaled Q to F_samples
    pre_calc_alloc_values.F_sam .= Matrix(pre_calc_alloc_values.Fqr.Q)
    rmul!(pre_calc_alloc_values.F_sam, sqrt(pre_calc_alloc_values.n-1))   # scale in-place
    return ()
end






function sample_gamma_sigma ()
        
    # one function here (sampling sig and gamma)
    # Sample Σ and γ#
    pre_calc_alloc_values.Xstar[1:pre_calc_alloc_values.n, (pre_calc_alloc_values.p+1):(pre_calc_alloc_values.p+pre_calc_alloc_values.K)] .= F_sam; # update Xstar with F_sam
    
    # use MNIW to sample γ Σ #
    invVγstar = cholesky(pre_calc_alloc_values.Xstar'pre_calc_alloc_values.Xstar); # ??? does this use the original invVγstar in precalc alloc function? 
    #mul!(μγstar, transpose(Xstar), Ystar); μγstar = invVγstar \ μγstar;
    pre_calc_alloc_values.μγstar .= invVγstar \ (pre_calc_alloc_values.Xstar'pre_calc_alloc_values.Ystar); # ??? this one does use original from precalc prealloc
    pre_calc_alloc_values.Y_Xm .= pre_calc_alloc_values.Ystar - pre_calc_alloc_values.Xstar * μγstar;      
    pre_calc_alloc_values.bstar .= [pre_calc_alloc_values.bΣ[ind] + 0.5 * (norm(pre_calc_alloc_values.Y_Xm[:, ind])^2) for ind in 1:pre_calc_alloc_values.q]; 
    initialization_values.Σ_sam .= [rand(InverseGamma(pre_calc_alloc_values.astar, pre_calc_alloc_values.bstar[ind]), 1)[1] for ind in 1:pre_calc_alloc_values.q];          # sample Σ
    randn!(pre_calc_alloc_values.u)  # Fills u with standard normal samples
    initialization_values.γ_sam .= (invVγstar.U \ reshape(pre_calc_alloc_values.u, (pre_calc_alloc_values.p + pre_calc_alloc_values.K), pre_calc_alloc_values.q)) * 
                    Diagonal(sqrt.(initialization_values.Σ_sam)) + pre_calc_alloc_values.μγstar;          # sample γ 
end





function run_mcmc(precalc_alloc::NamedTuple, initialization::NamedTuple, N_sam::Int, result_path::AbstractString)
    using DelimitedFiles, ProgressMeter, Random
    
    # Ensure the results directory exists
    if !isdir(result_path)
        mkpath(result_path)
    end

    # Initialize progress bar
    prog = Progress(N_sam, 1, "Running MCMC...", 50)
    
    # Seed RNG for reproducibility
    Random.seed!(11)

    # Precompute file paths
    f_file     = joinpath(result_path, "F_sam.csv")
    sigma_file = joinpath(result_path, "Σ_sam.csv")
    gamma_file = joinpath(result_path, "γ_sam.csv")

    # Optionally, write headers or initial zeros
    writedlm(f_file, fill(0.0, 1, precalc_alloc.K), ", ")
    writedlm(sigma_file, fill(0.0, 1, precalc_alloc.q), ", ")
    writedlm(gamma_file, fill(0.0, 1, precalc_alloc.q), ", ")

    # Loop over MCMC iterations
    for l in 1:N_sam
        # 1️⃣ Sample F
        F_result = sample_f(precalc_alloc, initialization, l)
        F_sam = F_result.F_sam

        # Append F_sam to CSV
        open(f_file, "a") do io
            writedlm(io, F_sam, ", ")
        end

        # 2️⃣ Sample γ and Σ
        sample_gamma_sigma(precalc_alloc, initialization, F_sam)

        # Append Σ_sam and γ_sam to CSV
        open(sigma_file, "a") do io
            writedlm(io, initialization_values.Σ_sam, ", ")
        end

        open(gamma_file, "a") do io
            writedlm(io, initialization_values.γ_sam, ", ")
        end

        # Update progress
        next!(prog)
    end
end



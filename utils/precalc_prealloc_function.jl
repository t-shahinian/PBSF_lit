## MAYBE LATER: RENAME INTO SOMETHING BETTER BECAUSE THIS HOLDS ALL THE UTILS FOR THE SIMULATION CLEANUP

## 1/4: ask gpt to compare the imported files and see if it runs as expected. as of right now im getting

## 1/5: everything runs and produces (reasonable) output
## issue: csv files don't quite match. they are close though (off by a couple of decimals)
## investigate that further
## i only copy pasted the post processing stuff. i need to fix it tomorrow

## an error for the run mcmc chunk, but i think it should be pretty close


# FOR DEBUGGING: go line by line and compare to original code


# NOTE: on the outside we need to define pre_calc_alloc_values,initialization_values as the functions (i.e. pre_calc_alloc_values=pre_calc_alloc_function(...))

#let precalc_alloc_function equal a variable called pre_calc_alloc_values so i can access it
#let initialization_function equal a variable called initialization_values so i can access it 


### function in julia which does pre-calculation and pre-allocation

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
    NN = BuildNN(coords_ord, m);

## Precalculation II ##
    # 2.C: Calculate Cholesky decompositions V_lambda and V_beta (?? i think it's just initializing)

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

    # Pre-allocation for F updates
    nsam = (n * q) + (K * n);
    Xtilde = spzeros(Float64, (q+K)*n, n*K);
    Ytilde = Array{Float64}(undef, nsam);
    F_sam = Array{Float64, 2}(undef, n, K);
    # tati commented out (don't think this is needed, but confirm ???) Fqr = qr(zeros(Float64, n, K)) 
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
    # tati commented out (don't think this is needed, but confirm ???) invVγstar = cholesky(sparse(1.0I, p+K, p+K)); # doesn't fine in-place update for this (changed this to 1.0I)
    u = Array{Float64}(undef, (p + K) * q);  # Pre-allocate space for random samples;
    Y_Xm = spzeros(n + p + K, q); # store the residual
    bstar = fill(0.0, q); astar = aΣ + 0.5 * (n);

# specifying what to return which will be used in future steps/functions
    return (invΣhalf=invΣhalf,N_sam=N_sam,F_m=F_m,F_sam=F_sam,NN=NN,μγstar=μγstar,u=u,astar=astar,bstar=bstar,Y_Xm=Y_Xm,bΣ=bΣ,Xstar=Xstar,Ystar=Ystar,v=v,D=D,I_A=I_A,nnIndx_row=nnIndx_row,nnIndx_col=nnIndx_col,Xtilde=Xtilde,Ytilde=Ytilde,A=A,Y_ord=Y_ord,X_ord=X_ord,coords_ord=coords_ord, K=K, n=n,p=p,q=q, inv_sqrt_Σ_diag=inv_sqrt_Σ_diag,invD_ele=invD_ele,invD=invD,ϕ_sam=ϕ_sam)
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




function results(result_path::AbstractString,pre_calc_alloc_values)
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
    return (gamma_file=gamma_file, sigma_file=sigma_file, F_file=F_file,result_path=result_path)
end


# function for sampling f, function for sampling gamma and sigma, and then a function that wraps all of this chunk together
function sample_f(pre_calc_alloc_values,initialization_values,l)
        # Build the matrix D_Sigma_o^{1/2} #
    # Compute inverse square root in-place
    @. pre_calc_alloc_values.inv_sqrt_Σ_diag = 1 / sqrt(initialization_values.Σ_sam)
    # Efficient broadcasting for invD_ele
    pre_calc_alloc_values.invD_ele .= repeat(pre_calc_alloc_values.inv_sqrt_Σ_diag, inner=pre_calc_alloc_values.n)
    # Update sparse matrices
    pre_calc_alloc_values.invD .= sparse(1:(pre_calc_alloc_values.n*pre_calc_alloc_values.q), 1:(pre_calc_alloc_values.n*pre_calc_alloc_values.q), pre_calc_alloc_values.invD_ele)
    pre_calc_alloc_values.invΣhalf .= sparse(1:pre_calc_alloc_values.q, 1:pre_calc_alloc_values.q, 1 ./ sqrt.(initialization_values.Σ_sam)) # ??? is this using the same inverse sigma as in preallocation? it doesn't seem like it
    
    if l == 1
        for k in 1:pre_calc_alloc_values.K
            getAD(pre_calc_alloc_values.coords_ord, pre_calc_alloc_values.NN.nnIndx, pre_calc_alloc_values.NN.nnDist, pre_calc_alloc_values.NN.nnIndxLU, pre_calc_alloc_values.ϕ_sam[k, l], 0.5, pre_calc_alloc_values.A[k], pre_calc_alloc_values.D[k]);
            pre_calc_alloc_values.I_A[k] .= sparse(pre_calc_alloc_values.nnIndx_row, pre_calc_alloc_values.nnIndx_col, vcat(-pre_calc_alloc_values.A[k], ones(pre_calc_alloc_values.n)));
        end
    end
    pre_calc_alloc_values.Ytilde .= vcat(pre_calc_alloc_values.invD * vec(pre_calc_alloc_values.Y_ord - pre_calc_alloc_values.X_ord * initialization_values.γ_sam[1:p, :]), zeros(pre_calc_alloc_values.K * pre_calc_alloc_values.n));
    pre_calc_alloc_values.Xtilde .= vcat(kron(pre_calc_alloc_values.invΣhalf * sparse(transpose(initialization_values.γ_sam[(pre_calc_alloc_values.p + 1):(pre_calc_alloc_values.p + pre_calc_alloc_values.K), :])), 
                            sparse(1:pre_calc_alloc_values.n, 1:pre_calc_alloc_values.n, ones(pre_calc_alloc_values.n))),
             blockdiag([Diagonal(1 ./ sqrt.(pre_calc_alloc_values.D[k])) * pre_calc_alloc_values.I_A[k] for k in 1:pre_calc_alloc_values.K]...));
       
    # use LSMR to generate sample of F # 
    randn!(pre_calc_alloc_values.v)  # Fills v with standard normal samples
    elapsed_time = @elapsed begin
    pre_calc_alloc_values.F_sam .= reshape(lsmr(pre_calc_alloc_values.Xtilde, collect(pre_calc_alloc_values.Ytilde) + pre_calc_alloc_values.v), :, pre_calc_alloc_values.K);
    end
    
    # Print the elapsed time
    #println("Iteration $l: Time taken for lsmr step = $elapsed_time seconds")
    
    pre_calc_alloc_values.F_m .= mean(pre_calc_alloc_values.F_sam, dims = 1); # ??? should dimensions be 1, or general?
    pre_calc_alloc_values.F_sam .-= pre_calc_alloc_values.F_m; 

    # Perform thin QR on the tall slice of F_sam in-place
    Fqr = qr!(pre_calc_alloc_values.F_sam)  # Note the ! for in-place modification

    # Assign scaled Q to F_samples
    pre_calc_alloc_values.F_sam .= Matrix(Fqr.Q)
    rmul!(pre_calc_alloc_values.F_sam, sqrt(pre_calc_alloc_values.n-1))   # scale in-place
    
end






function sample_gamma_sigma(pre_calc_alloc_values,initialization_values)
        
    # one function here (sampling sig and gamma)
    # Sample Σ and γ#
    pre_calc_alloc_values.Xstar[1:pre_calc_alloc_values.n, (pre_calc_alloc_values.p+1):(pre_calc_alloc_values.p+pre_calc_alloc_values.K)] .= pre_calc_alloc_values.F_sam; # update Xstar with F_sam
    
    # use MNIW to sample γ Σ #
    invVγstar = cholesky(pre_calc_alloc_values.Xstar'pre_calc_alloc_values.Xstar);  
    #mul!(μγstar, transpose(Xstar), Ystar); μγstar = invVγstar \ μγstar;
    pre_calc_alloc_values.μγstar .= invVγstar \ (pre_calc_alloc_values.Xstar'pre_calc_alloc_values.Ystar); # ??? this one does use original from precalc prealloc
    pre_calc_alloc_values.Y_Xm .= pre_calc_alloc_values.Ystar - pre_calc_alloc_values.Xstar * pre_calc_alloc_values.μγstar;      
    pre_calc_alloc_values.bstar .= [pre_calc_alloc_values.bΣ[ind] + 0.5 * (norm(pre_calc_alloc_values.Y_Xm[:, ind])^2) for ind in 1:pre_calc_alloc_values.q]; 
    initialization_values.Σ_sam .= [rand(InverseGamma(pre_calc_alloc_values.astar, pre_calc_alloc_values.bstar[ind]), 1)[1] for ind in 1:pre_calc_alloc_values.q];          # sample Σ
    randn!(pre_calc_alloc_values.u)  # Fills u with standard normal samples
    initialization_values.γ_sam .= (invVγstar.U \ reshape(pre_calc_alloc_values.u, (pre_calc_alloc_values.p + pre_calc_alloc_values.K), pre_calc_alloc_values.q)) * 
                    Diagonal(sqrt.(initialization_values.Σ_sam)) + pre_calc_alloc_values.μγstar;          # sample γ 
end





function run_mcmc(pre_calc_alloc_values::NamedTuple, initialization_values::NamedTuple, results_values::NamedTuple)    

    # Initialize progress bar
    prog = Progress(pre_calc_alloc_values.N_sam, 1, "Running MCMC...", 50)

    # Loop over MCMC iterations
    for l in 1:pre_calc_alloc_values.N_sam
        # Sample F
        sample_f(pre_calc_alloc_values, initialization_values, l)

        # Append F_sam to file
        open(results_values.F_file, "a") do io
            writedlm(io, pre_calc_alloc_values.F_sam, ", ")
        end

        # Sample γ and Σ
        sample_gamma_sigma(pre_calc_alloc_values, initialization_values)

        # Append Σ_sam to file
        open(results_values.sigma_file, "a") do io
            writedlm(io, initialization_values.Σ_sam, ", ")
        end

        # Append γ_sam to file
        open(results_values.gamma_file, "a") do io
            writedlm(io, initialization_values.γ_sam, ", ")
        end

        # Update progress
        next!(prog)
    end
    
end




function post_processing(pre_calc_alloc_values::NamedTuple,results_values::NamedTuple)
    γ_chain = Matrix{Float64}(CSV.read(results_values.gamma_file, DataFrame))
    ind_γ_sam = (1): (pre_calc_alloc_values.p + pre_calc_alloc_values.K) :((pre_calc_alloc_values.p + pre_calc_alloc_values.K) * pre_calc_alloc_values.N_sam);

    Σ_chain = Matrix{Float64}(CSV.read(results_values.sigma_file, DataFrame))
    ind_Σ_sam = (1): pre_calc_alloc_values.q :(pre_calc_alloc_values.q * pre_calc_alloc_values.N_sam);

    F_chain = Matrix{Float64}(CSV.read(results_values.F_file, DataFrame))
    ind_F_sam = (1): pre_calc_alloc_values.n :(pre_calc_alloc_values.n * pre_calc_alloc_values.N_sam);

    N_warmup = floor(Int, 0.25 * pre_calc_alloc_values.N_sam); # ??? i replcaed this. it used to be Integer(trunc) gpt said that could give issues


    # post-processing #
    # step 1 compute the posterior mean
    Λ_mean = Array{Float64, 2}(undef, pre_calc_alloc_values.K, pre_calc_alloc_values.q);
    for k in 1:pre_calc_alloc_values.K
        Λ_mean[k, :] = mean(γ_chain[((k + pre_calc_alloc_values.p)+ N_warmup*(pre_calc_alloc_values.p+pre_calc_alloc_values.K)):(pre_calc_alloc_values.p+pre_calc_alloc_values.K):(k+pre_calc_alloc_values.p+(pre_calc_alloc_values.p+pre_calc_alloc_values.K)*(pre_calc_alloc_values.N_sam-1)), :], dims = 1);
    end
    Λ_mean

    # step 2 label switching
    labels = zeros(Int, pre_calc_alloc_values.N_sam, pre_calc_alloc_values.K);

    # For each iteration and each k
    for i in 1:pre_calc_alloc_values.N_sam
        for k in 1:pre_calc_alloc_values.K
            # Extract current row sample and corresponding mean row
            row_sample = γ_chain[((k + pre_calc_alloc_values.p)+ (i-1)*(pre_calc_alloc_values.p+pre_calc_alloc_values.K)), :];
            row_mean = Λ_mean[k,:];
            
            # Compute inner product and determine sign
            inner_prod = dot(row_sample, row_mean)
            labels[i,k] = sign(inner_prod)
            
            # Switch the sign when the inner product is negative
            if labels[i,k] == -1 
                γ_chain[((k + pre_calc_alloc_values.p)+ (i-1)*(pre_calc_alloc_values.p+pre_calc_alloc_values.K)), :] = -γ_chain[((k + pre_calc_alloc_values.p)+ (i-1)*(pre_calc_alloc_values.p+pre_calc_alloc_values.K)), :];
                F_chain[(1+(i-1)*pre_calc_alloc_values.n):((i)*pre_calc_alloc_values.n), k] = -F_chain[(1+(i-1)*pre_calc_alloc_values.n):((i)*pre_calc_alloc_values.n), k]; # flip sign for iteration i for F
            end
        end
    end
    [countmap(labels[(N_warmup+1):pre_calc_alloc_values.N_sam, col]) for col in 1:pre_calc_alloc_values.K]
    return (Λ_mean=Λ_mean)
end



### function in julia which does pre-calculation and pre-allocation

## i added p q and n beacuse they are specified in the .jld file. unsure if we want to get rid of that ASK
function initialization(p,q,n,coords,X,Y,ϕ_sam,K; μβ=nothing,μΛ= nothing,VΛ=nothing,Vβ=nothing,aΣ=2,bΣ=nothing,m=10, N_sam = 20000)    # Assign defaults inside the function after knowing p, q, K
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
    F_ord = F0[ordered_indices, :]; # not needed here, used in model comparison step?? maybe delete

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
    invVγstar = cholesky(sparse(I, p+K, p+K)); # doesn't fine in-place update for this 
    u = Array{Float64}(undef, (p + K) * q);  # Pre-allocate space for random samples;
    Y_Xm = spzeros(n + p + K, q); # store the residual
    bstar = fill(0.0, q); astar = aΣ + 0.5 * (n);

# specifying what to return which will be used in future steps/functions
    return (Y_ord=Y_ord,X_ord=X_ord,coords_ord=coords_ord,F_ord=F_ord,NN=NN)
end

## note: i didn't include these variables (becasue they weren't called in future code, but if we edit this code/make it more general we may need themF)
    # inv_VΛ, inv_Vr, inv_Lr, inv_LΛ (these may be useful)
    # nIndx (these didn't seem to be useful//only need inside this function)
    # i stopped at nIndx (so continue searching after for what to put in return)

# storing returned values so we can use later. this may need to be updated depending on what we need (i put this in jupyter, but leaving here for this note)
# pre_calc_alloc =  initialization(p,q,n,coords,X,Y,ϕ_sam,K; μβ=fill(0.0, p, q),μΛ= fill(0.0, K, q),VΛ,Vβ,aΣ=2,bΣ=fill(1, q),m=10, N_sam = 20000)

## when i get home: clear the notebook. rerun oNLY my code (not the original precalculation and preallocaiotn), and see what i need to change to get the rest to run
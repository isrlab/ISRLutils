module ISRLutils

# -----------------------------------------------------------------------------------
# This file contains functions which are needed by other packages developed at ISRLab
# -----------------------------------------------------------------------------------

# Convert linear index to multi-dimensional index
# ===============================================
function sub2ind(i,N)
    d = length(N);
    if d == 1
        return i;
    else
        x = i-1;
        ii = zeros(Integer,d);
        for j in d:-1:1 
            A =  prod(N[1:(j-1)]);
            ii[j] = div(x,A) + 1;
            x = x % A;
        end
        return ii
    end
end

# Generate ND grid
# ================
"""
G = GenerateNDGrid(lb,ub,N)
lb, ub, N are d-dimensional vectors specifying the limits and number of points along each dimension.
G is d x prod(N) matrix with grid points.
"""
function GenerateNDGrid(lb,ub,N)
    d = length(N);
    g = [range(lb[i],ub[i],length=N[i]) for i in 1:d];
    nmax = prod(N);
    G = zeros(d,nmax);

    for i in 1:nmax
        ii = sub2ind(i,N);
        for j = 1:d
            G[j,i] = g[j][ii[j]];
        end
    end
    return G
end

# Hypercube corners
# =================
function cornerPoints(x1,x2)
	x1 = x1[:]; x2 = x2[:]; # Vectorize
	d = length(x1);
	if d!=length(x2)
		error("dimension mismatch in x1, x2")
	end

	X = [x1 x2]; # d x 2 matrix
	p = []; # Not the best way to do it.

	ii = zeros(Int,d);
	for i in 0:(2^d-1)
		digits!(ii,i,base=2);
		xx = zeros(Float64,d);
		for j in 1:d
			xx[j] = X[j,ii[j]+1];
		end
		push!(p,xx');
	end
	return vcat(p...);
end

end # module

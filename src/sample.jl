using MathOptInterface
const MOI = MathOptInterface

## Sample solver
# problem is a function that takes in θ and outputs, e.g., iterations (or the quantity of interest)
# N samples will be taken in the hyperbox  lb ≤ θ ≤ ub
function sample_solver(problem,lb,ub,N;sampling=:grid)
  nth = length(ub);

  if(sampling==:MC)
	θs = [lb+ rand(nth).*(ub-lb) for i in 1:N]
  elseif(sampling==:grid)
	n_grid = Int(floor(N^(1/nth)))
	θs = Iterators.product([LinRange(lb[i],ub[i],n_grid) for i in 1:nth]...)
  end
  return reduce(hcat,[vcat(collect(θ),problem(collect(θ))) for θ in θs]) 
end
## Compute Chebyshev centers
# part is a vector of regions 
# (requires fields Ath and bth, which define the polyhedron Ath θ* ≤ bth)
function compute_centers(part)
  N = length(part);
  if(N==0)
	println("The partition is empty");
	return
  end
  nth = size(part[1].Ath,1);
  solver = MOI.OptimizerWithAttributes(GLPK.Optimizer,
                                       MOI.Silent()=>true,
                                       "tol_bnd"=>1e-7)
  lib = DefaultLibrary{Float64}(solver);
  centers = Vector{Float64}[];
  radii = Float64[];
  for (k,p) in enumerate(part)
	print("\rComputing Chebyshev centers: $(round(k/N*100,digits=2))%          ");
	P = polyhedron(hrep(Array(p.Ath'),p.bth),lib);
	d = try dim(P)
	catch
	  0
	end
	if(d==nth)
	  c,r = chebyshevcenter(P,proper=false);
	  push!(centers,c);
	  push!(radii,r);
	else # Lower dimensional Region 
	  push!(centers,Float64[]);
	  push!(radii,0);
	end
  end
  return centers,radii
end


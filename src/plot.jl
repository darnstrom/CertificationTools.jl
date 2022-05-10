using PGFPlotsX, LaTeXStrings 
## Main functions 
function plot_partition(part::Vector;inds=[1,2],slice_vals=Float64[], key=p->getfield(p,:iter),alpha=1.0,axis=[-1.0,1.0,-1.0,1.0],solver=:GLPK)
  if(length(part)==0 || isempty(part[end].Ath)) error("The partition is empty") end
  poly_table,vert_table = transform_part_to_tables(part,inds,slice_vals,key;solver)
  pgfplot_partition(poly_table,vert_table,inds,axis)
end
function plot_samples(samples;inds=[1,2],axis=[-1.0,1.0,-1.0,1.0])
  pgfplot_samples(samples,inds=inds,axis=axis)
end
## Auxiliary function 
function transform_part_to_tables(part::Vector,inds::Vector{Int64},slice_vals::Vector{Float64},key::Function;solver=:GLPK,verbose=true)
  nth = size(part[1].Ath,1);
  if(isempty(slice_vals))
	slice_vals = zeros(nth-2)
  end
  if(solver==:GLPK)
	lib = DefaultLibrary{Float64}(()->GLPK.Optimizer());
  end
  
  slice_inds= setdiff(1:nth,inds);
  part_vrep = NamedTuple[];
  n_max_vert = 0;
  N = length(part)
  for (k,p) in enumerate(part)
	verbose && (mod(k,10)==0) && print("\rForming V-reps: $(round(100*k/N,digits=0))%        ");
	b = p.bth-p.Ath[slice_inds,:]'*slice_vals;
	A = Array(p.Ath[inds,:]');
	P = polyhedron(hrep(A,b),lib);
	if(!isempty(P))
	  try # Sometimes planar_contour gives empty... => Disregard
		xs,ys = Polyhedra.planar_contour(P) # Get vertices 
		push!(part_vrep,(xs=xs,ys=ys,z=key(p)))
		n_max_vert = max(n_max_vert,length(xs))
	  catch
	  end
	end
  end
  verbose && print("\rForming V-reps: 100%         \n");
  id_offset = 0;
  poly_table= zeros(Int64,0,n_max_vert+1); # ids for vertices in each polyhedron (add color id...)
  vert_table = zeros(0,3); # x-coordinate, y-coordinate, z-coordinate of vertex 

  for p in part_vrep 
	n_vert = length(p.xs)
	poly_row = [id_offset*ones(Int64,n_max_vert-n_vert); #Padding for pgfplot 
				collect((id_offset):(id_offset+n_vert-1));
				p.z]
	poly_table = [poly_table; poly_row'];
	vert_table = [vert_table;p.xs p.ys zeros(n_vert)]
	id_offset+=n_vert
  end


  return poly_table,vert_table 
end
## PGF partition
function pgfplot_partition(poly_table::Matrix{Int64},vert_table::Matrix{Float64},inds,
	axis = [-1.0,1.0,-1.0,1.0];facet_color="black")
  max_iter = maximum(poly_table[:,end]);
  min_iter = minimum(poly_table[:,end]);
  cbar_ticks = collect(min_iter:Int64(ceil((max_iter-min_iter)/6)):max_iter)
  cbar_ticks .+= Int64(floor((max_iter-cbar_ticks[end])/2));
  push!(PGFPlotsX.CUSTOM_PREAMBLE,raw"\usepgfplotslibrary{patchplots}")
  @pgf Axis({view = (0, 90),
			 xlabel = latexstring("\\theta_"*string(inds[1])),
			 ylabel = latexstring("\\theta_"*string(inds[2])),
			 xmin=axis[1],xmax=axis[2],
			 ymin=axis[3],ymax=axis[4],
			 point_meta_min=min_iter-0.5,
			 point_meta_max=max_iter+0.5,
			 xlabel_style={yshift="15pt"},
			 ylabel_style={yshift="-20pt"},
			 xtick=axis[1:2],
			 ytick=axis[2:4],
			 colorbar,
			 colorbar_horizontal,
			 colorbar_sampled,
			 colorbar_style =
			 {
			  samples=max_iter-min_iter+2,
			  xlabel="\\# iterations",
			  xticklabel_style={yshift="13pt"},
			  xtick_style={draw="none"},
			  xtick=cbar_ticks
			  },
			 },
			Plot3(
				  {
				   patch,
				   line_width="0.5pt",
				   faceted_color=facet_color,
				   "patch type" = "polygon",
				   "vertex count" = size(poly_table,2)-1,
				   "table/row sep" = "\\\\",
				   patch_table_with_point_meta = TableData(poly_table)
				   },
				  Table(:x => vert_table[:,1],
						:y => vert_table[:,2],
						:z => vert_table[:,3])))
end

## PGF samples
function pgfplot_samples(Samples;inds=[1,2],axis=[-1.0,1.0,-1.0,1.0])
  max_iter = maximum(Samples[end,:]);
  min_iter = minimum(Samples[end,:]);
  cbar_ticks = collect(min_iter:Int64(ceil((max_iter-min_iter)/6)):max_iter)
  cbar_ticks .+= Int64(floor((max_iter-cbar_ticks[end])/2));
@pgf Axis(
		  {
		   xmin=axis[1],xmax=axis[2],
		   ymin=axis[3],ymax=axis[4],
		   point_meta_min=min_iter-0.5,
		   point_meta_max=max_iter+0.5,
		   xlabel_style={yshift="15pt"},
		   ylabel_style={yshift="-20pt"},
		   xtick=axis[1:2],
		   ytick=axis[2:4],
		   colorbar,
		   colorbar_horizontal,
		   colorbar_sampled,
		   colorbar_style =
		   {
			samples=max_iter-min_iter+2,
			xlabel="\\# iterations",
			xticklabel_style={yshift="13pt"},
			xtick_style={draw="none"},
			xtick=cbar_ticks
			},
		   },
		  Plot(
			   {
				scatter,
				only_marks,
				point_meta=raw"\thisrow{key}",
				mark_size=1
				},
			   Table(
					x = Samples[1,:], 
					y = Samples[2,:],
					key= Samples[3,:],
					)
			   )
		  )
end

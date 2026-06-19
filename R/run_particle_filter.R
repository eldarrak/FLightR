# run_particle.filter.R
# functions used during the main run

#' Run Particle Filter
#' 
#' Main function of FLightR, it takes fully prepared object created by \code{\link{make.prerun.object}} and produces a result object that can be used for plotting etc.
#' @param all.out An object created by \code{\link{make.prerun.object}}.
#' @param threads An amount of threads to use for the particle-filter step.
#'   Default is -1. Value 1 forces sequential evaluation. Note that
#'   parallelism in \code{\link{make.prerun.object}} can be useful, but recent
#'   validation benchmarks found the particle-filter PSOCK path slower than
#'   \code{threads = 1} for the tested data/backend; benchmark your own data
#'   before using \code{threads > 1} here.
#' @param cpus another way to specify  \code{threads}
#' @param nParticles  total amount of particles to be used with the run. 10 000 (1e4) is recommended for the preliminary run and 1 000 000 (1e6) for the final
#' @param known.last Set to FALSE if your bird was not at a known place during last twilight in the data
#' @param precision.sd if \code{known.last} then what is the precision of this information. Will be used to resample particles proportionally to their distance from the known last point with probability \code{P = dnorm(0, precision.sd)}
#' @param behav.mask.low.value Probability value that will be used instead of 0 in the behavioural mask. If set to 1 behavioural mask will not be active anymore
#' @param k Kappa parameter from vonMises distribution. Default is NA, otherwise will generate particles in a direction of a previous transitions with kappa = k
#' @param plot Should function plot preliminary map in the end of the run?
#' @param cluster.type see help to package parallel for details
#' @param a minimum distance that is used in the movement model - left boundary for truncated normal distribution of distances moved between twilights. Default is 45 for as default grid has a minimum distance of 50 km.
#' @param b Maximum distance allowed to fly between two consecutive twilights
#' @param L how many consecutive particles to resample
#' @param adaptive.resampling Above what level of ESS resampling should be skipped
#' @param check.outliers switches ON the online outlier routine 
#' @param sink2file will write run details in a file instead of showing on the screen
#' @param add.jitter will add spatial jitter inside a grid cell for the median estimates
#' @param profile.phases logical; if \code{TRUE}, stores detailed particle-filter
#'   phase timings in the returned object. Intended for diagnostics.
#' @param profile.top.level logical; if \code{TRUE}, stores top-level timing for
#'   particle-filter post-processing phases in the returned object.
#' @param propagation.backend Movement proposal backend. \code{"auto"} currently
#'   uses the fast cached backend when possible and falls back to legacy on
#'   construction failure. \code{"partial_cached"} lazily caches reachable
#'   destinations only for visited grid cells and is recommended for larger
#'   spatial extents; local validation benchmarks observed about 20-25x faster
#'   particle-filter runtime for a full-length 1e6-particle run on one tested
#'   dataset and machine. Actual speedups depend on grid size, backend, particle
#'   count, hardware, and thread settings. \code{"cached"} eagerly builds
#'   candidate lists for all grid cells and can be fastest on moderate grids.
#'   \code{"legacy"} keeps the original per-proposal distance calculation for
#'   fallback/reproducibility.
#' @param verbose Controls particle-filter console output. The default
#'   \code{FALSE} is quiet and suppresses normal progress messages while keeping
#'   warnings and errors visible. Use \code{TRUE} for normal progress messages,
#'   or \code{"debug"} for detailed loop diagnostics.
#' @return FLightR object, containing output and extracted results. It is a list with the following elements 
#' 
#'    \item{Indices}{List with prior information and indices}
#'    \item{Spatial}{Spatial data - Grid, Mask, spatial likelihood}
#'    \item{Calibration}{all calibration parameters}
#'    \item{Data}{original data}
#'    \item{Results}{The main results object. Main components of it are
#'       \describe{
#'       \item{Quantiles}{dataframe containing results on locations. Each line corresponds to a twilight}
#'       \item{Movement.results}{dataframe containing all the movement results, Note - time at line n means time of the end of transition between n and n-1}
#'       \item{outliers}{id of twilights excluded by online outlier detection tool}
#'       \item{LL}{-Log likelihood}
#'       \item{Points.rle}{run length encoding object with posterior distribution for every twilight. Note that numbers of points correspond to line numbers in \code{$Spatial$Grid}}
#'       \item{Transitions.rle}{run length encoding object with all the transitions}
#'        }
#'   }
#'
#' @examples
#' File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
#' # to run example fast we will cut the real data file by 2013 Aug 20
#' Proc.data<-get.tags.data(File, end.date=as.POSIXct('2013-07-02', tz='GMT'))
#' Calibration.periods<-data.frame(
#'        calibration.start=NA,
#'        calibration.stop=as.POSIXct("2013-08-20", tz='GMT'),
#'        lon=5.43, lat=52.93) 
#'        #use c() also for the geographic coordinates, if you have more than one calibration location
#'        # (e. g.,  lon=c(5.43, 6.00), lat=c(52.93,52.94))
#' print(Calibration.periods)
#' 
#' # NB Below likelihood.correction is set to FALSE for fast run! 
#' # Leave it as default TRUE for real examples
#' Calibration<-make.calibration(Proc.data, Calibration.periods, likelihood.correction=FALSE)
#' 
#' Grid<-make.grid(left=0, bottom=50, right=10, top=56,
#'   distance.from.land.allowed.to.use=c(-Inf, Inf),
#'   distance.from.land.allowed.to.stay=c(-Inf, Inf))
#'
#' all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93),
#'                              Calibration=Calibration, threads=2)
#' # here we will run only 1e4 partilces for a very short track.
#' # One should use 1e6 particles for the full run.
#' Result<-run.particle.filter(all.in, threads=1,
#'            nParticles=1e3, known.last=TRUE,
#'            precision.sd=25, check.outliers=FALSE)
#'
#' @author Eldar Rakhimberdiev
#' @export
run.particle.filter<-function(all.out, cpus=NULL, threads=-1, nParticles=1e6, known.last=TRUE, precision.sd=25, behav.mask.low.value=0.00, k=NA, plot=TRUE, cluster.type="PSOCK", a=45, b=1500, L=90, adaptive.resampling=0.99, check.outliers=FALSE, sink2file=FALSE, add.jitter=FALSE, profile.phases=FALSE, profile.top.level=FALSE, propagation.backend=c("auto", "cached", "legacy", "partial_cached"), verbose=FALSE) {
   run_pf_start<-proc.time()[["elapsed"]]
   cl<-match.call()
   propagation.backend<-match.arg(propagation.backend)
   if (!is.null(cpus)) {
      warning("use threads instead of cpus! cpus will be supressed in the newer versions\n")
      threads<-cpus
   }
   
   	all.out$Results<-list()
    all.out$Results$SD<-vector(mode = "double")
    all.out$Results$LL<-vector(mode = "double")

   Threads<-1
   parallel=FALSE
   if (threads!=1) warning("Current PSOCK particle-filter execution with threads != 1 may be much slower than sequential execution; recent validation benchmarks found run.particle.filter(..., threads = 1) faster than threads = 4. Use threads = 1 for run.particle.filter() unless this exact dataset/backend has been benchmarked.", call.=FALSE)
   if (propagation.backend=="partial_cached" && threads!=1) warning("The experimental partial_cached propagation backend is sequential-first; with PSOCK threads > 1 each worker may build an independent cache and performance has not been validated.", call.=FALSE)

   if (threads!=1){
   parallel=TRUE
   Possible.threads=parallel::detectCores()
   if (threads<=0) Threads=max(Possible.threads+threads, 1)
   if (threads>0) Threads=min(Possible.threads,threads)
    pf_message("creating cluster with", Threads, "threads", verbose=verbose)
    hosts <- rep("localhost",Threads)
    mycl <- parallel::makeCluster(hosts, type=cluster.type)
    #mycl <- parallel::makeCluster(Threads, type=cluster.type)
    parallel::clusterSetRNGStream(mycl)
	parallel::clusterEvalQ(mycl, library("FLightR")) 
	pf_message('   Done\n', verbose=verbose)

    pf_core_start<-proc.time()[["elapsed"]]
    tryCatch(Res<- pf.run.parallel.SO.resample(in.Data=all.out, threads=Threads, nParticles=nParticles, known.last=known.last, precision.sd=precision.sd, behav.mask.low.value=behav.mask.low.value, k=k, parallel=parallel, plot=FALSE, existing.cluster=mycl, cluster.type=cluster.type, a=a, b=b, L=L, sink2file=sink2file, adaptive.resampling=adaptive.resampling, RStudio=FALSE, check.outliers=check.outliers, profile.phases=profile.phases, propagation.backend=propagation.backend, verbose=verbose), finally = parallel::stopCluster(mycl))
    pf_core_seconds<-proc.time()[["elapsed"]]-pf_core_start
	} else {	
	mycl=NA
    pf_core_start<-proc.time()[["elapsed"]]
    Res<- pf.run.parallel.SO.resample(in.Data=all.out, threads=Threads, nParticles=nParticles, known.last=known.last, precision.sd=precision.sd, behav.mask.low.value=behav.mask.low.value, k=k, parallel=parallel, plot=FALSE, existing.cluster=mycl, cluster.type=cluster.type, a=a, b=b, L=L, sink2file=sink2file, adaptive.resampling=adaptive.resampling, RStudio=FALSE, check.outliers=check.outliers, profile.phases=profile.phases, propagation.backend=propagation.backend, verbose=verbose)
    pf_core_seconds<-proc.time()[["elapsed"]]-pf_core_start
	}
    # Part 2. Creating matrix of results.
	all.out$Results<-list()
	all.out$Results$outliers <- Res$Results$outliers
	all.out$Results$tmp.results<-Res$Results$tmp.results
    # Part 2a. Estimating log likelihood
    get_ll_start<-proc.time()[["elapsed"]]
    LL<- get.LL.PF(all.out, Res$Points)
    get_ll_seconds<-proc.time()[["elapsed"]]-get_ll_start
    pf_message("+----------------------------------+\n", verbose=verbose)
    pf_message("|     estimated negative Log Likelihood is",  LL, "\n", verbose=verbose)
    pf_message("+----------------------------------+\n", verbose=verbose)
    #save(LL, file=paste(LL, "time", format(Sys.time(), "%H-%m"), ".RData"))
	# Part 2b comparing the likelihood with previous estimate
    # now we need to save it..
      

	  # Part 3. Updating proposal
      pf_message("estimating results object\n", verbose=verbose)
      all.out.old<-all.out
      coordinates_start<-proc.time()[["elapsed"]]
      all.out<- get.coordinates.PF(Res$Points, all.out, add.jitter=add.jitter, verbose=verbose)
      coordinates_seconds<-proc.time()[["elapsed"]]-coordinates_start
      movement_start<-proc.time()[["elapsed"]]
      Movement.parameters<- estimate.movement.parameters(Res$Trans, all.out, fixed.parameters=NA, a=a, b=b, parallel=parallel, existing.cluster=mycl, nParticles=nParticles, verbose=verbose)
      movement_seconds<-proc.time()[["elapsed"]]-movement_start
	  
	all.out$Results$Movement.results=Movement.parameters$Movement.results
	all.out$Results$Transitions.rle=Movement.parameters$Transitions.rle	
	
	  all.out$Results$LL<-LL
	  
	phase_profile<-if (isTRUE(profile.phases)) Res$Results$phase_profile else NULL
	partial_cache_diagnostics<-Res$Results$partial_cache_diagnostics
	final_assembly_start<-proc.time()[["elapsed"]]
	all.out$Results<-list(

        #Final.Means=cbind(all.out$Results$Final.Means[-1,],
		#time=all.out$Indices$Matrix.Index.Table$time),
		Quantiles=all.out$Results$Quantiles,
		Movement.results=all.out$Results$Movement.results,
		outliers=all.out$Results$outliers,
		LL=all.out$Results$LL,
		SD=all.out$Results$SD,
		Points.rle=all.out$Results$Points.rle[-1],
		Transitions.rle=all.out$Results$Transitions.rle,
		tmp.results=all.out$Results$tmp.results,
		partial_cache_diagnostics=partial_cache_diagnostics)
	if (isTRUE(profile.phases)) {
	  all.out$Results$phase_profile<-phase_profile
	}
	# [-1] was added in ver 0.3.5, and it makes schedules correct.
	# still there are problems with transitions - the first transition we have is from the starting point and it has no Dawn or Dusk! So one should remove that one if he/she wants to pair these transitions with the time
    rm(Res)
    #rm(All.results.mat)
    # plotting resuls
    if (plot) {
	   lazy.result.plot(all.out)
    }
    gc()
  
  all.out$Spatial$tmp<-NULL
  all.out$call<-cl
  if (is.null(all.out$Metadata)) all.out$Metadata<-list()
  all.out$Metadata$ParticleFilter<-list(
    nParticles=nParticles,
    known.last=known.last,
    precision.sd=precision.sd,
    threads=Threads,
    requested.threads=threads,
    propagation.backend=propagation.backend,
    a=a,
    b=b,
    L=L,
    adaptive.resampling=adaptive.resampling,
    check.outliers=check.outliers,
    add.jitter=add.jitter
  )
  all.out$Results$FLightRver<-utils::packageVersion("FLightR")
  if (isTRUE(profile.phases) || isTRUE(profile.top.level)) {
    final_assembly_seconds<-proc.time()[["elapsed"]]-final_assembly_start
    total_run_seconds<-proc.time()[["elapsed"]]-run_pf_start
    top_level_profile<-data.frame(
      phase=c("pf.run.parallel.SO.resample", "get.LL.PF", "get.coordinates.PF", "estimate.movement.parameters", "final_result_assembly_cleanup_plotting", "total_run.particle.filter"),
      elapsed_seconds=c(pf_core_seconds, get_ll_seconds, coordinates_seconds, movement_seconds, final_assembly_seconds, total_run_seconds),
      nParticles=nParticles,
      threads=Threads,
      backend=propagation.backend,
      n_twilights=nrow(all.out$Indices$Main.Index),
      LL=LL,
      run_label=Sys.getenv("FLIGHTR_RUN_LABEL", unset=NA_character_),
      stringsAsFactors=FALSE
    )
    all.out$Results$top_level_profile<-top_level_profile
  }
  pf_message("DONE!\n", verbose=verbose)
  return(all.out)
}

pf_message<-function(..., verbose=FALSE) {
  if (isTRUE(verbose) || identical(verbose, "debug")) message(...)
}

pf_debug<-function(..., verbose=FALSE) {
  if (identical(verbose, "debug")) message(...)
}

build.grid.movement.candidates<-function(Grid, a=45, b=500) {
  Grid.lonlat<-as.matrix(Grid[, c(1, 2), drop=FALSE])
  Grid.df<-as.data.frame(Grid.lonlat)
  names(Grid.df)<-c("lon", "lat")
  Grid.sf<-sf::st_as_sf(Grid.df, coords=c("lon", "lat"), crs=4326)
  n.grid<-nrow(Grid.lonlat)
  to<-vector("list", n.grid)
  distance<-vector("list", n.grid)
  bearing<-vector("list", n.grid)
  for (from in seq_len(n.grid)) {
    dists<-as.numeric(sf::st_distance(Grid.sf[from,], Grid.sf))/1000
    candidates<-which(is.finite(dists) & dists>=a & dists<=b)
    if (length(candidates)==0) {
      candidates<-from
    }
    to[[from]]<-as.integer(candidates)
    distance[[from]]<-as.numeric(dists[candidates])
    bearing[[from]]<-as.numeric(geosphere::bearing(Grid.lonlat[from, , drop=FALSE], Grid.lonlat[candidates, , drop=FALSE]))
  }
  list(to=to, distance=distance, bearing=bearing, a=a, b=b, grid_n=n.grid)
}


wrap.longitudes<-function(lon) {
  ((lon+180) %% 360)-180
}

candidate.bbox.indices<-function(grid_lonlat, point_lonlat, b) {
  lon<-wrap.longitudes(grid_lonlat[,1])
  lat<-grid_lonlat[,2]
  lon0<-wrap.longitudes(point_lonlat[1])
  lat0<-point_lonlat[2]
  lat_delta<-b/100
  lat_min<-max(-90, lat0-lat_delta)
  lat_max<-min(90, lat0+lat_delta)
  lat_ok<-lat>=lat_min & lat<=lat_max
  cos_lat<-cos(lat0*pi/180)
  if (!is.finite(cos_lat) || abs(cos_lat)<1e-6) {
    lon_ok<-rep(TRUE, length(lon))
  } else {
    lon_delta<-min(180, b/(100*abs(cos_lat)))
    if (!is.finite(lon_delta) || lon_delta>=180) {
      lon_ok<-rep(TRUE, length(lon))
    } else {
      lon_min<-lon0-lon_delta
      lon_max<-lon0+lon_delta
      if (lon_min < -180) {
        lon_ok<-lon>=wrap.longitudes(lon_min) | lon<=lon_max
      } else if (lon_max > 180) {
        lon_ok<-lon>=lon_min | lon<=wrap.longitudes(lon_max)
      } else {
        lon_ok<-lon>=lon_min & lon<=lon_max
      }
    }
  }
  which(lat_ok & lon_ok)
}

build.partial.movement.cache<-function(Grid, a=45, b=500, max_cache_mb=512) {
  cache<-new.env(parent=emptyenv())
  cache$grid_lonlat<-as.matrix(Grid[, c(1, 2), drop=FALSE])
  cache$grid_n<-nrow(cache$grid_lonlat)
  cache$a<-a
  cache$b<-b
  cache$entries<-new.env(parent=emptyenv())
  cache$order<-integer()
  cache$hits<-0L
  cache$misses<-0L
  cache$builds<-0L
  cache$evictions<-0L
  cache$visited_origins<-integer()
  cache$total_candidate_entries<-0L
  cache$approx_cache_bytes<-0
  cache$max_cache_bytes<-max_cache_mb*1024^2
  cache
}

build.partial.movement.entry<-function(cache, from) {
  from<-as.integer(from)
  point<-cache$grid_lonlat[from,]
  idx<-candidate.bbox.indices(cache$grid_lonlat, point, cache$b)
  if (length(idx)>0) {
    src.df<-data.frame(lon=point[1], lat=point[2])
    dst.df<-data.frame(lon=cache$grid_lonlat[idx,1], lat=cache$grid_lonlat[idx,2])
    dists<-as.numeric(sf::st_distance(
      sf::st_as_sf(src.df, coords=c("lon", "lat"), crs=4326),
      sf::st_as_sf(dst.df, coords=c("lon", "lat"), crs=4326)
    ))/1000
    keep<-which(is.finite(dists) & dists>=cache$a & dists<=cache$b)
    to<-as.integer(idx[keep])
    dists<-as.numeric(dists[keep])
  } else {
    to<-integer()
    dists<-numeric()
  }
  if (length(to)==0) {
    to<-from
    dists<-0
  }
  bearings<-as.numeric(geosphere::bearing(cache$grid_lonlat[from, , drop=FALSE], cache$grid_lonlat[to, , drop=FALSE]))
  bearings[!is.finite(bearings)]<-0
  entry<-list(from=from, to=as.integer(to), distance=as.numeric(dists), bearing=as.numeric(bearings), n_candidates=length(to))
  entry$object_size_bytes<-as.numeric(utils::object.size(entry))
  entry
}

get.partial.movement.entry<-function(cache, from) {
  from<-as.integer(from)
  key<-as.character(from)
  cache$visited_origins<-unique(c(cache$visited_origins, from))
  if (exists(key, envir=cache$entries, inherits=FALSE)) {
    cache$hits<-cache$hits+1L
    return(get(key, envir=cache$entries, inherits=FALSE))
  }
  cache$misses<-cache$misses+1L
  entry<-build.partial.movement.entry(cache, from)
  entry_size<-entry$object_size_bytes
  if ((cache$approx_cache_bytes+entry_size) <= cache$max_cache_bytes) {
    assign(key, entry, envir=cache$entries)
    cache$order<-c(cache$order, from)
    cache$builds<-cache$builds+1L
    cache$total_candidate_entries<-cache$total_candidate_entries+entry$n_candidates
    cache$approx_cache_bytes<-cache$approx_cache_bytes+entry_size
  }
  entry
}

partial.movement.cache.diagnostics<-function(cache) {
  if (is.null(cache)) return(NULL)
  keys<-ls(cache$entries, all.names=TRUE)
  counts<-if (length(keys)>0) vapply(keys, function(key) get(key, envir=cache$entries, inherits=FALSE)$n_candidates, integer(1)) else integer()
  fractions<-if (length(counts)>0 && cache$grid_n>0) counts/cache$grid_n else numeric()
  data.frame(
    propagation_backend="partial_cached",
    grid_size=cache$grid_n,
    cached_origins=length(keys),
    visited_origins=length(unique(cache$visited_origins)),
    cache_hits=cache$hits,
    cache_misses=cache$misses,
    cache_builds=cache$builds,
    cache_evictions=cache$evictions,
    total_candidate_entries=sum(counts),
    median_candidate_count=if (length(counts)>0) stats::median(counts) else NA_real_,
    max_candidate_count=if (length(counts)>0) max(counts) else NA_integer_,
    median_candidate_fraction_of_grid=if (length(fractions)>0) stats::median(fractions) else NA_real_,
    max_candidate_fraction_of_grid=if (length(fractions)>0) max(fractions) else NA_real_,
    approximate_cache_object_bytes=cache$approx_cache_bytes,
    stringsAsFactors=FALSE
  )
}

generate.points.dirs.legacy<-function(x , in.Data, Current.Proposal, a=45, b=500) {
  # this function is needed to generate new points - it works as from input point Index and biological proposal
  ################
  # x has 3 columns
  # Index
  # number
  # nMoving
 # here is the important change for the mask started in vesion 3.0
 # the new idea is that we wnat to set probabilities of stop to 0 in case the bird is flying over the restricted habitat..
  
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  
  #if (!is.null(in.Data$Spatial$tmp$dDistance))in.Data$Spatial$tmp$Distance<-attach.big.matrix(in.Data$Spatial$tmp$dDistance)
  #if (!is.null(in.Data$Spatial$tmp$dAzimuths))in.Data$Spatial$tmp$Azimuths<-attach.big.matrix(in.Data$Spatial$tmp$dAzimuths)
  
  if (x[[3]]>0) {
	# here is the addition of clever mask
	#if (in.Data$Spatial$Grid[x[[1]],3]==0) x[[3]]=x[[2]]
	# end of addition
	#Dists.distr<-in.Data$Spatial$tmp$Distance[x[[1]],]	
	#Dists.distr<- sp::spDists(in.Data$Spatial$Grid[x[[1]], c(1,2), drop=FALSE] ,  in.Data$Spatial$Grid[,c(1,2)], longlat=TRUE)
	Dists.distr <- sf::st_distance(
	  sf::st_as_sf(as.data.frame(in.Data$Spatial$Grid[x[[1]], c(1,2), drop=FALSE]),coords=c('lon','lat'),crs=4326),
	  sf::st_as_sf(as.data.frame(in.Data$Spatial$Grid[, c(1,2)]),coords=c('lon','lat'),crs=4326))/1000
	
    Dists.probs<-truncnorm::dtruncnorm(as.numeric(Dists.distr), a=a, b=b, Current.Proposal$M.mean, Current.Proposal$M.sd)
    ###
    #fields::image.plot(fields::as.image(Dists.probs, x=in.Data$Spatial$Grid, nrow=50, ncol=50))
    #
    if (Current.Proposal$Kappa>0) {
      #Angles.dir<-in.Data$Spatial$tmp$Azimuths[x[[1]],]
      Angles.dir<-geosphere::bearing(in.Data$Spatial$Grid[x[[1]], c(1,2), drop=FALSE], in.Data$Spatial$Grid[,c(1,2)])
	  
      Angles.probs<-as.numeric(suppressWarnings(circular::dvonmises(Angles.dir/180*pi, mu=Current.Proposal$Direction/180*pi, kappa=Current.Proposal$Kappa)))
      Angles.probs[is.na(Angles.probs)]<-0
    }
    else {
      Angles.probs<-0.1591549
    }
    # this is neede to catch an error 
    if (any(c(Angles.probs, Dists.probs)<0)) {
      warning("PF produced weird probs \n")
      tmp.out<-list(Angles.probs=Angles.probs, Dists.probs=Dists.probs, Current.Proposal=Current.Proposal, Data=x)
      tmpfile<-paste0(tempfile(), '.RData')
      warning(paste('PF produced weird probs find them saved in', tmpfile))
      save(tmp.out, file=tmpfile)
      Biol.proposal<-pmax.int(Dists.probs*Angles.probs, 0)
    } else {
      Biol.proposal<-Dists.probs*Angles.probs
    }
    pos.biol<-suppressWarnings(sample.int(length(Biol.proposal), size=x[[3]], replace=TRUE, prob=Biol.proposal))
    # ok, now we want to return more complicated stuff - final indices!
    return(resample(as.integer(c(pos.biol, rep(x[[1]],(x[[2]]-x[[3]]))))))
  } else {
    return(as.integer(rep(x[[1]],(x[[2]]))))
  }
}

generate.points.dirs.cached<-function(x, in.Data, Current.Proposal, a=45, b=500) {
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  from<-as.integer(x[[1]])
  n.total<-as.integer(x[[2]])
  n.moving<-as.integer(x[[3]])
  if (n.moving<=0) return(as.integer(rep(from, n.total)))

  candidates<-in.Data$Spatial$tmp$movement_candidates
  if (is.null(candidates) || candidates$grid_n!=nrow(in.Data$Spatial$Grid) ||
      candidates$a!=a || candidates$b!=b) {
    return(generate.points.dirs.legacy(x, in.Data, Current.Proposal, a=a, b=b))
  }

  to<-candidates$to[[from]]
  dists<-candidates$distance[[from]]
  bearings<-candidates$bearing[[from]]
  if (length(to)==0 || length(dists)!=length(to)) {
    return(generate.points.dirs.legacy(x, in.Data, Current.Proposal, a=a, b=b))
  }

  sd<-Current.Proposal$M.sd
  if (!is.finite(sd) || sd<=0) return(generate.points.dirs.legacy(x, in.Data, Current.Proposal, a=a, b=b))

  Biol.proposal<-exp(-0.5*((dists-Current.Proposal$M.mean)/sd)^2)
  if (is.finite(Current.Proposal$Kappa) && Current.Proposal$Kappa>0) {
    angle.diff<-(bearings-Current.Proposal$Direction)*pi/180
    Biol.proposal<-Biol.proposal*exp(Current.Proposal$Kappa*cos(angle.diff))
  }
  Biol.proposal[!is.finite(Biol.proposal)]<-0
  if (all(Biol.proposal<=0)) return(generate.points.dirs.legacy(x, in.Data, Current.Proposal, a=a, b=b))

  Full.proposal<-numeric(nrow(in.Data$Spatial$Grid))
  Full.proposal[to]<-Biol.proposal
  pos.biol<-suppressWarnings(sample.int(length(Full.proposal), size=n.moving, replace=TRUE, prob=Full.proposal))
  return(resample(as.integer(c(pos.biol, rep(from, n.total-n.moving)))))
}


generate.points.dirs.partial<-function(x, in.Data, Current.Proposal, a=45, b=500) {
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  from<-as.integer(x[[1]])
  n.total<-as.integer(x[[2]])
  n.moving<-as.integer(x[[3]])
  if (n.moving<=0) return(as.integer(rep(from, n.total)))

  cache<-in.Data$Spatial$tmp$partial_movement_cache
  if (is.null(cache) || cache$grid_n!=nrow(in.Data$Spatial$Grid) || cache$a!=a || cache$b!=b) {
    return(generate.points.dirs.legacy(x, in.Data, Current.Proposal, a=a, b=b))
  }
  entry<-get.partial.movement.entry(cache, from)
  to<-entry$to
  dists<-entry$distance
  bearings<-entry$bearing
  if (length(to)==0 || length(dists)!=length(to)) {
    return(generate.points.dirs.legacy(x, in.Data, Current.Proposal, a=a, b=b))
  }

  sd<-Current.Proposal$M.sd
  if (!is.finite(sd) || sd<=0) return(generate.points.dirs.legacy(x, in.Data, Current.Proposal, a=a, b=b))

  Biol.proposal<-exp(-0.5*((dists-Current.Proposal$M.mean)/sd)^2)
  if (is.finite(Current.Proposal$Kappa) && Current.Proposal$Kappa>0) {
    angle.diff<-(bearings-Current.Proposal$Direction)*pi/180
    Biol.proposal<-Biol.proposal*exp(Current.Proposal$Kappa*cos(angle.diff))
  }
  Biol.proposal[!is.finite(Biol.proposal)]<-0
  if (all(Biol.proposal<=0)) return(generate.points.dirs.legacy(x, in.Data, Current.Proposal, a=a, b=b))

  Full.proposal<-numeric(nrow(in.Data$Spatial$Grid))
  Full.proposal[to]<-Biol.proposal
  pos.biol<-suppressWarnings(sample.int(length(Full.proposal), size=n.moving, replace=TRUE, prob=Full.proposal))
  return(resample(as.integer(c(pos.biol, rep(from, n.total-n.moving)))))
}

generate.points.dirs<-function(x, in.Data, Current.Proposal, a=45, b=500) {
  if (!is.null(in.Data$Spatial$tmp$movement_candidates)) {
    generate.points.dirs.cached(x, in.Data, Current.Proposal, a=a, b=b)
  } else if (!is.null(in.Data$Spatial$tmp$partial_movement_cache)) {
    generate.points.dirs.partial(x, in.Data, Current.Proposal, a=a, b=b)
  } else {
    generate.points.dirs.legacy(x, in.Data, Current.Proposal, a=a, b=b)
  }
}

weight_stack_columns<-function(start_col, active_cols, max_columns) {
  ((start_col-1L+seq_len(active_cols)-1L) %% max_columns)+1L
}

weight_stack_matrix<-function(values, start_col, active_cols) {
  values[, weight_stack_columns(start_col, active_cols, ncol(values)), drop=FALSE]
}

weight_stack_last_column<-function(start_col, active_cols, max_columns) {
  ((start_col+active_cols-2L) %% max_columns)+1L
}

weight_stack_last<-function(values, start_col, active_cols) {
  values[, weight_stack_last_column(start_col, active_cols, ncol(values))]
}

pf.run.parallel.SO.resample<-function(in.Data, threads=2, nParticles=1e6, known.last=TRUE, precision.sd=25, behav.mask.low.value=0.01, k=NA, parallel=TRUE, plot=TRUE, existing.cluster=NA, cluster.type="PSOCK", a=45, b=500, sink2file=FALSE, L=25, adaptive.resampling=0.5, RStudio=FALSE, check.outliers=FALSE, profile.phases=FALSE, propagation.backend=c("auto", "cached", "legacy", "partial_cached"), verbose=FALSE) {
  ### to make algorhythm work in a fast mode w/o directional proposal use k=NA
  if (sink2file & !RStudio)  sink(file=paste("pf.run.parallel.SO.resample", format(Sys.time(), "%H-%m"), "txt", sep="."))
  if (sink2file & RStudio) sink()
  #### if save.rle=TRUE function will save a out.rle object in out.rle.RData file in working directory BUT to make it work There have to be out.rle column in Main.Index that will contain true for the twilights where results needed.
  ####
  
 if (any(in.Data$Spatial$Behav.mask!=1)) { 
	smart.filter=TRUE
	pf_message("smart filter is ON\n", verbose=verbose)
	} else {
	smart.filter=FALSE
	pf_message("smart filter is OFF\n", verbose=verbose)
	}
  propagation.backend<-match.arg(propagation.backend)
  if (propagation.backend=="partial_cached") {
    if (parallel) warning("The experimental partial_cached propagation backend is sequential-first; PSOCK workers will not share cache state and may be slower.", call.=FALSE)
    if (is.null(in.Data$Spatial$tmp)) in.Data$Spatial$tmp<-list()
    in.Data$Spatial$tmp$movement_candidates<-NULL
    in.Data$Spatial$tmp$partial_movement_cache<-build.partial.movement.cache(in.Data$Spatial$Grid, a=a, b=b)
    pf_message("partial cached propagation backend is ON\n", verbose=verbose)
  } else if (propagation.backend!="legacy") {
    movement_candidates<-try(build.grid.movement.candidates(in.Data$Spatial$Grid, a=a, b=b), silent=TRUE)
    if (inherits(movement_candidates, "try-error")) {
      if (propagation.backend=="cached") stop("Cached propagation candidate construction failed: ", conditionMessage(attr(movement_candidates, "condition")), call.=FALSE)
      warning("Cached propagation candidate construction failed; using legacy propagation backend.", call.=FALSE)
      propagation.backend<-"legacy"
    } else {
      if (is.null(in.Data$Spatial$tmp)) in.Data$Spatial$tmp<-list()
      in.Data$Spatial$tmp$partial_movement_cache<-NULL
      in.Data$Spatial$tmp$movement_candidates<-movement_candidates
      pf_message("cached propagation backend is ON\n", verbose=verbose)
    }
  } else {
    if (!is.null(in.Data$Spatial$tmp)) {
      in.Data$Spatial$tmp$movement_candidates<-NULL
      in.Data$Spatial$tmp$partial_movement_cache<-NULL
    }
    pf_message("legacy propagation backend is ON\n", verbose=verbose)
  }  # so, the idea here will be that we don't need to create these complicated Indexes..
  in.Data.short<-list(Indices=in.Data$Indices,  Spatial=in.Data$Spatial)
  in.Data.short$Spatial$Behav.mask<-NULL
  in.Data.short$Spatial$Phys.Mat<-NULL
  
  Parameters<-list(in.Data=in.Data.short, a=a, b=b)
  #if (!is.null(Parameters$in.Data$Spatial$tmp$dDistance)) Parameters$in.Data$Spatial$tmp$Distance<-attach.big.matrix(Parameters$in.Data$Spatial$tmp$dDistance)
  #if (!is.null(Parameters$in.Data$Spatial$tmp$dAzimuths)) Parameters$in.Data$Spatial$tmp$Azimuths<-attach.big.matrix(Parameters$in.Data$Spatial$tmp$dAzimuths)
  if (parallel) {
    if (length(existing.cluster)==1) {
      hosts <- rep("localhost",threads)
      mycl <- parallel::makeCluster(hosts, type=cluster.type)    
      #mycl <- parallel::makeCluster(threads, type=cluster.type)
      parallel::clusterSetRNGStream(mycl)
      ### we don' need to send all parameters to node. so keep it easy..
      # cleaning dataset
    }
    else {
      mycl<-existing.cluster
    }
    parallel::clusterExport(mycl, "Parameters", envir=environment())
	#if (!is.null(in.Data$Spatial$tmp$dDistance))  {parallel::clusterEvalQ(mycl, {Parameters$in.Data$Spatial$tmp$Distance <- attach.big.matrix(Parameters$in.Data$Spatial$tmp$dDistance);1})
	#}
	#if (!is.null(in.Data$Spatial$tmp$dAzimuths))  parallel::clusterEvalQ(mycl, { Parameters$in.Data$Spatial$tmp$Azimuths <- attach.big.matrix(Parameters$in.Data$Spatial$tmp$dAzimuths);1})
  }
  
  
  #########################################
  # part from the main function
  
  in.Data$Spatial$Behav.mask[in.Data$Spatial$Behav.mask<behav.mask.low.value]<-behav.mask.low.value
  in.Data$Spatial$Phys.Mat<-in.Data$Spatial$Phys.Mat*in.Data$Spatial$Behav.mask
  #=========================================
  # get value for density for angles
  if (!is.na(k)) {
    D.kappa<-as.numeric(suppressWarnings(circular::dvonmises(0, mu=0, kappa=k)))
  }
  else {
    D.kappa<-as.numeric(suppressWarnings(circular::dvonmises(0, mu=0, kappa=0)))
  }
  # stage 1. create burn-in set
  #Last.Particles<-rep(as.integer(in.Data$Spatial$start.point), nParticles)
  #Last.Last.Particles<-Last.Particles
  
  # here I'll save results stack.
  #All.results
  
  Results.stack<-as.matrix(rep(as.integer(in.Data$Spatial$start.point), nParticles))
  
  Weights.stack.max.cols<-L+2L
  Weights.stack<-matrix(NA_real_, nrow=nParticles, ncol=Weights.stack.max.cols)
  Weights.stack[,1]<-1/nParticles
  Weights.stack.active<-1L
  Weights.stack.start<-1L
  
  New.weights<-rep(1/nParticles, nParticles)
  #All.results<-NULL
  in.Data$outliers<-c()
	  #Trans<-vector(mode = "list", length = nrow(in.Data$Indices$Main.Index))
	  Trans<-vector(mode = "list")
  	if (check.outliers) {
	  in.Data$AB.distance<-c()
	  #in.Data$AC.distance<-c()
	  in.Data$AC.distance2<-c()
	  in.Data$Dif.ang<-c()
	  #in.Data$BC.mean<-c()
  }
  
  phase_names<-c("proposal_lookup", "propagate_particles", "physical_weights",
                 "directional_weights", "smart_filter", "outlier_check",
                 "cumulative_weight_product", "ESS_calculation", "resampling",
                 "accepted_particles_assignment", "results_stack_append",
                 "weights_stack_append", "point_rle_creation",
                 "transition_rle_creation", "stack_drop", "final_smoothing",
                 "final_stack_flush")
  phase_profile<-NULL
  profile_phase<-function(name, expr) {
    start<-proc.time()[["elapsed"]]
    value<-force(expr)
    phase_times[[name]]<<-phase_times[[name]]+(proc.time()[["elapsed"]]-start)
    value
  }
  
  propagate.particles<-function(Last.Particles, Current.Proposal, parallel=TRUE, Parameters, mycl) {
    Order.vector<-order(Last.Particles)
    Particles.rle<-bit::intrle(as.integer(Last.Particles[Order.vector]))
    if (is.null(Particles.rle)) 	Particles.rle<-rle(as.integer(Last.Particles[Order.vector]))
    Last.State<-cbind(Particles.rle$values, Particles.rle$length)
    Last.State<-cbind(Last.State,nMoving=stats::rbinom(dim(Last.State)[[1]], size=Last.State[,2], prob=Current.Proposal$Decision))
    Last.State.List<-split(Last.State, row(Last.State))
    nSeq<-nrow(Last.State)
    if (nSeq==1 | (!parallel)) {
      New.Points<-lapply(Last.State.List, FUN=function(x,Current.Proposal, Parameters) do.call(generate.points.dirs, c(x=list(x), Parameters, Current.Proposal=list(Current.Proposal))), Current.Proposal, Parameters)
      #New.Points<-lapply(Last.State.List, pf.par.internal, Current.Proposal)
    }
    else {
      pf_debug(" initial diversity is ", nSeq, verbose=verbose)
      #New.Points<-clusterApplyLB(mycl, Last.State.List,  pf.par.internal, Current.Proposal)
      New.Points<-parallel::clusterApply(mycl, Last.State.List,  pf.par.internal, Current.Proposal)
	  pf_debug("  Done\n", verbose=verbose)
    }
    #=======================================================
    New.Particles<-unlist(New.Points)[sort.list(Order.vector, na.last = NA, method =  "quick")]
    if (isTRUE(profile.phases)) {
      attr(New.Particles, "n_groups")<-nSeq
      attr(New.Particles, "n_moving_groups")<-sum(Last.State[,"nMoving"]>0)
    }
    return(New.Particles)
  }
  
  
  ## for ver 1.6 
  #Prev.Weights<-rep(1, nParticles) 
  ResampleCount<-0
  steps.from.last<-2
  total_length<-nrow(in.Data$Indices$Main.Index)
  if (isTRUE(profile.phases)) phase_profile<-vector("list", total_length+1L)
  for (Time.Period in 1:total_length) {
    if (isTRUE(profile.phases)) {
      iteration_start<-proc.time()[["elapsed"]]
      phase_times<-stats::setNames(rep(0, length(phase_names)), phase_names)
      n_unique_last_particles<-length(unique(Results.stack[,ncol(Results.stack)]))
      n_groups<-NA_integer_
      n_moving_groups<-NA_integer_
      n_unique_new_particles<-NA_integer_
      n_unique_point_rle_values<-NA_integer_
      n_unique_transition_rle_values<-NA_integer_
      did_resample<-FALSE
      ESS<-NA_real_
    }
    steps.from.last=steps.from.last+1

	pf_debug("\n\n##########################\n     Time.Period", Time.Period, "of", total_length, "\n", verbose=verbose)
    #cat("prep. data:")
    Current.Proposal<-if (isTRUE(profile.phases)) {
      profile_phase("proposal_lookup", in.Data$Indices$Matrix.Index.Table[in.Data$Indices$Main.Index$Biol.Prev[Time.Period],])
    } else {
      in.Data$Indices$Matrix.Index.Table[in.Data$Indices$Main.Index$Biol.Prev[Time.Period],]
    }
    #=======================================
    pf_debug("generating new particles", verbose=verbose)
    New.Particles<-if (isTRUE(profile.phases)) {
      profile_phase("propagate_particles", propagate.particles(Last.Particles=Results.stack[,ncol(Results.stack)], Current.Proposal=Current.Proposal, parallel=parallel, Parameters=Parameters, mycl=mycl))
    } else {
      propagate.particles(Last.Particles=Results.stack[,ncol(Results.stack)], Current.Proposal=Current.Proposal, parallel=parallel, Parameters=Parameters, mycl=mycl)
    }
    if (isTRUE(profile.phases)) {
      n_groups<-attr(New.Particles, "n_groups", exact=TRUE)
      n_moving_groups<-attr(New.Particles, "n_moving_groups", exact=TRUE)
      n_unique_new_particles<-length(unique(as.integer(New.Particles)))
    }
    
    #=====================================================
    # resampling step 
    # we need to estimate weights of resulting points and then resample them proportionally to weights..
    Current.Phys.Weights<-if (isTRUE(profile.phases)) {
      profile_phase("physical_weights", in.Data$Spatial$Phys.Mat[New.Particles,Time.Period])
    } else {
      in.Data$Spatial$Phys.Mat[New.Particles,Time.Period]
    }
    
    if (!is.na(k) & Time.Period >1) {
      get.directional.weights<-function(in.Data, Last.Last.Particles, Last.Particles, New.Particles, k, parallel, mycl, D.kappa) {
        #===================================
        #Prev.Dirs<-in.Data$Spatial$tmp$Azimuths[cbind(Last.Last.Particles, Last.Particles)]/180*pi
		
        Prev.Dirs<-apply(matrix(c(Last.Last.Particles, Last.Particles), ncol=2), 1, dir_fun, in.Data)/180*pi
		
        #New.Dirs<-in.Data$Spatial$tmp$Azimuths[cbind(Last.Particles,New.Particles)]/180*pi
        New.Dirs<-apply(matrix(c(Last.Particles, New.Particles), ncol=2), 1, dir_fun, in.Data)/180*pi
		
        FromTo=matrix(c(Prev.Dirs, New.Dirs), ncol=2)
        if (parallel) {Angles.probs<-parallel::parApply(mycl, FromTo, 1, flightr.dvonmises, mykap=k)
        } else {
          Angles.probs<-apply(FromTo,1, FUN=function(x, k) as.numeric(suppressWarnings(circular::dvonmises(x[[2]], mu=x[[1]], kappa=k))), k=k)
        }
        Angles.probs[is.na(Angles.probs)]<-D.kappa
        return(Angles.probs)
      }
      
      Angles.probs<-if (isTRUE(profile.phases)) {
        profile_phase("directional_weights", get.directional.weights(in.Data, Results.stack[,ncol(Results.stack)-1], Results.stack[,ncol(Results.stack)], New.Particles, k, parallel, mycl, D.kappa))
      } else {
        get.directional.weights(in.Data, Results.stack[,ncol(Results.stack)-1], Results.stack[,ncol(Results.stack)], New.Particles, k, parallel, mycl, D.kappa)
      }
      Current.Weights<- Current.Phys.Weights*Angles.probs
    }
    else {Current.Weights<-Current.Phys.Weights*D.kappa}

	if (isTRUE(profile.phases)) {
	  profile_phase("smart_filter", {
	    if (smart.filter) {
	      Index.Stable<-which(New.Particles == Results.stack[,ncol(Results.stack)])
	      Current.Weights[Index.Stable]<-Current.Weights[Index.Stable]*(in.Data$Spatial$Behav.mask[New.Particles][Index.Stable])
	    }
	    NULL
	  })
	} else if (smart.filter) {
	#=========================
	# here we should introduce clever mask..
	# this will be something like that
	Index.Stable<-which(New.Particles == Results.stack[,ncol(Results.stack)])
	Current.Weights[Index.Stable]<-Current.Weights[Index.Stable]*(in.Data$Spatial$Behav.mask[New.Particles][Index.Stable])
	}
	#=======================================================
	# now I want to add 3.3
	# the idea is following - we should compare distances from the last t-2 to t-1 and to t
	outlier_start<-if (isTRUE(profile.phases)) proc.time()[["elapsed"]] else NA_real_
	if (Time.Period>2 & check.outliers) {
	# here I want to try directiondl outliers..
	# so we need to pick up particles that moved..
	# 
	#=================================================================

	# AB.distance<-stats::weighted.mean(	sp::spDists(in.Data$Spatial$Grid[Results.stack[,(ncol(Results.stack)-1)], c(1,2), drop=FALSE], 
	#                                                in.Data$Spatial$Grid[Results.stack[,ncol(Results.stack)], c(1,2), drop=FALSE], 
	#                                                longlat=TRUE, diagonal=TRUE), 
	#                                    Weights.stack[,ncol(Weights.stack)])

	AB.distance <- stats::weighted.mean(sf::st_distance(
	  sf::st_as_sf(as.data.frame(in.Data$Spatial$Grid[Results.stack[,(ncol(Results.stack)-1)], c(1,2), drop=FALSE]),coords=c('lon','lat'),crs=4326), 
	  sf::st_as_sf(as.data.frame(in.Data$Spatial$Grid[Results.stack[,ncol(Results.stack)], c(1,2), drop=FALSE]),coords=c('lon','lat'),crs=4326), by_element=TRUE)/1000, 
	  weight_stack_last(Weights.stack, Weights.stack.start, Weights.stack.active)) |> as.numeric()
	
	# AC.distance2<-	stats::weighted.mean(sp::spDists(in.Data$Spatial$Grid[Results.stack[,(ncol(Results.stack)-1)], c(1,2), drop=FALSE], 
	#                                                 in.Data$Spatial$Grid[New.Particles, c(1,2), drop=FALSE], 
	#                                                 longlat=TRUE, 
	#                                                 diagonal=TRUE), 
	#                                     Weights.stack[,ncol(Weights.stack)]*Current.Weights)	
	AC.distance2 <- stats::weighted.mean(sf::st_distance(
	  sf::st_as_sf(as.data.frame(in.Data$Spatial$Grid[Results.stack[,(ncol(Results.stack)-1)], c(1,2), drop=FALSE]),coords=c('lon','lat'),crs=4326),
	  sf::st_as_sf(as.data.frame(in.Data$Spatial$Grid[New.Particles, c(1,2), drop=FALSE]),coords=c('lon','lat'),crs=4326),
	  by_element=TRUE
	  )/1000, 
	  weight_stack_last(Weights.stack, Weights.stack.start, Weights.stack.active)*Current.Weights) |> as.numeric()
	
	pf_debug("AB.distance:", round(AB.distance, 2), "\n", verbose=verbose)
	pf_debug("AC.distance2:", round(AC.distance2, 2), "\n", verbose=verbose)

	Dif.ang=180
	if (AB.distance>50) { # go for angles only if disctances are high!
	resample <- function(x, ...) x[sample.int(length(x), ...)]
	# the AB ones wil have folowing..

	BA.dir<-apply(Results.stack[,ncol(Results.stack):(ncol(Results.stack)-1), drop=FALSE], 1, dir_fun, in.Data)
	
	BA.moved<-which(!is.na(BA.dir))
	BA.mean<-circular::mean.circular(circular::circular(resample(BA.dir[BA.moved], replace=TRUE,prob=weight_stack_last(Weights.stack, Weights.stack.start, Weights.stack.active)[BA.moved]), units="degrees"), na.rm=TRUE)
	BC.dir<-apply(matrix(c(Results.stack[,ncol(Results.stack)], New.Particles), ncol=2), 1, dir_fun, in.Data)
	
	BC.moved<-which(!is.na(BC.dir))
	BC.mean<-circular::mean.circular(circular::circular(resample(BC.dir[BC.moved], replace=TRUE, prob=(weight_stack_last(Weights.stack, Weights.stack.start, Weights.stack.active)*Current.Weights)[BC.moved]), units="degrees"), na.rm=TRUE)
	dif.ang<-function(x,y) {
	# this function provides the minimum angle between two angles..
	y=y*pi/180
	x=x*pi/180
	z<-c(y-x, y-x+2*pi, y-x-2*pi)
	z[which.min(abs(z))]*180/pi
	}
	if (!is.null(BA.mean) & !is.null(BC.mean)) {
	Dif.ang<-dif.ang(BA.mean, BC.mean)
	pf_debug("anglular change", round(Dif.ang,2), "\n", verbose=verbose)
	}
	}	

#=================================
# now I want to temporarily save this ..
  in.Data$AB.distance<-c(in.Data$AB.distance, AB.distance)
  #in.Data$AC.distance<-c(in.Data$AC.distance, AC.distance)
  in.Data$AC.distance2<-c(in.Data$AC.distance2, AC.distance2)
  in.Data$Dif.ang<-c(in.Data$Dif.ang, Dif.ang)
  #in.Data$BC.mean<-c(in.Data$BC.mean, BC.mean)
	# outlier detection 10
	# oulier detection
	#if (AB.distance>(AC.distance2*1.5) | (abs(Dif.ang)<90 & AB.distance>50)) {
	
	#if ( steps.from.last>2 & ((AB.distance>50 & AB.distance>(AC.distance2*1.3))| (abs(Dif.ang)<90))) {
	if ( steps.from.last>2 & ((AB.distance>50 & AB.distance>(AC.distance2*1.3))| (AB.distance>AC.distance2 & abs(Dif.ang)<100))) {
	steps.from.last=0
	warning("outlier number", length(in.Data$outliers)+1, "detected! removing observations from ", Time.Period-1, "!\n")
	in.Data$outliers<-c(in.Data$outliers, Time.Period-1)
	# and now the main question - how should we remove the weights??
	# probably the easiest way would be to add 1/nParticles into the last column..
	Weights.stack[, weight_stack_last_column(Weights.stack.start, Weights.stack.active, Weights.stack.max.cols)]<-rep(1/nParticles, nParticles)
	# ok and now I have to change the resample to t-1
	}
	}
	if (isTRUE(profile.phases)) phase_times[["outlier_check"]]<-phase_times[["outlier_check"]]+(proc.time()[["elapsed"]]-outlier_start)
	
    #=====================================================
    # here is the change from 1.6 - Adding cumulative weights now..
	# and this came back at the version 3.2!
    #Current.Weights.with.Prev.mat<-cbind(Weights.stack, Current.Weights)
    
	#rowProds <- function(a) exp(rowMeans(log(a)))
    rowProds <- function(a) exp(rowSums(log(a)))
	
	#Current.Weights.with.Prev<-	pmax(rowProds(Current.Weights.with.Prev.mat), 1e-323)
	
	Current.Weights.with.Prev<-if (isTRUE(profile.phases)) {
	  profile_phase("cumulative_weight_product", {
	if (Time.Period != total_length) {
	   	# here is the place - I'll check weights without use of the last stage information! 
     	pmax(rowProds(weight_stack_matrix(Weights.stack, Weights.stack.start, Weights.stack.active)), 1e-323)
    } else {
	    pmax(rowProds(cbind(weight_stack_matrix(Weights.stack, Weights.stack.start, Weights.stack.active), Current.Weights)) , 1e-323)
	}
	  })
	} else {
	if (Time.Period != total_length) {
	   	# here is the place - I'll check weights without use of the last stage information! 
     	pmax(rowProds(weight_stack_matrix(Weights.stack, Weights.stack.start, Weights.stack.active)), 1e-323)
    } else {
	    pmax(rowProds(cbind(weight_stack_matrix(Weights.stack, Weights.stack.start, Weights.stack.active), Current.Weights)) , 1e-323)
	}
	}
	#
    #================================================================
    # ver 1.7. 
    # ADAPTIVE RESAMPLING
    if (isTRUE(profile.phases)) {
      profile_phase("ESS_calculation", {
    if (adaptive.resampling!=1) {
      ESS<-(sum(Current.Weights.with.Prev)^2)/sum((Current.Weights.with.Prev)^2)
      pf_debug("ESS is ", ESS, verbose=verbose)
if (is.na(ESS)) {
	ESS=1
	Current.Weights.with.Prev.mat<-cbind(weight_stack_matrix(Weights.stack, Weights.stack.start, Weights.stack.active), Current.Weights)
	tmpfile<-paste0(tempfile(), '.RData')
    save(Current.Weights.with.Prev, Current.Weights.with.Prev.mat, Current.Weights, file=tmpfile)
    warning(paste('find the unexpected weights saved in', tmpfile))
	}		
    } else {
      ESS<-1
    }
      NULL
      })
    } else {
    if (adaptive.resampling!=1) {
      ESS<-(sum(Current.Weights.with.Prev)^2)/sum((Current.Weights.with.Prev)^2)
      pf_debug("ESS is ", ESS, verbose=verbose)
if (is.na(ESS)) {
	ESS=1
	Current.Weights.with.Prev.mat<-cbind(weight_stack_matrix(Weights.stack, Weights.stack.start, Weights.stack.active), Current.Weights)
	tmpfile<-paste0(tempfile(), '.RData')
    save(Current.Weights.with.Prev, Current.Weights.with.Prev.mat, Current.Weights, file=tmpfile)
    warning(paste('find the unexpected weights saved in', tmpfile))
	}		
    } else {
      ESS<-1
    }
    }
    did_resample<-((ESS<(nParticles*adaptive.resampling) & Time.Period>1) | Time.Period == nrow(in.Data$Indices$Main.Index))
    if (isTRUE(profile.phases)) {
      profile_phase("resampling", {
    if (did_resample) {	
	  if (length(unique(Current.Weights.with.Prev))==1) Current.Weights.with.Prev<-rep(1, nParticles)
	  Rows<-try(sample.int(nParticles, replace = TRUE, prob = Current.Weights.with.Prev))
			#if (class(Rows)=="try-error") {
			#	print(Current.Weights.with.Prev) 
			#	print(Rows)
			#	stop("weird..\n")														 
		#}#
		#}#
      ResampleCount<-ResampleCount+1
      pf_debug(" - resampling", ResampleCount, "\n", verbose=verbose)
      Results.stack<-Results.stack[Rows,]
	  Weights.stack[,1]<-1/nParticles
	  Weights.stack.active<-1L
	  Weights.stack.start<-1L
    } else {
      Rows<-1:nParticles # no resampling
      pf_debug("\n", verbose=verbose)
    }
      NULL
      })
    } else if (did_resample) {	
	  if (length(unique(Current.Weights.with.Prev))==1) Current.Weights.with.Prev<-rep(1, nParticles)
	  Rows<-try(sample.int(nParticles, replace = TRUE, prob = Current.Weights.with.Prev))
      ResampleCount<-ResampleCount+1
      pf_debug(" - resampling", ResampleCount, "\n", verbose=verbose)
      Results.stack<-Results.stack[Rows,]
	  Weights.stack[,1]<-1/nParticles
	  Weights.stack.active<-1L
	  Weights.stack.start<-1L
    } else {
      Rows<-1:nParticles # no resampling
      pf_debug("\n", verbose=verbose)
    }
	New.weights<-Current.Weights[Rows]

    #=====================================================
    #               Smoothing
    #====================================================
    
    #======================================================
    # start OF SMOOTHING POINTS ESTIMATION
    
    #  ==================================
    # ver. 1.8 - now use one function for all particle propagations:
    #	cat("  smoothing new particles")
    #
    #	New.Particles.MH<-propagate.particles(Last.Particles= Results.stack[,ncol(Results.stack)], Current.Proposal=Current.Proposal, parallel=parallel, Parameters=Parameters, mycl=mycl)
    #
    #	Current.Phys.Weights.MH<-in.Data$Spatial$Phys.Mat[New.Particles.MH,Time.Period]
    #	
    #   if (!is.na(k) & Time.Period > 1) {
    #	  Angles.probs<-get.directional.weights(in.Data, Results.stack[,ncol(Results.stack)-1], Results.stack[,ncol(Results.stack)], New.Particles.MH, k, parallel, mycl, D.kappa)
    #      Current.Weights.MH<- Current.Phys.Weights.MH*Angles.probs
    #    }
    #    else {Current.Weights.MH<-Current.Phys.Weights.MH*D.kappa}
    
    #   end of smoothing estimation    
    #============================================
    #============================================
    
    # now I need to compare resampled results with new resuklts
    #according to MH
    #	v<-runif(nParticles)
    #	Accepted.Particles<-New.Particles.MH
    #	Ratio<-Current.Weights.MH/(Current.Weights[Rows])
    #	Ratio[is.na(Ratio)]<-1 # NA is caused by 0/0
    #	Eq<-!as.logical(floor(pmin(1, Ratio)-v)+1)
    #	Accepted.Particles[Eq]<-New.Particles[Rows][Eq]
    #	New.weights<-Current.Weights.MH
    #	New.weights[Eq]<-Current.Weights[Rows][Eq]
    #	cat("smoothed ", nParticles-sum(Eq), "Particles \n")
    
    #============================================
    # this line is needed instead of all MH stuff
    Accepted.Particles<-if (isTRUE(profile.phases)) {
      profile_phase("accepted_particles_assignment", New.Particles[Rows])
    } else {
      New.Particles[Rows]
    }
    
    #All.results<-paste(All.results, Accepted.Particles, sep=".")
    ## adding points to the results matrix
    if (isTRUE(profile.phases)) {
      profile_phase("results_stack_append", {
        Results.stack<-cbind(Results.stack, Accepted.Particles)
        NULL
      })
    } else {
      Results.stack<-cbind(Results.stack, Accepted.Particles)
    }
    #==============================================================
    # 1.6
    if (isTRUE(profile.phases)) {
      profile_phase("weights_stack_append", {
        Weights.stack.active<-Weights.stack.active+1L
        if (Weights.stack.active>Weights.stack.max.cols) stop("Weights.stack exceeded preallocated capacity")
        Weights.stack[, weight_stack_last_column(Weights.stack.start, Weights.stack.active, Weights.stack.max.cols)]<-New.weights
        NULL
      })
    } else {
      Weights.stack.active<-Weights.stack.active+1L
      if (Weights.stack.active>Weights.stack.max.cols) stop("Weights.stack exceeded preallocated capacity")
      Weights.stack[, weight_stack_last_column(Weights.stack.start, Weights.stack.active, Weights.stack.max.cols)]<-New.weights
    }
    
    #Prev.Weights<-(Prev.Weights[Rows]) * New.weights
    #Prev.Weights<-Prev.Weights/mean(Prev.Weights)
    # end of 1.6 addition
    #===================================================
    
    #Last.Last.Particles<-Last.Particles[Rows]
    #Last.Particles<-Accepted.Particles
    
    #	if (Time.Period==2) Particles.at.2<-Accepted.Particles
    #	if (Time.Period>2) {
    #		Particles.at.2<-Particles.at.2[Rows]
    pf_debug("   unique P in point leaving stack", length(unique(Results.stack[,1])), "\n", verbose=verbose)
    #	}
    
    if (Time.Period<=L) pf_debug("creating stack\n", verbose=verbose)
	if (Time.Period==L) 	Points<-vector(mode = "list")
    if (Time.Period>L) {
	# the new idea is that we could skip the saving all results and save just points and transitions - we are not outputting them anyways... THis will help avoiding the sort of All.results, that proved to be very slow..
      # save points
      #if (is.null(Points))  Points<-vector(mode = "list")
	  Rle<-if (isTRUE(profile.phases)) {
	    profile_phase("point_rle_creation", {
	      tmp_rle<-bit::intrle(sort.int(Results.stack[,1], method="quick"))
          if (is.null(tmp_rle)) tmp_rle<-rle(sort.int(Results.stack[,1], method="quick"))
          tmp_rle
	    })
	  } else {
	    tmp_rle<-bit::intrle(sort.int(Results.stack[,1], method="quick"))
        if (is.null(tmp_rle)) tmp_rle<-rle(sort.int(Results.stack[,1], method="quick"))
        tmp_rle
	  }
	  if (isTRUE(profile.phases)) n_unique_point_rle_values<-length(Rle$values)
	  Points[length(Points)+1]<-list(Rle)
      #  All.results<-paste(Results.stack[,1], sep=".")
      #} else {
      #  All.results<-paste(All.results, Results.stack[,1], sep=".")
      #}
      # save transitions
      Trans[[Time.Period-L]]<-if (isTRUE(profile.phases)) {
        profile_phase("transition_rle_creation", get.transition.rle(Results.stack[,1], Results.stack[,2], n_grid=nrow(in.Data$Spatial$Grid)))
      } else {
        get.transition.rle(Results.stack[,1], Results.stack[,2], n_grid=nrow(in.Data$Spatial$Grid))
      }
      if (isTRUE(profile.phases)) n_unique_transition_rle_values<-length(Trans[[Time.Period-L]]$values)
      # clean Results.stack
      if (isTRUE(profile.phases)) {
        profile_phase("stack_drop", {
          Results.stack<-Results.stack[,-1]
          Weights.stack.start<-(Weights.stack.start %% Weights.stack.max.cols)+1L
          Weights.stack.active<-Weights.stack.active-1L
          NULL
        })
      } else {
        Results.stack<-Results.stack[,-1]
        # clean Weights.stack
        Weights.stack.start<-(Weights.stack.start %% Weights.stack.max.cols)+1L
        Weights.stack.active<-Weights.stack.active-1L
      }
pf_debug("******************\n", verbose=verbose)
    }
    if (isTRUE(profile.phases)) {
      phase_profile[[Time.Period]]<-data.frame(
        Time.Period=Time.Period,
        nParticles=nParticles,
        total_length=total_length,
        n_unique_last_particles=n_unique_last_particles,
        n_groups=n_groups,
        n_moving_groups=n_moving_groups,
        n_unique_new_particles=n_unique_new_particles,
        ESS=ESS,
        did_resample=did_resample,
        ResampleCount=ResampleCount,
        Results.stack_ncol=ncol(Results.stack),
        Weights.stack_ncol=Weights.stack.active,
        n_unique_point_rle_values=n_unique_point_rle_values,
        n_unique_transition_rle_values=n_unique_transition_rle_values,
        n_outliers=length(in.Data$outliers),
        Results.stack_size_bytes=as.numeric(utils::object.size(Results.stack)),
        Weights.stack_size_bytes=as.numeric(utils::object.size(Weights.stack)),
        total_iteration_time=proc.time()[["elapsed"]]-iteration_start,
        as.list(phase_times),
        stringsAsFactors=FALSE
      )
    }
  }
  #####################################
  #save(All.results, file="All.results.usmoothed.RData")
  # now we need to add final point!
  
  final_phase_times<-stats::setNames(rep(0, length(phase_names)), phase_names)
  if (known.last) {
    final_start<-proc.time()[["elapsed"]]
    Results.stack<-pf.final.smoothing(in.Data, Results.stack, precision.sd=precision.sd, nParticles=nParticles, last.particles=Results.stack[,ncol(Results.stack)])
    if (isTRUE(profile.phases)) final_phase_times[["final_smoothing"]]<-proc.time()[["elapsed"]]-final_start
  }
  
  if (!is.list(Points)) Points<-vector(mode = "list")
  # and here we need to add a thing that will finish All.results and Trans from the points that are still in the stack
  pf_message("adding last points form the stack to the resutls\n", verbose=verbose)
  Length<-ncol(Results.stack)

  final_flush_start<-proc.time()[["elapsed"]]
  for (rest in 1:Length) {
    # save points
	  Rle<-bit::intrle(sort.int(Results.stack[,rest], method="quick"))
      if (is.null(Rle)) Rle<-rle(sort.int(Results.stack[,rest], method="quick"))
	  Points[length(Points)+1]<-list(Rle)

	#All.results<-paste(All.results, Results.stack[,rest], sep=".")
    if (rest<Length) {
      # save transitions
      #Trans[[Time.Period-L+rest]]<-get.transition.rle(Results.stack[,rest], Results.stack[,rest+1])
      Trans[[length(Trans)+1]]<-get.transition.rle(Results.stack[,rest], Results.stack[,rest+1], n_grid=nrow(in.Data$Spatial$Grid))
    }
  }
  if (isTRUE(profile.phases)) {
    final_phase_times[["final_stack_flush"]]<-proc.time()[["elapsed"]]-final_flush_start
    phase_profile[[total_length+1L]]<-data.frame(
      Time.Period=NA_integer_,
      nParticles=nParticles,
      total_length=total_length,
      n_unique_last_particles=NA_integer_,
      n_groups=NA_integer_,
      n_moving_groups=NA_integer_,
      n_unique_new_particles=NA_integer_,
      ESS=NA_real_,
      did_resample=NA,
      ResampleCount=ResampleCount,
      Results.stack_ncol=ncol(Results.stack),
      Weights.stack_ncol=Weights.stack.active,
      n_unique_point_rle_values=NA_integer_,
      n_unique_transition_rle_values=NA_integer_,
      n_outliers=length(in.Data$outliers),
      Results.stack_size_bytes=as.numeric(utils::object.size(Results.stack)),
      Weights.stack_size_bytes=as.numeric(utils::object.size(Weights.stack)),
      total_iteration_time=sum(final_phase_times),
      as.list(final_phase_times),
      stringsAsFactors=FALSE
    )
  }
  #if (parallel)   parallel::clusterEvalQ(mycl, rm(Parameters)) 
  #if (length(existing.cluster)==1) parallel::stopCluster(cl = mycl)
  if (sink2file) sink()
  tmp.results<-list(AB.distance=in.Data$AB.distance, AC.distance2=in.Data$AC.distance2, Dif.ang=in.Data$Dif.ang)

  if (isTRUE(profile.phases)) {
    phase_profile<-do.call(rbind, phase_profile[!vapply(phase_profile, is.null, logical(1))])
  } else {
    phase_profile<-NULL
  }

  partial_cache_diagnostics<-partial.movement.cache.diagnostics(in.Data$Spatial$tmp$partial_movement_cache)
  return(list(Points=Points, Trans=Trans, Results=list(outliers=in.Data$outliers, tmp.results=tmp.results, phase_profile=phase_profile, partial_cache_diagnostics=partial_cache_diagnostics)))
}


get.coordinates.PF<-function(Points, in.Data, add.jitter=FALSE, verbose=FALSE) {
  #library("aspace")
  # this function will extract point coordinates from the output matrix.. 
  # the question is do we need only mean and sd or also median and quantiles?
  # I will start from mean and SD
  #plot(c(min(in.Data$Spatial$Grid[,1]),max(in.Data$Spatial$Grid[,1])), c(min(in.Data$Spatial$Grid[,2]),max(in.Data$Spatial$Grid[,2])), type="n")
  #log <- capture.output({
  #  Means=aspace::calc_box(id=1,  points=in.Data$Spatial$Grid[inverse.rle(Points[[1]]),1:2])
  # 
  #for (i in 2:length(Points)) {
  #  Means[i,]=aspace::calc_box(id=i,  points=in.Data$Spatial$Grid[inverse.rle(Points[[i]]), 1:2])
    #plot_box(plotnew=FALSE, plotpoints=FALSE)
  #}
	#})
   if (length(in.Data$Results)==0) in.Data$Results<-list()
  #in.Data$Results$Final.Means<-Means
  
  #############
  # new part for medians
  pf_message("estimating quantiles for positions\n", verbose=verbose)
	
		# check whether Grid was over dateline:
	overdateline<-ifelse(attr(in.Data$Spatial$Grid, 'left')>	attr(in.Data$Spatial$Grid, 'right'), TRUE, FALSE)

	
	Quantiles<-c()
	CIntervals<-c()
	cur_Grid<-in.Data$Spatial$Grid
	if (overdateline) cur_Grid[cur_Grid[,1]<0,1]<-cur_Grid[cur_Grid[,1]<0,1]+360
	
	for (i in 1:length(Points)) {
	

	#Mode_cur<-	cur_Grid[Points[[i]]$values[which.max(Points[[i]]$lengths)],1]
	#cur_Grid[,1]<-ifelse(cur_Grid[,1]<Mode_cur-180, cur_Grid[,1]+360, cur_Grid[,1])
	Quantiles<-rbind(Quantiles, c(summary(cur_Grid[inverse.rle(Points[[i]]),2]), Mode=cur_Grid[Points[[i]]$values[which.max(Points[[i]]$lengths)],2], summary(cur_Grid[inverse.rle(Points[[i]]),1]), Mode=cur_Grid[Points[[i]]$values[which.max(Points[[i]]$lengths)],1]))
	CIntervals<-rbind(CIntervals, c(stats::quantile(cur_Grid[inverse.rle(Points[[i]]),2], probs = c(0.025, 0.975)), stats::quantile(cur_Grid[inverse.rle(Points[[i]]),1], probs = c(0.025, 0.975))))
	}
	Quantiles<-as.data.frame(Quantiles)
	names(Quantiles)[1:6]<-paste(names(Quantiles)[1:6], "lat", sep="")
	names(Quantiles)[8:13]<-paste(names(Quantiles)[8:13], "lon", sep="")

	###########
	# doing jitter first
	if (add.jitter) {
    in.Data$Results$Quantiles<-Quantiles
	pf_message("adding jitter to medians\n", verbose=verbose)
	jitter_coords<-get.coords.jitter(in.Data)
	if (!is.null(jitter_coords)) {
	Quantiles$MedianlonJ<-jitter_coords[,1]
	Quantiles$MedianlatJ<-jitter_coords[,2]
	}
	} else {
	Quantiles$MedianlonJ<-Quantiles$Medianlon
	Quantiles$MedianlatJ<-Quantiles$Medianlat
	}
	names(Quantiles)<-gsub("\\s","", names(Quantiles))
	names(Quantiles)<-gsub("1","F", names(Quantiles))
	names(Quantiles)<-gsub("3","T", names(Quantiles))
	
	pf_message("adding 95% credibility intervals to medians\n", verbose=verbose)
	Quantiles$LCI.lat<-CIntervals[,1]
	Quantiles$UCI.lat<-CIntervals[,2]
	Quantiles$LCI.lon<-CIntervals[,3]
	Quantiles$UCI.lon<-CIntervals[,4]
	  
	if (overdateline) {
	   Columns<-c(8:15, 19, 20)
       for (i in Columns) {
           Quantiles[Quantiles[,i]>180,i]<-Quantiles[Quantiles[,i]>180,i]-360
       }	   
	}
	if (nrow(Quantiles)==nrow(in.Data$Indices$Matrix.Index.Table)) {
	   Quantiles<-cbind(Quantiles, time=in.Data$Indices$Matrix.Index.Table$time)
	} else {
	   Quantiles<-cbind(Quantiles[-1,], time=in.Data$Indices$Matrix.Index.Table$time)
	}
	in.Data$Results$Quantiles<-Quantiles
	
    in.Data$Results$Points.rle<-Points
  return(in.Data)
}


estimate.movement.parameters<-function(Trans, in.Data, fixed.parameters=NA, a=45, b=500, parallel=FALSE, existing.cluster=NULL, estimatetruncnorm=FALSE, nParticles=1e6, verbose=FALSE) {
  mycl=existing.cluster
  
  # old name update.proposal.PF
  # main function for proposal update
  # from the version 1.5 this function will do estimation of distance and it's SD..

  #=========================================
  # here is another place to change...
  # ver. 1.8
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #Trans<-vector(mode = "list", length = (dim(output.matrix)[2]-1))
  #cat("   extracting indexes\n")
  #for (i in 1:(dim(output.matrix)[2]-1)) {
  #  Trans[[i]]<-get.transition.rle(output.matrix[,i], output.matrix[,i+1])
  #}
  pf_message("   estimating distances\n", verbose=verbose)
  #####   let's try to get distance distribution:
  Distances<-Trans
  dist.fun<-function(Movement_Points) {
   #in.Data$Spatial$tmp$Distance[Movement_Points[,1], Movement_Points[,2]]
    # sp::spDists(in.Data$Spatial$Grid[Movement_Points[,1], c(1,2), drop=FALSE],
    #             in.Data$Spatial$Grid[Movement_Points[,2], c(1,2), drop=FALSE],
    #             longlat=TRUE,
    #             diagonal=TRUE)
    
    as.numeric(sf::st_distance(sf::st_as_sf(as.data.frame(in.Data$Spatial$Grid[Movement_Points[,1], c(1,2), drop=FALSE]),coords=c('lon','lat'),crs=4326),
                sf::st_as_sf(as.data.frame(in.Data$Spatial$Grid[Movement_Points[,2], c(1,2), drop=FALSE]),coords=c('lon','lat'),crs=4326),
                by_element = TRUE)/1000)
  }
  for (i in 1:length(Trans)) {
    Movement_Points<-decode_transition_rle_values(Trans[[i]], n_grid=nrow(in.Data$Spatial$Grid))
    Distances[[i]]$values<-dist.fun(Movement_Points)
  }
  pf_message("   estimating directions\n", verbose=verbose)
  #ok, now we want to get directions
  
  Directions<-Trans
  for (i in 1:length(Trans)) {
    Movement_Points<-decode_transition_rle_values(Trans[[i]], n_grid=nrow(in.Data$Spatial$Grid))
    Directions[[i]]$values<-apply(Movement_Points, 1, FUN= dir_fun, in.Data)
  }
  pf_message("   estimating mean directions and kappas\n", verbose=verbose)
  
  Mean.Directions<-unlist(lapply(Directions, FUN=function(x) CircStats::circ.mean(circular::circular(inverse.rle(list(lengths=x$lengths[!is.na(x$values)], values=x$values[!is.na(x$values)]*pi/180))))*180/pi))
  #plot(Mean.Directions)
  pf_message("   estimating kappas\n", verbose=verbose) #CircStats::est.kappa
  Kappas<-unlist(lapply(Directions, FUN=function(x) CircStats::est.kappa(inverse.rle(list(lengths=x$lengths[!is.na(x$values)], values=x$values[!is.na(x$values)]*pi/180)))))
  
  
  # now we want to get mean distance
  pf_message("   estimating mean dists\n", verbose=verbose)
  #
  #cat("   estimating mean and SD to report dists SD\n")
  # attempt to go just for simple distance estimation...
  Mean2report<-unlist(lapply(Distances, FUN=function(x) mean(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))
  SD2report<-unlist(lapply(Distances, FUN=function(x) stats::sd(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))
  pf_message("   estimating probs of migration\n", verbose=verbose)
  # ok, now we want to get parameters for distances
  #
  Probability.of.migration<-unlist(lapply(Distances, FUN=function(x) sum(x$lengths[x$values!=0])))/nParticles
  
  ## 
    if (estimatetruncnorm) {
  Mean.and.Sigma<-lapply(Distances, FUN=function(x)  mu.sigma.truncnorm(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])), a=a, b=b))
  #}
  Mean.Dists<-sapply(Mean.and.Sigma, "[[", i=1)
  pf_message("   estimating dists SD\n", verbose=verbose)
  Mean.SD<-sapply(Mean.and.Sigma, "[[", i=2)
  } else {
  Mean.Dists<-Mean2report
  Mean.SD<-SD2report
  }

  #unlist(lapply(Distances, FUN=function(x) mean(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))
  #plot(Mean.Dists)
  pf_message("   estimating median dists\n", verbose=verbose)
  Median.Dists<-unlist(lapply(Distances, FUN=function(x) stats::median(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))

  pf_message("   creating output", verbose=verbose)
  ################
  # create new proposal from the posteriors
  #New.Matrix.Index.Table<-data.frame(Direction=Mean.Directions, M.mean=Mean.Dists, M.sd=Mean.SD, Decision=Probability.of.migration, Kappa=Kappas)
  ################
  
  Movement.results<-in.Data$Indices$Matrix.Index.Table
  #all.arrays.object<-in.Data
  
  #----------------
  # actually I do not want to update this - so I rather save it outside if possible..
  Movement.results$Decision<-Probability.of.migration
  Movement.results$Direction<-Mean.Directions
  Movement.results$Kappa<-Kappas
  Movement.results$M.mean<-Mean.Dists
  Movement.results$M.sd<-Mean.SD
  Movement.results$M.medians<-Median.Dists
  Movement.results$Mean2report<-Mean2report
  Movement.results$SD2report<-SD2report
  
  #all.arrays.object$Indices$Matrix.Index.Table<-New.Matrix.Index.Table
  
  #all.arrays.object$Indices$Main.Index$Biol.Prev=1:(nrow(all.arrays.object$Indices$Matrix.Index.Table))
    
  if (is.list(fixed.parameters)) {
    # this part can be moved upper not to estimate parametrs in case they are fixed.
    if ("M.sd" %in% names(fixed.parameters)) Movement.results$M.sd<-fixed.parameters$M.sd
    
    if ("M.mean" %in% names(fixed.parameters)) Movement.results$M.mean<-fixed.parameters$M.mean
    
    if ("Kappa" %in% names(fixed.parameters)) Movement.results$Kappa<-fixed.parameters$Kappa
    
    if ("Direction" %in% names(fixed.parameters)) Movement.results$Direction<-fixed.parameters$Direction
  }
  #all.arrays.object$Indices$Main.Index$Biol.Next=2:(nrow(all.arrays.object$Indices$Matrix.Index.Table))
  pf_message(" DONE!\n", verbose=verbose)
  Res<-list(Movement.results=Movement.results, Transitions.rle=Trans)
  return(Res)
}


transition_base<-function(n_grid) {
  if (length(n_grid)!=1 || is.na(n_grid) || n_grid<1) stop("n_grid must be a positive scalar")
  as.numeric(n_grid)+1
}

encode_transition<-function(from, to, n_grid) {
  base<-transition_base(n_grid)
  if (any(from<1 | to<1 | from>n_grid | to>n_grid, na.rm=TRUE)) stop("Transition indices must be between 1 and n_grid")
  as.numeric(from)*base+as.numeric(to)
}

decode_transition<-function(key, n_grid=NULL, transition_base=NULL, warn_legacy=FALSE) {
  if (is.null(transition_base)) {
    if (!is.null(n_grid) && n_grid>=1e5 && isTRUE(warn_legacy)) {
      warning("Transition base metadata is missing; decoding with legacy base 1e5. Large-grid legacy objects may decode incorrectly.", call.=FALSE)
    }
    transition_base<-1e5
  }
  Movement_Points<-cbind(from=as.numeric(key)%/%transition_base,
                         to=as.numeric(key)%%transition_base)
  storage.mode(Movement_Points)<-"integer"
  Movement_Points
}

decode_transition_rle_values<-function(transition_rle, values=transition_rle$values, n_grid=NULL) {
  base<-attr(transition_rle, "transition_base", exact=TRUE)
  if (is.null(base)) {
    attr_n_grid<-attr(transition_rle, "n_grid", exact=TRUE)
    if (!is.null(attr_n_grid)) base<-transition_base(attr_n_grid)
  }
  decode_transition(values, n_grid=n_grid, transition_base=base, warn_legacy=TRUE)
}

get.transition.rle=function(From, To, n_grid=max(c(From, To), na.rm=TRUE)) {
  base<-transition_base(n_grid)
  Res<-rle(sort.int(encode_transition(From, To, n_grid), method = "quick"))
  attr(Res, "transition_base")<-base
  attr(Res, "n_grid")<-as.integer(n_grid)
  Res
}


mu.sigma.truncnorm<-function(x, a=45, b=500) {
# this is used optionally
  if (length(unique(x))>1) {
    tr.norm<-function(prm) {
      sum(-log(truncnorm::dtruncnorm(as.numeric(x),a=a,b=b,mean=prm[1],sd=prm[2])))
    }
    Res=try(stats::optim(c(mean(x),stats::sd(x)), tr.norm, method="BFGS"))
    if (inherits(Res, "try-error")) {
	tmpfile<-paste0(tempfile(), '.RData')
      save(x, Res, file=tmpfile)
      Res$par<-c(mean(x), stats::sd(x))
     warning(paste('find the unexpected mu.sigma.truncnorm estimates saved in', tmpfile))
    }
    return(c(Res$par[1], Res$par[2]))
  } else {
    return(c(mean(x), stats::sd(x)))
  }
}


coords.aeqd.jitter <- function(coords, r, n)
{
# coords is a vector of leghth 2 c(lon, lat)
# r should be in meters
# made on te basis of this:
#"http://gis.stackexchange.com/questions/121489/1km-circles-around-lat-long-points-in-many-places-in-world"
 stopifnot(length(coords) == 2)
  
  	#p = sp::SpatialPoints(matrix(coords, ncol=2), proj4string=sp::CRS("+proj=longlat +datum=WGS84"))
  	p = sf::st_as_sf(as.data.frame(matrix(coords, ncol=2)), coords=c('V1','V2'),crs = 4326)  
  	
     # aeqd <- sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
     #                 p@coords[[2]], p@coords[[1]])
     aeqd <- sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                     sf::st_coordinates(p)[[2]], sf::st_coordinates(p)[[1]])  
    
    #projected <- sp::spTransform(p, sp::CRS(aeqd))
    projected <- sf::st_transform(p,sf::st_crs(aeqd))
    
    #buffered <- rgeos::gBuffer(projected, width=r, byid=TRUE)
    buffered <- sf::st_buffer(projected, dist = r)
    
     # lambert <- sprintf("+proj=laea +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
     #                 p@coords[[2]], p@coords[[1]])
     lambert <- sprintf("+proj=laea +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                     sf::st_coordinates(p)[[2]], sf::st_coordinates(p)[[1]])  
    
	#buffered_eqarea <- sp::spTransform(buffered, sp::CRS(lambert))
	buffered_eqarea <- sf::st_transform(buffered,sf::st_crs(lambert))
	
	
	#random_points<-sp::spsample(buffered_eqarea,n=n,type="random")
	random_points <- sf::st_sample(buffered_eqarea, size = n, type = "random")
	
	# if (is.null(random_points)) random_points<-sp::spsample(buffered_eqarea,n=n,type="random", iter=40) 
	if (is.null(random_points)){
	  # Loop until you have the desired number of samples
	  counter<-1
	  while (nrow(random_points) == 0 & counter <= 40) {
	    # Sample a subset of points
	    random_points <- sf::st_sample(buffered_eqarea, size = n, type = "random")
	    counter<-counter+1
	  }

	  }

	
	if (is.null(random_points)) random_points<-p
    #sp::spTransform(random_points, p@proj4string)
    sf::st_transform(random_points, crs = sf::st_crs(p))
    
}

# wrapper for jitter
get.coords.jitter<-function(in.Data) {
	Distance<-in.Data$Spatial$tmp$Distance
	 #if (is.null(Distance)) Distance=sp::spDists(in.Data$Spatial$Grid[,1:2], longlat=TRUE)
   if (is.null(Distance)) Distance=sf::st_distance(sf::st_as_sf(as.data.frame(in.Data$Spatial$Grid[,1:2]),coords=c('lon','lat'),crs=4326))/1000
	Distance<-as.numeric(Distance[,1])
	JitRadius<-min(Distance[Distance>0])/2*1000 # in meters
	#now I want to generate random poitns in the radius of this
	coords=cbind(in.Data$Results$Quantiles$Medianlon, in.Data$Results$Quantiles$Medianlat)
	tmp<-try(apply(coords, 1,  coords.aeqd.jitter, r=JitRadius, n=1 ))
	jitter_coords<-NULL
	if( !inherits(tmp, "try-error")) {
	#jitter_coords<-t(sapply(tmp, sp::coordinates))
	jitter_coords<-t(sapply(tmp, sf::st_coordinates))
	
	}
	return(jitter_coords)
}


get.LL.PF<-function(in.Data, data) {
  # needed to estimate log likelihood of the optimization
  L=0
  if (is.list(data)) {
  for (i in 1:(length(data)-1)) { 
      L=L+log(mean(in.Data$Spatial$Phys.Mat[inverse.rle(data[[i]]),i]))
	  }
  } else{ 
    for (i in 1:(dim(data)[2]-1)) {
    L=L+log(mean(in.Data$Spatial$Phys.Mat[data[,i],i]))
	}
  }
  #L=L/(dim(All.results.mat)[2]-1)
  #LL=-sum(log(L))
  LL=-L
  return(LL)
}


pf.par.internal<-function(x, Current.Proposal) {
  # this simple function is needed to save I/0 during parallel run. Being executed at a slave it creates inoput for the next function taking Current proposal from the slave and Parameters from the master
  new.Parameters<-c(x=list(x), get("Parameters"), Current.Proposal=list(Current.Proposal))
  Res<-do.call( generate.points.dirs, new.Parameters)
  return(Res)
}


flightr.dvonmises<-function(x, mykap) {
  return(as.numeric(suppressWarnings(circular::dvonmises(x[[2]], mu=x[[1]], kappa=mykap))))}

  
pf.final.smoothing<-function(in.Data, results.stack, precision.sd=25, nParticles=1e6, last.particles=NA) {
  # this function simply resamples final points proportionally to the distance to known finish.
  Final.point.real<-in.Data$Spatial$stop.point
  # now we want to get distances.. I'll not index it as we will do this only once..
  Final.points.modeled=last.particles
  # Weights<-stats::dnorm(  sp::spDists(in.Data$Spatial$Grid[Final.points.modeled, c(1,2), drop=FALSE],
  #                                     in.Data$Spatial$Grid[Final.point.real, c(1,2), drop=FALSE],
  #                                     longlat=TRUE),
  #                         mean=0,
  #                         sd=precision.sd)
  Weights<-stats::dnorm(as.numeric(sf::st_distance(
    sf::st_as_sf(as.data.frame(in.Data$Spatial$Grid[Final.points.modeled, c(1,2), drop=FALSE]),coords=c('lon','lat'),crs=4326),
    sf::st_as_sf(as.data.frame(in.Data$Spatial$Grid[Final.point.real, c(1,2), drop=FALSE]),coords=c('lon','lat'),crs=4326)
                                          )/1000),
                          mean=0,
                          sd=precision.sd)

  Rows<- try(suppressWarnings(sample.int(nParticles, replace = TRUE, prob = Weights/sum(Weights))))
  if (inherits(Rows ,'try-error')) {
    tmpfile<-paste0(tempfile(), '.RData')
	save(last.particles, Weights, file=tmpfile)
    warning(paste('final smoothing failed, error data saved to the', tmpfile))
    return(results.stack)
	} else {
    return(results.stack[Rows,])
	}
}

dir_fun<-function(x, in.Data) {
	  geosphere::bearing(in.Data$Spatial$Grid[x[[1]], c(1,2), drop=FALSE], in.Data$Spatial$Grid[x[[2]], c(1,2), drop=FALSE])
}
    
dist.fun<-function(x, Result, transition_rle=NULL) {
  Movement_Points<-decode_transition_rle_values(transition_rle, x, n_grid=nrow(Result$Spatial$Grid))
  # sp::spDists(Result$Spatial$Grid[Movement_Points[,1], c(1,2), drop=FALSE], 
  #             Result$Spatial$Grid[Movement_Points[,2], c(1,2), drop=FALSE],
  #             longlat=TRUE, 
  #             diagonal=TRUE)
  as.numeric(sf::st_distance(sf::st_as_sf(as.data.frame(Result$Spatial$Grid[Movement_Points[,1], c(1,2), drop=FALSE]),coords=c('lon','lat'),crs=4326), 
                  sf::st_as_sf(as.data.frame(Result$Spatial$Grid[Movement_Points[,2], c(1,2), drop=FALSE]),coords=c('lon','lat'),crs=4326),
                  by_element = TRUE)/1000)
}

lazy.result.plot<-function(Result) {
    graphics::plot(Meanlat~Meanlon, type="p", data=Result$Results$Quantiles, pch=3, col="blue", main="mean positions")
    graphics::points(Meanlat~Meanlon, type="p", data=Result$Results$Quantiles, pch=3, col="blue")
    graphics::lines(Meanlat~Meanlon, data=Result$Results$Quantiles, col="blue")
    maps::map('world2', add=TRUE)
}



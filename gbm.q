//general math settings
pi:acos -1

// generate n normal variables with mean m, standard deviation sd
rnorm: {[n;m;sd]
    u1: n?1f;
    u2: n?1f;
    m + sd * sqrt[-2*log u1] * cos 2*u2*pi};

/// parameters: s0,r,vol,T,timestep,N, monte-carlo method, default general solution path
/// usage example - gbm[100f;0.1;0.2;1;252;10;`]
.gbm.path:{[s0;mu;vol;T;timestep;N;method]
	dt:T % timestep;
	drift:dt * mu - vol*vol % 2;
	diffusion: vol * rnorm[timestep*N;0;1] * sqrt dt;
	euler_term:1 + (mu*dt) + vol * rnorm[timestep*N;0;1] * sqrt dt;
	
	if[not method in ``general`euler;'"No path options available"];

	//default general solution
	// path options
	if[method in ``general;
		path: (enlist N#s0),s0 *\ N cut exp drift + diffusion];
	if[method in `euler;
		path: (enlist N#s0),s0 *\ N cut euler_term];
	
	path}

// analytical mean and variance
.gbm.am:{[s0;mu;T] s0 * exp mu*T};                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
.gbm.av:{[s0;mu;vol;T](s0 * s0 * exp 2*mu*T) * (-1 +exp vol*vol*T)};

// summary for numerical and analytical values of path
.gbm.stat:{[path]

	/numerical mean, variance
	nm:avg last path;
	nv:var last path;

	flip `type`mean`variance!(`Numerical`Analytical;(nm;.gbm.am); (nv;.gbm.av))}

\
//test case:
s0:100f
mu:0.1
vol:0.2
T:1
timestep:252
N:10
.gbm.path[100f;0.1;0.2;1;252;10;`]
.gbm.am[100f;0.1;1]
.gbm.av[100f;0.1;0.2;1]
/

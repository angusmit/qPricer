// load required script
\l constant.q

// tracking table schema
.gbm.tab:([] method:`$() ;s0:`float$(); mu:`float$(); vol:`float$(); T:(); timestep:`int$(); N:`int$();Amean:`float$();Avar:`float$();Nmean:`float$();Nvar:`float$())

/// parameters: s0,mu,vol,T,timestep,N, monte-carlo method, default general solution path
/// usage example - .gbm.path[100f;0.1;0.2;1;252;10;`]
.gbm.path:{[s0;mu;vol;T;timestep;N;method]
	dt:T % timestep;
	drift:dt * mu - vol*vol % 2;
	diffusion: vol * .const.norm01[timestep*N] * sqrt dt;
	euler_term:1 + (mu*dt) + vol * .const.norm01[timestep*N] * sqrt dt;
	
	if[not method in ``general`euler;'"No path options available"];

	//default general solution
	// path options
	if[method in ``general; 
		method:`general;
		path: (enlist N#s0),s0 *\ N cut exp drift + diffusion];
	if[method in `euler;
		path: (enlist N#s0),s0 *\ N cut euler_term];

	// analytical mean and variance
	am:s0 * exp mu*T;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
	av:(s0 * s0 * exp 2*mu*T) * (-1 +exp vol*vol*T);
	nm:avg last path;
	nv:var last path;
	`.gbm.tab upsert (method;s0;mu;vol;T;timestep;N;am;av;nm;nv);
	path}

// analytical mean and variance
.gbm.a:{[s0;mu;vol;T]`mean`variance!(s0 * exp mu*T; (s0 * s0 * exp 2*mu*T) * (-1 +exp vol*vol*T))}

// numerical mean and variance
.gbm.n:{[path]`mean`variance!(avg last path;var last path)}

/
// test case:
s0:100f
mu:0.1
vol:0.2
T:1
timestep:252
N:10
.gbm.path[100f;0.1;0.2;1;252;10;`]
.gbm.am[100f;0.1;1]
.gbm.av[100f;0.1;0.2;1]
.gbm.tab:([] method:`$() ;s0:`float$(); mu:`float$(); vol:`float$(); T:(); timestep:`int$(); N:`int$();Amean:`float$();Avar:`float$();Nmean:`float$();Nvar:`float$())
\
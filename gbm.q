/// parameters: s0,r,vol,T,timestep,N

/maths setting
pi:acos -1


/normal rv N~[0,1]
/https://wiki.thalesians.com/index.php/Programming/Kdb/Labs/Option_pricing
/or kx stat.q
nrv: {$[x=2*n:x div 2;raze sqrt[-2*log n?1f]*/:(sin;cos)@\:(2*pi)*n?1f;-1_.z.s 1+x]}

/ normal cdf
normal_cdf:{abs(x>0)-(exp[-.5*x*x]%sqrt 2*pi)*t*.31938153+t*-.356563782+t*1.781477937+t*-1.821255978+1.330274429*t:1%1+.2316419*abs x};

gbm:{[s0;r;vol;T;timestep;N]

	s0:s0;
	r:r;
	vol:vol;
	T:T;
	timestep:timestep;
	N:N;
	
	dt:T % timestep;
	
	f:{[s] s * exp (dt * r - (vol*vol % 2)) + vol * (nrv 1)[0] * sqrt dt};	
	s_0:N#s0;
	s:s_0;
	do[timestep;s_0:{f x} each s_0;s:s,'s_0];f:{[s] s * exp (dt * r - (vol*vol % 2)) + vol * (nrv 1)[0] * sqrt dt};
	s
 }

/analytical mean
/a_mean:s0 * exp r*T

/analytical variance
/a_variance: (s0 * s0 * exp 2*r*T) * (-1 +exp vol*vol*T)


bs_call:{[S0; K; r; T; volatility; dividend]

	d1: (1 % volatility * sqrt T) * (log S0 % K) + T * r + 0.5 * volatility * volatility;

	d2: d1 - volatility * sqrt T;(S0 * normal_cdf[d1] * exp neg dividend * T) - K * normal_cdf[d2] * exp neg r*T}







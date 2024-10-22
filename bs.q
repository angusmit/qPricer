// general math settings
pi:acos -1
linspace:{[s;e;n] step:(1%n) *e-s; s+step* til n+1}

.bs.tab:([] strike:(); option:(); impliedvol:())

.bs.Call:{[S0;K;r;T;vol;dividend]

	normal_cdf:{abs(x>0)-(exp[-.5*x*x]%sqrt 2*pi)*t*.31938153+t*-.356563782+t*1.781477937+t*-1.821255978+1.330274429*t:1%1+.2316419*abs x};

	d1: (1 % vol * sqrt T) * (log S0 % K) + T * r + 0.5 * vol * vol;
	d2: d1 - vol * sqrt T;

	(S0 * normal_cdf[d1] * exp neg dividend * T) - K * normal_cdf[d2] * exp neg r*T}

.bs.Put:{[S0;K;r;T;vol]
	dividend:0;
	normal_cdf:{abs(x>0)-(exp[-.5*x*x]%sqrt 2*pi)*t*.31938153+t*-.356563782+t*1.781477937+t*-1.821255978+1.330274429*t:1%1+.2316419*abs x};

	d1: (1 % vol * sqrt T) * (log S0 % K) + T * r - dividend + 0.5 * vol * vol;
	d2: d1 - vol * sqrt T;
	
	(K * normal_cdf[neg[d2]] * exp neg r*T) - S0 * normal_cdf[neg[d1]] * exp neg dividend * T}

.bs.option:{[S0;K;r;T;vol;opt] $[opt in `p`P; .bs.Put[S0;K;r;T;vol] ; .bs.Call[S0;K;r;T;vol]]}

// IV for call option
.bs.IV:{[call_price; S0; K; r; T; dividend; tol; max_iter]
	sigma_res:`;
	sigma_min:0.001;
	sigma_max:3;
	if[tol~`; tol: 10 xexp -6];
	if[max_iter~`; max_iter:1000];
	do[max_iter; sigma_mid:0.5 * sigma_min + sigma_max;
		price_mid: .bs.Call[S0;K;r;T;sigma_mid;dividend];
			$[tol>abs(price_mid - call_price); sigma_res:sigma_mid; $[price_mid > call_price; sigma_max:sigma_mid; sigma_min:sigma_mid]]];
	if[sigma_res~`;sigma_res:sigma_mid];sigma_res}
	
	
\
def find_vol_by_bisection(call_price, S0, K, r, T, dividend, tol=1e-6, max_iter=1000):
    sigma_min, sigma_max = 0.001, 2
    for _ in range(max_iter):
        sigma_mid = (sigma_min + sigma_max) / 2
        price_mid = bs_call(S0, K, r, T, sigma_mid, dividend)
        if abs(price_mid - call_price) < tol:
            return sigma_mid
        elif price_mid > call_price:
            sigma_max = sigma_mid
        else:
            sigma_min = sigma_mid
    return sigma_mid
S0:100;K:90;r:0.5;T:1;vol:0.2;dividend:0


params
S0:s0:100f;mu:0.1;vol:0.2;T:1;timestep:252;N:10
K:90;r:0.5;dividend:0

optprice:.bs.Call[s0;K;r;T;vol]
.bs.IV[optprice; s0; K; r; T;`;`;`]
case with gbm
gbm params
	
path:.gbm.path[s0;mu;vol;T;timestep;N;`]

bscall params

NOP:avg (exp neg[r]* T) * {max 0,x} each neg[K] + last path
.bs.IV[NOP; s0; K; r; T;`;`]

//risk neutral price of simulated paths
NOP:(exp -0+neg[r]* T) * avg {?[0>x;0;x]}'[last path-K]
.bs.IV[NOP; s0; K; r; T; dividend;`;`]

//implied vol testing for bs call
linspace:{[s;e;n] step:(1%n) *e-s; s+step* til n+1}
kl:linspace[S0*0.5;S0*1.5;10]
cl:{.bs.Call[S0;x;r;T;vol;dividend]}each kl
ivl:{.bs.IV[x; S0; y; r; T; dividend;`;`]}'[cl;kl]
t:([] strike:kl; option:cl; impliedvol:ivl)
linspace:{[s;e;n] step:(1%n) *e-s; s+step* til n+1}
/




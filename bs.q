// general math settings
pi:acos -1

.bs.Call:{[S0;K;r;T;vol;dividend]

	normal_cdf:{abs(x>0)-(exp[-.5*x*x]%sqrt 2*pi)*t*.31938153+t*-.356563782+t*1.781477937+t*-1.821255978+1.330274429*t:1%1+.2316419*abs x};

	d1: (1 % vol * sqrt T) * (log S0 % K) + T * r + 0.5 * vol * vol;
	d2: d1 - vol * sqrt T;

	(S0 * normal_cdf[d1] * exp neg dividend * T) - K * normal_cdf[d2] * exp neg r*T}

.bs.IV:{[call_price; S0; K; r; T; dividend; tol; max_iter]
	sigma_res:`;
	sigma_min:0.001;
	sigma_max:2;
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

call_price:.bs.Call[s0;K;r;T;vol;dividend]
.bs.IV[call_price; s0; K; r; T; dividend;`;`]
/
\
case with gbm
gbm params
s0:100f;mu:0.1;vol:0.2;T:1;timestep:252;N:10	
path:.gbm.path[s0;mu;vol;T;timestep;N;`]

bscall params
K:90;r:0.5;dividend:0
NOP:avg (exp neg[r]* T) * {max 0,x} each neg[K] + last path
.bs.IV[NOP; s0; K; r; T; dividend;`;`]

//risk neutral price of simulated paths
NOP:(exp -0+neg[r]* T) * avg {?[0>x;0;x]}'[last path-K]
.bs.IV[NOP; s0; K; r; T; dividend;`;`]
/



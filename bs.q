// load required script
\l constant.q

.bs.tab:([] strike:(); option:(); impliedvol:());
.bs.greeks:([] insertTime:`timestamp$();optionType:`$();stock:`float$();strike:`float$(); rfRate:`float$(); maturity:`float$(); volatility:`float$();
	delta:`float$(); gamma:`float$(); theta:`float$(); rho:`float$(); vega:`float$(); vanna:`float$(); vomma:`float$(); zomma:`float$());

.bs.Call:{[S0;K;r;T;vol]

	d1: (1 % vol * sqrt T) * (log S0 % K) + T * r + 0.5 * vol * vol;
	d2: d1 - vol * sqrt T;
	(S0 * .const.normal_cdf[d1] * exp neg T) - K * .const.normal_cdf[d2] * exp neg r*T};

.bs.Put:{[S0;K;r;T;vol]

	d1: (1 % vol * sqrt T) * (log S0 % K) + T * r  + 0.5 * vol * vol;
	d2: d1 - vol * sqrt T;
	
	(K * .const.normal_cdf[neg[d2]] * exp neg r*T) - S0 * .const.normal_cdf[neg[d1]] * exp neg T};

// cdf not matching with python
.bs.option:{[S0;K;r;T;vol;opt]
  d1: (1 % vol * sqrt T) * (log S0 % K) + T * r + 0.5 * vol * vol;
  d2: d1 - vol * sqrt T;

  pdf_d1: .const.normal_pdf[d1];
  cdf_d1: .const.normal_cdf[d1];
  cdf_neg_d1:.const.normal_cdf[neg[d1]];
  cdf_d2:.const.normal_cdf[d2];
  cdf_neg_d2:.const.normal_cdf[neg[d2]];

  gamma: pdf_d1 % s0 * vol * sqrt T;
  vega: s0 * pdf_d1 * sqrt T;

  vanna: vega * 1 - d1 % vol * sqrt T;
  vomma: vega * d1 * d2 % vol;
  zomma: vega * (-1 + d1 * d2)  % vol;

  $[opt in `put;
	// Put option payoff
	[delta:neg[cdf_neg_d1]; 
	  theta: (neg[S0] * pdf_d1 * vol % 2 * sqrt T) + r * K * cdf_neg_d2 * exp neg r*T;
	  rho: 0.01 * neg[K] * T * cdf_neg_d2 * exp neg r*T;
	  payoff:(K * cdf_neg_d2  * exp neg r*T) - S0 * cdf_neg_d1 * exp neg dividend * T];

	// Call option payoff
	[opt:`call;delta: cdf_d1;	
	  theta: (neg[S0] * pdf_d1 * vol % 2 * sqrt T) - r * K * cdf_d2 * exp neg r*T;
	  rho: 0.01 * K * T * cdf_d2 * exp neg r*T ;
	  payoff:(S0 * cdf_d1 * exp neg dividend * T) - K * cdf_d2 * exp neg r*T]
    ]; 
	0N!gamma;
	`.bs.greeks insert (.z.p;opt; S0; K; r; T; vol; delta; gamma; theta; rho; vega; vanna; vomma; zomma);
	:payoff
  }

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
	
// GREEKS

/
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
S0:s0:100f;mu:0.1;vol:0.2;T:1f;timestep:252f;N:10
K:90f;r:0.5%100;dividend:0
.bs.option[S0;K;r;T;vol;`opt]
.bs.greeks

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
\

/
Delta: Measures Impact of a Change in the Price of Underlying- the speed of option pps change.

Gamma: Measures the Rate of Change of Delta- the acceleration of option pps change.

Theta: Measures Impact of a Change in Time Remaining- option pps change per day closer to expiration.

Vega: Measures Impact of a Change in Volatility- option pps change per 1 value change in Volatility.

Rho: Measures Impact of a Change in Interest Rates- option pps change per 1 value change in Interest Rates.
\


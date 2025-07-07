// load required script
\l constant.q
\l impliedVol.q

.bs.tab:([] strike:(); option:(); impliedvol:());
.bs.greekstab:([] insertTime:`timestamp$();optionType:`$();payoff:`float$();stock:`float$();strike:`float$(); rfRate:`float$(); maturity:`float$(); volatility:`float$();
	delta:`float$(); gamma:`float$(); theta:`float$(); rho:`float$(); vega:`float$(); vanna:`float$(); vomma:`float$(); zomma:`float$());

.bs.call:{[S0;K;r;T;vol]

	d1: (1 % vol * sqrt T) * (log S0 % K) + T * r + 0.5 * vol * vol;
	d2: d1 - vol * sqrt T;
	(S0 * .const.normal_cdf[d1]) - K * .const.normal_cdf[d2] * exp neg r*T};

.bs.put:{[S0;K;r;T;vol]

	d1: (1 % vol * sqrt T) * (log S0 % K) + T * r  + 0.5 * vol * vol;
	d2: d1 - vol * sqrt T;
	
	(K * .const.normal_cdf[neg[d2]] * exp neg r*T) - S0 * .const.normal_cdf[neg[d1]] * exp neg T};

// cdf not matching with python
.bs.option:{[S0;K;r;T;vol;opt]
  // Validate input assumptions
  if[S0 < 0;'"stock price must be positive"];
  if[K <= 0; '"K must be > 0"];
  if[vol <= 0; '"volatility must be > 0"];
  if[T < 0; '"T must be >= 0"];
  if[not opt in `call`put; '"option type must be `call or `put"];

  d1: (1 % vol * sqrt T) * (log S0 % K) + T * r + 0.5 * vol * vol;
  d2: d1 - vol * sqrt T;

  pdf_d1: .const.normal_pdf[d1];
  cdf_d1: .const.normal_cdf[d1];
  cdf_neg_d1:.const.normal_cdf[neg[d1]];
  cdf_d2:.const.normal_cdf[d2];
  cdf_neg_d2:.const.normal_cdf[neg[d2]];

  gamma: pdf_d1 % S0 * vol * sqrt T;
  vega: S0 * pdf_d1 * sqrt T;

  vanna: 0.01 * vega * 1 - d1 % vol * sqrt T;
  vomma: vega * d1 * d2 % vol;
  zomma: 0.001 * vega * (-1 + d1 * d2)  % vol;

  $[opt in `put;
	// Put option payoff
	[delta:neg[cdf_neg_d1]; 
	  // annual theta
	  theta: (neg[S0] * pdf_d1 * vol % 2 * sqrt T) + r * K * cdf_neg_d2 * exp neg r*T;
	  rho: neg[K] * T * cdf_neg_d2 * exp neg r*T;
	  payoff:.bs.put[S0;K;r;T;vol]
	  ];

	// Call option payoff
	[delta: cdf_d1;	
	  // annual theta
	  theta: (neg[S0] * pdf_d1 * vol % 2 * sqrt T) - r * K * cdf_d2 * exp neg r*T;
	  rho: K * T * cdf_d2 * exp neg r*T ;
	  // (S0 * normal_cdf[d1] * exp neg dividend * T) - K * normal_cdf[d2] * exp neg r*T
	  payoff:.bs.call[S0;K;r;T;vol]
    ]]; 
	`.bs.greekstab insert (.z.p;opt;payoff; S0; K; r; T; vol; delta; gamma; theta; rho; vega; vanna; vomma; zomma);
	:payoff
  }

//default 
/ S0:100f;K:90f;T:0.5;r:0.05; vol: 0.2;opt:`call
// edge cases
// At-The-Money (ATM): S = K
S0:100f;K:100f;T:0.5;r:0.05; vol: 0.2;opt:`call
/ expected 
S0:100f;K:90f;T:0.5;r:0.05; vol: 0.2;opt:`call
.bs.option[S0;K;r;T;vol;opt]
// S=0
S0:0f;K:90f;T:0.5;r:0.05; vol: 0.2;opt:`call
.bs.option[S0;K;r;T;vol;opt]
// Very Short Time to Maturity (T ≈ 0)
S0:100f;K:90f;T:0.00001;r:0.05; vol: 0.2;opt:`call
.bs.option[S0;K;r;T;vol;opt]
// Deep In-The-Money Call (S >> K)
S0:1000f;K:90f;T:0.5;r:0.05; vol: 0.2;opt:`call
.bs.option[S0;K;r;T;vol;opt]
// Deep Out-of-The-Money Put (S >> K)
S0:100f;K:90f;T:0.5;r:0.05; vol: 0.2;opt:`put
.bs.option[S0;K;r;T;vol;opt]
// Zero Volatility (σ → 0)
S0:100f;K:90f;T:0.5;r:0.05; vol: 0.00001;opt:`call
.bs.option[S0;K;r;T;vol;opt]
// Zero Risk-Free Rate (r = 0)
S0:100f;K:90f;T:0.5;r:0f; vol: 0.2;opt:`call
.bs.option[S0;K;r;T;vol;opt]
// Negative Interest Rate (r < 0)
S0:100f;K:90f;T:0.5;r:-0.05; vol: 0.2;opt:`call
.bs.option[S0;K;r;T;vol;opt]
// Very High Volatility (σ >> 1)
S0:100f;K:90f;T:0.5;r:0.05; vol: 5f;opt:`call
.bs.option[S0;K;r;T;vol;opt]

/
// S=0
// Very Short Time to Maturity (T ≈ 0)
// Deep In-The-Money Call (S >> K)
// Deep Out-of-The-Money Put (S >> K)
// Zero Volatility (σ → 0)
// Zero Risk-Free Rate (r = 0)
// Negative Interest Rate (r < 0)
// Very High Volatility (σ >> 1)

Summary Table

| Scenario         | Test Case                      | Expected Call Price Behavior     | Delta   | Gamma   | Vega    | Theta   | Rho       |
|------------------|--------------------------------|----------------------------------|---------|---------|---------|---------|-----------|
| S = 0            | S0=0                           | 0                                | 0       | 0       | 0       | 0       | 0         |
| T ≈ 0            | T=0.00001                      | Intrinsic value                  | 1/0     | High/0  | 0       | High/0  | 0         |
| Deep ITM Call    | S0=1000, K=90                  | S0 - K e^{-rT}                   | 1       | 0       | 0       | 0       | Small +   |
| Deep OTM Put     | S0=1000, K=90, opt=put         | 0                                | 0       | 0       | 0       | 0       | Small -   |
| σ → 0            | vol=0.00001                    | Discounted intrinsic value       | 1/0     | 0       | 0       | 0       | 0         |
| r = 0            | r=0                            | Slightly higher than r>0         | ≈default| ≈default| ≈default| ≈default| 0         |
| r < 0            | r=-0.05                        | Higher than r=0                  | ≈default| ≈default| ≈default| ≈default| Unusual   |
| σ >> 1           | vol=5                          | Approaches S0                    | ~0.5    | High    | High    | High    | ≈default  |

\

/
// testing area
params
S0:s0:100f;mu:0.1;vol:0.2;T:1f;timestep:252f;N:10
K:90f;r:0.5;dividend:0
.bs.option[S0;K;r;T;vol;`opt]
.bs.greekstab
optprice: .bs.call[S0;K;r;T;vol] 
.bs.IVCall[optprice; s0; K; r; T;`;`]


// edge cases

.bs.IV:[call_price; S0; K; r; T; dividend; tol; max_iter]
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

// GREEKS
/
Delta: Measures Impact of a Change in the Price of Underlying- the speed of option pps change.

Gamma: Measures the Rate of Change of Delta- the acceleration of option pps change.

Theta: Measures Impact of a Change in Time Remaining- option pps change per day closer to expiration.

Vega: Measures Impact of a Change in Volatility- option pps change per 1 value change in Volatility.

Rho: Measures Impact of a Change in Interest Rates- option pps change per 1 value change in Interest Rates.
\


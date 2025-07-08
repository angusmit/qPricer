\l constant.q
\l bs.q 

/reference: https://www.codearmo.com/python-tutorial/options-trading-greeks-black-scholes
/
ffd:{[f;x;dx]
    // Forward Finite difference approximation of derivative
    // f: function to differentiate
    // x: point at which to evaluate the derivative
    // dx: small change in x
    (f[x+dx] - f[x]) % dx
 };

bfd:{[f;x;dx]
    // Backward finite difference approximation of derivative
    // f: function to differentiate
    // x: point at which to evaluate the derivative
    // dx: small change in x
    (f[x] - f[x-dx]) % dx
 };

cfd:{[f;x;dx]
    // Central finite difference approximation of derivative
    // f: function to differentiate
    // x: point at which to evaluate the derivative
    // dx: small change in x
    (f[x+dx] - f[x-dx]) % 2*dx
 };

.const.linspace:{[start;end;n] step:(1%n) *end-start; start+step* til n+1};
.const.arange:{[start;end;n] add:n+; end-:n; add\[end>;start]};

dxs:1.5 1 0.5 0.2 0.1; 
nums:.const.arange[1;5;0.01];
//nums: .const.linspace[1;5;8];

ffd_l: last each ffd[f;;]/:[nums;dxs];
bfd_l: last each bfd[f;;]/:[nums;dxs];   
cfd_l: last each cfd[f;;]/:[nums;dxs];
\


/ S0:100f;K:100f;T:1.0;r:0.00; vol: 0.25;opt:`call;dx:`; method:`forward;
/ .fdm.delta[S0;K;r;T;vol;dx;method]
/ .fdm.delta[S0;K;r;T;vol;dx;`central]
/ .fdm.delta[S0;K;r;T;vol;dx;`backward]
.fdm.delta: {[S0;K;r;T;vol;dx;opt;method]
  if[not opt in `call`put; '"option type must be `call or `put"];
  if[not method in `forward`backward`central; '"method must be `forward, `backward or `central"];
  // Validate input assumptions
  if[S0 < 0; '"stock price must be positive"];
  if[K <= 0; '"K must be > 0"];
  if[vol <= 0; '"volatility must be > 0"];  
  if[T < 0; '"T must be >= 0"];
  // Select the appropriate payoff function based on the option type
  ?[opt=`call; payoffFunc:`.bs.call; payoffFunc:`.bs.put];
  ?[dx~`; ds:10 xexp -5;ds:dx];
  // Forward difference method for delta
  if[method=`forward; 
    //res: (.bs.call[S0+ds;K;r;T;vol] - .bs.call[S0;K;r;T;vol]) % ds];
    res: (sum @[payoffFunc ./: ((S0+ds;K;r;T;vol);(S0;K;r;T;vol));1;-1*] ) % ds];
  // Backward difference method for delta
  if[method=`backward; 
    //res:(.bs.call[S0;K;r;T;vol] - .bs.call[S0 - ds;K;r;T;vol]) % ds];
    res:(sum @[payoffFunc ./: ((S0;K;r;T;vol);(S0-ds;K;r;T;vol));1;-1*] ) % ds];
    
  // Central difference method for delta
  if[method=`central; 
    //res:(.bs.call[S0 + ds;K;r;T;vol] - .bs.call[S0 - ds;K;r;T;vol]) % 2 * ds];
    res:(sum @[payoffFunc ./: ((S0+ds;K;r;T;vol);(S0-ds;K;r;T;vol));1;-1*] ) % 2 * ds];
    :res

  };

.gk.delta:{[S0;K;r;T;vol;opt]
  d1: (1 % vol * sqrt T) * (log S0 % K) + T * r + 0.5 * vol * vol;
  d2: d1 - vol * sqrt T;
  cdf_d1: .const.normal_cdf[d1];
  cdf_neg_d1:.const.normal_cdf[neg[d1]];
  ?[opt=`call; 
    // Call, put analytical delta 
    res:cdf_d1; res:neg[cdf_neg_d1];
  ];
  :res
 }

/ delta testing
/ S0:100f;K:100f;T:1.0;r:0.00; vol: 0.25;opt:`call;dx:`; method:`forward;
//.bs.option[S0;K;r;T;vol;`call]
/ s:.const.linspace[0;250;10]
/ true:.fdm.delta[;K;r;T;vol;dx;`call;`backward]'[s]
/ test:.gk.delta[s;K;r;T;vol;`call]
/ .fdm.delta[S0;K;r;T;vol;dx;`call;`backward]
/ .fdm.delta[S0;K;r;T;vol;dx;`call;`forward]

/ GAMMA
.gk.gamma:{[S0;K;r;T;vol]
  d1: (1 % vol * sqrt T) * (log S0 % K) + T * r + 0.5 * vol * vol;
  d2: d1 - vol * sqrt T;
  pdf_d1: .const.normal_pdf[d1];
  :pdf_d1 % S0 * vol * sqrt T
  }

/ GAMMA -> Central finite difference method
.fdm.gamma: {[S0;K;r;T;vol;dx;opt;method]
  if[not opt in `call`put; '"option type must be `call or `put"];
  if[not method in `forward`backward`central; '"method must be `forward, `backward or `central"];
  // Validate input assumptions
  if[S0 < 0; '"stock price must be positive"];
  if[K <= 0; '"K must be > 0"];
  if[vol <= 0; '"volatility must be > 0"];  
  if[T < 0; '"T must be >= 0"];
  // Select the appropriate payoff function based on the option type
  ?[opt=`call; payoffFunc:`.bs.call; payoffFunc:`.bs.put];  
  ?[dx~`; ds:10 xexp -5;ds:dx];
  
  // Forward difference method for gamma
  if[method=`forward;
    res: (sum @[payoffFunc ./: ((S0+ds;K;r;T;vol);(S0;K;r;T;vol);(S0-ds;K;r;T;vol));1;-2*])%ds*ds];
    
  // Backward difference method for gamma
  if[method=`backward;
    res: (sum @[payoffFunc ./: ((S0+2*ds;K;r;T;vol);(S0+ds;K;r;T;vol);(S0;K;r;T;vol));1;-2*])%ds*ds];
    
  // Central difference method for gamma
  if[method=`central;
    res: (sum @[payoffFunc ./: ((S0;K;r;T;vol);(S0-ds;K;r;T;vol);(S0-2*ds;K;r;T;vol));1;-2*])%ds*ds];
    
    :res
  };

/ gamma testing
/ S0:100f;K:100f;T:1.0;r:0.00; vol: 0.25;opt:`call;dx:`; method:`forward;
/ s:.const.linspace[0;250;1000]
/ true:.fdm.gamma[;K;r;T;vol;dx;`call;`backward]'[s]
/ test:.gk.gamma[s;K;r;T;vol]
/ .fdm.gamma[S0;K;r;T;vol;dx;`call;`backward]
/ .fdm.gamma[S0;K;r;T;vol;dx;`call;`forward]

/ VEGA
.gk.vega:{[S0;K;r;T;vol]
  d1: (1 % vol * sqrt T) * (log S0 % K) + T * r + 0.5 * vol * vol;
  pdf_d1: .const.normal_pdf[d1];
  :S0 * pdf_d1 * sqrt T
  }

.fdm.vega: {[S0;K;r;T;vol;dx;opt;method]
  if[not opt in `call`put; '"option type must be `call or `put"];
  if[not method in `forward`backward`central; '"method must be `forward, `backward or `central"];
  // Validate input assumptions
  if[S0 < 0; '"stock price must be positive"];
  if[K <= 0; '"K must be > 0"];
  if[vol <= 0; '"volatility must be > 0"];  
  if[T < 0; '"T must be >= 0"];
  // Select the appropriate payoff function based on the option type
  ?[opt=`call; payoffFunc:`.bs.call; payoffFunc:`.bs.put];  
  ?[dx~`; dv:10 xexp -4;dv:dx];
  
  // Forward difference method for vega
  if[method=`forward;
    res: (sum payoffFunc ./: ((S0;K;r;T;vol+dv);(S0;K;r;T;vol)))%dv];
    
  // Backward difference method for vega
  if[method=`backward;
    res: (sum payoffFunc ./: ((S0;K;r;T;vol);(S0;K;r;T;vol-dv)))%dv];
    
  // Central difference method for vega
  if[method=`central;
    res: (sum payoffFunc ./: ((S0;K;r;T;vol+dv);(S0;K;r;T;vol-dv)))%2*dv];
    
    :res
  };

// theta
.gk.theta:{[S0;K;r;T;vol;opt]
  d1: (1 % vol * sqrt T) * (log S0 % K) + T * r + 0.5 * vol * vol;
  d2: d1 - vol * sqrt T;
  pdf_d1: .const.normal_pdf[d1];
  cdf_d1: .const.normal_cdf[d1];
  cdf_neg_d1:.const.normal_cdf[neg[d1]];
  cdf_d2:.const.normal_cdf[d2];
  cdf_neg_d2:.const.normal_cdf[neg[d2]];
  
  ?[opt=`call;
    // Call, put analytical theta
    res:(neg[S0] * pdf_d1 * vol % 2 * sqrt T) - r * K * cdf_d2 * exp neg r*T;
    res: (neg[S0] * pdf_d1 * vol % 2 * sqrt T) + r * K * cdf_neg_d2 * exp neg r*T;
  ];
  :res
  };

.fdm.theta: {[S0;K;r;T;vol;dx;opt;method]
  if[not opt in `call`put; '"option type must be `call or `put"];
  if[not method in `forward`backward`central; '"method must be `forward, `backward or `central"];
  // Validate input assumptions
  if[S0 < 0; '"stock price must be positive"];
  if[K <= 0; '"K must be > 0"];
  if[vol <= 0; '"volatility must be > 0"];  
  if[T < 0; '"T must be >= 0"];
  // Select the appropriate payoff function based on the option type
  ?[opt=`call; payoffFunc:`.bs.call; payoffFunc:`.bs.put];  
  ?[dx~`; dt:0.00001;dt:dx];
  
  // Forward difference method for theta
  if[method=`forward;
    
    res: (sum @[payoffFunc ./: ((S0;K;r;T+dt;vol);(S0;K;r;T;vol));1;-1*])%dt];
    
  // Backward difference method for theta
  if[method=`backward;
    
    res: (sum @[payoffFunc ./: ((S0;K;r;T;vol);(S0;K;r;T-dt;vol));1;-1*])%dt];
    
  // Central difference method for theta
  if[method=`central;
    res: (sum @[payoffFunc ./: ((S0;K;r;T+dt;vol);(S0;K;r;T-dt;vol));1;-1*])%2*dt];
    
    :neg[res]
  };

/ theta testing
/ S0:100f;K:100f;T:1.0;r:0.00; vol: 0.25;opt:`call;dx:0.0001; method:`forward;
/ t:.const.linspace[1;0.05;6]
/ true:.fdm.theta[S0;K;r;;vol;dx;`call;`backward]'[t]
/ test:.gk.theta[S0;K;r;t;vol;opt]

/rho
.gk.rho:{[S0;K;r;T;vol;opt]
  d1: (1 % vol * sqrt T) * (log S0 % K) + T * r + 0.5 * vol * vol;
  d2: d1 - vol * sqrt T;
  cdf_d2:.const.normal_cdf[d2];
  cdf_neg_d2:.const.normal_cdf[neg[d2]];
  
  ?[opt=`call;
    // Call, put analytical rho
    res: K * T * cdf_d2 * exp neg r*T;
    res: neg[K] * T * cdf_neg_d2 * exp neg r*T;
  ];
  :res
  };

  .fdm.rho: {[S0;K;r;T;vol;dx;opt;method]
  if[not opt in `call`put; '"option type must be `call or `put"];
  if[not method in `forward`backward`central; '"method must be `forward, `backward or `central"];
  // Validate input assumptions
  if[S0 < 0; '"stock price must be positive"];
  if[K <= 0; '"K must be > 0"];
  if[vol <= 0; '"volatility must be > 0"];  
  if[T < 0; '"T must be >= 0"];
  
  // Select the appropriate payoff function based on the option type 
  ?[opt=`call; payoffFunc:`.bs.call; payoffFunc:`.bs.put];
  ?[dx~`; dr:10 xexp -5;dr:dx];
  // Forward difference method for rho
  if[method=`forward; 
    res: (sum @[payoffFunc ./: ((S0;K;r+dr;T;vol);(S0;K;r;T;vol));1;-1*] ) % dr];
  // Backward difference method for rho
  if[method=`backward; 
    res:(sum @[payoffFunc ./: ((S0;K;r;T;vol);(S0;K;r-dr;T;vol));1;-1*] ) % dr];
  // Central difference method for rho
  if[method=`central; 
    res:(sum @[payoffFunc ./: ((S0;K;r+dr;T;vol);(S0;K;r-dr;T;vol));1;-1*] ) % 2 * dr];
    :res

  };

/ rho testing
/ S0:100f;K:100f;T:1.0;r:0.00; vol: 0.25;opt:`call;dx:0.0001; method:`forward;
/ rate:.const.linspace[0;1;10]
/ true:.fdm.rho[S0;K;;T;vol;dx;`call;`backward]'[rate]
/ test:.gk.rho[S0;K;rate;T;vol;opt]
/maths setting
.const.pi:acos -1;
.const.linspace:{[s;e;n] step:(1%n) *e-s; s+step* til n+1};

/normal rv N~[0,1]
.const.nrv: {$[x=2*n:x div 2;raze sqrt[-2*log n?1f]*/:(sin;cos)@\:(2*.const.pi)*n?1f;-1_.z.s 1+x]};

.const.normal_pdf:{ (1 % sqrt 2 * .const.pi) * exp neg[0.5 * x*x]};

/ old normal cdf
//.const.normal_cdf:{abs(x>0)-(exp[-.5*x*x]%sqrt 2*.const.pi)*t*.31938153+t*-.356563782+t*1.781477937+t*-1.821255978+1.330274429*t:1%1+.2316419*abs x};

/ normal cdf
// polynomial approximation from Abramowitz and Stegun, ensuring high accuracy (within approximately 7.5e-8 error) for all real-valued inputs
.const.normal_cdf:{[x]
    xa: abs x;
    t: 1 % 1 + 0.2316419*xa;
    poly: t * (0.31938153 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));
    normalden: (exp neg (xa*xa) % 2) % sqrt 2 * .const.pi;
    comp: normalden * poly;
    ?[x >= 0; 1 - comp; comp]
  };

// generate n normal variables with mean m, standard deviation sd
// reference: https://armantee.github.io/sampling-with-kdb-p1/
.const.norm: {[n;m;sd]
    u1: n?1f;
    u2: n?1f;
    m + sd * sqrt[-2*log u1] * cos 2*u2*.const.pi
  };

// normal variable with mu = 0, sd = 1;
.const.norm01:.const.norm[;0;1];

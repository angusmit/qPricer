// load required script
\l constant.q
\l bs.q

.iv.CallBisection:{[call_price; S0; K; r; T; tol; max_iter]
	sigma_res:`;
	sigma_min:0.001;
	sigma_max:2;
	if[tol~`; tol: 10 xexp -6];
	if[max_iter~`; max_iter:1000];
	do[max_iter; sigma_mid:0.5 * sigma_min + sigma_max;
		price_mid: .bs.call[S0;K;r;T;sigma_mid];
			$[tol>abs(price_mid - call_price); sigma_res:sigma_mid; $[price_mid > call_price; sigma_max:sigma_mid; sigma_min:sigma_mid]]];
	if[sigma_res~`;sigma_res:sigma_mid];sigma_res}


// testing area
/
params
S0:s0:100f;mu:0.1;vol:0.2;T:1f;timestep:252f;N:10
K:90f;r:0.5;dividend:0
.bs.option[S0;K;r;T;vol;`opt]
.bs.greekstab
optprice: .bs.call[S0;K;r;T;vol] 
.iv.CallBisection[optprice; s0; K; r; T;`;`]
\
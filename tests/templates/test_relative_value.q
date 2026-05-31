\l core/init.q
/ ============================================================================
/ test_relative_value.q - the relativeValue template (new shape). Research OS R6.
/ (1) known-answer PnL: a spread that spikes then reverts -> the fade is short at the
/     spike and captures the reversion; frictionless cost is zero. (2) the spread-
/     stationarity gate PASSES a mean-reverting series and FAILS a random walk (it bites).
/ (3) end-to-end .template.rv.run on a synthetic curve history returns pnl/validation/meta.
/ Synthetic only, no HDB.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
apx:{[a;b;tol] tol>abs a-b};

/ --- (1) known-answer PnL on a constructed spread (frictionless) ---
cfg:.cfg.templates.rv,`lookback`entryZ`notional`txnCostRate!(3;1.0;1f;0f);
spread:100 100 100 110 100 100 100f;               / a spike at index 3, then it reverts
d:.template.rv.__signalPnl[spread; 2020.01.01+til 7; cfg];
/ at the spike the trailing z-score = (110-103.333)/4.714 = 1.414 > 1 -> fade SHORT (pos=-1).
chk[(-1f)=d[`pos] 3; "the fade must be SHORT (-1) at the spike (z>entryZ)"];
chk[0f=d[`pos] 2; "no position before the spike crosses the band"];
/ next day the spread reverts 110->100 (dSpread=-10); short position -> gross = (-1)*(-10) = +10.
chk[apx[10f; d[`gross] 4; 1e-9]; "fading the spike captures +10 when it reverts"];
chk[apx[10f; sum d`pnl; 1e-9]; "frictionless total PnL = gross = +10"];
chk[all 0f=d`cost; "txnCostRate=0 -> zero cost via the fill model (net == gross)"];

/ --- a small txnCostRate makes cost positive on the position changes (routes through .exec.fill) ---
cfgC:.cfg.templates.rv,`lookback`entryZ`notional`txnCostRate!(3;1.0;1f;0.01);
dc:.template.rv.__signalPnl[spread; 2020.01.01+til 7; cfgC];
chk[(sum dc`cost)>0f; "a positive txnCostRate -> positive cost"];
chk[(sum dc`pnl)<sum d`pnl; "net of cost is below the frictionless PnL"];

/ --- (2) the stationarity gate bites ---
ar:{[n] e:n#0.2 -0.3 0.1 -0.15 0.25 -0.2; r:0f; out:(); i:0; while[i<n; r:(0.5*r)+e i; out,:r; i+:1]; out}[200];
rw:sums 200#0.1 -0.1 0.2 -0.2 0.15 -0.05;
sAr:.template.rv.stationarity ar;
sRw:.template.rv.stationarity rw;
chk[sAr`stationary; "a mean-reverting spread PASSES the stationarity gate"];
chk[not sRw`stationary; "a random-walk spread FAILS the stationarity gate (AR1 ~ 1)"];
chk[(sAr`ar1)<sRw`ar1; "the mean-reverting AR(1) coef is below the random walk's"];
chk[(sRw`ar1)>0.95; "the random walk's AR(1) coef is near 1"];

/ --- (3) end-to-end run on a synthetic curve history ---
dates:2020.01.01+til 6;
ch:([] asofDate:dates,dates;
    tenor:(6#0.0833),6#0.25;
    price:(100 101 99 110 100 100f),98 99 97 100 98 98f;
    contractYM:(12)#202003;
    expiry:(12)#2020.03.20);
out:.template.rv.run[`commodity`curveHistory`overrides!(`CRUDE;ch;`lookback`entryZ`txnCostRate!(3;1.0;0f))];
chk[(`pnl`validation`meta)~key out; "run must return pnl/validation/meta"];
chk[6=count out`pnl; "one PnL row per curve date"];
chk[(out[`meta;`shape])~`relativeValue; "meta records the relativeValue shape"];
chk[(out[`meta;`edgeSource])~`structural; "relativeValue claims a STRUCTURAL edge"];
chk[`stationarity in key out`validation; "validation carries the template-specific stationarity gate"];

-1 "test_relative_value: PASS";

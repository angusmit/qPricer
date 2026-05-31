\l core/init.q
/ ============================================================================
/ test_factor_pca.q - the factor capability (.factor.pca). Research OS R8.
/ A SYNTHETIC panel from KNOWN loadings + scores -> PCA recovers them within tolerance,
/ sign-fixed (front loading positive), and IDENTICALLY across two runs (deterministic).
/ Explained variance sums sensibly. Pure, no HDB.
/ ============================================================================
chk:{[b;msg] if[not b; '"FAIL: ",msg]};
apx:{[a;b;tol] tol>abs a-b};

cfg:.cfg.factor;
Llevel:.factor.__norm 1 1 1 1f;            / level: flat across maturities
Lslope:.factor.__norm -3 -1 1 3f;          / slope: orthogonal to level
nObs:120;
s1:0.02*sin 0.3*til nObs;                   / level scores (dominant)
s2:0.01*cos 0.5*til nObs;                   / slope scores (smaller)
X:(s1 *\: Llevel) + s2 *\: Lslope;          / T x M, an EXACT 2-factor panel

d:.factor.pca[X;2;cfg];

/ --- recovery (sign-fixed): PC1 ~ level, PC2 ~ slope, |cosine| ~ 1 ---
chk[apx[1f; abs Llevel wsum (d`loadings)[;0]; 1e-3]; "PC1 must recover the level loading"];
chk[apx[1f; abs Lslope wsum (d`loadings)[;1]; 1e-3]; "PC2 must recover the slope loading"];

/ --- the sign convention is enforced: front loading positive on each PC ---
chk[0f<(d`loadings)[0;0]; "PC1 front loading must be positive (sign convention)"];
chk[0f<(d`loadings)[0;1]; "PC2 front loading must be positive (sign convention)"];

/ --- explained variance: a 2-factor panel -> 2 PCs sum to ~1, ordered (PC1 >= PC2) ---
chk[apx[1f; sum d`explainedVar; 1e-6]; "explained variance of 2 PCs on a 2-factor panel sums to ~1"];
chk[(d`explainedVar)[0] >= (d`explainedVar)[1]; "PCs are ordered by explained variance"];
chk[all 0f<=d`explainedVar; "explained variance is non-negative"];

/ --- residuals ~ 0 for an exact 2-factor panel ---
chk[1e-6 > max max abs d`residuals; "the 2-factor reconstruction leaves ~0 residual"];

/ --- DETERMINISM: two runs are byte-identical (fixed init / tol / sign) ---
chk[(d`loadings)~(.factor.pca[X;2;cfg])`loadings; "loadings are deterministic across runs"];
chk[(d`scores)~(.factor.pca[X;2;cfg])`scores; "scores are deterministic across runs"];

/ --- shape sanity ---
chk[(4 2)~(count d`loadings),count first d`loadings; "loadings are M x k"];
chk[(nObs;2)~(count d`scores;count first d`scores); "scores are T x k"];

/ --- the capability conforms to its R2 contract ---
chk[.factor.conforms `curvePCA; "the registered curvePCA capability conforms to the factor contract"];

-1 "test_factor_pca: PASS";

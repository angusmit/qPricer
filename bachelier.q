\l constant.q
\l bs.q
// Arithmetic Brownian motion, or Bachelier Model
// dS=μSdt+σe**(μ(t−T))dW
/ S B t := S0(1 + σWt),
/ https://onlinelibrary.wiley.com/doi/epdf/10.1002/fut.22315
/ https://quant.stackexchange.com/questions/32863/bachelier-model-call-option-pricing-formula
/ https://veeenu.github.io/firm/QuantFinAndDer.pdf
/  St; r;sigma; T; t;
/ call/put: https://arxiv.org/pdf/2104.08686
/ Q-price https://download.ssrn.com/19/09/23/ssrn_id3458344_code2794557.pdf?response-content-disposition=inline&X-Amz-Security-Token=IQoJb3JpZ2luX2VjEKf%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJIMEYCIQCHJlrP3rTaYzxHjm6UPcpXYaRabVjTJ8g3gCUaGMYAmQIhANwZCghFvB9WtIM3Mta3XU%2FvWKAtUCXdgsDHjHDY7a7QKsYFCLD%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQBBoMMzA4NDc1MzAxMjU3IgyersoaStfdULcFwqkqmgX%2F%2Fcfm74WB19HudRzuCzabtUFPsoFLg8fWIxuGV1Tudti7ZQTlQXzXPXMlZN1y5cmkXWaMGNBKjGaZqIatDhtR6dOcCERpP4FAMHATLhWXSy%2FLmR4lkZQMF3uppLhaxPqMkCqStAREFjAiKo0%2BVfNxbYOhZUaSZs93O%2FgsYd5E4atDJCuDFRzRBfNtXW5srNsG0ZmECpyKviSydO2cQRzHTXHWmBhPEPn%2BPlJTb4o8I5TH7hzWfl0nfMuv1ykp28ayeXxIF%2ByXmx7aYrFcMdBwWgO6pewEmzBxlhDadHF3fcvllKFVfCHZBP9TihKHSFLuG4YZHfaHTTJ%2BNPQMRUTICDsCP8ovpxCFrb21UcKzbH8Wst5wANAr%2FTZdYWA%2F39SKU1O%2BnJn%2FhMsjST0cwPqxguZqZDiGjy7%2FAJ9UQZ9PmLdNzBUFnfYp%2BBA2uk1kuVkl7T5wEA6rIjq1vc9OXyZo%2Fi2GP5dpb%2FGJ%2B3ejxckw1%2FJUSQcRh6KMrlkCpsY5cQYbOrln1yu%2BqkrIFTrdlTMU6hly9ouKEm1LreE2E7Q4CK5%2BSDaelHDKSJx5A3OVGmO6uxTh88OW3CnPZ7Wug4Y1TjnT3ZVj3NQDZI5TYQEg4P%2BSElHmIyh%2FjZ3JKMMYmuc%2BPpwibIBHYsBQfe%2FUjRWm0%2B8BDaFUrFe%2FgJ2ZGP24NhAJuz6lKYdk7CaZE%2BU8B5rxc9XpwyKBSz0Ajqv3YvTlCKLC6FfyzvuXULk6NMU%2B4dcpL7Kt4GnlM7SGnigzNtGtfH8Gh%2BwEfXDql2hyVsTU9FjaOBAevLoq3yyr98NcP%2F%2FGPEVRKUocVk2rDYeKuF8hydZa4jXr%2FD0pQ6SJqZzhDy%2BVzBQKXcA9fD3x1NNmTOQD6MGTJmFGO%2F8wxdm7wwY6sAG4EN4Xb1dIaa2KtNhWLrmSmH8SvEr3raZRlb56wNWtySKZ5y3KGI18LxzFTToQDkkLv7j1ARgALEjOkCDg5d2AJhmFdrkvPiw7aA2GdrWm%2BAOMUx61ou4nNn1cMi3NbzE%2BlydHWT4VCi9aAAIKfQs4j%2BoUfffXBW%2BryMO7Oz8wSffdnp%2BVOl3mQVh4%2Fm00Mm4HDs4vfyRVz9kKrM1p2%2Bw37ykHumKwDwL1SpYOhOyd5Q%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20250709T223514Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAUPUUPRWE4G6LNEMC%2F20250709%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=6e975ca7dc46ac8dd685b44aadff315949a20e001776059e450bc74e2a03d3d2&abstractId=3428994

/ risk-neutral  
.bachelier.call:{[S;K;r;sigma; T; t]
    // Validate input assumptions
    if[sigma <= 0; '"volatility must be > 0"];
    if[T < t; '"T must be >= t"];
    if[t< 0; '"t must be >= 0"];

    sqrtdt: sqrt T - t;
    $[r=0;
        [
            D: (S - K) % sigma * sqrtdt; 
            // undiscounted call option payoff when r = 0
            res:((S - K) * .const.normal_cdf[D]) + sigma * sqrtdt * .const.normal_pdf[D];
        ];
        [
            Bt: (exp neg r * (T - t));
            D: (S - K * Bt ) % sigma * ((1 - Bt * Bt) % 2 * r);
            // undiscounted call option payoff when r != 0
            res:(S - K * Bt) * .const.normal_cdf[D] + sigma * ((1 - Bt * Bt) % 2 * r) * .const.normal_pdf[D];
        ]
    ];
    :res
 }
/ S:100f;K:90f;T:0.5;t:0.00001;r:0.05; sigma: 0.2; opt:`call
/ S:140f;K:90f;T:0.5;t:0.00001;r:0.05; sigma: 0.2; opt:`call
/ S:90+0.000000001;K:90f;T:1.0;t:0.00001;r:0.05; sigma: 0.2; opt:`call
/ / s: .const.linspace[-200;200;1000]
/ s: 20 * 1 + til 10
/ .bs.call[S;K;r;T;sigma]
/ .bs.call[s;K;r;T;sigma]
/ .bachelier.call[s;K;r;sigma; T; t]
/ .bachelier.call[S;K;r;sigma; T; t]

/ verify when F0 = K, C_n[K] = sigma * sqrt T % 2 * .const.pi
/ eq2 - https://arxiv.org/pdf/2104.08686 

/ .bs.call[90+0.000000001;K;r;T;sigma%S] ~ .bachelier.call[K+0.000001;K;r;sigma; T; t] ~ ( sigma % S ) * 0.4 * S * sqrt T
/ 
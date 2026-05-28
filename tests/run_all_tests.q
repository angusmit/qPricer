/ run_all_tests.q - grouped test runner (v0.14)
\l lib/init.q

/ --- Test suites by domain ---

.test.coreFiles:(
    "tests/core/test_european_call.q";
    "tests/core/test_european_put.q";
    "tests/core/test_put_call_parity.q";
    "tests/core/test_grid_convergence.q";
    "tests/core/test_input_validation.q";
    "tests/core/test_scenario_risk.q";
    "tests/core/test_crank_nicolson_call.q";
    "tests/core/test_crank_nicolson_put.q";
    "tests/core/test_crank_nicolson_vs_explicit.q");

.test.greeksFiles:(
    "tests/greeks/test_greeks_call.q";
    "tests/greeks/test_greeks_put.q";
    "tests/greeks/test_greeks_extended_products.q");

.test.americanFiles:(
    "tests/american/test_american_put.q";
    "tests/american/test_american_call.q";
    "tests/american/test_american_call_zero_dividend.q";
    "tests/american/test_american_call_dividend.q";
    "tests/american/test_early_exercise_boundary.q");

.test.barrierFiles:(
    "tests/barrier/test_barrier_up_and_out_call.q";
    "tests/barrier/test_barrier_down_and_out_put.q";
    "tests/barrier/test_barrier_validation.q";
    "tests/barrier/test_barrier_up_and_out_put.q";
    "tests/barrier/test_barrier_down_and_out_call.q";
    "tests/barrier/test_barrier_up_and_in_call.q";
    "tests/barrier/test_barrier_up_and_in_put.q";
    "tests/barrier/test_barrier_down_and_in_call.q";
    "tests/barrier/test_barrier_down_and_in_put.q";
    "tests/barrier/test_barrier_parity.q";
    "tests/barrier/test_barrier_full_validation.q");

.test.localVolFiles:(
    "tests/localvol/test_local_vol_flat_equivalence_call.q";
    "tests/localvol/test_local_vol_flat_equivalence_put.q";
    "tests/localvol/test_local_vol_validation.q";
    "tests/localvol/test_local_vol_skew_sanity.q";
    "tests/localvol/test_local_vol_time_dependence.q";
    "tests/localvol/test_local_vol_american_put.q";
    "tests/localvol/test_local_vol_american_call.q";
    "tests/localvol/test_local_vol_american_flat_equivalence_put.q";
    "tests/localvol/test_local_vol_american_flat_equivalence_call.q";
    "tests/localvol/test_local_vol_american_validation.q";
    "tests/localvol/test_local_vol_barrier_flat_equivalence.q";
    "tests/localvol/test_local_vol_barrier_up_and_out_call.q";
    "tests/localvol/test_local_vol_barrier_down_and_out_put.q";
    "tests/localvol/test_local_vol_barrier_knock_in_parity.q";
    "tests/localvol/test_local_vol_barrier_validation.q");

.test.marketFiles:(
    "tests/market/test_market_data_book.q";
    "tests/market/test_market_data_book_missing_symbol.q";
    "tests/market/test_vol_surface_lookup.q";
    "tests/market/test_vol_surface_validation.q";
    "tests/market/test_vol_surface_nearest_lookup.q";
    "tests/market/test_surface_market_data_pricing.q");

.test.portfolioFiles:(
    "tests/portfolio/test_portfolio_pricing.q";
    "tests/portfolio/test_portfolio_greeks.q";
    "tests/portfolio/test_portfolio_scenario_risk.q";
    "tests/portfolio/test_portfolio_multi_symbol_pricing.q";
    "tests/portfolio/test_portfolio_multi_symbol_scenario_risk.q";
    "tests/portfolio/test_portfolio_extended_products.q";
    "tests/portfolio/test_scenario_extended_products.q";
    "tests/portfolio/test_portfolio_monte_carlo_products.q";
    "tests/portfolio/test_portfolio_basket_products.q";
    "tests/portfolio/test_portfolio_lookback_products.q";
    "tests/portfolio/test_portfolio_heston_products.q";
    "tests/portfolio/test_calibrated_model_pricing.q";
    "tests/portfolio/test_model_comparison_pricing.q";
    "tests/portfolio/test_sabr_surface_pricing.q";
    "tests/portfolio/test_portfolio_merton_products.q";
    "tests/portfolio/test_portfolio_bates_products.q";
    "tests/portfolio/test_portfolio_var.q";
    "tests/portfolio/test_portfolio_historical_replay.q";
    "tests/portfolio/test_portfolio_limit_monitoring.q");

.test.impliedVolFiles:(
    "tests/impliedvol/test_implied_vol_call.q";
    "tests/impliedvol/test_implied_vol_put.q";
    "tests/impliedvol/test_option_chain_implied_vols.q";
    "tests/impliedvol/test_implied_vol_invalid_prices.q";
    "tests/impliedvol/test_option_chain_partial_failures.q");

.test.reportingFiles:(
    "tests/reporting/test_report.q";
    "tests/reporting/test_pnl_explain.q";
    "tests/reporting/test_batch.q";
    "tests/reporting/test_audit.q";
    "tests/reporting/test_regression.q";
    "tests/reporting/test_dashboard.q");

.test.infraFiles:(
    "tests/infra/test_result.q";
    "tests/infra/test_config.q";
    "tests/infra/test_timing.q";
    "tests/infra/test_testutil.q";
    "tests/infra/test_cache.q");

.test.stressFiles:(
    "tests/stress/test_stress.q";
    "tests/stress/test_perfdiag.q";
    "tests/stress/test_perfopt.q";
    "tests/stress/test_mc_diagnostics.q";
    "tests/stress/test_model_limit_stress.q");

.test.monteCarloFiles:(
    "tests/montecarlo/test_monte_carlo_paths.q";
    "tests/montecarlo/test_monte_carlo_european.q";
    "tests/montecarlo/test_asian_arithmetic_call.q";
    "tests/montecarlo/test_asian_arithmetic_put.q";
    "tests/montecarlo/test_asian_geometric_validation.q";
    "tests/montecarlo/test_monte_carlo_greeks.q";
    "tests/montecarlo/test_correlation_matrix.q";
    "tests/montecarlo/test_correlated_paths.q";
    "tests/montecarlo/test_basket_option_call.q";
    "tests/montecarlo/test_basket_option_put.q";
    "tests/montecarlo/test_basket_option_validation.q";
    "tests/montecarlo/test_basket_greeks.q";
    "tests/montecarlo/test_path_diagnostics.q";
    "tests/montecarlo/test_lookback_fixed_call.q";
    "tests/montecarlo/test_lookback_fixed_put.q";
    "tests/montecarlo/test_lookback_floating_call.q";
    "tests/montecarlo/test_lookback_floating_put.q";
    "tests/montecarlo/test_lookback_validation.q";
    "tests/montecarlo/test_lookback_greeks.q";
    "tests/montecarlo/test_antithetic_variance_reduction.q";
    "tests/montecarlo/test_moment_matching.q";
    "tests/montecarlo/test_control_variate_european.q";
    "tests/montecarlo/test_control_variate_geometric_asian.q";
    "tests/montecarlo/test_mc_convergence_table.q";
    "tests/montecarlo/test_mc_confidence_intervals.q";
    "tests/montecarlo/test_heston_paths.q";
    "tests/montecarlo/test_heston_european_call.q";
    "tests/montecarlo/test_heston_european_put.q";
    "tests/montecarlo/test_heston_black_scholes_limit.q";
    "tests/montecarlo/test_heston_validation.q";
    "tests/montecarlo/test_heston_greeks.q";
    "tests/montecarlo/test_merton_jump_paths.q";
    "tests/montecarlo/test_merton_european_call.q";
    "tests/montecarlo/test_merton_european_put.q";
    "tests/montecarlo/test_merton_black_scholes_limit.q";
    "tests/montecarlo/test_merton_validation.q";
    "tests/montecarlo/test_merton_greeks.q";
    "tests/montecarlo/test_bates_paths.q";
    "tests/montecarlo/test_bates_european_call.q";
    "tests/montecarlo/test_bates_european_put.q";
    "tests/montecarlo/test_bates_limit_cases.q";
    "tests/montecarlo/test_bates_validation.q";
    "tests/montecarlo/test_bates_greeks.q");

.test.calibrationFiles:(
    "tests/calibration/test_calibration_objective.q";
    "tests/calibration/test_black_scholes_vol_calibration.q";
    "tests/calibration/test_surface_calibration_diagnostics.q";
    "tests/calibration/test_heston_grid_calibration.q";
    "tests/calibration/test_calibration_validation.q";
    "tests/calibration/test_model_comparison.q";
    "tests/calibration/test_calibration_residuals.q";
    "tests/calibration/test_calibration_bucket_report.q";
    "tests/calibration/test_calibration_model_ranking.q";
    "tests/calibration/test_calibration_report_export.q";
    "tests/calibration/test_sabr_implied_vol.q";
    "tests/calibration/test_sabr_smile_generation.q";
    "tests/calibration/test_sabr_calibration.q";
    "tests/calibration/test_sabr_validation.q";
    "tests/calibration/test_sabr_model_comparison.q";
    "tests/calibration/test_merton_grid_calibration.q";
    "tests/calibration/test_merton_model_comparison.q";
    "tests/calibration/test_bates_grid_calibration.q";
    "tests/calibration/test_bates_model_comparison.q");

.test.modelcheckFiles:(
    "tests/modelcheck/test_model_limit_black_scholes.q";
    "tests/modelcheck/test_model_limit_heston_merton_bates.q";
    "tests/modelcheck/test_model_limit_local_vol_sabr.q";
    "tests/modelcheck/test_model_consistency_report.q";
    "tests/modelcheck/test_modelcheck_validation.q");

.test.riskFiles:(
    "tests/risk/test_pnl_distribution.q";
    "tests/risk/test_var_expected_shortfall.q";
    "tests/risk/test_var_report.q";
    "tests/risk/test_risk_contribution.q";
    "tests/risk/test_worst_loss_report.q";
    "tests/risk/test_var_validation.q";
    "tests/risk/test_historical_shock_table.q";
    "tests/risk/test_historical_replay_pricing.q";
    "tests/risk/test_historical_replay_distribution.q";
    "tests/risk/test_historical_var_report.q";
    "tests/risk/test_historical_worst_events.q";
    "tests/risk/test_historical_replay_validation.q";
    "tests/risk/test_historical_var_regression.q";
    "tests/risk/test_limit_table_validation.q";
    "tests/risk/test_var_limit_monitoring.q";
    "tests/risk/test_greek_limit_monitoring.q";
    "tests/risk/test_pnl_limit_monitoring.q";
    "tests/risk/test_limit_summary_report.q";
    "tests/risk/test_limit_monitoring_validation.q";
    "tests/risk/test_daily_risk_run.q";
    "tests/risk/test_daily_risk_run_validation.q");

.test.realdataFiles:(
    "tests/realdata/test_barchart_parser.q";
    "tests/realdata/test_barchart_multiday_replay.q";
    "tests/realdata/test_barchart_performance.q";
    "tests/realdata/test_barchart_model_pricing.q";
    "tests/realdata/test_barchart_model_vs_market.q");


.test.commodityFiles:(
    "tests/commodity/test_futures_curve.q";
    "tests/commodity/test_front_roll_backtest.q";
    "tests/commodity/test_calendar_spread_replay.q";
    "tests/commodity/test_black76_price.q";
    "tests/commodity/test_black76_implied_vol.q";
    "tests/commodity/test_spread_option.q";
    "tests/commodity/test_electricity_foundation.q";
    "tests/core/test_assetclass_registry.q";
    "tests/commodity/test_schwartz_transition_moments.q";
    "tests/commodity/test_schwartz_stationary_variance.q";
    "tests/commodity/test_schwartz_futures_price.q";
    "tests/commodity/test_schwartz_futures_curve_shape.q";
    "tests/commodity/test_schwartz_european_call.q";
    "tests/commodity/test_schwartz_european_put.q";
    "tests/commodity/test_schwartz_gbm_limit.q";
    "tests/commodity/test_schwartz_validation.q";
    "tests/commodity/test_schwartz2_factor_moments.q";
    "tests/commodity/test_schwartz2_futures_price.q";
    "tests/commodity/test_schwartz2_futures_curve_shape.q";
    "tests/commodity/test_schwartz2_european_call.q";
    "tests/commodity/test_schwartz2_european_put.q";
    "tests/commodity/test_schwartz2_correlation_effect.q";
    "tests/commodity/test_schwartz2_onefactor_limit.q";
    "tests/commodity/test_schwartz2_validation.q";
    "tests/commodity/test_mrjump_validation.q";
    "tests/commodity/test_mrjump_paths.q";
    "tests/commodity/test_mrjump_jump_diagnostics.q";
    "tests/commodity/test_mrjump_zero_jump_limit.q";
    "tests/commodity/test_mrjump_zero_jump_size.q";
    "tests/commodity/test_mrjump_option_call.q";
    "tests/commodity/test_mrjump_option_put.q";
    "tests/commodity/test_mrjump_scenario_paths.q";
    "tests/commodity/test_mrjump_lambda_zero_price_benchmark.q";
    "tests/commodity/test_mrjump_stationary_variance.q";
    "tests/commodity/test_mrjump_jump_intensity_sensitivity.q";
    "tests/commodity/test_mrjump_jump_mean_sensitivity.q";
    "tests/commodity/test_modelreport_price_all.q";
    "tests/commodity/test_modelreport_price_differences.q";
    "tests/commodity/test_modelreport_summary.q";
    "tests/commodity/test_modelreport_scenario_pnl.q";
    "tests/commodity/test_modelreport_validation.q";
    "tests/commodity/test_modelreport_ranking.q";
    "tests/commodity/test_modelreport_greeks_black76.q";
    "tests/commodity/test_modelreport_greeks_schwartz.q";
    "tests/commodity/test_modelreport_greeks_schwartz2.q";
    "tests/commodity/test_modelreport_greeks_mrjump.q";
    "tests/commodity/test_modelreport_greeks_all.q";
    "tests/commodity/test_modelreport_greeks_summary.q";
    "tests/commodity/test_modelreport_comparison_with_greeks.q";
    "tests/commodity/test_modelreport_primary_sensitivity.q";
    "tests/commodity/test_modelreport_price_disagreement.q";
    "tests/commodity/test_modelreport_greeks_disagreement.q";
    "tests/commodity/test_modelreport_scenario_disagreement.q";
    "tests/commodity/test_modelreport_disagreement_alerts.q";
    "tests/commodity/test_modelreport_comparison_risk.q";
    "tests/commodity/test_modelreport_disagreement_insufficient_models.q");

/ --- Combine all suites ---

.test.suites:(
    (`core;       .test.coreFiles);
    (`greeks;     .test.greeksFiles);
    (`american;   .test.americanFiles);
    (`barrier;    .test.barrierFiles);
    (`localvol;   .test.localVolFiles);
    (`market;     .test.marketFiles);
    (`portfolio;  .test.portfolioFiles);
    (`impliedvol; .test.impliedVolFiles);
    (`reporting;  .test.reportingFiles);
    (`infra;      .test.infraFiles);
    (`stress;     .test.stressFiles);
    (`montecarlo; .test.monteCarloFiles);
    (`calibration; .test.calibrationFiles);
    (`modelcheck; .test.modelcheckFiles);
    (`risk;       .test.riskFiles);
    (`realdata;   .test.realdataFiles);
    (`commodity;  .test.commodityFiles));

/ --- Runner ---

.test.pass:0;
.test.fail:0;

.test.runOne:{[testPath]
    -1 "  --- ",testPath," ---";
    codeLines:read0 hsym `$testPath;
    filteredLines:codeLines where not codeLines like "\\l *";
    codeBlock:"\n" sv filteredLines;
    testResult:@[{value x;`OK};codeBlock;{x}];
    if[testResult~`OK; .test.pass+:1];
    if[not testResult~`OK; -2 "    FAILED: ",$[10h=type testResult;testResult;string testResult]; .test.fail+:1];
 };

.test.runSuite:{[suitePair]
    suiteName:suitePair 0;
    testFiles:suitePair 1;
    -1 "";
    -1 "--- Suite: ",string[suiteName]," (",string[count testFiles]," tests) ---";
    .test.runOne each testFiles;
 };

.test.runAll:{[]
    -1 "=============================================================================";
    -1 " qFDM v",.qfdm.version," Test Suite";
    -1 "=============================================================================";
    .test.runSuite each .test.suites;
    -1 "";
    -1 "=============================================================================";
    -1 " Results: ",string[.test.pass]," passed, ",string[.test.fail]," failed";
    -1 "=============================================================================";
    if[.test.fail=0; -1 "All tests passed."];
    if[.test.fail>0; '"Some tests failed: ",string[.test.fail]," failures"];
 };

.test.runAll[];

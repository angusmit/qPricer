/ run_all_tests.q - grouped test runner (v0.14; v0.48 adds timing instrumentation)
/ Available test commands:
/   q tests/run_all_tests.q                    - full suite, prints timing report
/   q tests/run_fast_tests.q                   - fast tier (deterministic, no MC), < ~2 min
/   q tests/run_smoke_tests.q                  - very short core smoke tests, < 30 s
/   q tests/run_group.q <group> -q             - single group, e.g. strategy, core, commodity
\l core/init.q

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
    "tests/realdata/test_barchart_model_vs_market.q";
    "tests/realdata/test_parser_crude_filename.q";
    "tests/realdata/test_parser_crude_loadcontract.q";
    "tests/realdata/test_parser_crude_curve.q";
    "tests/realdata/test_parser_crude_curve_empty.q";
    "tests/realdata/test_hdb_ingest_query_equivalence.q");


.test.commodityFiles:(
    "tests/commodity/test_futures_curve.q";
    "tests/commodity/test_front_roll_backtest.q";
    "tests/commodity/test_calendar_spread_replay.q";
    "tests/commodity/test_black76_price.q";
    "tests/commodity/test_black76_implied_vol.q";
    "tests/commodity/test_spread_option.q";
    "tests/commodity/test_spread_option_margrabe.q";
    "tests/commodity/test_spread_option_mc.q";
    "tests/commodity/test_seasonality_factor.q";
    "tests/commodity/test_seasonality_overlay_baseline.q";
    "tests/commodity/test_seasonality_overlay_pattern.q";
    "tests/commodity/test_calibrate_curve_roundtrip.q";
    "tests/commodity/test_calibrate_curve_economic.q";
    "tests/commodity/test_calibrate_curve_guards.q";
    "tests/commodity/test_convenience_yield_series_tracking.q";
    "tests/commodity/test_convenience_yield_series_skip.q";
    "tests/commodity/test_kalman_filter_correctness.q";
    "tests/commodity/test_kalman_missing_observations.q";
    "tests/commodity/test_kalman_schwartz_smith_roundtrip.q";
    "tests/commodity/test_parser_futures.q";
    "tests/commodity/test_seasonality_fit.q";
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
    "tests/commodity/test_modelreport_disagreement_insufficient_models.q";
    "tests/commodity/test_modelreport_portfolio_positions.q";
    "tests/commodity/test_modelreport_run_position_risk.q";
    "tests/commodity/test_modelreport_run_portfolio_risk.q";
    "tests/commodity/test_modelreport_portfolio_alert_summary.q";
    "tests/commodity/test_modelreport_portfolio_exposure.q";
    "tests/commodity/test_modelreport_worst_offenders.q";
    "tests/commodity/test_modelreport_portfolio_dashboard.q";
    "tests/commodity/test_modelreport_portfolio_error_isolation.q";
    "tests/commodity/test_modelreport_limit_check.q";
    "tests/commodity/test_modelreport_portfolio_limits_ok.q";
    "tests/commodity/test_modelreport_portfolio_limits_warning.q";
    "tests/commodity/test_modelreport_portfolio_limits_breach.q";
    "tests/commodity/test_modelreport_limit_breach_report.q";
    "tests/commodity/test_modelreport_portfolio_risk_with_limits.q";
    "tests/commodity/test_modelreport_portfolio_limits_error_isolation.q";
    "tests/commodity/test_modelreport_limit_snapshot.q";
    "tests/commodity/test_modelreport_limit_history_append.q";
    "tests/commodity/test_modelreport_limit_history_summary.q";
    "tests/commodity/test_modelreport_limit_breach_trend.q";
    "tests/commodity/test_modelreport_repeated_breaches.q";
    "tests/commodity/test_modelreport_limit_history_dashboard.q";
    "tests/commodity/test_modelreport_limit_history_wrapper.q");

.test.allocFiles:(
    "tests/portfolio/test_alloc_weights.q";
    "tests/portfolio/test_alloc_backtest.q");

.test.regimeFiles:(
    "tests/regime/test_regime_label.q";
    "tests/regime/test_regime_breakdown.q";
    "tests/regime/test_regime_analogue_distance.q";
    "tests/regime/test_regime_analogue_nearest.q");

.test.govFiles:(
    "tests/gov/test_gov_deflated_sharpe.q";
    "tests/gov/test_gov_ledger.q";
    "tests/gov/test_gov_gates.q";
    "tests/gov/test_gov_zones.q";
    "tests/gov/test_gov_holdout_oneshot.q";
    "tests/gov/test_gov_verdict_failsafe.q";
    "tests/gov/test_gov_skeptic_annotation.q");

.test.registryFiles:(
    "tests/registry/test_registry_core.q";
    "tests/registry/test_contracts_conform.q";
    "tests/registry/test_existing_registered.q");

.test.cardsFiles:(
    "tests/cards/test_cards_core.q";
    "tests/cards/test_cards_validation_status.q";
    "tests/cards/test_cards_audit.q";
    "tests/cards/test_cards_gated_run.q";
    "tests/cards/test_cards_audit_dynamic.q");

.test.templatesFiles:(
    "tests/templates/test_template_registry.q";
    "tests/templates/test_directional_template_faithful.q";
    "tests/templates/test_relative_value.q");

.test.workflowFiles:(
    "tests/workflow/test_workflow_loop.q";
    "tests/workflow/test_workflow_bounded.q";
    "tests/workflow/test_workflow_runreplay.q");

.test.factorFiles:(
    "tests/factor/test_factor_pca.q";
    "tests/factor/test_factor_template.q");

.test.stateFiles:(
    "tests/state/test_asof_accessor.q";
    "tests/state/test_market_state.q");

.test.curveFiles:(
    "tests/curve/test_curve_engine.q";
    "tests/curve/test_curve_snapshot.q");

.test.rollFiles:(
    "tests/roll/test_roll_active.q";
    "tests/roll/test_roll_events_continuous.q");

.test.replayFiles:(
    "tests/replay/test_replay_loop.q";
    "tests/replay/test_replay_execution.q");

.test.evidenceFiles:(
    "tests/evidence/test_audit_passes.q";
    "tests/evidence/test_audit_bites.q";
    "tests/evidence/test_gated_run.q");

.test.seasonFiles:(enlist "tests/season/test_seasonality.q");
.test.carryFiles:(enlist "tests/carry/test_carry.q");
.test.attributionFiles:(enlist "tests/attribution/test_attribution.q");

.test.executionFiles:(
    "tests/execution/test_exec_fill.q";
    "tests/execution/test_exec_commodity_bt_legacy.q";
    "tests/execution/test_exec_realism.q";
    "tests/execution/test_exec_equity_legacy.q";
    "tests/execution/test_exec_equity_realism.q";
    "tests/execution/test_exec_equity_legcost.q");

.test.strategyFiles:(
    "tests/strategy/test_strategy_registry.q";
    "tests/strategy/test_strategy_path_synthetic.q";
    "tests/strategy/test_strategy_path_barchart_adapter.q";
    "tests/strategy/test_strategy_gamma_scalp_synthetic.q";
    "tests/strategy/test_strategy_gamma_scalp_band.q";
    "tests/strategy/test_strategy_gamma_scalp_reconciliation.q";
    "tests/strategy/test_strategy_gamma_scalp_barchart.q";
    "tests/strategy/test_strategy_hedge_helper_shared.q";
    "tests/strategy/test_strategy_short_variance_gate.q";
    "tests/strategy/test_strategy_short_variance_synthetic.q";
    "tests/strategy/test_strategy_short_variance_premium.q";
    "tests/strategy/test_strategy_short_variance_recon.q";
    "tests/strategy/test_strategy_calendar_roll_legs.q";
    "tests/strategy/test_strategy_calendar_roll_event.q";
    "tests/strategy/test_strategy_calendar_roll_accounting.q";
    "tests/strategy/test_strategy_calendar_roll_recon.q";
    "tests/strategy/test_strategy_calendar_roll_unhedged.q";
    "tests/strategy/test_strategy_portfolio_value_helper.q";
    "tests/strategy/test_strategy_risk_reversal.q";
    "tests/strategy/test_strategy_risk_reversal_gate.q";
    "tests/strategy/test_strategy_risk_reversal_accounting.q";
    "tests/strategy/test_strategy_model_disagreement.q";
    "tests/strategy/test_strategy_model_disagreement_gate.q";
    "tests/strategy/test_strategy_model_disagreement_accounting.q";
    "tests/strategy/test_strategy_delta_vega_hedge.q";
    "tests/strategy/test_strategy_delta_vega_hedge_accounting.q";
    "tests/strategy/test_strategy_delta_vega_hedge_neutrality.q";
    "tests/strategy/test_strategy_path_ensemble.q";
    "tests/strategy/test_strategy_portfolio_run_ensemble.q";
    "tests/strategy/test_strategy_portfolio_performance.q";
    "tests/strategy/test_strategy_portfolio_correlation.q";
    "tests/strategy/test_strategy_portfolio_dashboard.q";
    "tests/strategy/test_strategy_long_vol.q";
    "tests/strategy/test_strategy_long_vol_gate.q";
    "tests/strategy/test_strategy_long_vol_accounting.q";
    "tests/strategy/test_strategy_collar_tail_hedge.q";
    "tests/strategy/test_strategy_collar_tail_hedge_mode.q";
    "tests/strategy/test_strategy_collar_tail_hedge_accounting.q";
    "tests/strategy/test_strategy_put_ratio_backspread.q";
    "tests/strategy/test_strategy_put_ratio_backspread_accounting.q";
    "tests/strategy/test_strategy_iron_condor.q";
    "tests/strategy/test_strategy_iron_condor_payoff.q";
    "tests/strategy/test_strategy_iron_condor_accounting.q";
    "tests/strategy/test_strategy_barrier_hedge.q";
    "tests/strategy/test_strategy_barrier_hedge_knockout.q";
    "tests/strategy/test_strategy_barrier_hedge_accounting.q";
    "tests/strategy/test_strategy_jump_premium.q";
    "tests/strategy/test_strategy_jump_premium_gate.q";
    "tests/strategy/test_strategy_jump_premium_accounting.q";
    "tests/strategy/test_strategy_path_correlated_determinism.q";
    "tests/strategy/test_strategy_path_correlated_schema.q";
    "tests/strategy/test_strategy_dispersion.q";
    "tests/strategy/test_strategy_dispersion_correlation.q";
    "tests/strategy/test_strategy_dispersion_accounting.q";
    "tests/strategy/test_strategy_path_futures_curve_determinism.q";
    "tests/strategy/test_strategy_path_futures_curve_schema.q";
    "tests/strategy/test_strategy_power_spike.q";
    "tests/strategy/test_strategy_power_spike_gate.q";
    "tests/strategy/test_strategy_power_spike_accounting.q";
    "tests/strategy/test_strategy_commodity_calendar.q";
    "tests/strategy/test_strategy_commodity_calendar_roll.q";
    "tests/strategy/test_strategy_commodity_calendar_accounting.q";
    "tests/strategy/test_strategy_path_correlated_curves_schema.q";
    "tests/strategy/test_strategy_path_correlated_curves_determinism.q";
    "tests/strategy/test_strategy_spark_spread.q";
    "tests/strategy/test_strategy_spark_spread_accounting.q";
    "tests/strategy/test_strategy_crack_spread.q";
    "tests/strategy/test_strategy_crack_spread_accounting.q";
    "tests/strategy/test_commodity_signals_roll.q";
    "tests/strategy/test_commodity_signals_causality.q";
    "tests/strategy/test_strategy_convenience_carry.q";
    "tests/strategy/test_strategy_chi_reversion.q";
    "tests/strategy/test_strategy_real_curve_calendar_roll.q";
    "tests/strategy/test_commodity_score.q";
    "tests/strategy/test_strategy_time_series_momentum.q";
    "tests/strategy/test_strategy_two_timescale.q";
    "tests/strategy/test_strategy_storage_cash_carry.q";
    "tests/strategy/test_strategy_carry_momentum_combo.q";
    "tests/strategy/test_strategy_curve_relative_value.q";
    "tests/strategy/test_walk_forward_splits.q";
    "tests/strategy/test_walk_forward_aggregate.q";
    "tests/strategy/test_walk_forward_run.q";
    "tests/strategy/test_deseasonalize_roundtrip.q";
    "tests/strategy/test_commodity_signals_deseasonalize.q";
    "tests/strategy/test_cross_commodity.q";
    "tests/strategy/test_seasonal_calendar_spread.q");

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
    (`commodity;  .test.commodityFiles);
    (`execution;  .test.executionFiles);
    (`alloc;      .test.allocFiles);
    (`regime;     .test.regimeFiles);
    (`gov;        .test.govFiles);
    (`registry;   .test.registryFiles);
    (`cards;      .test.cardsFiles);
    (`templates;  .test.templatesFiles);
    (`workflow;   .test.workflowFiles);
    (`factor;     .test.factorFiles);
    (`state;      .test.stateFiles);
    (`curve;      .test.curveFiles);
    (`roll;       .test.rollFiles);
    (`replay;     .test.replayFiles);
    (`evidence;   .test.evidenceFiles);
    (`season;     .test.seasonFiles);
    (`carry;      .test.carryFiles);
    (`attribution; .test.attributionFiles);
    (`strategy;   .test.strategyFiles));

/ --- Runner with timing ---

.test.pass:0;
.test.fail:0;
.test.timings:();
.test.currentGroup:`unknown;

.test.runOne:{[testPath]
    -1 "  --- ",testPath," ---";
    codeLines:read0 hsym `$testPath;
    filteredLines:codeLines where not codeLines like "\\l *";
    codeBlock:"\n" sv filteredLines;
    startTime:.z.p;
    testResult:@[{value x;`OK};codeBlock;{x}];
    elapsedMs:(.z.p-startTime)%1000000;
    .test.timings,:enlist `groupName`testPath`elapsedMs`status!(.test.currentGroup;testPath;`float$elapsedMs;$[testResult~`OK;`OK;`FAIL]);
    if[testResult~`OK; .test.pass+:1];
    if[not testResult~`OK; -2 "    FAILED: ",$[10h=type testResult;testResult;string testResult]; .test.fail+:1];
 };

.test.runSuite:{[suitePair]
    suiteName:suitePair 0;
    testFiles:suitePair 1;
    -1 "";
    -1 "--- Suite: ",string[suiteName]," (",string[count testFiles]," tests) ---";
    .test.currentGroup:suiteName;
    .test.runOne each testFiles;
 };

.test.printTimingReport:{[]
    if[0=count .test.timings; :()];
    timingsTbl:.test.timings;
    -1 "";
    -1 "=============================================================================";
    -1 " Timing report";
    -1 "=============================================================================";
    totalMs:sum timingsTbl`elapsedMs;
    -1 " Total wall time: ",string[`long$totalMs]," ms (",string[`long$totalMs%1000]," s)";
    -1 "";
    -1 " Per-group subtotals:";
    groupSums:0!select totalMs:sum elapsedMs, testCount:count i by groupName from timingsTbl;
    groupSumsSorted:`totalMs xdesc groupSums;
    {[r] -1 "   ",(-12$string r`groupName)," ",(-8$string `long$r`totalMs)," ms (",string[r`testCount]," tests)"} each groupSumsSorted;
    -1 "";
    -1 " Slowest 20 tests:";
    slowestSorted:`elapsedMs xdesc timingsTbl;
    takeN:20&count slowestSorted;
    slowest:takeN#slowestSorted;
    {[r] -1 "   ",(-8$string `long$r`elapsedMs)," ms  ",r`testPath} each slowest;
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
    .test.printTimingReport[];
    if[.test.fail>0; '"Some tests failed: ",string[.test.fail]," failures"];
 };

/ Run unless a downstream runner pre-set .test.skipAutoRun (used by run_group.q,
/ run_fast_tests.q to reuse suite definitions without firing the full suite).
if[not `skipAutoRun in key `.test; .test.runAll[]];

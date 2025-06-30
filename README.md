# qPricer — Option Pricing in q/kdb+

A high-performance financial pricing engine written in `q`, designed for expressive modeling and real-time analytics in derivatives pricing.

Work in Progress — Bringing classic and exotic option pricing models to kdb+/q.

---

## Planned Features

- Greeks  
  Sensitivity measures: delta, gamma, vega, theta, rho — computed in pure q

- Jump Diffusion Models  
  Support for discontinuous asset price models such as Merton or Kou processes

- Asian Option Pricing  
  Implementation of geometric and arithmetic average pricing (analytical and simulation-based)

- Implied Volatility Solver  
  Numerical solver to infer market-implied volatility from observed prices

---

## Reference

For analytical solutions to geometric Asian options in Python:  
[Geometric Asian Option – Analytical Solution](https://github.com/JynxC98/quantitative_finance/blob/main/geometric-asian-option/python_method/analytical_solution.py)


## Technology Stack

- Core Logic: `q` (kdb+)  
- Future Extensions: Integration with Python using `embedPy` or `PyQ`  

---

## Contributing

The project is in early development. Star the repository, fork it, or open an issue to get involved. Contributions and feedback are welcome.

# MOJuMP.jl: Multi-objective Optimization extension for JuMP using AUGMECON in Julia

Install by executing in the REPL,

`Pkg.clone("git://github.com/UtkarshChemE/MOJuMP.jl.git")`

This package defines two function `tradeoff_table` and `eff_pareto`. Use Julia help `?` at the REPL to get more information about their use. The methods are based on the reference cited. Look for examples in `notebooks` folder. This is a WIP package.

**Disclaimer**: MOJuMP is *not* developed or maintained by the JuMP developers.

> **Warning**: Always define your objectives as minimization problem. Current implementation can only handle minimization problem. Convert your maximization objective to minimization by multiplying it with $-1$

# Reference
Mavrotas, George. "Effective implementation of the ε-constraint method in multi-objective mathematical programming problems." Applied mathematics and computation 213.2 (2009): 455-465.


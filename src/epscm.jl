function create_eps_model()
    # P = [:Lignite, :Oil, :Gas, :RES] # Power Gen Units
    # L = [:base, :middle, :peak]
    # DIRS = [:Cost, :CO2, :Endo]
    Annual_Demand = 64000.0
    df = [0.6, 0.3, 0.1]
    demand = df .* 64000.0
    pData = [31000 15000 22000 10000; 30 75 60 90; 1.44 0.72 0.45 0]
    PLset = collect(Iterators.product(1:4, 1:3))[[2,8,9]]

    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), 
            "NonConvex" => 1, "OutputFlag" => 0))
    @variables(m, begin
            x[1:4,1:3] ≥ 0
            Z[1:3]
    end)

    @constraints(m, begin
                objcost, sum(pData[2,:]'*x) == Z[1]
                objCO2, sum(pData[3,:]'*x) == Z[2]
                objEndo, sum(x[[1,4],:]) == -Z[3]
                xzero_set[p=1:4,l=1:3; (p,l) in PLset], x[p,l] == 0
                defcap, sum(x; dims=2) .≤ pData[1, :]
                defdem, sum(x; dims=1)[:] .≥ demand
     end)
    return m, Z
end

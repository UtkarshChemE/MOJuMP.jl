module MOJuMP

export tradeoff_table, eff_pareto

import JuMP, Sobol
import Combinatorics: permutations
import Base.Iterators: partition

"""
This is a little modification on `upper_bound` function to handle clean non-existance of upper bound
"""
get_upper_bound(v::JuMP.VariableRef) = JuMP.has_upper_bound(v) ? JuMP.upper_bound(v) : nothing

"""
Reset the upper bound to original provided in the model or delete it if original upper bound doesn't exists
"""
reset_upper_bound(v::JuMP.VariableRef, ogbound::Real) = JuMP.set_upper_bound(v, ogbound)
reset_upper_bound(v::JuMP.VariableRef, ogbound::Nothing) = JuMP.delete_upper_bound(v)

# """
# Break the series into N chunks and return either last one or everything other than that
# """
# function (collection, n, last=false)
#     nchunks = ceil(Int, length(collection)/n)
#     chunks = collect(partition(collection, n))
#     !last && return reduce(vcat, chunks[1:end-1])


"""
`tradeoff_table(m::JuMP.Model, OBJs::Vector{JuMP.VariableRef}, VARs::Union{Nothing,Vector{JuMP.VariableRef}}=nothing)`

Generate the tradeoff table for the model and objectives passed as variables
Objectives are always be minimized
"""
function tradeoff_table(m::JuMP.Model, OBJs::Vector{JuMP.VariableRef},
                        VARs::Union{Nothing,Vector{JuMP.VariableRef}}=nothing)
    #Error Checking
    @assert length(OBJs) ≥ 2 "Atleast two objectives need to specified"
    @assert all(JuMP.is_valid.(Ref(m), OBJs)) "A objective is not defined in the model"
    isnothing(VARs) || @assert all(JuMP.is_valid.(Ref(m), VARs)) "A variable to collect is not present in the model"

    #Create Arrays to store result
    ogbound = get_upper_bound.(OBJs)
    table_vals = zeros(factorial(length(OBJs)), length(OBJs))
    varsValues = isnothing(VARs) ? nothing : zeros(factorial(length(OBJs)), length(VARs))
    for (i,p) in enumerate(permutations(OBJs))
        for obj in p
            JuMP.@objective(m, Min, obj)
            JuMP.optimize!(m)
            JuMP.set_upper_bound(obj, JuMP.value(obj))
        end
        table_vals[i,:] .= JuMP.value.(OBJs)
        if !isnothing(varsValues)
            varsValues[i,:] .= JuMP.value.(VARs)
        end
        reset_upper_bound.(OBJs, ogbound)
    end
    return table_vals, varsValues
end     

"""
`eff_pareto(m::JuMP.Model, OBJs::Vector{JuMP.VariableRef}, VARs::Union{Nothing,Vector{JuMP.VariableRef}}=nothing; kwargs...)`

Find Pareto optimal solution using AUGMECON. The output points should be almost pareto optimal, if not the function will throw warning. If precomputed tradeoff table is not provided , the function will compute it for you.

# Arguments
- `samples::Int=20`: Number of sampling points on pareto surface
- `tOffTab::Matrix`: Precomputed tradeoff table of size N!×N, where `N=length(OBJs)`
- `tOffVarsValue::Matrix`: Precomputed value matrix for vars to be collected of size N!×M, where `N, M = length.((OBJs, VARs))`. Must be provided if both VARs and tOffTab are non-empty.
"""
function eff_pareto(m::JuMP.Model, OBJs::Vector{JuMP.VariableRef}, 
                    VARs::Union{Nothing,Vector{JuMP.VariableRef}}=nothing;
                    tOffTab::Union{Nothing,Matrix{Real}}=nothing, 
                    tOffVarsValue::Union{Nothing,Matrix{Real}}=nothing, samples::Int=20)
    N = length(OBJs)
    #Assertion Section
    @assert N ≥ 2 "Atleast two objectives need to specified"
    @assert all(JuMP.is_valid.(Ref(m), OBJs)) "A objective is not defined in the model"
    samples ≥ N || throw(DomainError(samples, "Atleast $(N) sample needs to be collected"))
    isnothing(tOffTab) || @assert all(size(tOffTab) .== (factorial(N), N)) "The trade-off table provided is not of appropiate size. It should be a $(factorial(N))×$(N) matrix. Providing trade-off table is not necessary"
    isnothing(VARs) || @assert all(JuMP.is_valid.(Ref(m), VARs)) "A variable to collect is not present in the model"
    if !isnothing(VARs) && !isnothing(tOffTab)
        @assert isnothing(tOffTab) == isnothing(tOffVarsValue) "Both objective value for tradeoff table and variables to be collect needs to be provided or neither should be provided"
    end

    ogbound = get_upper_bound.(OBJs)
    tOffTab, tOffVarsValue = isnothing(tOffTab) ? tradeoff_table(m, OBJs, VARs) : (tOffTab, tOffVarsValue)
    JuMP.set_lower_bound.(OBJs, minimum(tOffTab; dims=1)[:])
    JuMP.set_upper_bound.(OBJs, maximum(tOffTab; dims=1)[:])
    obj_range = reduce.(-, extrema(tOffTab; dims=1)) .|> abs |> vec
    #Create random Sobol sampler
    s = Sobol.SobolSeq(minimum(tOffTab[:, 2:end]; dims = 1),
        maximum(tOffTab[:, 2:end]; dims = 1))
    ϵ = zeros(samples, length(OBJs)-1)
    slacks = zeros(length(OBJs)-1) #Ideally this slacks should be zero
    nadirPoints = zeros(samples, length(OBJs))
    varsValues = isnothing(VARs) ? nothing : zeros(samples, length(VARs))
    
    # Always minimizing with respect to first objective. Set other objectives as epsilon constraints
    JuMP.@objective(m, Min, OBJs[1])
    for i = 1:(samples - samples÷N)
        ϵ[i, :] .= Sobol.next!(s)
        JuMP.set_upper_bound.(OBJs[2:end], ϵ[i, :])
        JuMP.optimize!(m)
        nadirPoints[i, :] .= JuMP.value.(OBJs)
        slacks .= (nadirPoints[i, 2:end] - ϵ[i, :]) ./ obj_range[2:end]
        any(x -> x ≥ 1e-3, slacks) && @warn "Slack not zero for eps=$(ϵ[i, :])"
        if !isnothing(varsValues)
            varsValues[i,:] .= JuMP.value.(VARs)
        end
    end

    # This is done to ensure uniform sampling in first objective's direction
    rOBJs = reverse(OBJs)
    tOffTab = reverse(tOffTab; dims=2)
    obj_range = reverse(obj_range)
    s = Sobol.SobolSeq(minimum(tOffTab[:, 2:end]; dims = 1),
        maximum(tOffTab[:, 2:end]; dims = 1))
    JuMP.@objective(m, Min, rOBJs[1])
    for i = (samples - samples÷N + 1):samples
        ϵ[i, :] .= Sobol.next!(s)
        JuMP.set_upper_bound.(rOBJs[2:end], ϵ[i, :])
        JuMP.optimize!(m)
        nadirPoints[i, :] .= JuMP.value.(OBJs)
        slacks .= (JuMP.value.(rOBJs)[2:end] - ϵ[i, :]) ./ obj_range[2:end]
        any(x -> x ≥ 1e-3, slacks) && @warn "Slack not zero for eps=$(ϵ[i, :])"
        if !isnothing(varsValues)
            varsValues[i,:] .= JuMP.value.(VARs)
        end
    end
    reset_upper_bound.(OBJs, ogbound)
    tOffTab = reverse(tOffTab; dims=2)
    nadirPoints = [nadirPoints; tOffTab]
    varsValues = [varsValues; tOffVarsValue]
    return nadirPoints, varsValues
end

end # module

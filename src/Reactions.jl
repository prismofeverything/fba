module Reactions

import JSON
import MathProgBase
import GLPKMathProgInterface

function readnetwork(path)
    JSON.parsefile("network/$path.json")
end

function selectkeys(dict, ks)
    Dict([pair for pair in collect(dict) if pair[1] in ks])
end

logicoperations = Dict(
    "is" => x -> isa(x[1], Bool) ? x[1] : x[1] > 0,
    "not" => x -> isa(x[1], Bool) ? !x[1] : x[1] == 0,
    ">" => x -> x[1] > x[2],
    "and" => all,
    "or" => any
)

function runlogic(logic, state, ops)
    if isa(logic, Number)
        logic
    elseif isa(logic, String)
        logic[1] == '"' ? logic[2:end-1] : get(state, logic, 0)
    elseif !isempty(logic)
        op = ops[logic[1]]
        args = map(sublogic -> runlogic(sublogic, state, ops), logic[2:end])
        op(args)
    else
        true
    end
end

function applylogic(logicmap, state, ops)
    Dict(key => runlogic(logic, state, ops) for (key, logic) in logicmap)
end

function moleculestate(network)
    keysets = [keys(reaction["reaction"]) for reaction in values(network["reactions"])]
    molecules = union(keysets...)
    state = Dict(zip(molecules, zeros(length(molecules))))
    initial = selectkeys(network["initial"], molecules)

    merge(state, initial)
end

function boundsvector(reactions, bounds, extreme)
    [get(bounds, reaction, extreme) for reaction in reactions]
end

function buildcolumn(molecules, reaction)
    map(m -> m in keys(reaction) ? reaction[m] : 0, molecules)
end

function buildstoichiometry(molecules, reactions, network)
    columns = [buildcolumn(molecules, network["reactions"][reaction]["reaction"]) for reaction in reactions]
    hcat(columns...)
end

function initialize(network, ops)
    state = moleculestate(network)
    molecules = sort(collect(keys(state)))
    reactions = sort(collect(keys(network["reactions"])))

    regulation = applylogic(network["regulation"], merge(state, network["initial"]), ops)
    active = filter(r -> runlogic(network["reactions"][r]["regulation"], regulation, ops), reactions)

    lower = boundsvector(active, network["lower"], -Inf)
    upper = boundsvector(active, network["upper"], Inf)
    stoichiometry = buildstoichiometry(molecules, active, network)

    Dict(
        "molecules" => molecules,
        "reactions" => reactions,
        "active" => active,
        "lower" => lower,
        "upper" => upper,
        "stoichiometry" => stoichiometry
    )
end

function fluxbalanceanalysis(network, maximize)
    conditions = initialize(network, logicoperations)
    solver = GLPKMathProgInterface.GLPKSolverLP(method=:Simplex, presolve=true)
    objective = [get(maximize, reaction, 0.0) for reaction in conditions["active"]]
    println(objective)

    solution = MathProgBase.linprog(
        objective,
        conditions["stoichiometry"],
        '=',
        0.0,
        conditions["lower"],
        conditions["upper"],
        solver
    )

    Dict(zip(conditions["active"], solution.sol))
end

end

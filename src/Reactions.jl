module Reactions

import JSON
import MathProgBase
import GLPKMathProgInterface

"Read a local network from the example network directory"
function readnetwork(name)
    JSON.parsefile("network/$name.json")
end

"Select a subset of a dictionary with the given keys"
function selectkeys(dict, ks)
    Dict([pair for pair in collect(dict) if pair[1] in ks])
end

"Default logical operations for evaluating regulation expressions"
logicoperations = Dict(
    "is" => x -> isa(x[1], Bool) ? x[1] : x[1] > 0,
    "not" => x -> isa(x[1], Bool) ? !x[1] : x[1] == 0,
    ">" => x -> x[1] > x[2],
    "and" => all,
    "or" => any
)

"Evaluate logical expressions governing regulation"
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

"Substitute the values of a dictionary representing logical expressions with their values given the supplied state"
function applylogic(logicmap, state, ops)
    Dict(key => runlogic(logic, state, ops) for (key, logic) in logicmap)
end

"Extract the initial state of the molecules in a network"
function moleculestate(network)
    keysets = [keys(reaction["reaction"]) for reaction in values(network["reactions"])]
    molecules = union(keysets...)
    state = Dict(zip(molecules, zeros(length(molecules))))
    initial = selectkeys(network["initial"], molecules)

    merge(state, initial)
end

"Build the bounds vector with the given default extreme if missing"
function boundsvector(reactions, bounds, extreme)
    [get(bounds, reaction, extreme) for reaction in reactions]
end

"Build one column of the stoichiometric matrix from the dictionary describing the molecules involved in the reaction"
function buildcolumn(molecules, reaction)
    map(m -> m in keys(reaction) ? reaction[m] : 0, molecules)
end

"Build the stoichiometric matrix from the given dictionary of reactions describing what molecules are consumed and produced"
function buildstoichiometry(molecules, reactions, active)
    columns = [buildcolumn(molecules, reactions[reaction]["reaction"]) for reaction in active]
    hcat(columns...)
end

"Initialize the parameters of the linear programming problem given the network description"
function initialize(network, ops)
    state = moleculestate(network)
    molecules = sort(collect(keys(state)))
    reactions = sort(collect(keys(network["reactions"])))

    regulation = applylogic(network["regulation"], merge(state, network["initial"]), ops)
    active = filter(r -> runlogic(network["reactions"][r]["regulation"], regulation, ops), reactions)

    lower = boundsvector(active, network["lower"], -Inf)
    upper = boundsvector(active, network["upper"], Inf)
    stoichiometry = buildstoichiometry(molecules, network["reactions"], active)

    Dict(
        "molecules" => molecules,
        "reactions" => reactions,
        "active" => active,
        "lower" => lower,
        "upper" => upper,
        "stoichiometry" => stoichiometry
    )
end

"Perform the linear programming for the flux balance analysis problem defined by the given network and maximization vector"
function fluxbalanceanalysis(network, objective)
    conditions = initialize(network, logicoperations)
    solver = GLPKMathProgInterface.GLPKSolverLP(method=:Simplex, presolve=true)
    objective = [get(objective, reaction, 0.0) for reaction in conditions["active"]]

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

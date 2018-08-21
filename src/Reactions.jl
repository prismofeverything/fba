module Reactions

import JSON

function readtsv(path)
    reactions = Dict{String, Dict{String, Float64}}()
    open(path) do file
        for line in enumerate(eachline(file))
            if line[1] > 1
                parts = split(line[2], "	")
                key = JSON.parse(parts[1])
                inner = JSON.parse(parts[2])
                reactions[key] = inner
            end
        end
    end
    reactions
end

function readjson(path)
    JSON.parsefile(path)
end

logicoperations = Dict(
    "is" => x -> isa(x[1], Bool) ? x[1] : x[1] > 0,
    "not" => x -> isa(x[1], Bool) ? !x[1] : x[1] == 0,
    ">" => x -> x[1] > x[2],
    "and" => all,
    "or" => any
)

function moleculestate(network)
    molecules = Dict{String, Float64}()

    for molecule in keys(network["initial"])
        molecules[molecule] = network["initial"][molecule]
    end

    for reaction in keys(network["reactions"])
        for molecule in keys(network["reactions"][reaction]["reaction"])
            if !(molecule in keys(molecules))
                molecules[molecule] = 0
            end
        end
    end

    molecules
end

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

function buildcolumn(molecules, reaction)
    map(m -> m in keys(reaction) ? reaction[m] : 0, molecules)
end

function buildrow(molecules, regulation, reaction, ops)
    allowed = runlogic(reaction["regulation"], regulation, ops)
    if allowed
        buildcolumn(molecules, reaction["reaction"])
    else
        zeros(Float64, length(molecules))
    end
end

function buildstoichiometry(molecules, reactions, network, state, ops)
    regulation = applylogic(network["regulation"], state, ops)
    hcat(map(r -> buildrow(molecules, regulation, network["reactions"][r], ops), reactions)...)
end

function initialize(network, ops)
    state = moleculestate(network)
    molecules = sort(collect(keys(state)))
    reactions = sort(collect(keys(network["reactions"])))
    stoichiometry = buildstoichiometry(molecules, reactions, network, state, ops)

    Dict(
        "molecules" => molecules,
        "reactions" => reactions,
        "stoichiometry" => stoichiometry
    )
end

end

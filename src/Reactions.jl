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

logicoperations = Dict(
    "if" => x -> x[1] > 0,
    "not" => x -> x[1] == 0,
    ">" => x -> x[1] > x[2],
    "and" => all,
    "or" => any
)

function runlogic(logic, state, ops)
    if typeof(logic) == String
        if logic[1] == '"'
            logic[2:end-1]
        else
            state[logic]
        end
    else
        op = ops[logic[1]]
        args = map(sublogic -> runlogic(sublogic, state, ops), logic[2:end])
        op(args)
    end
end

end

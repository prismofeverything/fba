module Reactions
import JSON

function readreactions(path)
    reactions = Dict{String, Dict{String, Float64}}()
    open(path) do file
        for line in enumerate(eachline(file))
            if line[1] > 1
                parts = line[2].split("	")
                inner = JSON.parse(parts[2])
                reactions[parts[1]] = JSON.parse(parts[2])
            end
        end
    end
    reactions
end

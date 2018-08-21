module Reactions

import JSON

function readreactions(path)
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

end

# Reactions

A minimal implementation of [Flux Balance Analysis](https://en.wikipedia.org/wiki/Flux_balance_analysis) with regulation in [Julia](https://julialang.org/).

## example

In your Julia repl, activate the project:

    julia> ]
    (v1.0) pkg> activate .
    (Reactions) pkg> [backspace]
    julia>

Now import the `Reactions` module and load one of the predefined networks (here, `threonine`):

    julia> import Reactions
    julia> network = Reactions.readnetwork("threonine")

    Dict{String,Any} with 5 entries:
      "regulation" => Dict{String,Any}()
      "initial"    => Dict{String,Any}()
      "lower"      => Dict{String,Any}("VthrC"=>0,"Bhms"=>0,"Bphs"=>0)
      "reactions"  => Dict{String,Any}("BATP"=>Dict{String,Any}("regulation"=>Any[],"reaction"=>Dict{String,Any}("ATP"=>1)),"VthrB"=>Dict{String,Any}("regulation"=>Any[],"reaction"=>Dict{String,Any}("ATP"=>-1,"H+"=>2,"ADP"=>1,"hms"=>-1,"phs"=>1)),"VthrC"=>Dict{String,Any}("regu…
      "upper"      => Dict{String,Any}("Bhms"=>10,"Bphs"=>0)

Then run the flux balance analysis on this network:

    julia> Reactions.fluxbalanceanalysis(threonine, Dict("Bthr" => 1))

    Dict{String,Float64} with 10 entries:
      "BATP"  => 10.0
      "VthrB" => 10.0
      "BH2O"  => 10.0
      "Bthr"  => -10.0
      "VthrC" => 10.0
      "BADP"  => -10.0
      "Bhms"  => 10.0
      "Bphs"  => 0.0
      "Bphos" => -10.0
      "BH+"   => -20.0

## initialization

This library works with networks defined by a dictionary containing several keys:

* `reactions` - a dictionary of each reaction in terms of molecules it consumes and produces, and also their activation with respect to regulation.
* `regulation` - a dictionary of reactions and their conditions for activation.
* `lower` - any lower bounds placed on the given reactions.
* `upper` - upper bounds on the reactions.
* `initial` - the initial state of the system, given to activate or deactivate various regulation proteins.

This dictionary can be instantiated directly in code, or it can be read from a json file:

    network = JSON.parsefile("/absolute/path/to/network.json")

If you want to instantiate one of the provided example networks, you can use the `readnetwork` function which looks in the local `network` directory:

    network = readnetwork("threonine")

Each key in the network dictionary is described below.

### reactions

Under the reactions key is a dictionary where each key is the name of a reaction, and each valued is a dictionary containing two keys: `reaction` and `regulation`. An example (all examples use JSON syntax):

    "R2a": {"reaction": {"B": -1, "ATP": 2, "NADH": 2, "C": 1}, "regulation": ["not", "RPb"]},

The `reaction` key contains a dictionary of molecule names to numbers representing how much of that molecule is consumed or produced in the reaction. So in the above example:

    {"B": -1, "ATP": 2, "NADH": 2, "C": 1}

This reaction consumes one `B` and produces two `ATP`, two `NADH` and one `C`. 

The `regulation` key provides a small expression that relates the state of the regulatory proteins to the activation of this reaction. The way these expressions work is that they are lists where the first element is a logical operation and the rest are arguments to that operation. So in the above case

    ["not", "RPb"]

`not` is the operation and `RPb` is the argument. This means this reaction is active when `RPb` is not active. The possible operations are 

* `is` - only active when the arguments are present (true)
* `not` - only active when the arguments are absent (false)
* `or` - performs a logical OR of two or more suboperations
* `and` - performs a logical AND of two or more suboperations

To use the `or` or `and` you can combine multiple conditions, like this:

    ["and", ["not", "RPb"], ["is", "RPO2"], ["not", "RPc1"]]

These nestings can be arbitrarily deep, though in practice they will probably not go beyond one or two levels. 

### regulation

Under the `regulation` key there is a dictionary of regulatory proteins to expressions governing their activation. These are much like the expressions above, but instead of being dependent on the state of regulatory proteins, they *define* the state of the regulatory proteins based on the values in `initial` representing the initial state of the system. For example:

    "RPO2": ["not", "Oxygen"]

says that `RPO2` is only active if there is no `Oxygen` in the system. 

In addition to the operations above, these can also be defined in terms of inequalities relative to state values, either fluxes or molecules:

    "RPh": [">", "Th", 0]

Means that `RPh` is active when the flux for the `Th` reaction is greater than zero.

### lower

This is a dictionary of reactions to their lower bounds. Any reaction left out of this dictionary is assumed to have a lower bound of `-Inf`.

    "R2a": -5,

### upper

Same as lower, but the upper bounds. Missing values are assumed to have an upper bound of `Inf`.

    "R3": 10.5,

### initial

The initial state of the molecules and fluxes in the system. These are read by the regulatory expressions to determine initial regulation:

    {"Oxygen": 1,
     "Carbon1": 1,
     "Th": 0,
     "R2b": 0}

## usage

Once you have a working network, you can solve the flux balance analysis problem directly with 

    julia> solution = Reactions.fluxbalanceanalysis(network, objective)
    Dict{String,Float64} with 10 entries:
      "BATP"  => 10.0
      "VthrB" => 10.0
      "BH2O"  => 10.0
      "Bthr"  => -10.0
      "VthrC" => 10.0
      "BADP"  => -10.0
      "Bhms"  => 10.0
      "Bphs"  => 0.0
      "Bphos" => -10.0
      "BH+"   => -20.0

The `network` argument is the network described above. The `objective` argument is a dictionary of reaction names to numbers representing the relative importance to the optimization problem. Any reactions omitted from this dictionary are assumed to be zero.

If you want to see what the stoichiometric matrix that represents the reactions is before solving the problem you can call `initialize` directly:

    julia> state = Reactions.initialize(network, Reactions.logicoperations)
    Dict{String,Array} with 6 entries:
      "active"        => ["BADP", "BATP", "BH+", "BH2O", "Bhms", "Bphos", "Bphs", "Bthr", "VthrB", "VthrC"]
      "molecules"     => ["ADP", "ATP", "H+", "H2O", "hms", "phos", "phs", "thr"]
      "lower"         => Real[-Inf, -Inf, -Inf, -Inf, 0, -Inf, 0, -Inf, -Inf, 0]
      "reactions"     => ["BADP", "BATP", "BH+", "BH2O", "Bhms", "Bphos", "Bphs", "Bthr", "VthrB", "VthrC"]
      "upper"         => Real[Inf, Inf, Inf, Inf, 10, Inf, 0, Inf, Inf, Inf]
      "stoichiometry" => [1 0 … 1 0; 0 1 … -1 0; … ; 0 0 … 1 -1; 0 0 … 0 1]

The matrix lives under the `stoichiometry` key:

    julia> state["stoichiometry"]
    8×10 Array{Int64,2}:
     1  0  0  0  0  0  0  0   1   0
     0  1  0  0  0  0  0  0  -1   0
     0  0  1  0  0  0  0  0   2   0
     0  0  0  1  0  0  0  0   0  -1
     0  0  0  0  1  0  0  0  -1   0
     0  0  0  0  0  1  0  0   0   1
     0  0  0  0  0  0  1  0   1  -1
     0  0  0  0  0  0  0  1   0   1

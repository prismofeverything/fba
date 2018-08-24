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

## usage


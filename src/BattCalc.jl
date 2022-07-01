module BattCalc

using Unitful, Parameters, LinearAlgebra, Measurements, PeriodicTable
export Cathode, F, R

include("Functions/Cathode.jl")
include("Functions/Anode.jl")
include("Functions/Stack.jl")

F,R = 96485.332123u"C/mol", 8.314462618

end

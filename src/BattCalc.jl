module BattCalc

using Unitful, Parameters, LinearAlgebra, Measurements, PeriodicTable
export Battery, Electrode, Parser, F, R, FormulaDict

include("Functions/Cathode.jl")
include("Functions/Anode.jl")
include("Functions/Stack.jl")

F,R = 96485.332123u"C/mol", 8.314462618

end

FormulaDict = Dict("Li"=>"lithium","Co"=>"cobalt", "O"=>"oxygen","Fe"=>"iron","Ni"=>"nickel","P"=>"phosphate",
                   "Mn"=>"manganese","Al"=>"aluminium","Si"=>"silicon","C" =>"carbon")

function Parser(List::String)
        
    rg1 = r"([[:alpha:]]+|[[:blank:]])"
    rg2 = r"([0-9]*[.]*[0-9]*[0-9]*[0-9])"

    chars = map(eachmatch(rg1, List)) do m
        m.match
    end

    nums = map(eachmatch(rg2, List)) do m
        parse(Float64, m.match)
    end
     
    return chars, nums
end


function Battery(Definition::String, Formula::String, Calc::String)

        chars, nums = Parser(Formula)
        CoatingThickness, Porousity, Area, ArealCap, Loading, Mass, Energy_Density, Capacity, NomVoltage, Specific_Cap = Electrode(chars, nums, Definition, Calc)

        return CoatingThickness, Porousity, Area, ArealCap, Loading, Mass, Energy_Density, Capacity, NomVoltage, Specific_Cap
end

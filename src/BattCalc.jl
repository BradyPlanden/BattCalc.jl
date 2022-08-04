module BattCalc


    using Unitful, Parameters, LinearAlgebra, Measurements, PeriodicTable
    export Pouch!, Anode, Cathode, Separator, Params
    export Electrode!, Parser, F, R, FormulaDict, SpecificCap, NomVoltage
    export MultiCell, PackStruct, Impedance, Electrolyte_

    include("BattCalcTypes.jl")


    F,R = 96485.332u"C/mol", 8.314463

    FormulaDict = Dict("Li"=>"lithium","Co"=>"cobalt", "O"=>"oxygen","Fe"=>"iron","Ni"=>"nickel","P"=>"phosphate",
                    "Mn"=>"manganese","Al"=>"aluminium","Si"=>"silicon","C" =>"carbon")

    SpecificCap = Dict("LCO"=> (175± 3.5)u"mA*hr/g", "NCM811" => (200±4.2)u"mA*hr/g", "NCM532" => (184±3.68)u"mA*hr/g", 
                    "LFP" => (157±3.14)u"mA*hr/g", "NCM622" => (190±4)u"mA*hr/g", "Graphite" => (265±5.1)u"mA*hr/g", "Li" => (3680±73.6)u"mA*hr/g")

    NomVoltage = Dict("LCO"=> (3.98±0.01)u"V", "NCM811" => (3.84±0.01)u"V", "NCM532" => (3.87±0.01)u"V", 
                    "LFP" => (3.37±0.01)u"V", "NCM622" => (3.85±0.01)u"V", "Graphite" => (0.17±0.01)u"V", "Li" => (0±0.0)u"V")

    CrystalDensity = Dict("NCM811"=> (4.95±0.1)u"g/cm^3", "Graphite" => (2.26±0.45)u"g/cm^3")




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


    function Electrode!(chars, nums, ElectrodeName, Calc, ElectrodeDef, Cell)
        """

        This function constructs the electrode structure. 

        """

        if ElectrodeDef == "Pos"
            Electrode = Cell.Pos
        elseif ElectrodeDef == "Neg"
            Electrode = Cell.Neg
        end

        if ElectrodeName == "Graphite"
            nums .= nums./maximum(nums)
        end

        #Sum Elements
        MolecularMass = 0.0u"g/mol"
        ComponentDensity = 0.0u"g/cm^3"

            if ElectrodeDef == "Pos" || ElectrodeName == "Graphite"
                for j in 1:length(chars)
                    MolecularMass += ustrip(elements[FormulaDict[chars[j]]].atomic_mass)*(u"g/mol")*nums[j]
                    #if chars[j] != "Li" && chars[j] != "O"
                    #    ComponentDensity += nums[j]*(elements[FormulaDict[chars[j]]].density)
                    #end
                end
            Electrode.Loading = (1-Electrode.Porousity)*(Electrode.CoatingThickness/1e4u"μm/cm")*(CrystalDensity[ElectrodeName]*1e3u"mg/g")
            elseif ElectrodeName == "Li"
                for j in 1:length(chars)
                    MolecularMass += ustrip(elements[FormulaDict[chars[j]]].atomic_mass)*(u"g/mol")*nums[j]
                    ComponentDensity += nums[j]*(elements[FormulaDict[chars[j]]].density)
                end
                Electrode.Loading = (1-Electrode.Porousity)*(Electrode.CoatingThickness/1e4u"μm/cm")*(ComponentDensity*1e3u"mg/g")
            end

        #Calculations
        
        Electrode.Area = Cell.Width*Cell.Length

            if ElectrodeDef == "Pos"
                Electrode.CCMass = (elements["aluminium"].density*1000u"mg/g"*Electrode.CollectorThickness/1e4u"μm/cm")*Electrode.Area
                Electrode.AMMass = uconvert(u"g",(Electrode.Area*Electrode.Loading)/1e3u"mg/g")
                Electrode.TotalMass = Electrode.CCMass+Electrode.AMMass
            else 
                Electrode.CCMass = (elements["copper"].density*1000u"mg/g"*Electrode.CollectorThickness/1e4u"μm/cm")*Electrode.Area
                Electrode.AMMass = uconvert(u"g",(Electrode.Area*Electrode.Loading)/1e3u"mg/g")
                Electrode.TotalMass = Electrode.CCMass+Electrode.AMMass
            end

            if Calc == "Theor"
                Electrode.ArealCap = (Electrode.Loading*uconvert(u"mA*hr/g",(elements[FormulaDict[chars[1]]].shells[end]*F/MolecularMass)))/1000u"mg/g"
            else
                Electrode.ArealCap = (Electrode.Loading*SpecificCap[ElectrodeName])/1000u"mg/g"
            end

        Electrode.Capacity = Electrode.Area*Electrode.ArealCap
        Electrode.Energy = uconvert(u"W*hr",(Electrode.Capacity/(1000u"mA/A")*NomVoltage[ElectrodeName]))
        Electrode.Ω = Impedance(Electrode) 

    end


    #kwargs: Number of Layers, area, composition, capacity, etc.
    function Pouch!(Cell, CathodeName::String,AnodeName::String, CathodeFormula::String, AnodeFormula::String, Calc::String, StorageMech::String, Layers::Int)
    
        #Cathode
        chars, nums = Parser(CathodeFormula)
        ElectrodeDef = "Pos"
        Electrode!(chars, nums, CathodeName, Calc, ElectrodeDef, Cell)

        #Anode 
        chars, nums = Parser(AnodeFormula)
        ElectrodeDef = "Neg"
        Electrode!(chars, nums, AnodeName, Calc, ElectrodeDef, Cell)

        #Cell Capacity
        if StorageMech == "Intercalation"
            Cell.Capacity = 2*min(Cell.Pos.Capacity,Cell.Neg.Capacity)*Layers #Scale capacities from double-coated collectors
        elseif StorageMech == "Deposition"
            Cell.Capacity = 2*Cell.Pos.Capacity*Layers
        end

        #Separator
        Cell.Sep.Area = Cell.Pos.Area * Layers

        #Electrolyte
        Cell.Electrolyte.ρₑ = 0.091u"(L*g)/(mol*cm^3)"*Cell.Electrolyte.cₑ+1.1u"g/cm^3" # Linear Correlation to Salt Concentration
        Cell.Electrolyte.Mass = Cell.Electrolyte.VolumeRatio*Cell.Capacity/1000u"mA*hr/A*hr"

        #Full Cell
        Cell.Thickness = uconvert(u"mm",(2*Cell.Pos.CoatingThickness+Cell.Pos.CollectorThickness+2*Cell.Sep.Thickness+2*Cell.Neg.CoatingThickness+Cell.Neg.CollectorThickness)*Layers+2*Cell.CasingThickness)
        Cell.Area = Cell.Pos.Area * Layers
        
        Negative_Tab = (Cell.CasingThickness*(1u"cm^2")*elements["copper"].density*1000u"mg/g")
        Positive_Tab = (Cell.CasingThickness*(1u"cm^2")*elements["aluminium"].density*1000u"mg/g")
        Cell.EChemMass = (2*Cell.Pos.AMMass+2*Cell.Neg.AMMass+2*uconvert(u"g",(Cell.Sep.Area*Cell.Sep.Loading*Cell.Sep.Thickness/1e4u"μm/cm")/1e3u"mg/g")+Cell.Pos.CCMass+Cell.Neg.CCMass)*Layers
        Cell.Mass = Cell.EChemMass + Negative_Tab + Positive_Tab + Cell.Electrolyte.Mass + 2 * (Cell.Pos.Area + Cell.Thickness*(Cell.Length + Cell.Width))* Cell.CasingThickness * elements["aluminium"].density*1000u"mg/g"

        Cell.Energy = uconvert(u"W*hr",(NomVoltage[CathodeName]-NomVoltage[AnodeName])*Cell.Capacity/1000u"mA/A")
        Cell.Energy_Density = uconvert(u"W*hr/kg",(Cell.Energy/Cell.Mass)*(1000u"g/kg"))
        Cell.VolDensity = uconvert(u"W*hr/L",(Cell.Energy/(Cell.Thickness/1e4u"μm/cm"*Cell.Pos.Area)))
        Cell.Ω = Cell.Pos.Ω/Layers + Cell.Neg.Ω/Layers

        return nothing

    end
    


    function MultiCell(Pack, N, M, Pₐ)
            # Energy_Density
            # NonActive_Mass

            α = N*M
            Pack.Capacity = Pack.Cell.Capacity*M
            Pack.Vₚ_Max = Pack.Cell.Voltage_Max*N
            Pack.Vₚ_Min = Pack.Cell.Voltage_Min*N
            Pack.Vₚ_Nom = Pack.Cell.Voltage_Nom*N
            Pack.Iₚ = uconvert(u"A",Pₐ/Pack.Vₚ_Nom)


            Pack.Mass = Pack.Cell.Mass*α
            Pack.Volume = Pack.Cell.Volume*α

            Ωₛ = Pack.Cell.Ω*N+(167e-5u"Ω"*N*2)
            Pack.Ω = 1/(1/Ωₛ*M)
            Pack.Ẇₕ = uconvert(u"kW",(Pack.Iₚ^2)*Pack.Ω)

            Pack.η = uconvert(NoUnits,(Pₐ-Pack.Ẇₕ)/(Pₐ))
            Pack.Energy = uconvert(u"kW*hr",Pack.Cell.Energy*α)
            Pack.Energy_Density = uconvert(u"W*hr/kg",Pack.Energy/Pack.Mass)
            Pack.Usable_Eₚ = uconvert(u"kW*hr",Pack.Energy*Pack.η)
            Usable_Eₖ = Pack.Usable_Eₚ/Pack.Mass
    end


    function Impedance(Electrode)

        Rp = uconvert(u"Ω",7.4*(Electrode.CoatingThickness/1e4u"μm/cm")/((Electrode.Porousity)*Electrode.Area*0.0089u"S/cm"))
        Rct = (119.98u"Ω*cm^2"*ustrip(Electrode.CoatingThickness)^(-0.775))/Electrode.Area

        if Rct/Rp <= 0.2 #Transport Limited
            L = sqrt(Rct*Rp)
        elseif Rct/Rp >= 0.66 # Kinetically Limited
            L = Rct + Rp/3
        else
            L = sqrt(Rct/Rp)*coth(1/sqrt(Rct/Rp))*Rp
        end

        return L

    end

        
end #module

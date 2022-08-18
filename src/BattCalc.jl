module BattCalc


    using Unitful, Parameters, LinearAlgebra, Measurements, PeriodicTable
    export Pouch!, Anode, Cathode, Separator, ElectrolyteStruct
    export Electrode!, Parser, F, R, FormulaDict, SpecificCap, NomVoltage
    export CellStruct, Module!, ModuleStruct, Pack!, PackStruct, Impedance

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
        Electrode.Ω = Impedance(Electrode, Cell.Storage) 

    end

    function Pouch!(Cell, CathodeName::String,AnodeName::String, CathodeFormula::String, AnodeFormula::String, Calc::String, Layers::Int)
    
        #Cathode
        chars, nums = Parser(CathodeFormula)
        ElectrodeDef = "Pos"
        Electrode!(chars, nums, CathodeName, Calc, ElectrodeDef, Cell)

        #Anode 
        chars, nums = Parser(AnodeFormula)
        ElectrodeDef = "Neg"
        Electrode!(chars, nums, AnodeName, Calc, ElectrodeDef, Cell)

        #Cell Capacity
        if Cell.Storage == "Intercalation"
            Cell.Capacity = 2*min(Cell.Pos.Capacity,Cell.Neg.Capacity)*Layers #Scale capacities from double-coated collectors
        elseif Cell.Storage == "Deposition"
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

    function Module!(Module, N, M, Pₐ)

            α = N * M
            Module.Capacity = Module.Cell.Capacity * M
            Module.Vₚ_Max = Module.Cell.Vₚ_Max * N
            Module.Vₚ_Min = Module.Cell.Vₚ_Min * N
            Module.Vₚ_Nom = Module.Cell.Vₚ_Nom * N   
            Module.Iₚ = uconvert(u"A",Pₐ/Module.Vₚ_Nom)

            Module.Mass = (Module.Cell.Mass * α) + Module.NonActive_Mass
            Module.Width = Module.Cell.Width * M
            Module.Length = Module.Cell.Thickness * N
            Module.Height = Module.Cell.Length
            Module.Volume = Module.Width * Module.Length * Module.Height

            Ωₛ = Module.Cell.Ω * N + (167e-5u"Ω" * N * 2)
            Module.Ω = 1 / (1 / Ωₛ * M)
            Module.Ẇₕ = uconvert(u"kW",(Module.Iₚ^2) * Module.Ω)

            Module.η = uconvert(NoUnits,(Pₐ-Module.Ẇₕ)/(Pₐ))
            Module.Energy = uconvert(u"kW*hr",Module.Cell.Energy * α)
            Module.Energy_Density = uconvert(u"W*hr/kg",Module.Energy/Module.Mass)
            Module.Usable_Eₚ = uconvert(u"kW*hr",Module.Energy * Module.η)
            Module.Usable_Density = uconvert(u"W*hr/kg",Module.Usable_Eₚ/Module.Mass)
            
            return nothing
    end

    function Pack!(Pack, N, M, Pₐ)

        α = N * M
        Pack.Capacity = Pack.Module.Capacity * M
        Pack.Vₚ_Max = Pack.Module.Vₚ_Max * N
        Pack.Vₚ_Min = Pack.Module.Vₚ_Min * N
        Pack.Vₚ_Nom = Pack.Module.Vₚ_Nom * N
        Pack.Iₚ = uconvert(u"A",Pₐ/Pack.Vₚ_Nom)

        Pack.Volume = uconvert(u"L",Pack.Module.Volume * α)
        Pack.Height = Pack.Module.Height
        Pack.Length = Pack.Module.Width * N
        Pack.Width = Pack.Module.Length
        Pack.CasingVolume = uconvert(u"L",(Pack.Height+Pack.CasingThickness*2)*(Pack.Width+Pack.CasingThickness*2)*(Pack.Length+Pack.CasingThickness*2) - Pack.Volume)
        Pack.CasingMass = Pack.CasingVolume * Pack.CasingDensity
        Pack.Mass = Pack.Module.Mass * α + Pack.CasingMass

        Ωₛ = Pack.Module.Ω*N+(167e-5u"Ω" * N * 2)
        Pack.Ω = 1 / (1 / Ωₛ * M)
        Pack.Ẇₕ = uconvert(u"kW",(Pack.Iₚ^2) * Pack.Ω)

        Pack.η = uconvert(NoUnits,(Pₐ-Pack.Ẇₕ)/(Pₐ))
        Pack.Energy = uconvert(u"kW*hr",Pack.Module.Energy * α)
        Pack.Energy_Density = uconvert(u"W*hr/kg",Pack.Energy/Pack.Mass)
        Pack.Usable_Eₚ = uconvert(u"kW*hr",Pack.Energy * Pack.η)
        Pack.Usable_Density = uconvert(u"W*hr/kg",Pack.Usable_Eₚ/Pack.Mass)
            
        return nothing
end


    function Impedance(Electrode, Storage)
        
        if Storage == "Deposition"
            Rp = uconvert(u"Ω",7.4*(Electrode.CoatingThickness/1e4u"μm/cm")/((Electrode.Porousity)*Electrode.Area*0.0089u"S/cm"))
            Rct = (119.98u"Ω*cm^2"*ustrip(Electrode.CoatingThickness)^(-0.775))/Electrode.Area
        else
            Rp = uconvert(u"Ω",7*(Electrode.CoatingThickness/1e4u"μm/cm")/((Electrode.Porousity)*Electrode.Area*0.015u"S/cm"))
            Rct = (50u"Ω*cm^2"*ustrip(Electrode.CoatingThickness)^(-0.775))/Electrode.Area
        end

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

module BattCalc


using Unitful, Parameters, LinearAlgebra, Measurements, PeriodicTable
export Battery, Electrode, Parser, F, R, FormulaDict, SpecificCap, NomVoltage
export Cell


F,R = 96485.332123u"C/mol", 8.314462618
FormulaDict = Dict("Li"=>"lithium","Co"=>"cobalt", "O"=>"oxygen","Fe"=>"iron","Ni"=>"nickel","P"=>"phosphate",
                   "Mn"=>"manganese","Al"=>"aluminium","Si"=>"silicon","C" =>"carbon")

SpecificCap = Dict("LCO"=> (175± 3.44)u"mA*hr/g", "NCM811" => (220±4.1)u"mA*hr/g", "NCM532" => (184±3.68)u"mA*hr/g", 
                    "LFP" => (157±3.14)u"mA*hr/g", "NCM622" => (190±4)u"mA*hr/g", "Graphite" => (285±4)u"mA*hr/g", "Li" => (3680±10)u"mA*hr/g")

NomVoltage = Dict("LCO"=> (3.98±0.01)u"V", "NCM811" => (3.84±0.01)u"V", "NCM532" => (3.87±0.01)u"V", 
                    "LFP" => (3.37±0.01)u"V", "NCM622" => (3.85±0.01)u"V", "Graphite" => (0.17±0.01)u"V", "Li" => (0±0.0)u"V")



@with_kw mutable struct Cathode
    Porousity = (0.61±0.0122)u"cm^3/cm^3"
    CoatingThickness = (45.0±0.9)u"μm"
    CollectorThickness = (13±0.26)u"μm"
    Area = (85±1.7)u"cm^2"
end

@with_kw mutable struct Anode
    Porousity = (0.2±0.004)u"cm^3/cm^3"
    CoatingThickness = (65±2)u"μm"
    CollectorThickness = (8±0.16)u"μm"
    Area = (85±1.7)u"cm^2"
end

@with_kw mutable struct Separator
    Porousity = (0.61±0.0122)u"cm^3/cm^3"
    Thickness = (14.0±0.28)u"μm"
    Area = (85±1.7)u"cm^2"
    Loading = (1.01±0.0202)u"mg/cm^3"
end


@with_kw mutable struct Params
    Neg::Anode
    Pos::Cathode
    Sep::Separator
end

Cell = Params(Anode(),Cathode(),Separator())


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



function Battery(CathodeName::String,AnodeName::String, CathodeFormula::String, AnodeFormula::String, Calc::String, StorageMech::String)
        
    
    #Cathode
    chars, nums = Parser(CathodeFormula)
    ElectrodeDef = "Pos"
    PosArealCap, PosLoading, PosMass, PosEnergy, PosCapacity = Electrode(chars, nums, CathodeName, Calc, ElectrodeDef, Cell)

    #Anode 
    chars, nums = Parser(AnodeFormula)
    ElectrodeDef = "Neg"
    NegArealCap, NegLoading, NegMass, NegEnergy, NegCapacity = Electrode(chars, nums, AnodeName, Calc, ElectrodeDef, Cell)


    #Full Cell
    TotalMass = PosMass+NegMass+2*uconvert(u"g",(Cell.Sep.Area*Cell.Sep.Loading*Cell.Sep.Thickness/1e4u"μm/cm")/1e3u"mg/g")
    if StorageMech == "Intercalation"
        TotaCapacity = min(PosCapacity,NegCapacity)
    elseif StorageMech == "Deposition"
        TotaCapacity = PosCapacity
    end
    TotalEnergy = uconvert(u"W*hr",(NomVoltage[CathodeName]-NomVoltage[AnodeName])*TotaCapacity/1000u"mA/A")
    TotalEnergyDensity = (2*TotalEnergy/TotalMass)*(1000u"g/kg")
    TotalThickness= (2*Cell.Pos.CoatingThickness+Cell.Pos.CollectorThickness+2*Cell.Sep.Thickness+2*Cell.Neg.CoatingThickness+Cell.Neg.CollectorThickness)
    TotalVolDensity = uconvert(u"W*hr/L",(2*TotalEnergy/(TotalThickness/1e4u"μm/cm"*Cell.Pos.Area)))

    return (PosArealCap, PosLoading, PosMass, PosEnergy, PosCapacity, NegArealCap, NegLoading, NegMass, NegEnergy, NegCapacity, TotalMass, TotaCapacity, TotalEnergy, TotalEnergyDensity, TotalThickness, TotalVolDensity)

end



function Electrode(chars, nums, ElectrodeName, Calc, ElectrodeDef, Cell)
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
    ComponentDensity = 0.0u"cm^3/g"
    if ElectrodeDef == "Pos" || ElectrodeName == "Graphite"
        for j in 1:length(chars)
            MolecularMass += ustrip(elements[FormulaDict[chars[j]]].atomic_mass)*(u"g/mol")*nums[j]
            if chars[j] != "Li" && chars[j] != "O"
                ComponentDensity += nums[j]/(elements[FormulaDict[chars[j]]].density)
            end
            
        end
    elseif ElectrodeName == "Li"
        for j in 1:length(chars)
            MolecularMass += ustrip(elements[FormulaDict[chars[j]]].atomic_mass)*(u"g/mol")*nums[j]
            ComponentDensity += nums[j]/(elements[FormulaDict[chars[j]]].density)
        end
    end


    #Calculations
    Loading = (1-Electrode.Porousity)*(Electrode.CoatingThickness/1e4u"μm/cm")/(ComponentDensity/1e3u"mg/g")#

    if ElectrodeDef == "Pos"
        Mass = uconvert(u"g",(2*Electrode.Area*Loading)/1e3u"mg/g"+(elements["aluminium"].density*1000u"mg/g"*Electrode.CollectorThickness/1e4u"μm/cm")*Electrode.Area)
    else 
        Mass = uconvert(u"g",(2*Electrode.Area*Loading)/1e3u"mg/g"+(elements["copper"].density*1000u"mg/g"*Electrode.CollectorThickness/1e4u"μm/cm")*Electrode.Area)
    end

    if Calc == "Theor"
        ArealCap = (Loading*uconvert(u"mA*hr/g",(elements[FormulaDict[chars[1]]].shells[end]*F/MolecularMass)))/1000u"mg/g"
    else
        ArealCap = (Loading*SpecificCap[ElectrodeName])/1000u"mg/g"
    end

    Capacity = Electrode.Area*ArealCap
    Energy = uconvert(u"W*hr",(Capacity/(1000u"mA/A")*NomVoltage[ElectrodeName]))

    return ArealCap, Loading, Mass, Energy, Capacity

end


    #Data:
    #Molecular 
    #Areal Capacity
    #Tortuosity
    #Lithium Concentration (max / min)
    #Loading 
    #Area
    #Specific Resistivity


    #To Do:
    #Energy Density (Nom. Voltage * Capcity / Mass)
    #Cathode Resistivity (Specific Capacity * Area)
    #Total Mass (Specific Loading * Area)
    #Total Energy (Nom. Voltage * Capacity)
    #Volumetric Energy Density (Nom. Voltage * Capacity / Volume)
    #Total Volume (Target Capacity / Areal Capacity)
    #Predicted Cycle Life (Coulombic efficiency until SOH target)
    #Max Power (Polarisation limit? Specific Resistance of Cathode)
    #Power Density (Power / Mass)

    #1Li:1Co:2.0O
    #1Li:0.6Ni:0.2Co:0.2Mn:2O

        
end #module

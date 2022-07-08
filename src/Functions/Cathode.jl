function Cathode(ActiveMaterial::String)
    """

    This function constructs the cathode structure. 

    """
    # For a given ActiveMaterial, construct the corresponding cathode structure from theoretical molecular values
    # To add: composition summation, i.e. NCM811 => 0.8Ni:0.1Co:Mn will automate to 0.8Ni:0.1Co:0.1Mn since (Ni+Co+Mn) = 1.0

    # Checks: ActiveMaterial is found in dict, Loading is within global limits, 


    #parsing input
    FormulaDict = Dict("Li"=>"lithium","Co"=>"cobalt", "O"=>"oxygen","Fe"=>"iron","Ni"=>"nickel","P"=>"phosphate",
                        "Mn"=>"manganese","Al"=>"aluminium","Si"=>"silicon","C" =>"carbon")
    rg1 = r"([[:alpha:]]+|[[:blank:]])"
    rg2 = r"([0-9]*[.]*[0-9]*[0-9]*[0-9])"

    chars = map(eachmatch(rg1, ActiveMaterial)) do m
        m.match
    end

    numbers = map(eachmatch(rg2, ActiveMaterial)) do m
        parse(Float64, m.match)
    end


    MolecularMass = 0.0u"g/mol"
    ComponentDensity = 0.0u"cm^3/g"
    for j in 1:length(chars)
        MolecularMass += ustrip(elements[FormulaDict[chars[j]]].atomic_mass)*(u"g/mol")*numbers[j]
        if chars[j] != "Li" && chars[j] != "O"
        ComponentDensity += numbers[j]/(elements[FormulaDict[chars[j]]].density)
        end
    end
    Porousity = (0.60±0.012)u"cm^3/cm^3"
    CoatingThickness = (50.0±2)u"μm"
    CollectorThickness = (13±0.26)u"μm"
    Loading = (1-Porousity)*(CoatingThickness/1e4u"μm/cm")/(ComponentDensity/1e3u"mg/g")#
    
    Resistivity = (15±0.75)u"Ω/cm^2"
    Area = (85±1.7)u"cm^2"
    Electrode_Mass = uconvert(u"g",(Area*Loading)/1e3u"mg/g"+(elements["aluminium"].density*1000u"mg/g"*CollectorThickness/1e4u"μm/cm")*Area)

    
    
    #Experimental
    LCO = 180u"mA*hr/g"
    Qt = (205±4.1)u"mA*hr/g" #NCM811
    NCM532 = 200u"mA*hr/g"
    LFP = 166u"mA*hr/g"

    LCO_V = 3.8
    NMC811_V = 3.7
    Electrode_V = (3.87±0.01)u"V" #NMC622_V
    LFP = 3.4


    #Theoretical
    #Qt = uconvert(u"mA*hr/g",(elements[FormulaDict[chars[1]]].shells[end]*F/MolecularMass))*0.7

    ArealCap = (Loading*Qt)/1000u"mg/g"
    Electrode_Capacity = Area*ArealCap
    Electrode_Energy_Density = uconvert(u"W*hr/kg",(Electrode_Capacity/(1000u"mA/A")*Electrode_V)/(Electrode_Mass/(1000u"g/kg")))


    return MolecularMass, Loading, CoatingThickness, ArealCap, Resistivity, Qt, Electrode_Capacity, Electrode_Mass, ComponentDensity, Electrode_Energy_Density


    #Data:
    #Molecular 
    #Areal Capacity
    #Tortuosity
    #Lithium Concentration (max / min)
    #Loading 
    #Area
    #Specific Resistivity


    #Calculations:
    #Energy Density (Nom. Voltage * Capcity / Mass)
    #Cathode Resistivity (Specific Capacity * Area)
    #Total Mass (Specific Loading * Area)
    #Total Energy (Nom. Voltage * Capacity)
    #Volumetric Energy Density (Nom. Voltage * Capacity / Volume)
    #Total Volume (Target Capacity / Areal Capacity)
    #Predicted Cycle Life (Coulombic efficiency until SOH target)
    #Max Power (Polarisation limit? Specific Resistance of Cathode)
    #Power Density (Power / Mass)


end


function Anode(Type::String, Loading, Capacity)
    """

    This function constructs the anode structure. 

    """


    if ActiveMaterial == "Graphite"
        #MolecularMass = Li+C+2*Ou"g/mol"
        Loading = (2±0.1)u"mg/cm^2"
        CoatingThickness = (40±2)u"μm"
        ArealCap = (3±0.15)u"mA*hr/cm^2"
        Resistivity = (15±0.75)u"Ω/cm^2"
    end

    return Loading, CoatingThickness, ArealCap, Resistivity

    #Data:
    #Molecular 
    #Areal Capacity
    #Tortuosity
    #Lithium Concentration (max / min)
    #Loading 
    #Area
    #Specific Resistivity


    #Calculations:
    #Energy Density (Nom. Voltage * Capcity / Mass)
    #Cathode Resistivity (Specific Capacity * Area)
    #Total Mass (Specific Loading * Area)
    #Total Energy (Nom. Voltage * Capacity)
    #Volumetric Energy Density (Nom. Voltage * Capacity / Volume)
    #Total Volume (Target Capacity / Areal Capacity)
    #Predicted Cycle Life (Coulombic efficiency until SOH target)
    #Max Power (Polarisation limit? Specific Resistance of Cathode)
    #Power Density (Power / Mass)


end
    
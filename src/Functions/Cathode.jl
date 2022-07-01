function Cathode(ActiveMaterial::String)
    """

    This function constructs the cathode structure. 

    """
    # For a given ActiveMaterial, construct the corresponding cathode structure from theoretical molecular values
    # To add: composition summation, i.e. NCM811 => 0.8Ni:0.1Co:Mn will automate to 0.8Ni:0.1Co:0.1Mn since (Ni+Co+Mn) = 1.0

    # Checks: ActiveMaterial is found in dict, Loading is within global limits, 


    #parsing input
    FormulaDict = Dict("Li"=>"lithium","Co"=>"cobalt", "O"=>"oxygen","Fe"=>"iron","Ni"=>"nickel","P"=>"phosphate","Mn"=>"manganese","Al"=>"aluminium", "Si"=>"silicon")
    rg1 = r"([[:alpha:]]+|[[:blank:]])"
    rg2 = r"([0-9]*[.]*[0-9]*[0-9]*[0-9])"

    chars = map(eachmatch(rg1, ActiveMaterial)) do m
        m.match
    end

    numbers = map(eachmatch(rg2, ActiveMaterial)) do m
        parse(Float64, m.match)
    end


    MolecularMass = 0.0u"g/mol"
    for j in 1:length(chars)
        MolecularMass += ustrip(elements[FormulaDict[chars[j]]].atomic_mass)*(u"g/mol")*numbers[j]
    end

    Loading = (2±0.1)u"mg/cm^2"
    Thickness = (40±2)u"μm"
    ArealCap = (3±0.15)u"mA*hr/cm^2"
    Resistivity = (15±0.75)u"Ω/cm^2"
    
    
    Qt = uconvert(u"mA*hr/g",(elements[FormulaDict[chars[1]]].shells[end]*F/MolecularMass))

    return MolecularMass, Loading, Thickness, ArealCap, Resistivity, Qt


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
        Thickness = (40±2)u"μm"
        ArealCap = (3±0.15)u"mA*hr/cm^2"
        Resistivity = (15±0.75)u"Ω/cm^2"
    end

    return Loading, Thickness, ArealCap, Resistivity

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
    
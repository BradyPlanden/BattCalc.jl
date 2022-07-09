function Electrode(chars, nums, Definition, Calc)
    """

    This function constructs the cathode structure. 

    """

    #Experimental Data
    Specific_Cap = Dict("LCO"=> (172± 3.44)u"mA*hr/g", "NCM811" => (222±4.1)u"mA*hr/g", "NCM532" => (184±3.68)u"mA*hr/g", "LFP" => (157±3.14)u"mA*hr/g")
    NomVoltage = Dict("LCO"=> (3.98±0.01)u"V", "NCM811" => (3.84±0.01)u"V", "NCM532" => (3.87±0.01)u"V", "LFP" => (3.37±0.01)u"V")
    Porousity = (0.60±0.012)u"cm^3/cm^3"
    CoatingThickness = (50.0±2)u"μm"
    CollectorThickness = (13±0.26)u"μm"
    Area = (85±1.7)u"cm^2"
    Resistivity = (15±0.75)u"Ω/cm^2"


    MolecularMass = 0.0u"g/mol"
    ComponentDensity = 0.0u"cm^3/g"
    for j in 1:length(chars)
        MolecularMass += ustrip(elements[FormulaDict[chars[j]]].atomic_mass)*(u"g/mol")*nums[j]
        if chars[j] != "Li" && chars[j] != "O"
        ComponentDensity += nums[j]/(elements[FormulaDict[chars[j]]].density)
        end
    end


    #Theoretical Specific Capacity
    if Calc == "Theor"
        Specific_Cap = uconvert(u"mA*hr/g",(elements[FormulaDict[chars[1]]].shells[end]*F/MolecularMass))*0.7
    end

    Loading = (1-Porousity)*(CoatingThickness/1e4u"μm/cm")/(ComponentDensity/1e3u"mg/g")#
    Mass = uconvert(u"g",(Area*Loading)/1e3u"mg/g"+(elements["aluminium"].density*1000u"mg/g"*CollectorThickness/1e4u"μm/cm")*Area)
    ArealCap = (Loading*Specific_Cap[Definition])/1000u"mg/g"
    Capacity = Area*ArealCap
    Energy_Density = uconvert(u"W*hr/kg",(Capacity/(1000u"mA/A")*NomVoltage[Definition])/(Mass/(1000u"g/kg")))


    return CoatingThickness, Porousity, Area, ArealCap, Loading, Mass, Energy_Density, Capacity, NomVoltage[Definition], Specific_Cap[Definition]

end




    # For a given ActiveMaterial, construct the corresponding cathode structure from theoretical molecular values
    # To add: composition summation, i.e. NCM811 => 0.8Ni:0.1Co:Mn will automate to 0.8Ni:0.1Co:0.1Mn since (Ni+Co+Mn) = 1.0

    # Checks: ActiveMaterial is found in dict, Loading is within global limits, 


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
        
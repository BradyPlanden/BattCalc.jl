using Parameters

@with_kw mutable struct Cathode
    Porousity = (0.61±0.0122)u"cm^3/cm^3"
    CoatingThickness = (45.0±0.9)u"μm"
    CollectorThickness = (13±0.26)u"μm"
    Area = (85±1.7)u"cm^2"
    Loading = (0.0±0.0)u"mg/cm^2"
    TotalMass = (0.0±0.0)u"g"
    ArealCap = (0.0±0.0)u"mA*hr/cm^2"
    Capacity =  (0.0±0.0)u"mA*hr"
    Energy = (0.0±0.0)u"W*hr"
    CCMass = (0.0±0.0)u"g"
    AMMass = (0.0±0.0)u"g"
end

@with_kw mutable struct Anode
    Porousity = (0.2±0.004)u"cm^3/cm^3"
    CoatingThickness = (65±2)u"μm"
    CollectorThickness = (8±0.16)u"μm"
    Area = (85±1.7)u"cm^2"
    Loading = (0.0±0.0)u"mg/cm^2"
    TotalMass = (0.0±0.0)u"g"
    ArealCap = (0.0±0.0)u"mA*hr/cm^2"
    Capacity =  (0.0±0.0)u"mA*hr"
    Energy = (0.0±0.0)u"W*hr"
    CCMass = (0.0±0.0)u"g"
    AMMass = (0.0±0.0)u"g"
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
    Layers::Int64 = 20
    Energy = (0±0)u"hr*W"
    Voltage_Max = (4.2±0.0)u"V"
    Voltage_Min = (3.0±0.0)u"V"
    Voltage_Nom = (3.8±0.0)u"V"
    Capacity = (5±0.00)u"A*hr"
    Ωᵢ = (0.005±0.00)u"Ω"
    Ωᵪ = (0.005±0.00)u"Ω"
    Mass = (50±0.0)u"g"
    Volume = (0±0)u"L"
    Energy_Density = (0±0)u"W*hr/kg"
    Thickness = (0±0)u"W*hr/kg"
    VolDensity = (0±0)u"W*hr/L"
    Area = (0±0)u"cm^2"
end

@with_kw mutable struct PackStruct
    Cell::Params
    Series::Int64 = 10
    Parallel::Int64 = 5
    Energy = (0±0)u"hr*W"
    Vₚ_Max = (0.0±0.0)u"V"
    Vₚ_Min = (0.0±0.0)u"V"
    Vₚ_Nom = (0.0±0.00)u"V"
    Capacity = (0.0±0.00)u"A*hr"
    Ω = (0.00±0.000)u"Ω"
    Mass = (0±0.0)u"g"
    Volume = (0±0)u"L"
    Energy_Density = (0±0)u"W*hr/kg"
    Thickness = (0±0)u"W*hr/kg"
    VolDensity = (0±0)u"W*hr/L"
    Ẇₕ = (0±0)u"W"
    η = (0±0)
    Usable_Eₚ = (0±0)u"kW*hr"
    Iₚ = (0±0)u"A"
end
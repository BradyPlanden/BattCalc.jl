using Parameters

@with_kw mutable struct Cathode
    Porousity = (0.61±0.0122)u"cm^3/cm^3"
    CoatingThickness = (45.0±0.9)u"μm"
    CollectorThickness = (13±0.26)u"μm"
    Area = (0.0±0.0)u"cm^2"
    Loading = (0.0±0.0)u"mg/cm^2"
    TotalMass = (0.0±0.0)u"g"
    ArealCap = (0.0±0.0)u"mA*hr/cm^2"
    Capacity =  (0.0±0.0)u"mA*hr"
    Energy = (0.0±0.0)u"W*hr"
    CCMass = (0.0±0.0)u"g"
    AMMass = (0.0±0.0)u"g"
    Ω = (0.0±0.0)u"Ω"
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
    Ω = (0.0±0.0)u"Ω"
end

@with_kw mutable struct Separator
    Porousity = (0.61±0.0122)u"cm^3/cm^3"
    Thickness = (14.0±0.28)u"μm"
    Area = (85±1.7)u"cm^2"
    Loading = (1.01±0.0202)u"mg/cm^3"
end

@with_kw mutable struct ElectrolyteStruct
    cₑ = (1.1±0.022)u"mol/L"
    ρₑ = (0.0±0.0)u"g/cm^3"
    VolumeRatio = (1.5±0.054)u"g/A*hr"
    Mass = (0±0.0)u"g"
end

@with_kw mutable struct CellStruct
    Neg::Anode
    Pos::Cathode
    Sep::Separator
    Electrolyte::ElectrolyteStruct
    Storage = "Deposition"
    Layers::Int64 = 20
    Energy = (0±0)u"hr*W"
    Vₚ_Max = (4.2±0.0)u"V"
    Vₚ_Min = (3.0±0.0)u"V"
    Vₚ_Nom = (3.8±0.0)u"V"
    Capacity = (5±0.00)u"A*hr"
    Ω = (0.005±0.00)u"Ω"
    EChemMass = (50±0.0)u"g"
    Mass = (50±0.0)u"g"
    Volume = (0±0)u"L"
    Energy_Density = (0±0)u"W*hr/kg"
    Thickness = (0±0)u"W*hr/kg"
    VolDensity = (0±0)u"W*hr/L"
    Area = (0±0)u"cm^2"
    Width = (4±0.08)u"cm"
    Length = (12.5±0.25)u"cm"
    CasingThickness = (0.021±0.00042)u"cm"
end

@with_kw mutable struct ModuleStruct
    Cell::CellStruct
    Series::Int64 = 20
    Parallel::Int64 = 2
    Energy = (0±0)u"hr*W"
    Vₚ_Max = (0.0±0.0)u"V"
    Vₚ_Min = (0.0±0.0)u"V"
    Vₚ_Nom = (0.0±0.00)u"V"
    Capacity = (0.0±0.00)u"A*hr"
    Ω = (0.00±0.000)u"Ω"
    NonActive_Mass = (1.0±0.02)u"kg"
    Mass = (0±0.0)u"g"
    Width = (0.0±0.0)u"m"
    Length = (0.0±0.0)u"m"
    Height = (0.0±0.0)u"m"
    Volume = (0±0)u"L"
    Energy_Density = (0±0)u"W*hr/kg"
    Thickness = (0±0)u"W*hr/kg"
    VolDensity = (0±0)u"W*hr/L"
    Ẇₕ = (0±0)u"W"
    η = (0±0)
    Usable_Eₚ = (0±0)u"kW*hr"
    Usable_Density = (0±0)u"W*hr/kg"
    Iₚ = (0±0)u"A"
end

@with_kw mutable struct PackStruct
    Module::ModuleStruct
    Series::Int64 = 6
    Parallel::Int64 = 1
    Energy = (0±0)u"hr*W"
    Vₚ_Max = (0.0±0.0)u"V"
    Vₚ_Min = (0.0±0.0)u"V"
    Vₚ_Nom = (0.0±0.00)u"V"
    Capacity = (0.0±0.00)u"A*hr"
    Ω = (0.00±0.000)u"Ω"
    Mass = (0±0.0)u"kg"
    Width = (0.0±0.0)u"m"
    Length = (0.0±0.0)u"m"
    Height = (0.0±0.0)u"m"
    Volume = (0±0)u"L"
    Energy_Density = (0±0)u"W*hr/kg"
    Thickness = (0±0)u"W*hr/kg"
    VolDensity = (0±0)u"W*hr/L"
    CasingVolume = (0±0)u"L"
    CasingDensity = (4.0±0.08)u"kg/L"
    CasingThickness = (7.24±0.1448)u"mm"
    CasingMass = (0±0.0)u"kg"
    Ẇₕ = (0±0)u"W"
    η = (0±0)
    Usable_Eₚ = (0±0)u"kW*hr"
    Usable_Density = (0±0)u"W*hr/kg"
    Iₚ = (0±0)u"A"
end
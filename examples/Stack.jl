using BattCalc, Measurements, Unitful, Plots, UnitfulRecipes

Pack = PackStruct(Cell=Params(Neg=Anode(),Pos=Cathode(),Sep=Separator()))
Stacks = 10
Pack.Cell.Pos.Porousity = (0.61±0.0)u"cm^3/cm^3"
Pack.Cell.Pos.CoatingThickness = (45.0±0.0)u"μm"
Pack.Cell.Pos.CollectorThickness = (13±0.0)u"μm"
Pack.Cell.Pos.Area = (85±0.0)u"cm^2"

Pack.Cell.Neg.Porousity = (0.2±0.00)u"cm^3/cm^3"
Pack.Cell.Neg.CoatingThickness = (65±0.0)u"μm"
Pack.Cell.Neg.CollectorThickness = (8±0.0)u"μm"
Pack.Cell.Neg.Area = (85±0.0)u"cm^2"

Pack.Cell.Sep.Porousity = (0.61±0.0)u"cm^3/cm^3"
Pack.Cell.Sep.Thickness = (14.0±0.0)u"μm"
Pack.Cell.Sep.Area = (85±0.0)u"cm^2"
Pack.Cell.Sep.Loading = (1.01±0.0)u"mg/cm^3"

function plotsf()
    Rng = 0.01:0.01:0.99
    Thickness = (ones(length(Rng)^2,size(Rng,2)))*u"μm"
    Density =  (ones(length(Rng)^2,size(Rng,2)).±0.00)*u"W*hr/kg"
    Impedance_lp =  (ones(length(Rng)^2,size(Rng,2)).±0.00)*u"Ω"
    Porousity =  (ones(length(Rng)^2,size(Rng,2)).±0.00)*u"cm^3/cm^3"
    k=1
    for i ∈ Rng
        for j ∈ Rng
            Pack.Cell.Pos.Porousity = (i±0.02*i)u"cm^3/cm^3"
            Pack.Cell.Neg.Porousity = (j±0.02*j)u"cm^3/cm^3"
            Pouch!(Pack.Cell,"NCM811","Graphite","1Li:0.8Ni:0.1Co:0.1Mn:2O","1.0Li6.0C", "Exper", "Intercalation", Stacks)
            MultiCell(Pack,100,2,40u"kW")
            Porousity[k] =  Measurements.value(Pack.Cell.Pos.Porousity)
            #Thickness[k] =  Measurements.value(Pack.Cell.Pos.CoatingThickness)
            Density[k] = Pack.Energy_Density
            Impedance_lp[k] = Pack.Cell.Ω
            k+=1
        end
    end
    return Thickness, Density, Impedance_lp, Porousity
end

Thickness, Density, Impedance_lp, Porousity = plotsf()

# nominal, ± σ:
nominal = Measurements.value.(Density)
nom_plus_std = nominal .- Measurements.uncertainty.(Density)
nom_minus_std = nominal .+ Measurements.uncertainty.(Density)

plot(Porousity, nom_minus_std, fillrange=nom_plus_std, fillalpha=0.2, linealpha = 0.3, legend = false,size=(1280,720))
plot!(Porousity,nominal)







# #Intercalation
# @with_kw mutable struct Cathode
#     Porousity = (0.61±0.0122)u"cm^3/cm^3"
#     CoatingThickness = (45.0±0.9)u"μm"
#     CollectorThickness = (13±0.26)u"μm"
#     Area = (85±1.7)u"cm^2"
# end

# @with_kw mutable struct Anode
#     Porousity = (0.2±0.004)u"cm^3/cm^3"
#     CoatingThickness = (65±2)u"μm"
#     CollectorThickness = (8±0.16)u"μm"
#     Area = (85±1.7)u"cm^2"
# end

# @with_kw mutable struct Separator
#     Porousity = (0.61±0.0122)u"cm^3/cm^3"
#     Thickness = (14.0±0.28)u"μm"
#     Area = (85±1.7)u"cm^2"
#     Loading = (1.01±0.0202)u"mg/cm^3"
# end



# #Deposition - AF
# @with_kw mutable struct Cathode
#     Porousity = (0.61±0.0122)u"cm^3/cm^3"
#     CoatingThickness = (45.0±0.9)u"μm"
#     CollectorThickness = (13±0.26)u"μm"
#     Area = (85±1.7)u"cm^2"
# end

# @with_kw mutable struct Anode
#     Porousity = (0.0±0.0)u"cm^3/cm^3"
#     CoatingThickness = (17.0±0.34)u"μm"
#     CollectorThickness = (8±0.16)u"μm"
#     Area = (85±1.7)u"cm^2"
# end

# @with_kw mutable struct Separator
#     Porousity = (0.61±0.0122)u"cm^3/cm^3"
#     Thickness = (14.0±0.28)u"μm"
#     Area = (85±1.7)u"cm^2"
#     Loading = (1.01±0.0202)u"mg/cm^3"
# end


# #Deposition - Excess Lithium
# @with_kw mutable struct Cathode
#     Porousity = (0.61±0.0122)u"cm^3/cm^3"
#     CoatingThickness = (45.0±0.9)u"μm"
#     CollectorThickness = (13±0.26)u"μm"
#     Area = (85±1.7)u"cm^2"
# end

# @with_kw mutable struct Anode
#     Porousity = (0.0±0.0)u"cm^3/cm^3"
#     CoatingThickness = (34.0±0.68)u"μm"
#     CollectorThickness = (8±0.16)u"μm"
#     Area = (85±1.7)u"cm^2"
# end

# @with_kw mutable struct Separator
#     Porousity = (0.61±0.0122)u"cm^3/cm^3"
#     Thickness = (14.0±0.28)u"μm"
#     Area = (85±1.7)u"cm^2"
#     Loading = (1.01±0.0202)u"mg/cm^3"
# end
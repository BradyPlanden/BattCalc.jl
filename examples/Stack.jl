using BattCalc, Measurements, Unitful, Plots, UnitfulRecipes

# Setup Simulation
Pack = PackStruct(Cell=Params(Neg=Anode(),Pos=Cathode(),Sep=Separator()))
Stacks = 10
Rng = 0.01:0.02:0.99
Pack.Cell.Pos.Porousity = (0.61±0.0)u"cm^3/cm^3"
Pack.Cell.Pos.CoatingThickness = (45.0±0.0)u"μm"
Pack.Cell.Pos.CollectorThickness = (13±0.0)u"μm"
Pack.Cell.Pos.Area = (85±0.0)u"cm^2"

Pack.Cell.Neg.Porousity = (0.2±0.00)u"cm^3/cm^3"
Pack.Cell.Neg.CoatingThickness = (205±0.0)u"μm"
Pack.Cell.Neg.CollectorThickness = (8±0.0)u"μm"
Pack.Cell.Neg.Area = (85±0.0)u"cm^2"

Pack.Cell.Sep.Porousity = (0.61±0.0)u"cm^3/cm^3"
Pack.Cell.Sep.Thickness = (14.0±0.0)u"μm"
Pack.Cell.Sep.Area = (85±0.0)u"cm^2"
Pack.Cell.Sep.Loading = (1.01±0.0)u"mg/cm^3"

function Sweep(Rng)
    Thickness = (ones(length(Rng),size(Rng,2)))*u"μm"
    Density =  (ones(length(Rng),size(Rng,2)).±0.00)*u"W*hr/kg"
    Impedance_lp =  (ones(length(Rng),size(Rng,2)).±0.00)*u"mΩ"
    Porousity =  (ones(length(Rng),size(Rng,2)).±0.00)*u"cm^3/cm^3"
    Loading = (ones(length(Rng),size(Rng,2)).±0.00)*u"mg/cm^2"
    Areal_Cap = (ones(length(Rng),size(Rng,2)).±0.00)*u"mA*hr/cm^2"
    k=1

    for i ∈ Rng
        Pack.Cell.Pos.Porousity = (i±0.02*i)u"cm^3/cm^3"
        Pack.Cell.Neg.Porousity = (i±0.02*i)u"cm^3/cm^3"
        Pouch!(Pack.Cell,"NCM811","Graphite","1Li:0.8Ni:0.1Co:0.1Mn:2O","1.0Li6.0C", "Exper", "Intercalation", Stacks)
        MultiCell(Pack,100,2,40u"kW")
        Porousity[k] =  Measurements.value(Pack.Cell.Pos.Porousity)
        Density[k] = Pack.Energy_Density
        Impedance_lp[k] = Pack.Cell.Ω
        Loading[k] = Pack.Cell.Pos.Loading
        Areal_Cap[k] = Pack.Cell.Pos.ArealCap
        k+=1
    end

    return Thickness, Density, uconvert.(u"mΩ",Impedance_lp), Porousity, Loading, Areal_Cap

end

Thickness, Density, Impedance_lp, Porousity, Loading, Areal_Cap = Sweep(Rng)

# strip, nominal, ± σ:
Densityw = ustrip(Measurements.value.(Density))
nominal = Measurements.value.(Density)
Areal_Capw = ustrip(Measurements.value.(Areal_Cap))
Impedance_lpw = ustrip(Measurements.value.(Impedance_lp))
nom_plus_std = nominal .- Measurements.uncertainty.(Density)
nom_minus_std = nominal .+ Measurements.uncertainty.(Density)


gr()
plot(Measurements.value.(Areal_Cap), 
    Measurements.value.(Impedance_lp), 
    axis=:l, 
    y_foreground_color_text=:blue,
    y_guidefontcolor = :blue,
    color=:blue,
    bottom_margin=5Plots.mm, 
    left_margin = 5Plots.mm, 
    right_margin = 25Plots.mm, 
    ylabel = "Impedance", 
    xlabel = "Positive Electrode Areal Capacity",
    size=(1280,720),
    framestyle=:box,
    linestyle = :dash,
    ylims = (0,125)
    )

plot!(Measurements.value.(Areal_Cap), 
     Measurements.value.(nominal), 
     axis=:r, 
     color=:black,
     bottom_margin=5Plots.mm, 
     left_margin = 5Plots.mm, 
     right_margin = 25Plots.mm, 
     ylabel = "Energy Density", 
     xlabel = "Positive Electrode Areal Capacity",
     size=(1280,720),
     framestyle=:box,
     #markershape = :auto
     )
     plot!(twinx(),
          Measurements.value.(Areal_Cap), 
          Measurements.value.(Impedance_lp),
          label = "")
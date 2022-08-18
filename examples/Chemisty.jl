using BattCalc, Measurements, Unitful, Plots

Pack = PackStruct(Module=ModuleStruct(Cell=CellStruct(Neg=Anode(),Pos=Cathode(),Sep=Separator(), Electrolyte=ElectrolyteStruct())))
Stacks = 19

#Pack.Module.Cell.Storage = "Intercalation" 
Pack.Module.Cell.Storage = "Deposition"

Pack.Module.Cell.Length = (12.5±0.25)u"cm"
Pack.Module.Cell.Width = (4±0.08)u"cm"
Pack.Module.Cell.Pos.Porousity = (0.335±0.0122)u"cm^3/cm^3"
Pack.Module.Cell.Pos.CoatingThickness = (45.0±0.9)u"μm"

    if Pack.Module.Cell.Storage == "Deposition"
        Pack.Module.Cell.Neg.Porousity = (0.02±0.0004)u"cm^3/cm^3"
        Pack.Module.Cell.Neg.CoatingThickness = (15.5±0.31)u"μm" #(35.5±0.71)
        Pack.Module.Cell.Electrolyte.VolumeRatio = (2±0.04)u"g/A*hr"
        Pouch!(Pack.Module.Cell,"NCM811","Li","0.833Li:0.8Ni:0.08Co:0.08Mn:0.04Al:2O","1.0Li", "Exper", Stacks)
        Module!(Pack.Module, 24, 2, (80/6)u"kW")
        Pack!(Pack, 6, 1, 80u"kW")
    else
        Pack.Module.Cell.Neg.Porousity = (0.25±0.004)u"cm^3/cm^3"
        Pack.Module.Cell.Neg.CoatingThickness = (72.5±1.3)u"μm"
        Pack.Module.Cell.Electrolyte.VolumeRatio = (2±0.04)u"g/A*hr"
        Pouch!(Pack.Module.Cell,"NCM811","Graphite","0.833Li:0.8Ni:0.08Co:0.08Mn:0.04Al:2O","0.9Li0.1Si6.0C", "Exper", Stacks)
        Module!(Pack.Module, 24, 2, (80/6)u"kW")
        Pack!(Pack, 6, 1, 80u"kW")
    end

   
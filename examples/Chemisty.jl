using BattCalc, Measurements, Unitful, Plots

Pack = PackStruct(Cell=Params(Neg=Anode(),Pos=Cathode(),Sep=Separator(), Electrolyte=Electrolyte_()))
Stacks = 17

Storage = "Intercalation" 
#Storage = "Deposition"

Pack.Cell.Length = (12.5±0.25)u"cm"
Pack.Cell.Width = (4±0.08)u"cm"
Pack.Cell.Pos.Porousity = (0.335±0.0122)u"cm^3/cm^3"
Pack.Cell.Pos.CoatingThickness = (45.0±0.9)u"μm"

    if Storage == "Deposition"
        Pack.Cell.Neg.Porousity = (0.02±0.0004)u"cm^3/cm^3"
        Pack.Cell.Neg.CoatingThickness = (31.5±0.62)u"μm"
        Pack.Cell.Electrolyte.VolumeRatio = (2±0.04)u"g/A*hr"
        Pouch!(Pack.Cell,"NCM811","Li","0.833Li:0.8Ni:0.08Co:0.08Mn:0.04Al:2O","1.0Li", "Exper", Storage, Stacks)
    else
        Pack.Cell.Neg.Porousity = (0.25±0.004)u"cm^3/cm^3"
        Pack.Cell.Neg.CoatingThickness = (72.5±1.3)u"μm"
        Pack.Cell.Electrolyte.VolumeRatio = (2±0.04)u"g/A*hr"
        Pouch!(Pack.Cell,"NCM811","Graphite","0.833Li:0.8Ni:0.08Co:0.08Mn:0.04Al:2O","0.9Li0.1Si6.0C", "Exper", Storage, Stacks)
    end

   
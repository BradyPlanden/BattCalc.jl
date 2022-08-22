using BattCalc, Measurements, Unitful

Pack = PackStruct(Module=ModuleStruct(Cell=CellStruct(Neg=Anode(),Pos=Cathode(),Sep=Separator(), Electrolyte=ElectrolyteStruct())))
Cell_Layers = 12

#Set Cell Chemistry
Pack.Module.Cell.Pos.Name = "NCM811"
Pack.Module.Cell.Pos.Formula = "0.833Li:0.8Ni:0.08Co:0.08Mn:0.04Al:2O"
Pack.Module.Cell.Neg.Name = "Li"
Pack.Module.Cell.Neg.Formula = "1.0Li"

#Select Cell Storage Mechanism
#Pack.Module.Cell.Storage = "Intercalation" 
Pack.Module.Cell.Storage = "Deposition"

#Set Cell Dimensions and Overwrite Defaults
Pack.Module.Cell.Length = (12.5±0.25)u"cm"
Pack.Module.Cell.Width = (4±0.08)u"cm"
Pack.Module.Cell.Pos.Porousity = (0.335±0.0067)u"cm^3/cm^3"
Pack.Module.Cell.Pos.CoatingThickness = (72.5±1.45)u"μm"

#Form Pack Depending on Storage Mechanism
if Pack.Module.Cell.Storage == "Deposition"
    Pack.Module.Cell.Neg.Porousity = (0.02±0.0004)u"cm^3/cm^3"
    Pack.Module.Cell.Neg.CoatingThickness = (43.5±0.87)u"μm" #(23.5±0.47)(43.5±0.94)
    Pack.Module.Cell.Electrolyte.VolumeRatio = (1.6±0.032)u"g/A*hr"
    Pack.Module.Cell.Vₚ_Offset = (0.09)u"V" #Nom Offset for high-voltage capability
    Pack.Module.Cell.Vₚ_Max = (4.5)u"V"

    #Pack Formation
    Pack!(Pack, 6, 1, 22, 2, 80u"kW", "Exper", Cell_Layers)

else
    Pack.Module.Cell.Neg.Porousity = (0.25±0.005)u"cm^3/cm^3"
    Pack.Module.Cell.Neg.CoatingThickness = (85.2±1.704)u"μm"
    Pack.Module.Cell.Sep.Thickness = (12±0.24)u"μm"
    Pack.Module.Cell.Electrolyte.VolumeRatio = (1.6±0.032)u"g/A*hr"
    Pack.Module.Cell.Vₚ_Offset = (0.0)u"V" #Nom Offset for high-voltage capability

    #Pack Formation
    Pack!(Pack, 6, 1, 24, 2, 80u"kW", "Exper", Cell_Layers)

end
OpenDatabase("buoyancy9000.bov")
AddPlot("Contour", "buoyancy")
#AddOperator("Contour");
#cnt = ContourAttributes()
#cnt.contourPercent = 50
#SetOperatorOptions(cnt)
DrawPlots()
e = ExportDBAttributes()
e.db_type = "VTK"
e.filename = "abhishek"
e.variables = ("hardyglobal")
ExportDatabase(e)

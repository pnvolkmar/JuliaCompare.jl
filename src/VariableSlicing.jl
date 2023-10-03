using Revise
import JuliaCompare as J
import JuliaCompare: db_files, Canada
import PromulaDBA as P
using DataFrames, DataFramesMeta

J.greet()

CODE_FOLDER = raw"\\Pink\c\2020CanadaCatalpa\Engine"
DATA_FOLDER = raw"\\Pink\c\2020CanadaCatalpa\2020Model\Ref22"
DATA_FOLDER2 = raw"C:\2020CanadaCatalpa\2020Model\Ref22_AgCH4_170"

vars = J.list_vars(CODE_FOLDER, DATA_FOLDER, db_files);
loc1 = J.Loc(vars, DATA_FOLDER);
loc2 = J.Loc(vars, DATA_FOLDER2);

totpol_ref = J.var("TotPol", loc1);
totpol_ag = J.var("TotPol", loc2);

totpol_diff = J.diff(totpol_ref, totpol_ag; name1="ref", name2="ag")
eccs = J.plot_diff(totpol_diff)

subset!(totpol_diff, :Poll => x -> x .== "Methane")
subset!(totpol_diff, :Area => x -> x .== "California")
subset!(totpol_diff, :ECC => x -> x .== "Animal Production")
subset!(totpol_diff, :Year => x -> x .>= 2020)
subset!(totpol_diff, :Year => x -> x .<= 2040)

df = totpol_diff
using Makie, CairoMakie
using Plots, Colors
colors = distinguishable_colors(13)
# colors = palette(:viridis, 13)

fig = Figure()
ax = Axis(fig[1, 1], title="California Ag Emissions")

# Plot
ag = lines!(ax, df.Year, df.ag)
ref = lines!(ax, df.Year, df.ref)

# Legend
Legend(fig[1, 2],
  [ag, ref],
  ["Ag Update", "Reference"],
  "Scenarios")

fig

using Plots, Colors
colors = distinguishable_colors(13)
# colors = palette(:viridis, 13)

fig = Figure()
ax = Axis(fig[1, 1], title="test")

# Plot
barplot!(ax, df.Year, df.Diff,
  stack=levelcode.(df.Area),
  color=colors[levelcode.(df.Area)])

# Legend
labels = levels(df.Area)
elements = [PolyElement(polycolor=colors[i]) for i in 1:length(labels)]
title = "Groups"

Legend(fig[1, 2], elements, labels, title)

fig

import JuliaCompare as J
import JuliaCompare: db_files, Canada
import PromulaDBA as P
using DataFrames, DataFramesMeta
using Makie, CairoMakie, Plots, Colors

J.greet()

BASE_FOLDER = raw"\\Coral\c\2020CanadaAshPlus"
BASE_FOLDER = raw"\\Scarlet\c\2020CanadaBanyanPlus"
SCENARIO = "CER_UnitLimit"

CODE_FOLDER = joinpath(BASE_FOLDER, "Engine")
DATA_FOLDER = joinpath(BASE_FOLDER, "2020Model", SCENARIO)

vars = J.list_vars(CODE_FOLDER, DATA_FOLDER, db_files);
loc1 = J.Loc(vars, DATA_FOLDER);
# loc2 = J.Loc(vars, DATA_FOLDER2);

psoecc = J.var("PSoECC", loc1)
@subset!(psoecc, :Area .== "Alberta", occursin.("SAGD", :ECC))#, :Year .> 2019)

Plots.plot(psoecc.Year, psoecc.PSoECC)

cggen = J.var("CgGen", loc1)
@subset!(cggen, :Area .== "Alberta", occursin.("SAGD", :ECC))#, :Year .> 2019)
cggen = groupby(cggen, :Year)
cggen = @combine cggen :CgGen = sum(:CgGen)
Plots.plot(cggen.Year, cggen.CgGen)

vunpcfmax = J.var("vUnPCFMax", loc1)
@subset!(vunpcfmax, :Unit .!= "")

unname = J.var("UnName", loc1)

@rsubset(unname, (:Unit == "AB00001300201NG") | (:Unit == "AB06100000120NG"))

unpcfmax = J.var("UnPCFMax", loc1)
@subset!(unpcfmax, :Unit .!= "")

unpcf = J.var("UnPCF", loc1)
@subset!(unpcf, :Unit .!= "")

df = J.join_vars(unpcf, unpcfmax, vunpcfmax)

@subset(df, :UnPCF .>= :UnPCFMax)
@subset(df, :UnPCF .> :UnPCFMax)

unlimited = J.var("UnLimited", loc1)
@subset!(unlimited, :Unit .!= "")

df = J.join_vars(df, unlimited)

probs = @subset(df, :UnPCF .> :UnPCFMax, :UnLimited .== 1, :Year .> 2034 - 1984)
length(unique(probs.Unit))

unmustrun = J.var("UnMustRun", loc1)
@subset!(unmustrun, :Unit .!= "")

df = J.join_vars(df, unmustrun)
probs = @subset(df, :UnPCF .> :UnPCFMax * 1.001, :UnLimited .== 1, :Year .> 2034 - 1984, :UnMustRun .== 0)
length(unique(probs.Unit))

code = "AB00001300201NG"
code = "AB06100000120NG"
uneaf = J.var("UnEAF", loc1)
@subset(uneaf, :Unit .== code)

unegc = J.var("unegc", loc1)
egc = @subset(unegc, :Unit .== code)

ungc = J.var("ungc", loc1)
gc = @subset(ungc, :Unit .== code)

segc = J.join_vars(egc, gc)
@rtransform! segc :S_UnEGC = ifelse(:UnGC == 0, 0, :UnEGC / :UnGC)

using Statistics
gdf = groupby(segc, [:Unit, :Year])
s_unegc = @combine gdf :S_UnEGC = mean(:S_UnEGC)

# plotting

plot_unit = function (code::String; unit_index_year::Bool = false)
  name = unname[unname.Unit.==code, :UnName]
  name = only(name)

  df2 = outerjoin(
    @subset(unpcf, :Unit .== code),
    @subset(unpcfmax, :Unit .== code),
    on=[:Unit, :Year]
  )

  df2 = outerjoin(df2,
    @subset(uneaf, :Unit .== code, :Month .== "Winter"),
    on=[:Unit, :Year]
  )

  df2 = J.join_vars(df2, s_unegc)

  if unit_index_year == true
    df2.Year .+= 1984
  end
  @subset!(df2, :Year .< 2051, :Year .>= 2025)

  colors = distinguishable_colors(4)
  # colors = palette(:viridis, 13)

  fig = Figure()
  ax = Axis(fig[1, 1], title=name)

  # Plot
  ag = lines!(ax, df2.Year, df2.UnPCF, linewidth=3)
  ref1 = lines!(ax, df2.Year, df2.UnPCFMax, linewidth=3)
  ref2 = lines!(ax, df2.Year, df2.UnEAF, linestyle=:dot, linewidth=3)
  ref3 = lines!(ax, df2.Year, df2.S_UnEGC, linestyle=:dash, linewidth=3)

  # Legend
  Legend(fig[1, 2],
    [ag, ref1, ref2, ref3],
    ["UnPCF", "UnPCFMax", "UnEAF", "S_UnEGC"],
    "Variables")

  fig
end

plot_unit(code; unit_index_year=true)

totpol_ag = J.var("TotPol", loc2);

totpol_diff = J.diff(totpol_ref, totpol_ag; name1="ref", name2="ag")
eccs = J.plot_diff(totpol_diff)

subset!(totpol_diff, :Poll => x -> x .== "Methane")
subset!(totpol_diff, :Area => x -> x .== "California")
subset!(totpol_diff, :ECC => x -> x .== "Animal Production")
subset!(totpol_diff, :Year => x -> x .>= 2020)
subset!(totpol_diff, :Year => x -> x .<= 2040)

df = totpol_diff
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

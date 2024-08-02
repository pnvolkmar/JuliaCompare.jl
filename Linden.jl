import JuliaCompare as J
import PromulaDBA as P
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

BASE_FOLDER = raw"\\Pink\c\2020CanadaLinden"
BASE_FOLDER = raw"\\Silver\c\2020CanadaAshPlus"
# BASE_FOLDER2 = raw"C:\2020Git\2020BetaJulia\2020Model\Ref23\database.hdf5"
# BASE_FOLDER2 = raw"\\Pink\c\2020CanadaMyrtle"
# BASE_FOLDER2 = raw"\\Pink\c\24.02.17 2020CanadaAspen"
SCENARIO1 = ""
SCENARIO2 = "Ref23A_170_A"

CODE_FOLDER = joinpath(BASE_FOLDER, "Engine")
DATA_FOLDER1 = joinpath(BASE_FOLDER, "2020Model")
DATA_FOLDER2 = joinpath(BASE_FOLDER, "2020Model", SCENARIO2)
DATA_FOLDER_REF23 = joinpath(BASE_FOLDER, "2020Model", "Ref23")
E2020_Folder = joinpath(BASE_FOLDER, "2020Model")

vars = J.list_vars(CODE_FOLDER, E2020_Folder, db_files);
loc1 = J.Loc(vars, DATA_FOLDER1);
Ref23 = J.Loc(vars, DATA_FOLDER_REF23);

# fnames = CSV.read("mce.csv", DataFrame; types=[String, String])
# db files to be unzipped
# dbas = sort(unique(fnames.db))
# J.unzip_dbas(DATA_FOLDER1, dbas)
# J.unzip_dbas(DATA_FOLDER2, dbas)

# names = fnames.db .* "/" .* fnames.var
# P.data(joinpath(DATA_FOLDER1, fnames.db[1] * ".dba"), string(fnames.var[1]))
# P.data(joinpath(DATA_FOLDER1, "SpInput.dba"), "xRPPImportsROW")
# vars1 = [P.data(joinpath(DATA_FOLDER1, fnames.db[i] * ".dba"), fnames.var[i]) for i in 1:nrow(fnames)];
# vars2 = [P.data(joinpath(DATA_FOLDER2, fnames.db[i] * ".dba"), fnames.var[i]) for i in 1:nrow(fnames)];

xGoalPol = J.var("xGoalPol", loc1)
@rsubset! xGoalPol :xGoalPol != 0
@rsubset! xGoalPol :Market == "Market(161)"

plot(xGoalPol.Year, xGoalPol.xGoalPol)

xAdjust = J.var("xRPPAdjustArea", loc1)
vAdjust = J.var("vRPPAdjustmentsArea", loc1)
vStock = J.var("vRPPStockChangeArea", loc1)
vTransfer = J.var("vRPPTransfersArea", loc1)
df = J.var("xRPPProdArea", loc1)

rename!(vAdjust, Dict("vArea" => "Area"))
unique(vAdjust.Area)
unique(xAdjust.Area)
@rsubset! xAdjust :Area in canada
df = J.diff(xAdjust, vAdjust)
@rsubset df :Diff != 0
@rsubset df :old != 0

J.plot_diff(@rsubset df :Fuel == "Diesel"; dim="Area")


J.plot_sets(@rsubset xAdjust :Area == "AB"; dim="Fuel")
J.plot_sets(@rsubset vAdjust :vArea == "AB"; dim="Fuel")
J.plot_sets(vStock; dim="Fuel")
J.plot_sets(vTransfer; dim="Fuel")


totpol = J.var("TotPol", loc1)
@subset! totpol :Poll .== "CO2"
totpol0 = deepcopy(totpol)
totpol0.TotPol .= 0

df = J.diff(totpol, totpol0; name1="ref", name2="zero");
# diffs = J.comparedata(vars1, vars2, fnames)

anmap = J.var("anmap", loc1)
@subset! anmap :ANMap .== 1
select!(anmap, Not(:ANMap))

@subset! df :Year .>= 2020
leftjoin!(df, anmap, on=:Area)
@subset! df :Nation .== "CN"

classes = J.var("eccclmap", loc1)
@subset! classes :ECCCLMap .== 1
select!(classes, Not(:ECCCLMap))
df = leftjoin(df, classes, on=:ECC)
@subset! df :Class .== "Com"

J.plot_diff(@subset df :Area .== "AB"; dim="ECC")

EuDem = J.var("COutput/EuDem", loc1)
@rsubset! EuDem :Area in ["ON", "AB"]
@rsubset! EuDem :EC in ["Offices", "Health", "Education"]
@rsubset! EuDem :Year >= 2020
J.plot_sets(EuDem; dim="EC")
J.plot_sets(EuDem; dim="FuelEP")
J.plot_sets(EuDem; dim="Enduse")

EuDem = J.var("IOutput/EuDem", loc1)
@rsubset! EuDem :EC == "PulpPaperMills"
leftjoin!(EuDem, anmap, on=:Area)
@rsubset! EuDem :Nation == "CN"
select!(EuDem, Not(:Nation))

J.plot_sets(EuDem; dim="Area")
dropmissing!(EuDem)
dropmissing!(EuDem, :Year)

unega = J.var("ungc", loc1)
uncogen = J.var("uncogen", loc1)
unplant = J.var("unplant", loc1)

df = rightjoin(uncogen, unega, on=:Unit)
@rsubset! df :UnCogen == 1
df = rightjoin(unplant, df, on=:Unit)

# @rsubset! df :UnCogen == 1 :UnPlant == "BiomassCCS"

unsector = J.var("unsector", loc1)
df = rightjoin(unsector, df, on=:Unit)

@rsubset! df :Year > 2020
@rsubset! df :UnSector == "PulpPaperMills"
J.plot_sets(df; dim="UnPlant")

# 25.	Air Passenger - emissions higher in NZ than Ref23
df = J.diff("TotPol", Ref23, loc1; name1="Ref23", name2="NZ_CT")

# df = @by(df, [:ECC, :Area, :Year], :Ref23 = sum(:Ref23), :NZ_CT = sum(:NZ_CT))
@rsubset! df :Area in canada :Poll == "CO2"
df = @by(df, [:ECC, :Year], :Ref23 = sum(:Ref23), :NZ_CT = sum(:NZ_CT))
@rsubset! df :Year == 2050
@. df.Diff = df.Ref23 - df.NZ_CT
@rsubset df :Diff < 0
@rsubset df :ECC == "AirPassenger"

ref23 = J.var("TotPol", Ref23)
@rsubset! ref23 :Area in canada :Poll == "CO2" :Year == 2050
ref23 = @by(ref23, [:ECC], :TotPol = sum(:TotPol))
@rsubset ref23 :ECC == "AirPassenger"

ref23a = J.var("TotPol", loc2)
@rsubset! ref23a :Area in canada :Poll == "CO2" :Year == 2050
ref23a = @by(ref23a, [:ECC], :TotPol = sum(:TotPol))
@rsubset ref23a :ECC == "AirPassenger"

nz = J.var("totpol", loc1)
@rsubset! nz :Area in canada :Poll == "CO2" :Year == 2050
# nz = @by(nz, [:ECC], :TotPol = sum(:TotPol))
J.plot_sets(@rsubset nz :ECC == "AirPassenger"; dim="Area")
@rsubset ref23a :ECC == "AirPassenger"
@rsubset ref23 :ECC == "AirPassenger"

totdemand = J.var("totdemand", loc1)
df = @rsubset totdemand :Year >= 2020 :ECC == "AirPassenger" :Area in canada
J.plot_sets(df; dim="Fuel", num=10)

totdemand = J.var("totdemand", loc2)
df = @rsubset totdemand :Year >= 2020 :ECC == "AirPassenger" :Area in canada
J.plot_sets(df; dim="Fuel", num=10)

v = J.var("voaprod", loc1)
x = J.var("xoaprod", loc1)
o = J.var("oaprod", loc1)

dims = names(x)[1:end-1]
rename!(v, [dims; "OAProd"])
rename!(x, [dims; "OAProd"])
d = J.diff(v, x)
@rsubset d :Diff != 0 :Area in canada

unretire = J.var("UnRetire", Ref23)
uncogen = J.var("uncogen", Ref23)
unsector = J.var("unsector", Ref23)
unplant = J.var("unplant", Ref23)
unonline = J.var("unonline", Ref23)


@rsubset! unretire :Unit != ""
leftjoin!(unretire, uncogen, on=:Unit)
leftjoin!(unretire, unsector, on=:Unit)
leftjoin!(unretire, unplant, on=:Unit)
leftjoin!(unretire, unonline, on=:Unit)

@rsubset! unretire :UnCogen == 1 :UnSector == "PulpPaperMills"
@rsubset! unretire :UnPlant in ["BiomassCCS"; "Biomass"]
@rsubset! unretire :Year == 2030
@rsubset! unretire :UnRetire > 2023
@rsubset unretire :UnRetire != 2200

J.var("CD", Ref23)


using CSV
using XLSX
myfile = raw"C:\Users\Owner\SSI Dropbox\2020 Data\Linden\Petroleum products_input_data - Ref24.xlsx"

function get_data(sheetname)
  prod = XLSX.readdata(myfile, sheetname, "A:F")
  prod = DataFrame(prod[2:end, :], prod[1, :])
  dropmissing!(prod)
end

prod = get_data("vRPPAProd")
exports = get_data("vRPPExports")
imports = get_data("vRPPImports")
regional = get_data(4)
stock = get_data(6)
product = get_data(7)
other = get_data(8)
xRPPDemandArea = J.var("xRPPDemandArea", loc1)
xRPPProdArea = J.var("xRPPProdArea", loc1)
xRPPImportsArea = J.var("xRPPImportsArea", loc1)
xRPPExportsArea = J.var("xRPPExportsArea", loc1)
xRPPIntracountryArea = J.var("xRPPIntracountryArea", loc1)
xRPPAdjustArea = J.var("xRPPAdjustArea", loc1)
xTotDemand = J.var("xTotDemand", loc1)
xTotDemand_LPG = @rsubset xTotDemand :ECC == "Petroleum" :Fuel == "LPG"
select!(xTotDemand_LPG, Not(:ECC))
@rsubset! xTotDemand_LPG :Area ∈ Canada :Year >= 1990 :Year <= 2050
@rsubset! xTotDemand_LPG :Year == 2021
xTotDemand_LPG.xTotDemand .*= 1.054615

xRPP = J.join_vars(xRPPDemandArea, xRPPProdArea, xRPPImportsArea)
xRPP = J.join_vars(xRPP, xRPPExportsArea, xRPPIntracountryArea)
xRPP = J.join_vars(xRPP, xRPPAdjustArea)
@rsubset! xRPP :Fuel == "LPG" :Year >= 1990 :Year <= 2025
@rsubset! xRPP :Area ∈ Canada

@rsubset xRPP :Year == 2021

J.plot_sets(xRPPAdjustArea, dim="Area", title="xRPPAdjustArea", num=13)
J.plot_sets(xRPPAdjustArea, dim="Fuel", title="xRPPAdjustArea", num=13)
xRPP = J.join_vars(xRPPDemandArea, xRPPExportsArea)
xRPP = stack(xRPP, [:xRPPDemandArea, :xRPPExportsArea])
xRPP.value = Float64.(xRPP.value)
@rsubset! xRPP :value != 0 :Area ∈ Canada :Year <= 2021

J.plot_sets(xRPP, dim="variable", title="Canadian LHS")

xRPPR = J.join_vars(xRPPProdArea, xRPPImportsArea, xRPPIntracountryArea, xRPPAdjustArea)
xRPPR = stack(xRPPR, [:xRPPProdArea, :xRPPImportsArea, :xRPPIntracountryArea, :xRPPAdjustArea])
xRPPR.value = Float64.(xRPPR.value)
@rsubset! xRPPR :value != 0 :Area ∈ Canada :Year <= 2021

J.plot_sets(xRPPR, dim="variable", title="Canadian RHS")

vRPPAdjustmentsArea = J.var("vRPPAdjustmentsArea", loc1)
vRPPStockChangeArea = J.var("vRPPStockChangeArea", loc1)
vRPPTransfersArea = J.var("vRPPTransfersArea", loc1)

vRPPAdjustmentsArea.Variable .= "vRPPAdjustmentsArea"
vRPPStockChangeArea.Variable .= "vRPPStockChangeArea"
vRPPTransfersArea.Variable .= "vRPPTransfersArea"

select!(vRPPAdjustmentsArea, :Variable, Not(:Variable))
select!(vRPPStockChangeArea, :Variable, Not(:Variable))
select!(vRPPTransfersArea, :Variable, Not(:Variable))

rename!(vRPPAdjustmentsArea, :vRPPAdjustmentsArea => :Data)
rename!(vRPPStockChangeArea, :vRPPStockChangeArea => :Data)
rename!(vRPPTransfersArea, :vRPPTransfersArea => :Data)

Adjustments = vcat(vRPPAdjustmentsArea, vRPPStockChangeArea, vRPPTransfersArea)
@rsubset! Adjustments :Data != 0

J.plot_sets(Adjustments, dim="Variable", title="Adjustments")

J.plot_sets(xRPPAdjustArea, dim="Area", num=13, title="xRPPAdjustArea")


df = vcat(prod, exports, imports, regional, stock, product, other)
df.Data = Float64.(df.Data)
df.Variable = String.(df.Variable)
df.Year = Int.(df.Year)

J.plot_sets(df, dim="Variable", title="Canadian Data Breakdown")


import JuliaCompare as J
import PromulaDBA as P
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

BASE_FOLDER = raw"\\Blue\c\2020CanadaWalnut"
BASE_FOLDER2 = raw"\\Blue\c\23.11.15 2020CanadaAsh Blue"
BASE_FOLDER3 = raw"\\Coral\c\2020CanadaWalnut"
SCENARIO1 = "Ref24"
SCENARIO2 = "Ref23"
SCENARIO3 = "EfficiencyPrograms-DEESwTest"
SCENARIO4 = "OGRef"

CODE_FOLDER1 = joinpath(BASE_FOLDER, "Engine")
CODE_FOLDER2 = joinpath(BASE_FOLDER2, "Engine")
CODE_FOLDER3 = joinpath(BASE_FOLDER3, "Engine")
DATA_FOLDER1 = joinpath(BASE_FOLDER, "2020Model", SCENARIO1)
DATA_FOLDER2 = joinpath(BASE_FOLDER2, "2020Model", SCENARIO2)
DATA_FOLDER3 = joinpath(BASE_FOLDER3, "2020Model")
DATA_FOLDER4 = joinpath(BASE_FOLDER3, "2020Model", SCENARIO4)
E2020_Folder = joinpath(BASE_FOLDER, "2020Model")

vars = J.list_vars(CODE_FOLDER1, DATA_FOLDER1, db_files);
vars2 = J.list_vars(CODE_FOLDER2, DATA_FOLDER2, db_files);
vars3 = J.list_vars(CODE_FOLDER3, DATA_FOLDER3, db_files);
loc1 = J.Loc(vars, DATA_FOLDER1);
loc2 = J.Loc(vars2, DATA_FOLDER2);
loc3 = J.Loc(vars3, DATA_FOLDER3);
loc4 = J.Loc(vars3, DATA_FOLDER4);

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

CERSM = J.var("ICalDB/CERSM", loc2)
Enduse = J.var("IInput/EUKey", loc2)
@rsubset! CERSM :Enduse == "1"
select!(CERSM, Not(:Enduse))
@rsubset! CERSM :Area == "AB" :EC ∈ ["SAGDOilSands", "CSSOilSands"] :Year >= 2020

CERSM24 = J.var("ICalDB/CERSM", loc1)
@rsubset! CERSM24 :Enduse == "Heat"
select!(CERSM24, Not(:Enduse))
@rsubset! CERSM24 :Area == "AB" :EC ∈ ["SAGDOilSands", "CSSOilSands"] :Year >= 2020

df = J.diff(CERSM, CERSM24; name1="Ref23", name2="Ref24")
s = @rsubset df :EC == "SAGDOilSands"
c = @rsubset df :EC == "CSSOilSands"
fig = Figure()
ax = Axis(fig[1, 1]; xlabel="CERSM", ylabel="CERSM")
lines!(ax, s.Year, s.Ref23; color="Blue")
lines!(ax, s.Year, s.Ref24; color="Black")
lines!(ax, c.Year, c.Ref23; color="Blue")
lines!(ax, c.Year, c.Ref24; color="Black")
display(fig)

base = s[s.Year.==2022, :Ref23]
base = base[1]
@transform! s @byrow :Adjust = :Ref23 ./ base

base = c[c.Year.==2022, :Ref23]
base = base[1]
@transform! c @byrow :Adjust = :Ref23 ./ base


EuDemand_24 = J.var("EuDemand", loc3)
EuDemand_23 = J.var("EuDemand", loc2)
EuDemand_d = J.diff(EuDemand_24, EuDemand_23; name1="Ref24", name2="Ref23")
@rsubset! EuDemand_d :ECC ∈ ["SAGDOilSands", "CSSOilSands"]
@rsubset! EuDemand_d :Area == "AB"
@rsubset! EuDemand_d :Fuel == "NaturalGas"
@rsubset! EuDemand_d :Year > 2020 :Year <= 2050
J.plot_sets(EuDemand_d, dim="ECC", num=2, title="Ref24-Ref23 EuDemand for OilSands")
df = @by(EuDemand_d, by = ["ECC", "Year"], :Ref24 = sum(:Ref24), :Ref23 = sum(:Ref23))
df.Diff = df.Ref24 - df.Ref23
@rsubset df :ECC == "SAGDOilSands"
@rsubset df :ECC == "CSSOilSands"



f! = function (df)
  n = names(df)
  if in("ECC", n)
    @rsubset! df :ECC ∈ ["SAGDOilSands", "CSSOilSands"]
  end
  if in("Area", n)
    @rsubset! df :Area == "AB"
  end
  if in("Fuel", n)
    @rsubset! df :Fuel == "NaturalGas"
  end
  if in("Year", n)
    @rsubset! df :Year > 2020 :Year <= 2050
  end
end

FPF_24 = J.var("FPF", loc3)
FPF_23 = J.var("FPF", loc2)
FPF_d = J.diff(FPF_24, FPF_23)
@rsubset! FPF_d :ES == "Industrial" :Fuel == "NaturalGas"
@rsubset! FPF_d :Area == "AB"
@rsubset! FPF_d :Year > 2020 :Year <= 2050
J.plot_sets(FPF_d, dim="Area", num=1, title="Ref24-Ref23 FPF for OilSands")
df = @by(FPF_d, by = ["ECC", "Year"], :new = sum(:new), :old = sum(:old))
df.Diff = df.new - df.old

DEE_24 = J.var("IOutput/DEE", loc1)
DEE_23 = J.var("IOutput/DEE", loc2)
DEE_New = J.var("IOutput/DEE", loc3)
DEMM = J.var("ICalDB/DEMM", loc3)
DEE_OG = J.var("IOutput/DEE", loc4)
DEE_x = J.var("IInput/xDEE", loc3)
Enduse = J.var("IInput/EUKey", loc2)
DEE_23 = J.join_vars(DEE_23, Enduse)
select!(DEE_23, Not(:Enduse))
DEE_23.Enduse = DEE_23.EUKey
select!(DEE_23, names(DEE_24))
DEE_dNew = J.diff(DEE_New, DEE_23; name1="Ref24_New", name2="Ref23")
DEE_d = J.diff(DEE_24, DEE_23; name1="Ref24", name2="Ref23")
@rsubset! DEE_d :EC ∈ ["SAGDOilSands", "CSSOilSands"]
@rsubset! DEE_d :Area == "AB" :Tech == "Gas" :Enduse == "Heat"
@rsubset! DEE_d :Year > 2020 :Year <= 2050
@rsubset! DEE_dNew :EC ∈ ["SAGDOilSands", "CSSOilSands"]
@rsubset! DEE_dNew :Area == "AB" :Tech == "Gas" :Enduse == "Heat"
@rsubset! DEE_dNew :Year > 2020 :Year <= 2050
J.plot_sets(DEE_d, dim="EC", title="Ref24-Ref23 DEE for Gas Tech for OilSands")
J.plot_sets(DEE_dNew, dim="EC", title="Ref24New-Ref23 DEE for Gas Tech for OilSands")
df = @by(DEE_d, by = ["ECC", "Year"], :new = sum(:new), :old = sum(:old))
df.Diff = df.new - df.old

df = J.join_vars(DEE_d, DEE_New)
df = J.join_vars(df, DEE_x)
df = J.join_vars(df, DEMM)
df = J.join_vars(df, DEE_OG)

PEE_24 = J.var("IOutput/PEE", loc1)
PEE_23 = J.var("IOutput/PEE", loc2)
PEE_23 = J.join_vars(PEE_23, Enduse)
select!(PEE_23, Not(:Enduse))
PEE_23.Enduse = PEE_23.EUKey
select!(PEE_23, names(PEE_24))
PEE_d = J.diff(PEE_24, PEE_23; name1="Ref24", name2="Ref23")
@rsubset! PEE_d :EC ∈ ["SAGDOilSands", "CSSOilSands"]
@rsubset! PEE_d :Area == "AB" :Tech == "Gas" :Enduse == "Heat"
@rsubset! PEE_d :Year > 2020 :Year <= 2050

@transform! PEE_d @byrow :PDiff = :Ref23 / :Ref24
@rsubset PEE_d :EC == "SAGDOilSands"

DEESw_24 = J.var("IInput/DEESw", loc1)
DEESw_23 = J.var("IInput/DEESw", loc3)
DEESw_d = J.diff(DEESw_24, DEESw_23; name1="Ref24", name2="Sw")
@rsubset! DEESw_d :EC ∈ ["SAGDOilSands", "CSSOilSands"]
@rsubset! DEESw_d :Area == "AB" :Tech == "Gas" :Enduse == "Heat"
@rsubset! DEESw_d :Year > 2020 :Year <= 2050

@transform! DEESw_d @byrow :PDiff = :Ref23 / :Ref24
@rsubset DEESw_d :EC == "SAGDOilSands"


ECFP_24 = J.var("IOutput/ECFP", loc1)
ECFP_23 = J.var("IOutput/ECFP", loc2)
ECFP_23 = J.join_vars(ECFP_23, Enduse)
select!(ECFP_23, Not(:Enduse))
ECFP_23.Enduse = ECFP_23.EUKey
select!(ECFP_23, names(ECFP_24))
ECFP_d = J.diff(ECFP_24, ECFP_23; name1="Ref24", name2="Ref23")
@rsubset! ECFP_d :EC ∈ ["SAGDOilSands", "CSSOilSands"]
@rsubset! ECFP_d :Area == "AB" :Tech == "Gas" :Enduse == "Heat"
@rsubset! ECFP_d :Year > 2020 :Year <= 2030
J.plot_sets(ECFP_d, dim="EC", title="Ref24-Ref23 ECFP for Gas Tech for OilSands")
df = @by(ECFP_d, by = ["ECC", "Year"], :new = sum(:new), :old = sum(:old))
df.Diff = df.new - df.old

PCostTech_24 = J.var("IOutput/PCostTech", loc1)
PCostTech_23 = J.var("IOutput/PCostTech", loc2)
PCostTech_d = J.diff(PCostTech_24, PCostTech_23; name1="Ref24", name2="Ref23")
@rsubset! PCostTech_d :EC ∈ ["SAGDOilSands", "CSSOilSands"]
@rsubset! PCostTech_d :Area == "AB" :Tech == "Gas"
@rsubset! PCostTech_d :Year > 2020 :Year <= 2030
J.plot_sets(PCostTech_d, dim="EC", title="Ref24-Ref23 PCostTech for Gas Tech for OilSands")
df = @by(PCostTech_d, by = ["ECC", "Year"], :new = sum(:new), :old = sum(:old))
df.Diff = df.new - df.old

DEMMM_24 = J.var("IOutput/DEMMM", loc1)
DEMMM_23 = J.var("IOutput/DEMMM", loc2)
DEMMM_23 = J.join_vars(DEMMM_23, Enduse)
select!(DEMMM_23, Not(:Enduse))
DEMMM_23.Enduse = DEMMM_23.EUKey
select!(DEMMM_23, names(DEMMM_24))
DEMMM_d = J.diff(DEMMM_24, DEMMM_23; name1="Ref24", name2="Ref23")
@rsubset! DEMMM_d :EC ∈ ["SAGDOilSands", "CSSOilSands"]
@rsubset! DEMMM_d :Area == "AB" :Tech == "Gas" :Enduse == "Heat"
@rsubset! DEMMM_d :Year > 2020 :Year <= 2030
J.plot_sets(DEMMM_d, dim="EC", title="Ref24-Ref23 DEMMM for Gas Tech for OilSands")
df = @by(DEMMM_d, by = ["ECC", "Year"], :new = sum(:new), :old = sum(:old))
df.Diff = df.new - df.old

FPCFSTech_24 = J.var("IOutput/FPCFSTech", loc1)
FPCFSTech_23 = J.var("IOutput/FPCFSTech", loc2)
FPCFSTech_d = J.diff(FPCFSTech_24, FPCFSTech_23; name1="Ref24", name2="Ref23")
@rsubset! FPCFSTech_d :EC ∈ ["SAGDOilSands", "CSSOilSands"]
@rsubset! FPCFSTech_d :Area == "AB" :Tech == "Gas"
@rsubset! FPCFSTech_d :Year > 2020 :Year <= 2030
J.plot_sets(FPCFSTech_d, dim="EC", title="Ref24-Ref23 FPCFSTech for Gas Tech for OilSands")
df = @by(FPCFSTech_d, by = ["ECC", "Year"], :new = sum(:new), :old = sum(:old))
df.Diff = df.new - df.old

NewPrice_23 = J.join_vars(ECFP_23, PCostTech_23, FPCFSTech_23)
NewPrice_24 = J.join_vars(ECFP_24, PCostTech_24, FPCFSTech_24)
@rsubset! NewPrice_23 :Area == "AB" :Year >= 2020 :Year <= 2030 :Enduse == "Heat" :Tech == "Gas"
@rsubset! NewPrice_24 :Area == "AB" :Year >= 2020 :Year <= 2030 :Enduse == "Heat" :Tech == "Gas"
@rsubset! NewPrice_23 :EC ∈ ["SAGDOilSands", "CSSOilSands"]
@rsubset! NewPrice_24 :EC ∈ ["SAGDOilSands", "CSSOilSands"]
@transform! NewPrice_23 @byrow :NewPrice = :ECFP - :PCostTech - :FPCFSTech
@transform! NewPrice_24 @byrow :NewPrice = :ECFP - :PCostTech - :FPCFSTech
@select! NewPrice_23 :EC :Year :NewPrice
@select! NewPrice_24 :EC :Year :NewPrice
NewPrice_d = J.diff(NewPrice_24, NewPrice_23; name1="Ref24", name2="Ref23")

Dev = J.var("Dev", loc2)
OGName = J.var("OGName", loc2)
OGArea = J.var("OGArea", loc2)
OGProcess = J.var("OGProcess", loc2)
Dev = J.join_vars(Dev, OGName, OGArea, OGProcess)
Dev.Year .+= 1984
@rsubset! Dev :OGProcess ∈ ["SAGDOilSands", "CSSOilSands"]
@rsubset! Dev :OGArea == "AB"
@rsubset! Dev :Year > 2019 :Year <= 2050
J.plot_sets(Dev; dim="OGProcess", title="Development of OG Resources")

# Row │ ECC           Year   new      old      Diff      
# │ String        Int64  Float64  Float64  Float64   
# ─────┼────────────────────────────────────────────────── 
# 11 │ SAGDOilSands   2050  437.513  582.633  -145.12   
# 12 │ CSSOilSands    2050  147.894  202.205   -54.3105 

TotDemand_24 = J.var("TotDemand", loc1)
TotDemand_23 = J.var("TotDemand", loc2)
TotDemand_23.Year .+= 1984
TotDemand_d = J.diff(TotDemand_23, TotDemand_24)
@rsubset! TotDemand_d :ECC ∈ ["SAGDOilSands", "CSSOilSands"] :Fuel == "NaturalGas"
@rsubset! TotDemand_d :Area ∈ Canada
@rsubset! TotDemand_d :Year > 2020 :Year <= 2050
J.plot_sets(TotDemand_d, dim="Area", num=1)
df = @by(TotDemand_d, by = ["ECC", "Year"], :new = sum(:new), :old = sum(:old))
df.Diff = df.new - df.old

# Row │ ECC           Year   new      old      Diff      
# │ String        Int64  Float64  Float64  Float64   
# ─────┼────────────────────────────────────────────────── 
# 11 │ SAGDOilSands   2050  437.513  582.633  -145.12   
# 12 │ CSSOilSands    2050  147.894  202.205   -54.3105 

df2 = J.diff("IOutput/SqDmd", loc1, loc2)
ref24 = J.var("IOutput/SqDmd", loc1)
ref23 = J.var("IOutput/SqDmd", loc2)
ref23.Year .+= 1984
df = J.diff(ref23, ref24)
@rsubset! df :Area ∈ Canada :EC ∈ ["SAGDOilSands", "CSSOilSands"]
@rsubset df :Year == 2050 :Diff != 0





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



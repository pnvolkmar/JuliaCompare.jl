import JuliaCompare as J
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

BASE_FOLDER = raw"\\Mint\c\2020CanadaAshPlus"
# BASE_FOLDER2 = raw"C:\2020Git\2020BetaJulia\2020Model\Ref23\database.hdf5"
# BASE_FOLDER2 = raw"\\Pink\c\2020CanadaMyrtle"
# BASE_FOLDER2 = raw"\\Pink\c\24.02.17 2020CanadaAspen"
SCENARIO1 = "NZ_CT_500"
SCENARIO2 = "Ref23A_170_A"

CODE_FOLDER = joinpath(BASE_FOLDER, "Engine")
DATA_FOLDER1 = joinpath(BASE_FOLDER, "2020Model", SCENARIO1)
DATA_FOLDER2 = joinpath(BASE_FOLDER, "2020Model", SCENARIO2)
DATA_FOLDER_REF23 = joinpath(BASE_FOLDER, "2020Model", "Ref23")
E2020_Folder = joinpath(BASE_FOLDER, "2020Model")

vars = J.list_vars(CODE_FOLDER, DATA_FOLDER1, db_files);
loc1 = J.Loc(vars, DATA_FOLDER1);
loc2 = J.Loc(vars, DATA_FOLDER2);
Ref23 = J.Loc(vars, DATA_FOLDER_REF23);

# fnames = CSV.read("mce.csv", DataFrame; types=[String, String])
# db files to be unzipped
# dbas = sort(unique(fnames.db))
# J.unzip_dbas(DATA_FOLDER1, dbas)
# J.unzip_dbas(DATA_FOLDER2, dbas)

# names = fnames.db .* "/" .* fnames.var
import PromulaDBA as P
# P.data(joinpath(DATA_FOLDER1, fnames.db[1] * ".dba"), string(fnames.var[1]))
# P.data(joinpath(DATA_FOLDER1, "SpInput.dba"), "xRPPImportsROW")
# vars1 = [P.data(joinpath(DATA_FOLDER1, fnames.db[i] * ".dba"), fnames.var[i]) for i in 1:nrow(fnames)];
# vars2 = [P.data(joinpath(DATA_FOLDER2, fnames.db[i] * ".dba"), fnames.var[i]) for i in 1:nrow(fnames)];

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

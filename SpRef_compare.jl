import JuliaCompare as J
import JuliaCompare: db_files, Canada
using DataFrames

BASE_FOLDER = raw"\\Pink\c\2020CanadaLinden"
SCENARIO1 = "Calib4_10p"
SCENARIO2 = "Calib4_AB"

CODE_FOLDER = joinpath(BASE_FOLDER, "Engine")
DATA_FOLDER1 = joinpath(BASE_FOLDER, "2020Model", SCENARIO1)
DATA_FOLDER2 = joinpath(BASE_FOLDER, "2020Model", SCENARIO2)
E2020_Folder = joinpath(BASE_FOLDER, "2020Model")

vars = J.list_vars(CODE_FOLDER, DATA_FOLDER1, db_files);
loc1 = J.Loc(vars, DATA_FOLDER1);
loc2 = J.Loc(vars, DATA_FOLDER2);

using CSV, DataFrames, DataFramesMeta
fnames = CSV.read("SpRef_filenames.csv", DataFrame; types=[String, String])

# db files to be unzipped
dbas = sort(unique(fnames.db))
J.unzip_dbas(DATA_FOLDER1, dbas)
J.unzip_dbas(DATA_FOLDER2, dbas)

# names = fnames.db .* "/" .* fnames.var
import PromulaDBA as P
# P.data(joinpath(DATA_FOLDER1, fnames.db[1] * ".dba"), string(fnames.var[1]))
# P.data(joinpath(DATA_FOLDER1, "SpInput.dba"), "xRPPImportsROW")
vars1 = [P.data(joinpath(DATA_FOLDER1, fnames.db[i] * ".dba"), fnames.var[i]) for i in 1:nrow(fnames)];
vars2 = [P.data(joinpath(DATA_FOLDER2, fnames.db[i] * ".dba"), fnames.var[i]) for i in 1:nrow(fnames)];

J.var("RfArea", loc1)

diffs = J.comparedata(vars1, vars2, fnames)
@subset! diffs :Diff .!= 0
@orderby diffs -:Diff

df = J.diff("RfCapEffective", loc1, loc2; name1=SCENARIO1, name2=SCENARIO2)
@subset! df :Diff .!= 0
@subset df :Calib4_AB .== 0
@subset df :RfUnit .== "RfUnit(4)"
xrfcap = J.diff("xRfCap", loc1, loc2; name1=SCENARIO1, name2=SCENARIO2)
@subset! xrfcap :Diff .!= 0

rfoor = J.diff("RfOOR", loc1, loc2; name1=SCENARIO1, name2=SCENARIO2)
@subset! rfoor :Diff .!= 0
@orderby rfoor :Year :RfUnit
@subset rfoor :RfUnit .== "RfUnit(1)" :Fuel .== "Asphalt"

xRPPProdArea = J.diff("xRPPProdArea", loc1, loc2; name1=SCENARIO1, name2=SCENARIO2)
@subset! xRPPProdArea :Diff .!= 0
RPPProdArea = J.diff("RPPProdArea", loc1, loc2; name1=SCENARIO1, name2=SCENARIO2)
@subset! RPPProdArea :Diff .> 0.1 .|| :Diff .< -0.1
df = @by(RPPProdArea, [:Area, :Year], :Calib4_10p = sum(:Calib4_10p), :Calib4_AB = sum(:Calib4_AB))
df.Diff = df.Calib4_10p - df.Calib4_AB
@subset df :Area .== "ON"

xRfProd = J.diff("xRfProd", loc1, loc2; name1=SCENARIO1, name2=SCENARIO2)
@subset! xRfProd :Diff .!= 0
RfProd = J.diff("RfProd", loc1, loc2; name1=SCENARIO1, name2=SCENARIO2)
@subset! RfProd :Diff .> 0.1 .|| :Diff .< -0.1
df = @by(xRfProd, [:RfUnit, :Year], :Calib4_10p = sum(:Calib4_10p), :Calib4_AB = sum(:Calib4_AB))
df.Diff = df.Calib4_10p - df.Calib4_AB
@subset df :RfUnit .== "RfUnit(1)"

import JuliaCompare as J
import JuliaCompare: db_files, Canada
using DataFrames

BASE_FOLDER = raw"\\Pink\c\2020CanadaLinden"
SCENARIO1 = "Calib4"
SCENARIO2 = "Calib4_rppprd"

CODE_FOLDER = joinpath(BASE_FOLDER, "Engine")
DATA_FOLDER1 = joinpath(BASE_FOLDER, "2020Model", SCENARIO1)
DATA_FOLDER2 = joinpath(BASE_FOLDER, "2020Model", SCENARIO2)

vars = J.list_vars(CODE_FOLDER, DATA_FOLDER1, db_files);
loc1 = J.Loc(vars, DATA_FOLDER1);
loc2 = J.Loc(vars, DATA_FOLDER2);

using CSV, DataFrames, DataFramesMeta
fnames = CSV.read("SpRef_filenames.csv", DataFrame; types=[String, String])

# db files to be unzipped
sort(unique(fnames.db))
names = fnames.db .* "/" .* fnames.var
import PromulaDBA as P
P.data(joinpath(DATA_FOLDER1, fnames.db[1] * ".dba"), string(fnames.var[1]))
P.data(joinpath(DATA_FOLDER1, "SpInput.dba"), "xRPPImportsROW")
vars1 = [P.data(joinpath(DATA_FOLDER1, fnames.db[i] * ".dba"), fnames.var[i]) for i in 1:nrow(fnames)];
vars2 = [P.data(joinpath(DATA_FOLDER2, fnames.db[i] * ".dba"), fnames.var[i]) for i in 1:nrow(fnames)];

diffs = J.comparedata(vars1, vars2, fnames)
@subset diffs :Diff .!= 0

using Revise
import JuliaCompare as J
import PromulaDBA as P
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

BASE_FOLDER = raw"\\Pink\c\2020CanadaSpruce"
BASE_FOLDERT = raw"\\Pink\c\2020CanadaTanoak"
SCENARIO1 = "Calib"
SCENARIO2 = "Calib2"

CODE_FOLDER = joinpath(BASE_FOLDER, "Engine")
DATA_FOLDERS1 = joinpath(BASE_FOLDER, "2020Model", SCENARIO1)
DATA_FOLDERS2 = joinpath(BASE_FOLDER, "2020Model", SCENARIO2)
DATA_FOLDERT1 = joinpath(BASE_FOLDERT, "2020Model", SCENARIO1)
DATA_FOLDERT2 = joinpath(BASE_FOLDERT, "2020Model", SCENARIO2)
E2020_Folder = joinpath(BASE_FOLDER, "2020Model")
HDF5_pathT1 = joinpath(DATA_FOLDERT1, "database.hdf5")
HDF5_pathT2 = joinpath(DATA_FOLDERT2, "database.hdf5")

vars = J.list_vars(CODE_FOLDER, DATA_FOLDER1, db_files);
vars_j = J.list_vars(HDF5_path)
SpruceC1 = J.Loc_p(vars, DATA_FOLDERS1, "SpruceCalib");
SpruceC2 = J.Loc_p(vars, DATA_FOLDERS2, "SpruceCalib2");
TanoakC1 = J.Loc_j(vars_j, HDF5_pathT1, "TanoakCalib");
TanoakC2 = J.Loc_j(vars_j, HDF5_pathT2, "TanoakCalib2");

# SaEC=ElecDmd-CgEC+PSoECC 
SaEC = J.diff("SOutput/SaEC", SpruceC2, TanoakC2); J.f_al(SaEC) # SpruceC2 - TanoakC2
ElecDmd = J.diff("SOutput/ElecDmd", SpruceC2, TanoakC2); J.f_al(ElecDmd) # SpruceC2 - TanoakC2
CgEC = J.diff("SOutput/CgEC", SpruceC2, TanoakC2); J.f_al(CgEC) # SpruceC2 - TanoakC2
PSoECC = J.diff("SOutput/PSoECC", SpruceC2, TanoakC2); J.f_al(PSoECC) # SpruceC2 - TanoakC2
# CgGen[fuel,ecc,area] = CgGen[fuel,ecc,area] + UnEGA[unit]
CgGen = J.diff("SOutput/CgGen", SpruceC2, TanoakC2); J.f_al(CgGen) # SpruceC2 - TanoakC2
UnEGA = J.diff("UnEGA", SpruceC2, TanoakC2) # SpruceC2 - TanoakC2
FlPlnMap = J.diff("EGInput/FlPlnMap", SpruceC2, TanoakC2) # Matches
@rsubset FlPlnMap :Diff != 0

using Revise
import JuliaCompare as J
import PromulaDBA as P
import SmallModel as M 
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

BASE_FOLDER = raw"\\Silver\c\2020CanadaSpruce"
BASE_FOLDER2 = raw"\\Silver\c\2020CanadaTanoak"
SCENARIO1 = ""
SCENARIO2 = ""

CODE_FOLDER = joinpath(BASE_FOLDER, "Engine")
DATA_FOLDER1 = joinpath(BASE_FOLDER, "2020Model")
DATA_FOLDER2 = joinpath(BASE_FOLDER2, "2020Model")

DATA_FOLDER1 = SCENARIO1 == "" ? DATA_FOLDER1 : joinpath(DATA_FOLDER1,SCENARIO1)
DATA_FOLDER2 = SCENARIO2 == "" ? DATA_FOLDER2 : joinpath(DATA_FOLDER2,SCENARIO1)
E2020_Folder = joinpath(BASE_FOLDER, "2020Model")
HDF5_path = joinpath(DATA_FOLDER2, "database.hdf5")

vars = J.list_vars(CODE_FOLDER, DATA_FOLDER1, db_files);
vars_j = J.list_vars(HDF5_path)
loc1 = J.Loc_p(vars, DATA_FOLDER1, "Spruce");
loc2 = J.Loc_j(vars_j, HDF5_path, "Tanoak");

EuFPol = J.diff("EuFPol", loc1, loc2)
@rsubset EuFPol :FuelEP == "Biomass"
@rsubset! EuFPol :FuelEP == "Biomass"
select!(EuFPol, Not([:FuelEP]))
@rsubset! EuFPol :Diff != 0

J.plot_diff(EuFPol; dim="ECC", num=10, title="Biomass diffs by ECC")

# UtilityGen appears to have a lot of historic issues. Let's tackle that first.

@rsubset! EuFPol :ECC == "UtilityGen"

J.plot_diff(@rsubset EuFPol :Year > 2025; dim="Poll", num=10, title="Biomass diffs by Poll")

# issues across Polls with steady ratios historically concentrated in NOX and COX
@rsubset! EuFPol :Poll == "NOX"
J.plot_diff(EuFPol; dim="Area", num=10, title="Biomass diffs by Area")

@rsubset EuFPol :Year == 1996

# issues not solely in ON, but ON is by far the worst offender

gd = groupby(EuFPol, :Year) 
gd = combine(gd, :Diff .=> sum)
gd = sort!(gd, :Diff_sum, rev=true)

# 1996 is the worst year

UnPol = J.diff("EGOutput/UnPol", loc1, loc2)

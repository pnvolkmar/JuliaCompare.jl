using Revise
import JuliaCompare as J
import PromulaDBA as P
import SmallModel as M 
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

BASE_FOLDER = raw"\\Pink\c\2020CanadaSpruce"
BASE_FOLDER2 = raw"\\Pink\c\2020CanadaTanoak"
SCENARIO1 = "Ref24"
SCENARIO2 = "Ref24"

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

@time EuFPol = J.diff("EuFPol", loc1, loc2)
@time J.var("EuFPol", loc2)
using Profile, ProfileView
ProfileView.@profview EuFPol = J.diff("EuFPol", loc1, loc2)

# Clear any previous profiling data
Profile.clear()

# Profile with default settings
@profile EuFPol = J.diff("EuFPol", loc1, loc2)

# View results in the terminal
Profile.print()


@rsubset! EuFPol :Area ∈ Canada
J.plot_diff(EuFPol; dim="ECC", num=10, title="EuFPol diffs by ECC")
# J.plot_diff(@rsubset EuFPol :FuelEP == "Biomass"; dim="ECC", num=10, title="Biomass diffs by ECC")


@rsubset! EuFPol :ECC == "UtilityGen"
# Almost all Natural Gas in the future and Coal in the past
J.plot_diff(EuFPol; dim="FuelEP", num=10, title="EuFPol diffs by FuelEP") 
# All CO2 in the future and SOX/NOX in the past
J.plot_diff(EuFPol; dim="Poll", num=10, title="EuFPol diffs by Poll") 
# SK, ON, AB are the biggest contributors, though many are present
J.plot_diff(EuFPol; dim="Area", num=10, title="EuFPol diffs by Area") 

TotPol_j = J.var("TotPol", loc2)

TotPol = J.diff("TotPol", loc1, loc2)
@rsubset! TotPol :Area ∈ Canada
J.plot_diff(TotPol; dim="ECC", num=10, title="TotPol diffs by ECC")

@rsubset! TotPol :ECC == "UtilityGen"
J.plot_diff(TotPol; dim="FuelEP", num=10, title="Biomass diffs by FuelEP") # Almost all Natural Gas
J.plot_diff(TotPol; dim="Poll", num=10, title="TotPol diffs by Poll") # All CO2
J.plot_diff(TotPol; dim="Area", num=10, title="TotPol diffs by Area") # Mostly CA a little in Mtn

# UtilityGen appears to have a lot of historic issues. Let's tackle that first.

UnPolGross_p = P.data(joinpath(DATA_FOLDER1,"EGOutput3.dba"),"UnPolGross")
db = loc2.HDF5_path
UnPolGross = M.ReadDisk(db, "EGOutput/UnPolGross")
UnArea = M.ReadDisk(db,"EGInput/UnArea")
UnCode = M.ReadDisk(db,"EGInput/UnCode")
UnCode_p = P.data(joinpath(DATA_FOLDER1,"EGInput.dba"), "UnCode")
ECC = M.ReadDisk(db,"SInput/ECC")
# We have an issue with UnCode
Codes = DataFrame(Spruce = UnCode_p, Tanoak = UnCode, Agree = UnCode .== UnCode_p)
@rsubset! Codes :Spruce != "" && :Tanoak != "Null"
@rsubset Codes :Agree == false
Codes = J.reconcile_codes(Codes)
@rsubset Codes :MatchAfterTransform == false

@rsubset Codes :Spruce == "NL_Cg_ECC34_OGSteam"

UnPlant_j = M.ReadDisk(db,"EGInput/UnPlant")
UnPlant_p = P.data(joinpath(DATA_FOLDER1,"EGInput.dba"), "UnPlant")
Plants = DataFrame(Spruce = UnPlant_p, Tanoak = UnPlant_j, Agree = UnPlant_j .== UnPlant_p)
@rsubset! Plants :Spruce != "" && :Tanoak != "Null"
@rsubset Plants :Agree == false
# All pants Agree
Codes.Plant = Plants.Spruce

@rsubset Codes :MatchAfterTransform == false
@rsubset Codes :MatchAfterTransform == false :Plant == "OGCC"

FuelEP = M.ReadDisk(db,"E2020DB/FuelEP")
Poll = M.ReadDisk(db, "E2020DB/Poll")
import SmallModel: Select, Yr
fuelep = Select(FuelEP, "NaturalGas")
poll = Select(Poll, "CO2")
year = Yr(2050)
df = DataFrame(UnCode = UnCode, UnArea = UnArea, UnPlant = UnPlant_j,
  Spruce = UnPolGross_p[:,fuelep,poll,year], Tanoak = UnPolGross[:,fuelep,poll,year])

df.Diff = df.Spruce - df.Tanoak
@rsubset! df :UnArea == "SK" abs(:Diff) > 1e-7
sort(df, :UnPlant)

UnEGA = J.diff("UnEGA", loc1, loc2)
UnArea
UnCode_p
UnArea = DataFrame(Unit = UnCode_p, UnArea = UnArea)
AddPlant = DataFrame(Unit = UnCode_p, UnPlant = UnPlant_p)
UnEGA = DataFrames.leftjoin(UnEGA, UnArea, on = :Unit)
UnEGA = DataFrames.leftjoin(UnEGA, AddPlant, on = :Unit)

UnRetire = M.ReadDisk(db,"EGInput/UnRetire", year)
AddRetire = DataFrame(Unit = UnCode_p, UnRetire = UnRetire)
UnEGA = DataFrames.leftjoin(UnEGA, AddRetire, on = :Unit)

UnOnline = M.ReadDisk(db,"EGInput/UnOnLine")
AddOnline = DataFrame(Unit = UnCode_p, UnOnline = UnOnline)
UnEGA = DataFrames.leftjoin(UnEGA, AddOnline, on = :Unit)

issues = @rsubset UnEGA :UnArea == "SK" :Year == 2039 abs(:Diff) > 1e-5
issue_codes = unique(issues.Unit)

UnGC = J.diff("UnGC", loc1, loc2)
df = @rsubset UnGC :Unit ∈ issue_codes abs(:Diff) > 0.01;
sort(df, [:Year,:Diff])

UnGCCR = J.diff("UnGCCR", loc1, loc2)
df = @rsubset UnGCCR :Unit ∈ issue_codes abs(:Diff) > 0.01;
sort(df, [:Year,:Diff])

UnGCCE = J.diff("UnGCCE", loc1, loc2)
df = @rsubset UnGCCE :Unit ∈ issue_codes abs(:Diff) > 0.01;
sort(df, [:Year,:Diff])

xUnGCCR = J.diff("xUnGCCR", loc1, loc2)
df = @rsubset xUnGCCR :Unit ∈ issue_codes abs(:Diff) > 0.01;
sort(df, [:Year,:Diff])

UnGCCI = J.diff("UnGCCI", loc1, loc2)
df = @rsubset UnGCCI :Unit ∈ issue_codes abs(:Diff) > 0.01;
sort(df, [:Year,:Diff])

xUnGCCI = J.diff("xUnGCCI", loc1, loc2)
df = @rsubset xUnGCCI :Unit ∈ issue_codes abs(:Diff) > 0.01;
sort(df, [:Year,:Diff])

CgGCCI = J.diff("CgGCCI", loc1, loc2)
df = @rsubset CgGCCI :Unit ∈ issue_codes abs(:Diff) > 0.01;
sort(df, [:Year,:Diff])

max(HDGCCI[plant,node,genco,area],IPGCCI[plant,node,genco,area],RnGCCI[plant,node,genco,area])

HDGCCI = J.diff("HDGCCI", loc1, loc2)
df = @rsubset HDGCCI :Plant == "SK" :Node == "SK" :Area == "SK" :Year == 2023;
sort(df, [:Year,:Diff])


@rsubset xUnGCCR :Unit == "SK_New_SolarPV"  abs(:Diff) > 0.01

xUnGCCI = J.diff("xUnGCCI", loc1, loc2)
@rsubset xUnGCCR :Diff == :Spruce :Diff != 0
@rsubset xUnGCCI :Unit == "SK_New_SolarPV" :Year > 2020 :Year <2040

GAProd = J.var("GAProd", loc2)
@rsubset GAProd isnan(:Value)

import SmallModel: ReadDisk
CurTime::Float64 = ReadDisk(db,"CInput/CurTime")[1] # Year for capital costs [tv]
YrDCC::Int = Int(CurTime)

xungccr = ReadDisk(db, "EGInput/xUnGCCR");
xungccr_df = ReadDisk(DataFrame, db, "EGInput/xUnGCCR");
xungccr_df[2124,:]
xungccr[2124,40]

wtf = J.var("EGInput/xUnGCCR", loc2)
@rsubset wtf :Value == 19.264
wtf2024 = @rsubset wtf :Year == 2024
findall(wtf2024.Value .== 19.264)
wtf2024[2124,:]
UnCode = ReadDisk(db, "EGInput/UnCode")
UnCode[2124]

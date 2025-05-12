using Revise
import JuliaCompare as J
import PromulaDBA as P
import SmallModel as M 
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta, TidierData

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

EuFPol = J.diff("EuFPol", loc1, loc2)
@rsubset! EuFPol :Area ∈ Canada
EuFPol = @mutate(EuFPol, PDiff = case_when(Diff == 0 => 0,
Spruce != 0 => Diff/Spruce,
Spruce == 0 => Diff/Tanoak))

issues = @rsubset EuFPol abs(:PDiff) > .2 :Diff > .1
sort(@by(issues, :ECC, :Count = length(:Diff)), :Count, rev= true)


ifelse.(EUFPol)
J.plot_diff(EuFPol; dim="ECC", num=10, title="EuFPol diffs by ECC")
# J.plot_diff(@rsubset EuFPol :FuelEP == "Biomass"; dim="ECC", num=10, title="Biomass diffs by ECC")
@rsubset! EuFPol :ECC != "ForeignPassenger" :ECC != "ForeignFreight" :ECC != "UtilityGen" :ECC != "ResidentialOffRoad"

# This is an issue with Spruce and I'll ignore it. 
@rsubset! EuFPol :ECC == "ResidentialOffRoad"
# Gasoline and Ethanol
J.plot_diff(EuFPol; dim="FuelEP", num=10, title="EuFPol diffs by FuelEP") 
# CO2 and COX, mainly
J.plot_diff(EuFPol; dim="Poll", num=10, title="EuFPol diffs by Poll") 
# QC and ON are the biggest contributors, though many are present
J.plot_diff(EuFPol; dim="Area", num=10, title="EuFPol diffs by Area") 
@by(EuFPol, [:Year], :Diff = sum(abs.(:Diff)))

dimension_filters= Dict{Symbol,Any}(:EC => "Passenger", :ECC => "Passenger")
EuFPol_p = @rsubset EuFPol :ECC == "Passenger"
# Gasoline and Diesel
J.plot_diff(EuFPol_p; dim="FuelEP", num=10, title="EuFPol_p diffs by FuelEP") 
# CO2 a little COX
J.plot_diff(EuFPol_p; dim="Poll", num=10, title="EuFPol_p diffs by Poll") 
# ON biggest contributors, though many are present
J.plot_diff(EuFPol_p; dim="Area", num=10, title="EuFPol_p diffs by Area") 

push!(dimension_filters, :FuelEP => ["Gasoline", "Diesel"])
push!(dimension_filters, :Poll => ["CO2", "COX"])
push!(dimension_filters, :Area => "ON")
push!(dimension_filters, :Year => "2030")

# # This is an issue with Spruce and I'll ignore it. 
# @rsubset! EuFPol :ECC == "ForeignPassenger"
# # Jet fuel
# J.plot_diff(EuFPol; dim="FuelEP", num=10, title="EuFPol diffs by FuelEP") 
# # All CO2 
# J.plot_diff(EuFPol; dim="Poll", num=10, title="EuFPol diffs by Poll") 
# # SK, ON, AB are the biggest contributors, though many are present
# J.plot_diff(EuFPol; dim="Area", num=10, title="EuFPol diffs by Area") 
# @by(EuFPol, [:Year], :Diff = sum(abs.(:Diff)))

Polute = J.diff_fast("Polute", loc1, loc2; sec = 'T', 
dimension_filters = dimension_filters)

# OffRoad 
Polute_p = J.var("TCalDB/Polute", loc1)
Polute_p = J.f_on(Polute_p)

Polute_j = J.var("TCalDB/Polute", loc2)
Polute_j = J.f_on(Polute_j)

Polute = J.diff(Polute_p, Polute_j; name1 = "Spruce", name2 = "Tanoak")
@rsubset Polute :Diff != 0

# 2×10 DataFrame
#  Row │ Enduse    FuelEP   Tech     EC                  Poll    Area    Year   Spruce    Tanoak    Diff           
#      │ String    String   String   String              String  String  Int64  Float64?  Float64?  Float64
# ─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────
#    1 │ Carriage  Ethanol  OffRoad  ResidentialOffRoad  COX     ON       2020   12222.9   12222.9     0.000191788
#    2 │ Carriage  Ethanol  OffRoad  ResidentialOffRoad  COX     ON       2021   13027.8   13185.5  -157.682

DmdFEPTech = J.diff("TOutput/DmdFEPTech", loc1, loc2)
DmdFEPTech = J.f_on(DmdFEPTech)
import SmallModel: ReadDisk

POCA = P.data(joinpath(DATA_FOLDER1,"TOutput2.dba"),"POCA");
POCA_j = M.ReadDisk(loc2.HDF5_path, "TOutput/POCA");
sets = M.ReadSets(loc2.HDF5_path, "TOutput/POCA")
dimension_filters= Dict(:EC => "ResidentialOffRoad", :FuelEP => "Ethanol", :Tech => "OffRoad", :Area => "ON", :Year => ["2020","2021"])
POCA_p, set_p = J.subset_array(POCA, sets, dimension_filters)
POCA_dfp = J.to_tidy_dataframe(POCA_p, set_p)
POCA_j, set_j = J.subset_array(POCA_j, sets, dimension_filters)
POCA_dfj = J.to_tidy_dataframe(POCA_j, set_j)
POCA_df = J.diff(POCA_dfp, POCA_dfj; name1 = "Spruce", name2 = "Tanoak")
@rsubset! POCA_df :Diff != 0
sort(POCA_df, [:Year, :Diff])

POCX = P.data(joinpath(BASE_FOLDER,"2020Model","Calib3","TInput.dba"),"POCX");
POCX = P.data(joinpath(DATA_FOLDER1,"TInput.dba"),"POCX");
POCX = P.data(joinpath(BASE_FOLDER,"2020Model","Calib4","TInput.dba"),"POCX");
POCX_j = M.ReadDisk(loc2.HDF5_path, "TInput/POCX");
sets = M.ReadSets(loc2.HDF5_path, "TInput/POCX")
dimension_filters= Dict(:EC => "ResidentialOffRoad", :FuelEP => ["Ethanol", "Gasoline"], 
  :Tech => "OffRoad", :Area => "ON", :Year => ["2020","2021"], :Poll => ["COX", "CO2"])
POCX_p, set_p = J.subset_array(POCX, sets, dimension_filters)
POCX_dfp = J.to_tidy_dataframe(POCX_p, set_p)
POCX_j, set_j = J.subset_array(POCX_j, sets, dimension_filters)
POCX_dfj = J.to_tidy_dataframe(POCX_j, set_j)
POCX_df = J.diff(POCX_dfp, POCX_dfj; name1 = "Spruce", name2 = "Tanoak")
@rsubset! POCX_df :Diff != 0
sort(POCX_df, [:Year, :Diff])

xTrEnFPol = P.data(joinpath(BASE_FOLDER,"2020Model","Calib3","TInput.dba"),"xTrEnFPol");
xTrEnFPol = P.data(joinpath(DATA_FOLDER1,"TInput.dba"),"xTrEnFPol");
xTrEnFPol = P.data(joinpath(BASE_FOLDER,"2020Model","Calib4","TInput.dba"),"xTrEnFPol");
xTrEnFPol_j = M.ReadDisk(loc2.HDF5_path, "TInput/xTrEnFPol");
sets = M.ReadSets(loc2.HDF5_path, "TInput/xTrEnFPol")
dimension_filters= Dict(:EC => "ResidentialOffRoad", :FuelEP => ["Ethanol", "Gasoline"], 
  :Tech => "OffRoad", :Area => "ON", :Year => ["2020","2021"], :Poll => ["COX", "CO2"])
xTrEnFPol_p, set_p = J.subset_array(xTrEnFPol, sets, dimension_filters)
xTrEnFPol_dfp = J.to_tidy_dataframe(xTrEnFPol_p, set_p)
xTrEnFPol_j, set_j = J.subset_array(xTrEnFPol_j, sets, dimension_filters)
xTrEnFPol_dfj = J.to_tidy_dataframe(xTrEnFPol_j, set_j)
xTrEnFPol_df = J.diff(xTrEnFPol_dfp, xTrEnFPol_dfj; name1 = "Spruce", name2 = "Tanoak")
@rsubset! xTrEnFPol_df :Diff != 0
sort(xTrEnFPol_df, [:Year, :Diff])

xOREnFPol = P.data(joinpath(BASE_FOLDER,"2020Model","Calib3","SInput.dba"),"xOREnFPol");
xOREnFPol = P.data(joinpath(BASE_FOLDER,"2020Model","Calib4","SInput.dba"),"xOREnFPol");
xOREnFPol = P.data(joinpath(DATA_FOLDER1,"SInput.dba"),"xOREnFPol");
xOREnFPol_j = M.ReadDisk(loc2.HDF5_path, "SInput/xOREnFPol");
sets = M.ReadSets(loc2.HDF5_path, "SInput/xOREnFPol")
dimension_filters= Dict(:EC => "ResidentialOffRoad", :ECC => "ResidentialOffRoad", :FuelEP => ["Ethanol", "Gasoline"], 
  :Tech => "OffRoad", :Area => "ON", :Year => ["2020","2021"], :Poll => ["COX", "CO2"])
xOREnFPol_p, set_p = J.subset_array(xOREnFPol, sets, dimension_filters)
xOREnFPol_dfp = J.to_tidy_dataframe(xOREnFPol_p, set_p)
xOREnFPol_j, set_j = J.subset_array(xOREnFPol_j, sets, dimension_filters)
xOREnFPol_dfj = J.to_tidy_dataframe(xOREnFPol_j, set_j)
xOREnFPol_df = J.diff(xOREnFPol_dfp, xOREnFPol_dfj; name1 = "Spruce", name2 = "Tanoak")
@rsubset! xOREnFPol_df :Diff != 0
sort(xOREnFPol_df, [:Year, :Diff])

push!(dimension_filters, :Poll => "COX")

POEM = P.data(joinpath(DATA_FOLDER1,"TOutput2.dba"),"POEM");
POEM_j = M.ReadDisk(loc2.HDF5_path, "TOutput/POEM");
sets = M.ReadSets(loc2.HDF5_path, "TOutput/POEM")
POEM_p, set_p = J.subset_array(POEM, sets, dimension_filters)
POEM_dfp = J.to_tidy_dataframe(POEM_p, set_p)
POEM_j, set_j = J.subset_array(POEM_j, sets, dimension_filters)
POEM_dfj = J.to_tidy_dataframe(POEM_j, set_j)
POEM_df = J.diff(POEM_dfp, POEM_dfj; name1 = "Spruce", name2 = "Tanoak")
@rsubset! POEM_df :Diff != 0
sort(POEM_df, [:Year, :Diff])

POEMA = P.data(joinpath(DATA_FOLDER1,"TOutput2.dba"),"POEMA");
POEMA_j = M.ReadDisk(loc2.HDF5_path, "TOutput/POEMA");
sets = M.ReadSets(loc2.HDF5_path, "TOutput/POEMA")
POEMA_p, set_p = J.subset_array(POEMA, sets, dimension_filters)
POEMA_dfp = J.to_tidy_dataframe(POEMA_p, set_p)
POEMA_j, set_j = J.subset_array(POEMA_j, sets, dimension_filters)
POEMA_dfj = J.to_tidy_dataframe(POEMA_j, set_j)
POEMA_df = J.diff(POEMA_dfp, POEMA_dfj; name1 = "Spruce", name2 = "Tanoak")
@rsubset! POEMA_df :Diff != 0
sort(POEMA_df, [:Year, :Diff])

DER = P.data(joinpath(DATA_FOLDER1,"TOutput.dba"),"DER");
DER_j = M.ReadDisk(loc2.HDF5_path, "TOutput/DER");
sets = M.ReadSets(loc2.HDF5_path, "TOutput/DER")
DER_p, set_p = J.subset_array(DER, sets, dimension_filters)
DER_dfp = J.to_tidy_dataframe(DER_p, set_p)
DER_j, set_j = J.subset_array(DER_j, sets, dimension_filters)
DER_dfj = J.to_tidy_dataframe(DER_j, set_j)
DER_df = J.diff(DER_dfp, DER_dfj; name1 = "Spruce", name2 = "Tanoak")
@rsubset! DER_df :Diff != 0
sort(DER_df, [:Year, :Diff])

RM = J.diff("TOutput/RM", loc1, loc2)
RM = J.f_on(RM)

StockAdjustment = J.diff("TInput/StockAdjustment", loc1, loc2)
StockAdjustment = J.f_on(StockAdjustment)

POCA = J.diff("TOutput/POCA", loc1, loc2)
POCA = J.f_on(POCA)
# UtilityGen

# This is an issue with xUnGCCR that I have to get to the bottom of. 
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
UnPolGross = ReadDisk(db, "EGOutput/UnPolGross")
UnArea = M.ReadDisk(db,"EGInput/UnArea")
UnCode = M.ReadDisk(db,"EGInput/UnCode")
M.WriteDisk(loc2.HDF5_path, "EGInput/Unit", UnCode)
M.WriteDisk(loc2.HDF5_path, "EGOutput/Unit", UnCode)
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
sort(issues, :Diff)
issue_codes = unique(issues.Unit)

UnGC = J.diff("UnGC", loc1, loc2)
df = @rsubset UnGC :Year ∈ [2024,2025] :Unit ∈ issue_codes abs(:Diff) != 0;
sort(df, [:Year, :Diff])

UnGCCR = J.diff("UnGCCR", loc1, loc2)
df = @rsubset! UnGCCR :Year ∈ [2024,2025] :Unit ∈ issue_codes abs(:Diff) != 0;
sort(df, :Diff)

UnGCCE = J.diff("UnGCCE", loc1, loc2)
df = @rsubset! UnGCCE :Year ∈ [2024,2025] :Unit ∈ issue_codes abs(:Diff) != 0;
sort(df, [:Year, :Diff])
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

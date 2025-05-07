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
dimension_filters = Dict{Symbol,Any}()
EuFPol = J.diff_fast("EuFPol", loc1, loc2; dimension_filters)
@rsubset! EuFPol :Area ∈ Canada
J.plot_diff(EuFPol; dim="ECC", num=10, title="EuFPol diffs by ECC")
# J.plot_diff(@rsubset EuFPol :FuelEP == "Biomass"; dim="ECC", num=10, title="Biomass diffs by ECC")


@rsubset! EuFPol :ECC == "Passenger"
# Almost all Gasoline and Diesel
J.plot_diff(EuFPol; dim="FuelEP", num=10, title="EuFPol diffs by FuelEP") 
push!(dimension_filters, :FuelEP => ["Gasoline","Diesel"])

# CO2 mainly with a little COX
J.plot_diff(EuFPol; dim="Poll", num=10, title="EuFPol diffs by Poll") 
push!(dimension_filters, :Poll => ["CO2","COX"])
# ON, AB, BC are the biggest contributors, though many are present
J.plot_diff(EuFPol; dim="Area", num=10, title="EuFPol diffs by Area")
push!(dimension_filters, :Area => "ON")

# ON, AB, BC are the biggest contributors, though many are present
J.plot_diff(EuFPol; dim="Tech", num=10, title="EuFPol diffs by Tech")
push!(dimension_filters, :Area => "ON")

Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters)
sec = 'T'
Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters, sec)
J.plot_diff(Polute; dim="Tech", num=10, title="Polute diffs by Tech")
push!(dimension_filters, :Tech => ["LDVDiesel","LDVGasoline"])
push!(dimension_filters, :Year => 2020:2050)

POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters, sec)
i = findall(vars.Variable .== "POCA" .&& first.(vars.Database) .== 'T')
i = i[1]
vars[i,:Database] = "TOutput2"
dimension_filters[:Year] = "2040"
POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters, sec)
@rsubset POCA :Diff != 0

POCX = J.diff_fast("POCX", loc1, loc2; dimension_filters, sec)
@rsubset POCX :Diff != 0

DmdFEPTech = J.diff_fast("DmdFEPTech", loc1, loc2; dimension_filters, sec)
@rsubset DmdFEPTech isapprox(:Diff, 0)

dim_f = Dict(:Area => "ON") 
@btime POCX = J.diff_fast("POCX", loc1, loc2; dimension_filters = dim_f, sec)
@rsubset POCX :EC == "Passenger" :Poll == "COX" !isapprox(:Diff, 0, atol=1e-3)
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

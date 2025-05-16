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
dimension_filters = Dict{Symbol,Any}()

      # PEDC[ecc,area,year] = ((xPE[ecc,area,year]*xInflation[area,year]-
      #   FPTaxF[fuel,es,area,year]*EEConv/1000*xInflation[area,year])/(1+FPSMF[fuel,es,area,year])-
      #   PPUC[area,year-1]+ExportsPE[area,year-1])/xInflation[area,year]

sec = 'E'
push!(dimension_filters, :Area => "ON")
push!(dimension_filters, :ECC => "SingleFamilyDetached")
push!(dimension_filters, :Year => ["2024", "2025"])
PEDC = J.diff_fast("PEDC", loc1, loc2; dimension_filters, sec)
xPE = J.diff_fast("xPE", loc1, loc2; dimension_filters, sec)
xInflation = J.diff_fast("xInflation", loc1, loc2; dimension_filters, sec)
FPTaxF = J.diff_fast("FPTaxF", loc1, loc2; dimension_filters, sec)
@rsubset FPTaxF abs(:Diff) >= 0.001
EEConv = J.diff("SInput/EEConv", loc1, loc2)

M.ReadDisk(db,"SInput/EEConv")
EEConv::Float64 = M.ReadDisk(db,"SInput/EEConv")[1] # Electric Energy Conversion (Btu/KWh)
FPSMF = J.diff_fast("FPSMF", loc1, loc2; dimension_filters, sec)
@rsubset FPSMF abs(:Diff) >= 0.001
push!(dimension_filters, :Year => string.(2023:2025))
PPUC = J.diff_fast("PPUC", loc1, loc2; dimension_filters, sec) # causing the primary diff
@rsubset PPUC abs(:Diff) >= 0.001
ExportsPE = J.diff_fast("ExportsPE", loc1, loc2; dimension_filters, sec) # bad

# @. @finite_math PPUC[Areas] = PUCT[Areas]/TSales[Areas]*1000
PUCT = J.diff_fast("PUCT", loc1, loc2; dimension_filters, sec) # causing the primary diff
@rsubset PUCT abs(:Diff) >= 0.001
TSales = J.diff_fast("TSales", loc1, loc2; dimension_filters, sec) # causing the primary diff
@rsubset TSales abs(:Diff) >= 0.001

#       PUCT[area] = PUCTBI[area]+PUCTSM[area]
PUCTBI = J.diff_fast("PUCTBI", loc1, loc2; dimension_filters, sec) # causing the primary diff
@rsubset PUCTBI abs(:Diff) >= 0.001
PUCTSM = J.diff_fast("PUCTSM", loc1, loc2; dimension_filters, sec) # causing the primary diff
@rsubset PUCTSM abs(:Diff) >= 0.001
PPCT = J.diff_fast("PPCT", loc1, loc2; dimension_filters, sec) # causing the primary diff
@rsubset PPCT abs(:Diff) >= 0.001

# PPCT[ppset,area] = sum(PPEGTM[area,timep,month]*PPVCost[ppset,month,area] for
#   month in Months, timep in TimePs)/1000
push!(dimension_filters, :PPSet => "BasePP")
PPVCost = J.diff_fast("EOutput/PPVCost", loc1, loc2; dimension_filters, sec) # bad
PPEGTM = J.diff_fast("PPEGTM", loc1, loc2; dimension_filters, sec) # causing the primary diff
@by(PPEGTM, :Year, :Diff = sum(:Diff))
push!(dimension_filters, :TimeP => "TimeP6")
push!(dimension_filters, :Month => "Winter")
push!(dimension_filters, :Year => "2050")
CnEG = J.diff_fast("CnEG", loc1, loc2; dimension_filters, sec) # causing the primary diff
@rsubset CnEG :Diff != 0

dimension_filters[:Year] = ["2024","2025"]
Capacity = J.diff_fast("Capacity", loc1, loc2; dimension_filters, sec) # causing the primary diff
@rsubset Capacity :Spruce != 0 || :Tanoak != 0
xCapacity = J.diff_fast("xCapacity", loc1, loc2; dimension_filters, sec) # causing the primary diff
@rsubset xCapacity :Spruce != 0 || :Tanoak != 0
df_20s = Dict{Symbol,Any}()
df_20s[:Year] = ["2023", "2024", "2025"]
xCapacity = J.diff_fast("xCapacity", loc1, loc2; dimension_filters=df_20s, sec) # causing the primary diff
@rsubset xCapacity (:Spruce != 0 || :Tanoak != 0)
@rsubset xCapacity :Diff != 0

######################
# Still need PUCTBI solved
######################
DCCost = J.diff_fast("DCCost", loc1, loc2; dimension_filters, sec) # a little off
DVCost = J.diff_fast("DVCost", loc1, loc2; dimension_filters, sec) # 

# DVCost[area] =
#   sum(EGBI[area,genco,plant]*UVCost[area,genco,plant] for plant in Plants, genco in GenCos) /
#   1000*ExchangeRate[area]
ExchangeRate = J.diff_fast("ExchangeRate", loc1, loc2; dimension_filters, sec) # 
EGBI = J.diff_fast("EGBI", loc1, loc2; dimension_filters, sec) # 
@rsubset EGBI :Diff != 0
@by(EGBI, :Year, :Diff = sum(:Diff))
UVCost = J.diff_fast("UVCost", loc1, loc2; dimension_filters, sec) # 
@rsubset UVCost :Diff != 0 :Plant == "OGCT"
@by(UVCost, :Year, :Diff = sum(:Diff))
UECOST = J.diff_fast("UECOST", loc1, loc2; dimension_filters, sec) # 
@rsubset UECOST :Diff != 0 :Plant == "OGCT"
@by(UECOST, :Year, :Diff = sum(:Diff))

######################
# Still need TSales solved
######################

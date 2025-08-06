using Revise
import JuliaCompare as J
import PromulaDBA as P
import SmallModel as M 
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

################################################################################
# Inputs for file locations ####################################################
################################################################################
BASE_FOLDER = raw"\\Pink\c\2020CanadaSpruce"
BASE_FOLDER2 = raw"\\Pink\c\2020CanadaTanoak"
SCENARIO1 = ""
SCENARIO2 = ""
################################################################################

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
tan_base = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaTanoak\\2020Model\\Base\\database.hdf5", "Tan_Base");
spr_base = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaSpruce\\2020Model\\Base", "Spr_Base");
tan_ogref = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaTanoak\\2020Model\\OGRef\\database.hdf5", "Tan_OGRef");
spr_ogref = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaSpruce\\2020Model\\OGRef", "Spr_OGRef");

################################################################################
# Initialize variables for analysis ############################################
################################################################################
sec = 'S'
dimension_filters = Dict{Symbol,Any}()

################################################################################
# Analysis of variables ########################################################
################################################################################
DevSw = J.diff_fast("DevSw", loc1, loc2; dimension_filters, sec) # 
@rsubset! DevSw :Diff != 0
using Statistics
@by(DevSw, :OGUnit, :begin = minimum(:Year), :end = maximum(:Year), :Spruce = mean(:Spruce), :Tanoak = mean(:Tanoak))
# 7×5 DataFrame
#  Row │ OGUnit               begin  end    Spruce   Tanoak  
#      │ String               Int64  Int64  Float64  Float32
# ─────┼─────────────────────────────────────────────────────
#    1 │ NL_HeavyOil_0001      2021   2021      2.0      0.0
#    2 │ NS_Gas_0001           2023   2046      0.0      2.0
#    3 │ SK_ConvGas_0001       2023   2048      0.0      2.0
#    4 │ YT_Gas_0001           2023   2039      0.0      2.0
#    5 │ ON_Gas_0001           2024   2024      0.0      2.0
#    6 │ NB_Gas_0001           2025   2050      0.0      2.0
#    7 │ SK_OS_Upgrader_0001   2027   2029      2.0      0.0
push!(dimension_filters,:OGUnit => "SK_OS_Upgrader_0001")
push!(dimension_filters,:Year => string.(2026:2030))
DevROIRef = J.diff_fast("DevROI", spr_ogref, tan_ogref; dimension_filters, sec) # 
OGRefNameDB = M.ReadDisk(loc2.HDF5_path,"E2020DB/OGRefNameDB") #  Oil/Gas Reference Case Name

push!(dimension_filters, :OGUnit => "AB_OS_SAGD_0001")
push!(dimension_filters, :Year => "2050")
pop!(dimension_filters, :Fuel)
OGFPrice = J.diff_fast("OGFPrice", loc1, loc2; dimension_filters, sec) 
J.add_pdiff!(OGFPrice)
@rsubset! OGFPrice :Diff != 0
sort!(OGFPrice, :PDiff)
OGPolCost = J.diff_fast("OGPolCost", loc1, loc2; dimension_filters, sec) 
J.add_pdiff!(OGPolCost)
# OGPolCosts[ogunit] = sum(OGFUse[ogunit,fuel]*OGPolPrice[ogunit,fuel] for fuel in Fuels)+OGFuCosts[ogunit]
OGFUse = J.diff_fast("OGFUse", loc1, loc2; dimension_filters, sec) 
J.add_pdiff!(OGFUse)
@rsubset OGFUse :Spruce != 0 || :Tanoak != 0
push!(dimension_filters, :Year => "2024")
sec = 'I'
DemRqPrior = J.diff_fast("DemRq", loc1, loc2; dimension_filters, sec) 

# DemRq[fuel,ecc,area] = sum(DmdRq[enduse,tech,ec,area]*
#   DmFrac[enduse,fuel,tech,ec,area] for enduse in Enduses,tech in Techs)
DmdRq = J.diff_fast("DmdRq", loc1, loc2; dimension_filters, sec) # Explains difference
push!(dimension_filters, :Enduse => "Heat", :Tech => "Gas")
DmFrac = J.diff_fast("DmFrac", loc1, loc2; dimension_filters, sec) # Matches

# @finite_math DmdRq[enduse,tech,ec,area] = DSt[enduse,ec,area]/PEE[enduse,tech,ec,area]/
#   DEE[enduse,tech,ec,area]*AMSF[enduse,tech,ec,area]*CERSM[enduse,ec,area]*ECUF[ecc,area]*
#   RPEI[enduse,tech,ec,area]/1.0e6
DSt = J.diff_fast("DSt", loc1, loc2; dimension_filters, sec) # match
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters, sec) # close
DEE = J.diff_fast("DEE", loc1, loc2; dimension_filters, sec) # close
AMSF = J.diff_fast("AMSF", loc1, loc2; dimension_filters, sec) # match
CERSM = J.diff_fast("CERSM", loc1, loc2; dimension_filters, sec) # close
ECUF = J.diff_fast("ECUF", loc1, loc2; dimension_filters, sec) # match
RPEI = J.diff_fast("RPEI", loc1, loc2; dimension_filters, sec) # match

ECC = pop!(dimension_filters,:ECC)
ECUF_test = J.diff_fast("ECUF", loc1, loc2; dimension_filters, sec) # match
ECUF_s2 = ECUF_test[1,:Spruce]
push!(dimension_filters, :ECC => ECC)

J.add_pdiff(DmdRq)
J.add_pdiff(PEE)
J.add_pdiff(DEE)
J.add_pdiff(CERSM)

DSt_s = DSt[1,:Spruce]
PEE_s = PEE[1,:Spruce]
DEE_s = DEE[1,:Spruce]
AMSF_s = AMSF[1,:Spruce]
CERSM_s = CERSM[1,:Spruce]
ECUF_s = ECUF[1,:Spruce]
RPEI_s = RPEI[1,:Spruce]

DSt_j = DSt[1,:Tanoak]
PEE_j = PEE[1,:Tanoak]
DEE_j = DEE[1,:Tanoak]
AMSF_j = AMSF[1,:Tanoak]
CERSM_j = CERSM[1,:Tanoak]
ECUF_j = ECUF[1,:Tanoak]
RPEI_j = RPEI[1,:Tanoak]


DmdRq_s = DSt_s/PEE_s/DEE_s*AMSF_s*CERSM_s*ECUF_s*RPEI_s/1.0e6
DmdRq_s = DSt_s/PEE_s/DEE_s*AMSF_s*CERSM_s*ECUF_s2*RPEI_s/1.0e6
  DmdRq# =DSt/PEE/DEE*AMSF*CERSM*ECUF*RPEI/1E6

DmdRq_j = DSt_j/PEE_j/DEE_j*AMSF_j*CERSM_j*ECUF_j*RPEI_j/1.0e6


# OGPolPrice
OGPolPrice = J.diff_fast("OGPolPrice", loc1, loc2; dimension_filters, sec) 

#       OGPolPrice[ogunit,fuel] = FPECCCPNet[fuel,ecc,area]+FPECCCFSNet[fuel,ecc,area]+FPECCOGEC[fuel,ecc,area]
push!(dimension_filters, :ECC => "SAGDOilSands", :EC => "SAGDOilSands", :Fuel => "NaturalGas", :Area => "AB")
FPECCCPNet = J.diff_fast("FPECCCPNet", loc1, loc2; dimension_filters, sec) 
FPECCCFSNet = J.diff_fast("FPECCCFSNet", loc1, loc2; dimension_filters, sec) 
FPECCOGEC = J.diff_fast("FPECCOGEC", loc1, loc2; dimension_filters, sec) 

@rsubset! OGPolPrice abs(:Diff) > 1e-5
J.add_pdiff!(OGPolPrice)
# VnCosts[ogunit]+FuCosts[ogunit]+FlCosts[ogunit]
VnCosts = J.diff_fast("VnCosts", loc1, loc2; dimension_filters, sec) 
@rsubset VnCosts abs(:Spruce) + abs(:Tanoak) > 0
FuCosts = J.diff_fast("FuCosts", loc1, loc2; dimension_filters, sec) 
FlCosts = J.diff_fast("FlCosts", loc1, loc2; dimension_filters, sec) 


OGArea = M.ReadDisk(loc2.HDF5_path, "SpInput/OGArea")
OGUnit = M.ReadDisk(loc2.HDF5_path, "SpInput/OGUnit")
Area_df = DataFrame(OGUnit = OGUnit, OGArea = OGArea)
leftjoin!(DevSw, Area_df, on = :OGUnit)
@rsubset DevSw :OGArea == "CA"
DevSw_b = J.diff_fast("DevSw", spr_base, tan_base; dimension_filters, sec) # 
@rsubset! DevSw_b :Diff != 0

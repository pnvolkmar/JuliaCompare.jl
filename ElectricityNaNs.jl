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
sec = 'I'
dimension_filters = Dict{Symbol,Any}()
push!(dimension_filters, :Area => Canada, :Year => string.(2023:2024))
################################################################################
# Analysis of variables ########################################################
################################################################################
EGFA = J.diff_fast("EGFA", loc1, loc2; dimension_filters, sec) # 
J.plot_diff(EGFA, dim="Area")
@rsubset! EGFA isnan(:Diff)
unique(EGFA.Fuel)
UnEGA = J.diff_fast("UnEGA", loc1, loc2; dimension_filters, sec) # 
@rsubset! UnEGA isnan(:Diff)
Fuel = M.ReadDisk(loc2.HDF5_path,"E2020DB/FuelKey")


AGFr = J.diff_fast("AGFr", loc1, loc2; dimension_filters, sec) # good
@rsubset! AGFr isnan(:Diff)
UnRCR = J.diff_fast("UnRCR", loc1, loc2; dimension_filters, sec) # good
@rsubset! UnRCR isnan(:Diff)


J.plot_diff(EuFPol, dim="FuelEP")
J.add_pdiff!(EuFPol)
EuFPol_pd = select(EuFPol, :FuelEP, :ECC, :Poll, :Area, :Year, :Spruce, :Tanoak, :PDiff => :Diff)
J.plot_diff(EuFPol_pd, dim="ECC")
@by(EuFPol, [:Year, :ECC], :PDiff = sum(:PDiff))
sort!(EuFPol, :PDiff)
@rsubset EuFPol :PDiff > 100

i = findall(vars.Variable .== "POCA" .&& first.(vars.Database) .== 'I')
i = i[1]
vars[i,:Database] = "IOutput2"

pop!(dimension_filters, :Poll)
POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters, sec) # Issue
J.add_pdiff!(POCA)
@rsubset! POCA abs(:PDiff) > 20

on = copy(dimension_filters)
push!(on, :Enduse => ["Heat", "Motors", "OthSub"])
push!(on, :FuelEP => ["HFO", "NaturalGas", "StillGas"])
push!(on, :EC => ["Petroleum"])
push!(on, :Area => ["ON"])

nb = copy(dimension_filters)
push!(nb, :EC => ["OtherNonferrous"])
push!(nb, :Area => ["NB"])
push!(nb, :Poll => ["SOX","PMT","PM25","PM10","BC"])

nl = copy(dimension_filters)
push!(nl, :Enduse => ["Heat", "Motors", "OthSub"])
push!(nl, :FuelEP => ["HFO", "LFO", "LPG"])
push!(nl, :EC => ["IronOreMining"])
push!(nl, :Area => ["NL"])

PolSw = J.diff_fast("PolSw", loc1, loc2; dimension_filters, sec) # matches, 2.0 everywhere
PolSw = J.diff_fast("PolSw", loc1, loc2; dimension_filters=on, sec) # matches, 2.0 everywhere
PolSw = J.diff_fast("PolSw", loc1, loc2; dimension_filters=nl, sec) # matches, 2.0 everywhere

POCX = J.diff_fast("POCX", loc1, loc2; dimension_filters, sec) # matches, everywhere
@rsubset POCX abs(:Diff) > 1e-1
POCX = J.diff_fast("POCX", loc1, loc2; dimension_filters=on, sec) # matches, everywhere
POCX = J.diff_fast("POCX", loc1, loc2; dimension_filters=nl, sec) # matches, everywhere

RM = J.diff_fast("RM", loc1, loc2; dimension_filters, sec) # issue!
@rsubset RM :Tanoak == 0
RM = J.diff_fast("RM", loc1, loc2; dimension_filters=nb, sec) # issue!
RM = J.diff_fast("RM", loc1, loc2; dimension_filters=nl, sec) # issue!

# @. RM = (1.0-RP)*(1.0-VR)*xRM
xRM = J.diff_fast("xRM", loc1, loc2; dimension_filters=nb, sec) # match
@rsubset xRM abs(:Diff) > 1e-1
RP = J.diff_fast("RP", loc1, loc2; dimension_filters, sec) # issue
@rsubset RP :Tanoak == 1

i = findall(vars.Variable .== "VR" .&& first.(vars.Database) .== 'I')
i = i[1]
vars[i,:Database] = "IOutput2"

VR = J.diff_fast("VR", loc1, loc2; dimension_filters=nb, sec) # match
@rsubset VR abs(:Diff) > 1e-1

          # RP[fuelep,ec,poll,area] = max(RPolicy[ec,poll,area],
          #   (1-minimum(POCS[:,fuelep,ec,poll,area])/
          #   max(maximum(POCX[:,fuelep,ec,poll,area]),0.000001)))*
          #   PCCov[fuelep,ec,poll,area]

RPolicy = J.diff_fast("RPolicy", loc1, loc2; dimension_filters=on, sec = sec)
RPolicy = J.diff_fast("RPolicy", loc1, loc2; dimension_filters=nl, sec = sec) # this is addressed

RPolicy = J.diff_fast("RPolicy", loc1, loc2; dimension_filters=on, sec = 'S')
@rsubset! RPolicy abs(:Diff) > 0
sort(RPolicy, :Diff)
# RPolECC = J.diff_fast("RPolECC", loc1, loc2; dimension_filters=nl, sec = sec) # this is addressed

push!(on, :Market => "Market180")
push!(on, :Year => string.(2020:2027))
xGoalPol = J.diff_fast("xGoalPol", loc1, loc2; dimension_filters=on, sec = sec)
AreaMarket = J.diff_fast("AreaMarket", loc1, loc2; dimension_filters=on, sec = sec)
ECCMarket = J.diff_fast("ECCMarket", loc1, loc2; dimension_filters=on, sec = sec)
@rsubset ECCMarket :Diff != 0
CapTrade = J.diff_fast("CapTrade", loc1, loc2; dimension_filters=on, sec = sec)
ECoverage = J.diff_fast("ECoverage", loc1, loc2; dimension_filters=on, sec = sec)
@rsubset ECoverage :Diff != 0
PCOVMarket = J.diff_fast("PCOVMarket", loc1, loc2; dimension_filters=on, sec = sec)
@rsubset PCOVMarket :Diff != 0
PollMarket = J.diff_fast("PollMarket", loc1, loc2; dimension_filters=on, sec = sec)
@rsubset PollMarket :Diff != 0

push!(nl, :Market => "Market91")
push!(nl, :Year => string.(2020:2027))
xGoalPol = J.diff_fast("xGoalPol", loc1, loc2; dimension_filters=nl, sec = sec)
AreaMarket = J.diff_fast("AreaMarket", loc1, loc2; dimension_filters=nl, sec = sec)
ECCMarket = J.diff_fast("ECCMarket", loc1, loc2; dimension_filters=nl, sec = sec)
@rsubset ECCMarket :Diff != 0
CapTrade = J.diff_fast("CapTrade", loc1, loc2; dimension_filters=nl, sec = sec)
ECoverage = J.diff_fast("ECoverage", loc1, loc2; dimension_filters=nl, sec = sec)
@rsubset ECoverage :Diff != 0
PCOVMarket = J.diff_fast("PCOVMarket", loc1, loc2; dimension_filters=nl, sec = sec)
@rsubset PCOVMarket :Diff != 0
PollMarket = J.diff_fast("PollMarket", loc1, loc2; dimension_filters=nl, sec = sec)
@rsubset PollMarket :Diff != 0

# @finite_math POCA[enduse,fuelep,tech,ec,poll,area] = POEM[enduse,fuelep,tech,ec,poll,area]/
#   DER[enduse,tech,ec,area]*1E6*RM[tech,ec,poll,area]

GoalPol = J.diff_fast("GoalPol", loc1, loc2; dimension_filters=nl, sec = sec) # match
GoalPol = J.diff_fast("GoalPol", loc1, loc2; dimension_filters=on, sec = sec) # match

GoalPolSw = J.diff_fast("GoalPolSw", loc1, loc2; dimension_filters=nl, sec = sec)
GoalPolSw = J.diff_fast("GoalPolSw", loc1, loc2; dimension_filters=on, sec = sec)


i = findall(vars.Variable .== "POEM" .&& first.(vars.Database) .== 'T')
i = i[1]
vars[i,:Database] = "TOutput2"
POEM = J.diff_fast("POEM", loc1, loc2; dimension_filters, sec) # close
DER = J.diff_fast("DER", loc1, loc2; dimension_filters, sec) # close
RM = J.diff_fast("RM", loc1, loc2; dimension_filters, sec) # Issue

# @. RM = (1.0-RP)*(1.0-VR)*xRM
RP = J.diff_fast("RP", loc1, loc2; dimension_filters, sec) # close
VR = J.diff_fast("VR", loc1, loc2; dimension_filters, sec) # close
xRM = J.diff_fast("xRM", loc1, loc2; dimension_filters, sec) # close

################################################################################
# Initialize variables for analysis ############################################
################################################################################
sec = 'E'
dimension_filters = Dict{Symbol,Any}()
push!(dimension_filters, :Year => "2022", :Poll => "BC")
push!(dimension_filters, :ECC => "UtilityGen", :FuelEP => "Diesel", :Area => "YT")
push!(dimension_filters, :Tech => "LDVGasoline", :EC => "Passenger")
################################################################################
# Analysis of variables ########################################################
################################################################################
EuFPol = J.diff_fast("EuFPol", loc1, loc2; dimension_filters, sec) # 
@rsubset! EuFPol :Diff != 0
J.add_pdiff!(EuFPol)
sort!(EuFPol, :PDiff)
@rsubset EuFPol :PDiff > 100


i = findall(vars.Variable .== "UnPolGross" .&& first.(vars.Database) .== 'E')
i = i[1]
vars[i,:Database] = "EGOutput3"
UnPolGross = J.diff_fast("UnPolGross", loc1, loc2; dimension_filters, sec) # 
@rsubset! UnPolGross :Diff != 0
UnArea = J.var("UnArea", loc1)
leftjoin!(UnPolGross,UnArea,on = :Unit)
@rsubset! UnPolGross :UnArea == "YT"
push!(dimension_filters, :Unit => UnPolGross.Unit)

# UnPolGross[unit,fuelep,poll] = UnDmd[unit,fuelep]*UnPOCA[unit,fuelep,poll]
UnDmd = J.diff_fast("UnDmd", loc1, loc2; dimension_filters, sec) # Tanoak is zero
UnPOCA = J.diff_fast("UnPOCA", loc1, loc2; dimension_filters, sec) # Matches

# UnDmd[unit,fuelep] = max(UnEGGross[unit]*UnHRt[unit]/1e6*UnFlFr[unit,fuelep],0)
UnEGGross = J.diff_fast("UnEGGross", loc1, loc2; dimension_filters, sec) # Tanoak is zero
UnHRt = J.diff_fast("UnHRt", loc1, loc2; dimension_filters, sec) # matches
UnFlFr = J.diff_fast("UnFlFr", loc1, loc2; dimension_filters, sec) # matches
UnEGA = J.diff_fast("UnEGA", loc1, loc2; dimension_filters, sec) # Tnaoak is zero
UnGC = J.diff_fast("UnGC", loc1, loc2; dimension_filters, sec) # matches
UnOnline = J.diff_fast("UnOnline", loc1, loc2; dimension_filters, sec) # matches
UnRetire = J.diff_fast("UnRetire", loc1, loc2; dimension_filters, sec) # Tanoak retires in 2019

################################################################################
# Initialize variables for analysis ############################################
################################################################################
sec = 'I'
dimension_filters = Dict{Symbol,Any}()
push!(dimension_filters, :Area => Canada)
push!(dimension_filters, :Year => string.(2000:2020), :Poll => "VOC")
push!(dimension_filters, :ECC => "Petroleum", :Area => "NS")
push!(dimension_filters, :Tech => "LDVGasoline", :EC => "Passenger")
################################################################################
# Analysis of variables ########################################################
################################################################################
EuFPol = J.diff_fast("EuFPol", loc1, loc2; dimension_filters, sec) # 
@rsubset! EuFPol :Diff != 0
J.add_pdiff!(EuFPol)

sort!(EuFPol, :PDiff)
J.plot_diff(EuFPol, dim = :ECC)

MEPol = J.diff_fast("MEPol", loc1, loc2; dimension_filters, sec)
@rsubset! MEPol :Diff != 0
J.add_pdiff!(MEPol)

sort!(MEPol, :PDiff)
J.plot_diff(MEPol, dim = :ECC)
# 153×8 DataFrame
#  Row │ ECC                   Poll    Area    Year   Spruce      Tanoak         Diff            PDiff   
#      │ String                String  String  Int64  Float64     Float32        Float64         Float64
# ─────┼─────────────────────────────────────────────────────────────────────────────────────────────────
#    1 │ Passenger             VOC     QC       2022     0.0         0.00168759    -0.00168759   -100.0
#  152 │ NGDistribution        VOC     NT       2022     3.26778     0.369255       2.89852       784.97
#  153 │ Petroleum             VOC     NS       2022  1246.41        0.0         1246.41          Inf

MEReduce = J.diff_fast("MEReduce", loc1, loc2; dimension_filters, sec) # Fine

# MEPol[ecc,poll,area] = (MEPOCA[ecc,poll,area]*MEDriver[ecc,area]*MEPolSwitch[ecc,poll,area])+
#   (xMEPol[ecc,poll,area]*(1-MEPolSwitch[ecc,poll,area]))
MEPOCA = J.diff_fast("MEPOCA", loc1, loc2; dimension_filters, sec) # missing in Tanoak
MEDriver = J.diff_fast("MEDriver", loc1, loc2; dimension_filters, sec) # missing in Tanoak
Driver = J.diff_fast("Driver", loc1, loc2; dimension_filters, sec) # missing in Tanoak
xDriver = J.diff_fast("xDriver", loc1, loc2; dimension_filters, sec) # matches
DriverMultiplier = J.diff_fast("DriverMultiplier", loc1, loc2; dimension_filters, sec) # matches
GRPAdj = J.diff_fast("GRPAdj", loc1, loc2; dimension_filters, sec) # matches
MEPolSwitch = J.diff_fast("MEPolSwitch", loc1, loc2; dimension_filters, sec) # matches, 1
xMEPol = J.diff_fast("xMEPol", loc1, loc2; dimension_filters, sec) # matches

# MEPOCA[ecc,poll,area] = MEPOCX[ecc,poll,area]*MERM[ecc,poll,area]*MEPOCM[ecc,poll,area]
MEPOCX = J.diff_fast("MEPOCX", loc1, loc2; dimension_filters, sec = 'M') # missing in Tanoak
yr = pop!(dimension_filters, :Poll)
MEPOCX = J.diff_fast("MEPOCX", loc1, loc2; dimension_filters, sec = 'M') # missing in Tanoak
@rsubset! MEPOCX abs(:Diff) > 1e-1
# issue starts in 2019, only in NS, Petroleum, VOC

i = findall(vars.Variable .== "UnPolGross" .&& first.(vars.Database) .== 'E')
i = i[1]
vars[i,:Database] = "EGOutput3"
UnPolGross = J.diff_fast("UnPolGross", loc1, loc2; dimension_filters, sec) # 
@rsubset! UnPolGross :Diff != 0
UnArea = J.var("UnArea", loc1)
leftjoin!(UnPolGross,UnArea,on = :Unit)
@rsubset! UnPolGross :UnArea == "YT"
push!(dimension_filters, :Unit => UnPolGross.Unit)

# UnPolGross[unit,fuelep,poll] = UnDmd[unit,fuelep]*UnPOCA[unit,fuelep,poll]
UnDmd = J.diff_fast("UnDmd", loc1, loc2; dimension_filters, sec) # Tanoak is zero
UnPOCA = J.diff_fast("UnPOCA", loc1, loc2; dimension_filters, sec) # Matches

# UnDmd[unit,fuelep] = max(UnEGGross[unit]*UnHRt[unit]/1e6*UnFlFr[unit,fuelep],0)
UnEGGross = J.diff_fast("UnEGGross", loc1, loc2; dimension_filters, sec) # Tanoak is zero
UnHRt = J.diff_fast("UnHRt", loc1, loc2; dimension_filters, sec) # matches
UnFlFr = J.diff_fast("UnFlFr", loc1, loc2; dimension_filters, sec) # matches
UnEGA = J.diff_fast("UnEGA", loc1, loc2; dimension_filters, sec) # Tnaoak is zero
UnGC = J.diff_fast("UnGC", loc1, loc2; dimension_filters, sec) # matches
UnOnline = J.diff_fast("UnOnline", loc1, loc2; dimension_filters, sec) # matches
UnRetire = J.diff_fast("UnRetire", loc1, loc2; dimension_filters, sec) # Tanoak retires in 2019

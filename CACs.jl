using Revise
import JuliaCompare as J
import PromulaDBA as P
import SmallModel as M 
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

################################################################################
# Inputs for file locations ####################################################
################################################################################
DATA_FOLDER1 = raw"\\Pink\c\2020CanadaSpruce\2020Model"
DATA_FOLDER2 = raw"\\Pink\c\2020CanadaTanoak\2020Model"
################################################################################

vars = J.list_vars(DATA_FOLDER1, db_files);
vars_j = J.list_vars(joinpath(DATA_FOLDER2,"database.hdf5"))
loc1 = J.Loc_p(vars, DATA_FOLDER1, "Spruce");
loc2 = J.Loc_j(vars_j, joinpath(DATA_FOLDER2, "database.hdf5"), "Tanoak");
tan_base = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaTanoak\\2020Model\\Base\\database.hdf5", "Tan_Base");
spr_base = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaSpruce\\2020Model\\Base", "Spr_Base");
tan_ogref = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaTanoak\\2020Model\\OGRef\\database.hdf5", "Tan_OGRef");
spr_ogref = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaSpruce\\2020Model\\OGRef", "Spr_OGRef");

################################################################################
# Initialize filters and sec for analysis ######################################
################################################################################
sec = 'T'
dimension_filters = Dict{Symbol,Any}()
push!(dimension_filters, :Area => Canada, :Year => string.(1986:2050))
offroad = copy(dimension_filters)
push!(offroad, :ECC => ["CommercialOffRoad", "ResidentialOffRoad"], :EC => ["CommercialOffRoad", "ResidentialOffRoad"])
push!(offroad, :Year => string.(2040:2043))
push!(offroad, :FuelEP => "Gasoline", :Tech => "OffRoad")
push!(offroad, :Poll => ["COX","VOC"])
push!(offroad, :Area => "ON") # on of many
nox = copy(dimension_filters)
push!(nox, :Poll => "NOX", :Year => string.(2020:2026))
push!(nox, :ECC => ["Freight"], :EC => ["Freight"])
push!(nox, :Area => ["BC"])
push!(nox, :FuelEP => "Diesel", :Tech => "TrainDiesel")
ugen = copy(dimension_filters)
push!(ugen, :ECC => "UtilityGen")
push!(ugen, :Year => string.(2022:2023))
push!(ugen, :Poll => ["SOX", "NOX"])
push!(ugen, :FuelEP => ["Coal", "HFO"])
push!(ugen, :Area => "NS")
push!(ugen, :Market => ["Market16", "Market96"])
################################################################################
# Analysis of variables: Overview ########################################################
################################################################################
EuFPol = J.diff_fast("EuFPol", loc1, loc2; dimension_filters, sec) # 
@rsubset! EuFPol :Diff != 0 :ECC != "ForeignPassenger" :ECC != "ForeignFreight" :Poll != "CO2"
J.plot_diff(EuFPol, dim="ECC")

################################################################################
# Analysis of variables: UtilityGen ########################################################
################################################################################

EuFPol_u = J.diff_fast("EuFPol", loc1, loc2; dimension_filters=ugen, sec) # 
@rsubset! EuFPol_u :Diff != 0 :ECC != "ForeignPassenger" :ECC != "ForeignFreight" :Poll != "CO2"
J.plot_diff(EuFPol_u, dim="Area")

i = findall(vars.Variable .== "UnPolGross" .&& first.(vars.Database) .== 'E')
i = i[1]
vars[i,:Database] = "EGOutput3"
i = findall(vars.Variable .== "EUVR" .&& first.(vars.Database) .== 'E')
i = i[1]
vars[i,:Database] = "EGOutput3"
i = findall(vars.Variable .== "EURP" .&& first.(vars.Database) .== 'E')
i = i[1]
vars[i,:Database] = "EGOutput3"

UnPolGross = J.diff_fast("UnPolGross", loc1, loc2; dimension_filters=ugen, sec='E') # 
@rsubset! UnPolGross :Diff != 0
UnArea = J.var("UnArea", loc1)
leftjoin!(UnPolGross,UnArea,on = :Unit)
@rsubset! UnPolGross :UnArea == "NS"
push!(ugen, :Unit => unique(UnPolGross.Unit))

# UnPolGross[unit,fuelep,poll] = UnDmd[unit,fuelep]*UnPOCA[unit,fuelep,poll]
UnDmd = J.diff_fast("UnDmd", loc1, loc2; dimension_filters=ugen, sec='E') # Matches
UnPOCA = J.diff_fast("UnPOCA", loc1, loc2; dimension_filters=ugen, sec='E') # zeroed out for Spruce
UnPOCX = J.diff_fast("UnPOCX", loc1, loc2; dimension_filters=ugen, sec='E') # Matches
UnRM = J.diff_fast("UnRM", loc1, loc2; dimension_filters=ugen, sec='E') # zeroed out Spruce
#         UnRM[unit,fuelep,poll] = (1-UnRP[unit,fuelep,poll])*
#  (1-EUVR[fuelep,plant,poll,area]) 
UnRP = J.diff_fast("UnRP", loc1, loc2; dimension_filters=ugen, sec='E') # switches to 1 in Spruce
EUVR = J.diff_fast("EUVR", loc1, loc2; dimension_filters=ugen, sec='E') # matches
# UnRP[unit,fuelep,poll] = max(EURP[fuelep,plant,poll,area],xUnRP[unit,fuelep,poll])
xUnRP = J.diff_fast("xUnRP", loc1, loc2; dimension_filters=ugen, sec='E') # matches
EURP = J.diff_fast("EURP", loc1, loc2; dimension_filters=ugen, sec='E') # switches to 1 in Spruce
# EURP[fuelep,plant,poll,area] = RPolicy[ecc,poll,area]*ECoverage[ecc,poll,pcov,area]
RPolicy = J.diff_fast("RPolicy", loc1, loc2; dimension_filters=ugen, sec='S') # switches to 1 in Spruce
ECoverage = J.diff_fast("ECoverage", loc1, loc2; dimension_filters=ugen, sec='S') # switches to 1 in Spruce

xGoalPol = J.diff_fast("xGoalPol", loc1, loc2; dimension_filters=on, sec = sec)
AreaMarket = J.diff_fast("AreaMarket", loc1, loc2; dimension_filters=ugen, sec = sec)
@rsubset! AreaMarket :Spruce == 1
ECCMarket = J.diff_fast("ECCMarket", loc1, loc2; dimension_filters=ugen, sec = sec)
@rsubset! ECCMarket :Spruce == 1
markets = intersect(unique(ECCMarket.Market), unique(AreaMarket.Market))
PollMarket = J.diff_fast("PollMarket", loc1, loc2; dimension_filters=ugen, sec = sec)
@rsubset! PollMarket :Spruce == 1
markets = intersect(markets, unique(PollMarket.Market))
GrossTot = J.diff_fast("GrossTot", loc1, loc2; dimension_filters=ugen, sec = sec)
CapTrade = J.diff_fast("CapTrade", loc1, loc2; dimension_filters=ugen, sec = sec)
GoalPol = J.diff_fast("GoalPol", loc1, loc2; dimension_filters=ugen, sec = sec)


# UnDmd[unit,fuelep] = max(UnEGGross[unit]*UnHRt[unit]/1e6*UnFlFr[unit,fuelep],0)
UnEGGross = J.diff_fast("UnEGGross", loc1, loc2; dimension_filters=ugen, sec='E') # Tanoak is zero
UnHRt = J.diff_fast("UnHRt", loc1, loc2; dimension_filters=ugen, sec='E') # matches
UnFlFr = J.diff_fast("UnFlFr", loc1, loc2; dimension_filters=ugen, sec='E') # matches
UnEGA = J.diff_fast("UnEGA", loc1, loc2; dimension_filters=ugen, sec='E') # Tnaoak is zero
UnGC = J.diff_fast("UnGC", loc1, loc2; dimension_filters=ugen, sec='E') # matches
UnOnline = J.diff_fast("UnOnline", loc1, loc2; dimension_filters=ugen, sec='E') # matches
UnRetire = J.diff_fast("UnRetire", loc1, loc2; dimension_filters=ugen, sec='E') # Tanoak retires in 2019

################################################################################
# Analysis of variables: OffRoad ########################################################
################################################################################


EuFPol_o1 = @rsubset EuFPol :ECC ∈ ["CommercialOffRoad", "ResidentialOffRoad"]
@rsubset! EuFPol_o1 :Year ∈ (2040:2043)

EuFPol_o = J.diff_fast("EuFPol", loc1, loc2; dimension_filters=offroad, sec) # 
@rsubset! EuFPol_o :Diff != 0 :Poll != "CO2"
J.plot_diff(EuFPol_o, dim="Poll"; title= "Offroad")

Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters=offroad, sec = sec) # this is fine
J.add_pdiff!(Polute)
@rsubset! Polute abs(:PDiff) != 0
J.plot_diff(Polute, dim = "Area")

# Everything is off by about 10 percent starting in 2023. MarineLight, Diesel is the biggest offender

i = findall(vars.Variable .== "POCA" .&& first.(vars.Database) .== 'T')
i = i[1]
vars[i,:Database] = "TOutput2"

POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters=offroad, sec) # Issue
J.add_pdiff!(POCA)
@rsubset! POCA abs(:PDiff) != 0

on = copy(dimension_filters)
push!(on, :Enduse => ["Heat", "Motors", "OthSub"])
push!(on, :FuelEP => ["HFO", "NaturalGas", "StillGas"])
push!(on, :EC => ["Petroleum"])
push!(on, :Area => ["ON"])

nl = copy(dimension_filters)
push!(nl, :Enduse => ["Heat", "Motors", "OthSub"])
push!(nl, :FuelEP => ["HFO", "LFO", "LPG"])
push!(nl, :EC => ["IronOreMining"])
push!(nl, :Area => ["NL"])

i = findall(vars.Variable .== "PolSw" .&& first.(vars.Database) .== 'T')
i = i[1]
vars[i,:Database] = "TInput2"

PolSw = J.diff_fast("PolSw", loc1, loc2; dimension_filters=offroad, sec) # matches, 2.0 everywhere
PolSw = J.diff_fast("PolSw", loc1, loc2; dimension_filters=nox, sec) # matches, 2.0 everywhere
PolSw = J.diff_fast("PolSw", loc1, loc2; dimension_filters=nl, sec) # matches, 2.0 everywhere

POCX = J.diff_fast("POCX", loc1, loc2; dimension_filters=offroad, sec) # does not match
POCX = J.diff_fast("POCX", loc1, loc2; dimension_filters=nox, sec) # does not match
@rsubset POCX abs(:Diff) > 1e-1
POCX = J.diff_fast("POCX", loc1, loc2; dimension_filters=on, sec) # matches, everywhere
POCX = J.diff_fast("POCX", loc1, loc2; dimension_filters=nl, sec) # matches, everywhere


i = findall(vars.Variable .== "TrMEPX" .&& first.(vars.Database) .== 'T')
i = i[1]
vars[i,:Database] = "TInput2"
TrMEPX = J.diff_fast("TrMEPX", loc1, loc2; dimension_filters=nox, sec) # does not match

RM = J.diff_fast("RM", loc1, loc2; dimension_filters, sec) # issue!
@rsubset RM abs(:Diff) > 1e-1
RM = J.diff_fast("RM", loc1, loc2; dimension_filters=on, sec) # issue!
RM = J.diff_fast("RM", loc1, loc2; dimension_filters=nl, sec) # issue!

# @. RM = (1.0-RP)*(1.0-VR)*xRM
xRM = J.diff_fast("xRM", loc1, loc2; dimension_filters, sec) # match
@rsubset xRM abs(:Diff) > 1e-1
RP = J.diff_fast("RP", loc1, loc2; dimension_filters, sec) # issue
@rsubset RP abs(:Diff) > 1e-1

i = findall(vars.Variable .== "VR" .&& first.(vars.Database) .== 'I')
i = i[1]
vars[i,:Database] = "IOutput2"

VR = J.diff_fast("VR", loc1, loc2; dimension_filters, sec) # match
@rsubset VR abs(:Diff) > 1e-1

          # RP[fuelep,ec,poll,area] = max(RPolicy[ec,poll,area],
          #   (1-minimum(POCS[:,fuelep,ec,poll,area])/
          #   max(maximum(POCX[:,fuelep,ec,poll,area]),0.000001)))*
          #   PCCov[fuelep,ec,poll,area]

RPolicy = J.diff_fast("RPolicy", loc1, loc2; dimension_filters=on, sec = sec)
RPolicy = J.diff_fast("RPolicy", loc1, loc2; dimension_filters=nl, sec = sec) # this is addressed

RPolicy = J.diff_fast("RPolicy", loc1, loc2; dimension_filters=nb, sec = 'S') # this is addressed
RPolicy = J.diff_fast("RPolicy", loc1, loc2; dimension_filters=on, sec = 'S')
@rsubset! RPolicy abs(:Diff) > 0
sort(RPolicy, :Diff)
# RPolECC = J.diff_fast("RPolECC", loc1, loc2; dimension_filters=nl, sec = sec) # this is addressed

push!(on, :Market => "Market180")
push!(on, :Year => string.(2020:2027))
xGoalPol = J.diff_fast("xGoalPol", loc1, loc2; dimension_filters=on, sec = sec)
AreaMarket = J.diff_fast("AreaMarket", loc1, loc2; dimension_filters=nb, sec = sec)
ECCMarket = J.diff_fast("ECCMarket", loc1, loc2; dimension_filters=nb, sec = sec)
@rsubset ECCMarket :Diff != 0
CapTrade = J.diff_fast("CapTrade", loc1, loc2; dimension_filters=nb, sec = sec)
ECoverage = J.diff_fast("ECoverage", loc1, loc2; dimension_filters=on, sec = sec)
@rsubset ECoverage :Diff != 0
PCOVMarket = J.diff_fast("PCOVMarket", loc1, loc2; dimension_filters=on, sec = sec)
@rsubset PCOVMarket :Diff != 0
PollMarket = J.diff_fast("PollMarket", loc1, loc2; dimension_filters=nb, sec = sec)
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
PCOVMarket = J.diff_fast("PCOVMarket", loc1, loc2; dimension_filters=nb, sec = sec)
PCCov = J.diff_fast("PCCov", loc1, loc2; dimension_filters=nb, sec = sec)
@rsubset PCCov :Tanoak != 1
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
RM = J.diff_fast("RM", loc1, loc2; dimension_filters=nb, sec) # Issue

# @. RM = (1.0-RP)*(1.0-VR)*xRM
RP = J.diff_fast("RP", loc1, loc2; dimension_filters, sec) # close
@rsubset RP :Tanoak == 1
VR = J.diff_fast("VR", loc1, loc2; dimension_filters, sec) # close
xRM = J.diff_fast("xRM", loc1, loc2; dimension_filters, sec) # close

nb = copy(dimension_filters)
push!(nb,:Area => "NB")
push!(nb,:EC => "OtherNonferrous")
push!(nb,:ECC => "OtherNonferrous")
push!(nb, :Poll => ["SOX","PMT","PM25","PM10","BC"])

RPolicy = J.diff_fast("RPolicy", loc1, loc2; dimension_filters=nb, sec = sec) # this is addressed
RPolicy = J.diff_fast("RPolicy", loc1, loc2; dimension_filters=nb, sec = 'S') # this has NaNs
@rsubset RPolicy isnan(:Diff)

GrossTot = J.diff_fast("GrossTot", loc1, loc2; dimension_filters=nb, sec = sec) # this is addressed
@rsubset GrossTot isnan(:Diff)
push!(nb, :Market => ["Market83","Market100","Market101","Market102","Market103"])

# GrossTot(Mkt)=sum(ECC,Poll,Area)(GrossPolPrior(ECC,Poll,Area)*PolConv(Poll))
push!(nb, :Year => string.(2022:2024))
GrossPol = J.diff_fast("GrossPol", loc1, loc2; dimension_filters, sec = sec) # this is addressed
@rsubset GrossPol isnan(:Diff)

# @finite_math GrossPol[ecc,poll,area] = GrossPol[ecc,poll,area]+
#                                      MEPol[ecc,poll,area]/MERM[ecc,poll,area]*
#                                      ECoverage[ecc,poll,pcov,area]
MEPol = J.diff_fast("MEPol", loc1, loc2; dimension_filters=nb, sec = sec) # this is fine
@rsubset MEPol isnan(:Diff)
MERM = J.diff_fast("MERM", loc1, loc2; dimension_filters=nb, sec = sec) # this is fine
@rsubset MERM isnan(:Diff)

    
# for ec in ECs
#   ecc = Select(ECC,EC[ec])
#   for area in Areas, poll in Polls
#     for fuelep in FuelEPs
#       PollTemp[fuelep,ec,poll,area] = sum(Polute[enduse,fuelep,ec,poll,area] for enduse in Enduses)
#     end
#     GrossPol[ecc,poll,area] = sum((PollTemp[fuelep,ec,poll,area]+
#       CgFPol[fuelep,ecc,poll,area])/RM[fuelep,ec,poll,area]*PCCov[fuelep,ec,poll,area]*
#       (1-ZeroFr[fuelep,poll,area]) for fuelep in FuelEPs)+
#       sum(FsPol[fuel,ec,poll,area]*ECoverage[ec,poll,noncombustion,area] for fuel in Fuels)
#   end
# end
Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters=nb, sec = sec) # this is fine
@rsubset Polute isnan(:Diff)
CgFPol = J.diff_fast("CgFPol", loc1, loc2; dimension_filters=nb, sec = sec) # this is fine
@rsubset CgFPol isnan(:Diff)
RM = J.diff_fast("RM", loc1, loc2; dimension_filters=nb, sec = sec) # this is fine
@rsubset RM isnan(:Diff)
PCCov = J.diff_fast("PCCov", loc1, loc2; dimension_filters=nb, sec = sec) # this is fine
@rsubset PCCov isnan(:Diff)
ZeroFr = J.diff_fast("ZeroFr", loc1, loc2; dimension_filters=nb, sec = sec) # this is fine
@rsubset ZeroFr isnan(:Diff)
FsPol = J.diff_fast("FsPol", loc1, loc2; dimension_filters=nb, sec = sec) # this is fine
@rsubset FsPol isnan(:Diff)
ECoverage = J.diff_fast("ECoverage", loc1, loc2; dimension_filters=nb, sec = sec) # this is fine
@rsubset ECoverage isnan(:Diff)

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

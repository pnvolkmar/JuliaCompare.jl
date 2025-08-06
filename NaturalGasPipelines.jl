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
sec = 'C'
dimension_filters = Dict{Symbol,Any}()
push!(dimension_filters, :Area => Canada, :Year => string.(1986:2050))

ngp = copy(dimension_filters)
push!(ngp, :ECC => "NGPipeline", :EC => "NGPipeline")
push!(ngp, :Poll => "CO2", :Year => string.(2000:2050))
push!(ngp, :Area => ["AB","SK", "BC", "ON"])
push!(ngp, :Fuel => "NaturalGas", :FuelEP => "NaturalGas")
push!(ngp, :Enduse => "OthSub")
push!(ngp, :Tech => ["Gas","Electric"])
push!(ngp, :ES => ["Commercial"])

################################################################################
# Analysis of variables: Overview ########################################################
################################################################################

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters=ngp, sec) # 
@rsubset! TotPol abs(:Diff) >= .01
J.add_pdiff!(TotPol)
J.plot_diff(TotPol, dim="Poll")
J.plot_diff(TotPol, dim="Area"; title = "Differences in NGPipeline TotPol: CO2")

    # TotPol[ecc,poll,area] = EnPol[ecc,poll,area]+NcPol[ecc,poll,area]+
    #   MEPol[ecc,poll,area]+VnPol[ecc,poll,area]+FlPol[ecc,poll,area]+
    #   FuPol[ecc,poll,area]+ORMEPol[ecc,poll,area]+SqPol[ecc,poll,area]
EnPol = J.diff_fast("EnPol", loc1, loc2; dimension_filters=ngp, sec) # Problem
NcPol = J.diff_fast("NcPol", loc1, loc2; dimension_filters=ngp, sec) # 
MEPol = J.diff_fast("MEPol", loc1, loc2; dimension_filters=ngp, sec) # 
VnPol = J.diff_fast("VnPol", loc1, loc2; dimension_filters=ngp, sec) # 
FlPol = J.diff_fast("FlPol", loc1, loc2; dimension_filters=ngp, sec) # 
FuPol = J.diff_fast("FuPol", loc1, loc2; dimension_filters=ngp, sec) # 
ORMEPol = J.diff_fast("ORMEPol", loc1, loc2; dimension_filters=ngp, sec) # 
SqPol = J.diff_fast("SqPol", loc1, loc2; dimension_filters=ngp, sec) # 

J.plot_diff(EnPol, dim = "Area")

EuFPol = J.diff_fast("EuFPol", loc1, loc2; dimension_filters=ngp, sec) # Problem
@rsubset! EuFPol abs(:Diff) >= 1
J.plot_diff(EuFPol, dim = "FuelEP"; title = "EuFPol NGPipelines CO2, BC,ON,QC, by FuelEP")
pop!(ngp, :EC)
Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters=ngp, sec = 'C') # this is fine
EC_I = unique(Polute.EC)
J.add_pdiff!(Polute)
@rsubset! Polute abs(:PDiff) != 0
J.plot_diff(Polute, dim = "Area")
unique(Polute.Enduse)

# Polute[enduse,fuelep,ec,poll,area] = EuDem[enduse,fuelep,ec,area]*
#     POCA[enduse,fuelep,ec,poll,area]+xPolute[enduse,fuelep,ec,poll,area]

i = findall(vars.Variable .== "POCA" .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters=ngp, sec) # 
@rsubset POCA :Diff != 0
EuDem = J.diff_fast("EuDem", loc1, loc2; dimension_filters=ngp, sec) # Problem
@rsubset! EuDem :Diff != 0
J.add_pdiff!(EuDem)
@rsubset EuDem abs(:PDiff) >0.1
sort(EuDem, :PDiff)
xPolute = J.diff_fast("xPolute", loc1, loc2; dimension_filters=ngp, sec) # 
J.add_pdiff!(EuDem)
@rsubset! EuDem abs(:PDiff) != 0

# EuDemF[enduse,fuel,ecc,area] = sum(Dmd[enduse,tech,ec,area]*
# DmFrac[enduse,fuel,tech,ec,area] for tech in Techs)
Dmd = J.diff_fast("Dmd", loc1, loc2; dimension_filters=ngp, sec) # 
J.add_pdiff!(Dmd)
sort(Dmd, :PDiff)
@rsubset! Dmd :PDiff != 0
unique(Dmd.Tech)
@rsubset! Dmd :Diff != 0 :Tech == "OffRoad"
J.add_pdiff!(Dmd)
DmFrac = J.diff_fast("DmFrac", loc1, loc2; dimension_filters=ngp, sec) # 
@rsubset! DmFrac :Diff != 0
J.add_pdiff!(DmFrac)

# @. Dmd = Dmd-EE
EE = J.diff_fast("EE", loc1, loc2; dimension_filters=ngp, sec) # 
@rsubset EE :Diff != 0
# @. Dmd = Dmd-DSMEU
DSMEU = J.diff_fast("DSMEU", loc1, loc2; dimension_filters=ngp, sec) # 
@rsubset! EE :Diff != 0
# Dmd[enduse,tech,ec,area] = Dmd[enduse,tech,ec,area]+
#                            SqDmd[tech,ec,area]*SqEUTechMap[enduse,tech]
SqDmd = J.diff_fast("SqDmd", loc1, loc2; dimension_filters=ngp, sec) # 
@rsubset SqDmd :Diff != 0
SqEUTechMap = J.diff_fast("SqEUTechMap", loc1, loc2; dimension_filters=ngp, sec) # 
@rsubset SqEUTechMap :Diff != 0
#
# @finite_math Dmd[enduse,tech,ec,area] = DER[enduse,tech,ec,area]*
# UMS[enduse,tech,ec,area]*CERSM[enduse,ec,area]*
# CUF[enduse,tech,ec,area]*WCUF[ec,area]*RPEI[enduse,tech,ec,area]/1.0e6*
# (TSLoad[enduse,ec,area]*
# (DDay[enduse,area]/DDayNorm[enduse,area])^DDCoefficient[enduse,ec,area]+
# (1.0-TSLoad[enduse,ec,area]))
DER = J.diff_fast("DER", loc1, loc2; dimension_filters=ngp, sec) # Problem
J.add_pdiff!(DER)
UMS = J.diff_fast("UMS", loc1, loc2; dimension_filters=ngp, sec)
CERSM = J.diff_fast("CERSM", loc1, loc2; dimension_filters=ngp, sec)
CUF = J.diff_fast("CUF", loc1, loc2; dimension_filters=ngp, sec)
WCUF = J.diff_fast("WCUF", loc1, loc2; dimension_filters=ngp, sec)
RPEI = J.diff_fast("RPEI", loc1, loc2; dimension_filters=ngp, sec)
TSLoad = J.diff_fast("TSLoad", loc1, loc2; dimension_filters=ngp, sec)
DDay = J.diff_fast("DDay", loc1, loc2; dimension_filters=ngp, sec)
DDayNorm = J.diff_fast("DDayNorm", loc1, loc2; dimension_filters=ngp, sec)
DDCoefficient = J.diff_fast("DDCoefficient", loc1, loc2; dimension_filters=ngp, sec)
EEImpact = J.diff_fast("EEImpact", loc1, loc2; dimension_filters=ngp, sec)
EESat = J.diff_fast("EESat", loc1, loc2; dimension_filters=ngp, sec)
xDmd = J.diff_fast("xDmd", loc1, loc2; dimension_filters=ngp, sec)

DERRRExo = J.diff_fast("DERRRExo", loc1, loc2; dimension_filters=ngp, sec)
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters=ngp, sec)
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters=ngp, sec)
i = findall(vars.Variable .∈ Ref(["DERV", "DERAV"]) .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
DERV = J.diff_fast("DERV", loc1, loc2; dimension_filters=ngp, sec) # Problem
select!(DERV, Not([:Enduse, :Tech]))
J.add_pdiff!(DERV)
@rsubset DERV :Vintage == "Vintage 1"
StockAdjustment = J.diff_fast("StockAdjustment", loc1, loc2; dimension_filters=ngp, sec) # fine
DERAV = J.diff_fast("DERAV", loc1, loc2; dimension_filters=ngp, sec)
DERA   = J.diff_fast("DERA", loc1, loc2; dimension_filters=ngp, sec) # problem
J.add_pdiff!(DERA)

PEEBeforeStd = J.diff_fast("PEEBeforeStd", loc1, loc2; dimension_filters=ngp, sec) # doesn't match
J.add_pdiff!(PEEBeforeStd)

PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters=ngp, sec) # doesn't match
J.add_pdiff!(PEE)
@rsubset PEE :PDiff != 0
PEESw = J.diff_fast("PEESw", loc1, loc2; dimension_filters=ngp, sec) # 1

PEEPrice = J.diff_fast("PEEPrice", loc1, loc2; dimension_filters=ngp, sec) # issue - fixed
J.add_pdiff!(PEEPrice)
@rsubset PEEPrice :PDiff != 0
PEM = J.diff_fast("PEM", loc1, loc2; dimension_filters=ngp, sec) # matches
J.add_pdiff!(PEM)
PEMM = J.diff_fast("PEMM", loc1, loc2; dimension_filters=ngp, sec) # matches
# PEEBeforeStd matches PEM = PEM*PEMM which means in Spruce, Electric and Oil are getting the else values
PEPM = J.diff_fast("PEPM", loc1, loc2; dimension_filters=ngp, sec) # 1
PFPN = J.diff_fast("PFPN", loc1, loc2; dimension_filters=ngp, sec) # matches
PFTC = J.diff_fast("PFTC", loc1, loc2; dimension_filters=ngp, sec) # matches
PCCR = J.diff_fast("PCCR", loc1, loc2; dimension_filters=ngp, sec) # matches
PCCRB = J.diff_fast("PCCR", spr_base, tan_base; dimension_filters=ngp, sec) # matches
PCTC = J.diff_fast("PCTC", loc1, loc2; dimension_filters=ngp, sec) # matches


H2PipelineMultiplier = J.diff_fast("H2PipelineMultiplier", loc1, loc2; dimension_filters=ngp, sec) # matches
@rsubset! H2PipelineMultiplier :Diff != 0

ngp_initial = copy(ngp)
push!(ngp_initial, :Year => "1985")

# PEM[enduse,ec,area] = maximum(PEEA[enduse,tech,ec,area,InitialYear]*
#   PEMX[enduse,tech,ec,area] for tech in Techs)

PEEA = J.diff_fast("PEEA", loc1, loc2; dimension_filters=ngp_initial, sec) # issue - fixed
J.add_pdiff!(PEEA)
PEMX = J.diff_fast("PEMX", loc1, loc2; dimension_filters=ngp, sec) # 

# @finite_math PEEA[enduse,tech,ec,area,InitialYear] = 
#   FPCI[enduse,tech,ec,area]*DSt[enduse,ec,area,InitialYear]/
#   PER[enduse,tech,ec,area,InitialYear]

FPCI = J.diff_fast("FPCI", loc1, loc2; dimension_filters=ngp_initial, sec) # unknown; temp var
Dst = J.diff_fast("Dst", loc1, loc2; dimension_filters=ngp_initial, sec) # matches
PER = J.diff_fast("PER", loc1, loc2; dimension_filters=ngp_initial, sec) # an issue
J.add_pdiff(PER)

# @finite_math FPCI[enduse,tech,ec,area] = 
#   sum(PCLV[age,ecc,area,InitialYear] for age in Ages)*
#   FPC[enduse,tech,ec,area]/TFPC[enduse,ec,area]
PCLV = J.diff_fast("PCLV", loc1, loc2; dimension_filters=ngp_initial, sec) # fine when summed
J.add_pdiff(PCLV)
sum(PCLV.Diff)
FPC = J.diff_fast("FPC", loc1, loc2; dimension_filters=ngp_initial, sec) # temp var
TFPC = J.diff_fast("TFPC", loc1, loc2; dimension_filters=ngp_initial, sec) # temp var, total of FPC

# @finite_math FPC[enduse,tech,ec,area] = 
#   PER[enduse,tech,ec,area,InitialYear]/DSt[enduse,ec,area,InitialYear]*
#   PDif[enduse,tech,ec,area]
PER = J.diff_fast("PER", loc1, loc2; dimension_filters=ngp_initial, sec) # issue - fixed
J.add_pdiff!(PER)
PDif = J.diff_fast("PDif", loc1, loc2; dimension_filters=ngp_initial, sec) # Matches

# PER[enduse,tech,ec,area,InitialYear] = 
#   DER[enduse,tech,ec,area,InitialYear]*DEEA[enduse,tech,ec,area,InitialYear]
DER = J.diff_fast("DER", loc1, loc2; dimension_filters=ngp_initial, sec) # fixed
J.add_pdiff!(DER)
DEEA = J.diff_fast("DEEA", loc1, loc2; dimension_filters=ngp_initial, sec) # issue - fixed
J.add_pdiff(DEEA)

DEE = J.diff_fast("DEE", loc1, loc2; dimension_filters=ngp_initial, sec) # issue - fixed

# @finite_math DEE[enduse,tech,ec,area,InitialYear] = 
#   DEM[enduse,tech,ec,area]*DEMM[enduse,tech,ec,area,InitialYear]*
#   (1/(1+(ECFP[enduse,tech,ec,area,InitialYear]/
#   Inflation[area,InitialYear]*DEPM[enduse,tech,ec,area,InitialYear]/
#   DFPN[enduse,tech,ec,area])^DFTC[enduse,tech,ec,area,InitialYear]))

# DEE[enduse,tech,ec,area,InitialYear] = 
#   max(DEE[enduse,tech,ec,area,InitialYear],
#     DEStd[enduse,tech,ec,area,InitialYear],
#     DEStdP[enduse,tech,ec,area,InitialYear])

DEStd = J.diff_fast("DEStd", loc1, loc2; dimension_filters=ngp_initial, sec) # match
DEStdP = J.diff_fast("DEStdP", loc1, loc2; dimension_filters=ngp_initial, sec) # match

DEM = J.diff_fast("DEM", loc1, loc2; dimension_filters=ngp_initial, sec) # match
DEMM = J.diff_fast("DEMM", loc1, loc2; dimension_filters=ngp_initial, sec) # missing in Spruce? it's zero
ECFP = J.diff_fast("ECFP", loc1, loc2; dimension_filters=ngp_initial, sec) # doesn't match - fixed
J.add_pdiff!(ECFP)
Inflation = J.diff_fast("Inflation", loc1, loc2; dimension_filters=ngp_initial, sec) # match
DEPM = J.diff_fast("DEPM", loc1, loc2; dimension_filters=ngp_initial, sec) # match
DFPN = J.diff_fast("DFPN", loc1, loc2; dimension_filters=ngp_initial, sec) # doesn't match - fixed
J.add_pdiff!(DFPN)
DFTC = J.diff_fast("DFTC", loc1, loc2; dimension_filters=ngp_initial, sec) # doesn't match - fixed
J.add_pdiff!(DFTC)

# @finite_math ECFP[enduse,tech,ec,area,curtime] = 
#   sum((FPEC[fuel,ec,area,curtime]+FPCFSNet[fuel,ec,area,curtime])*
#       DmdFuelTechPrior[enduse,fuel,tech,ec,area] for fuel in Fuels)/
#   sum(DmdFuelTechPrior[enduse,fuel,tech,ec,area] for fuel in Fuels)+
#   PCostTech[tech,ec,area,curtime]
# Note that curtime in this instance would be InitialYear
FPEC = J.diff_fast("FPEC", loc1, loc2; dimension_filters=ngp_initial, sec) # v off
FPCFSNet = J.diff_fast("FPCFSNet", loc1, loc2; dimension_filters=ngp_initial, sec) # zeros
DmdFuelTechPrior = J.diff_fast("DmdFuelTechPrior", loc1, loc2; dimension_filters=ngp_initial, sec) # temp var
PCostTech = J.diff_fast("PCostTech", loc1, loc2; dimension_filters=ngp_initial, sec) # zeros

# FPEC[fuel,ec,area,curtime] = FPF[fuel,es,area,curtime]
FPF = J.diff_fast("FPF", loc1, loc2; dimension_filters=ngp_initial, sec) # v off
push!(ngp_initial, :ES => ["Commercial", "Gas"])
pop!(ngp_initial, :ES)
@rsubset FPF :Diff != 0

pop!(ngp_initial, :ES)
ECESMap = J.diff_fast("ECESMap", loc1, loc2; dimension_filters=ngp_initial, sec) # v off

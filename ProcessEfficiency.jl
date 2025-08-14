using Revise
import JuliaCompare as J
import PromulaDBA as P
import SmallModel as M 
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

################################################################################
# Inputs for file locations ####################################################
################################################################################
DATA_FOLDER1 = raw"\\Pink\c\2020CanadaPine\2020Model\Ref25"
DATA_FOLDER2 = raw"\\Pink\c\2020CanadaRedwood\2020Model\Ref25"
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
sec = 'R'
filter = Dict{Symbol,Any}()
push!(filter, :Area => Canada, :Year => string.(1986:2050))
EC = J.var("ECKey", loc1; sec)
pe = copy(filter)

push!(pe, :ECC => EC.EC, :EC => EC.EC)
push!(pe, :Poll => "CO2", :Year => string.(2020:2050))
push!(pe, :Area => ["AB","SK", "BC", "ON"])
push!(pe, :Fuel => "NaturalGas", :FuelEP => "NaturalGas")
push!(pe, :Enduse => "OthSub")
push!(pe, :Tech => ["Gas","Electric"])
push!(pe, :ES => ["Commercial"])

################################################################################
# Analysis of variables: Overview ########################################################
################################################################################

PEE = J.var("PEE", locs; filter=pe, sec, diff = true, pdiff = true) # 
issues = @rsubset PEE abs(:Pine_pdiff_Redwood) > 1
@rsubset issues :Year ∈ 2023:2028 :Area == "NS"
@rsubset PEE :Year ∈ 2023:2028 :Area == "NS" :Enduse == "AC" :EC == "SingleFamilyAttached" :Tech == "Electric"
@rsubset PEE :Year ∈ 2023:2028 :Area == "MB" :Enduse == "AC" :EC == "SingleFamilyAttached" :Tech == "Electric"
exp_grow = @rsubset PEE :Year > 2020 :Area ∈ ["MB", "ON"] :Enduse == "AC"
J.plot_lines(exp_grow, locs; title = "Sum of PEE for residential ECCs in MB and ON's AC Enduse", units = "\$/Btu")
vars[vars.Variable .== "PEE",:]
J.plot_sets(PEE; col = :Pine_pdiff_Redwood, dim = :Area)
@rsubset PEE :Year ∈ 2002:2008 :Area == "YT" :Poll == "VOC"
    # PEE[ecc,poll,area] = EnPol[ecc,poll,area]+NcPol[ecc,poll,area]+
    #   MEPol[ecc,poll,area]+VnPol[ecc,poll,area]+FlPol[ecc,poll,area]+
    #   FuPol[ecc,poll,area]+ORMEPol[ecc,poll,area]+SqPol[ecc,poll,area]
EnPol = J.var("EnPol", locs; filter=pe, sec) # Problem
NcPol = J.var("NcPol", locs; filter=pe, sec) # 
MEPol = J.var("MEPol", locs; filter=pe, sec) # 
VnPol = J.var("VnPol", locs; filter=pe, sec) # 
FlPol = J.var("FlPol", locs; filter=pe, sec) # 
FuPol = J.var("FuPol", locs; filter=pe, sec) # 
ORMEPol = J.var("ORMEPol", locs; filter=pe, sec) # 
SqPol = J.var("SqPol", locs; filter=pe, sec) # 

J.plot_diff(EnPol, dim = "Area")

EuFPol = J.var("EuFPol", locs; filter=pe, sec) # Problem
@rsubset! EuFPol abs(:Diff) >= 1
J.plot_diff(EuFPol, dim = "FuelEP"; title = "EuFPol NGPipelines CO2, BC,ON,QC, by FuelEP")
pop!(pe, :EC)
Polute = J.var("Polute", locs; filter=pe, sec = 'C') # this is fine
EC_I = unique(Polute.EC)
J.add_pdiff!(Polute)
@rsubset! Polute abs(:PDiff) != 0
J.plot_diff(Polute, dim = "Area")
unique(Polute.Enduse)

# Polute[enduse,fuelep,ec,poll,area] = EuDem[enduse,fuelep,ec,area]*
#     POCA[enduse,fuelep,ec,poll,area]+xPolute[enduse,fuelep,ec,poll,area]

i = findall(vars.Variable .== "POCA" .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
POCA = J.var("POCA", locs; filter=pe, sec) # 
@rsubset POCA :Diff != 0
EuDem = J.var("EuDem", locs; filter=pe, sec) # Problem
@rsubset! EuDem :Diff != 0
J.add_pdiff!(EuDem)
@rsubset EuDem abs(:PDiff) >0.1
sort(EuDem, :PDiff)
xPolute = J.var("xPolute", locs; filter=pe, sec) # 
J.add_pdiff!(EuDem)
@rsubset! EuDem abs(:PDiff) != 0

# EuDemF[enduse,fuel,ecc,area] = sum(Dmd[enduse,tech,ec,area]*
# DmFrac[enduse,fuel,tech,ec,area] for tech in Techs)
Dmd = J.var("Dmd", locs; filter=pe, sec) # 
J.add_pdiff!(Dmd)
sort(Dmd, :PDiff)
@rsubset! Dmd :PDiff != 0
unique(Dmd.Tech)
@rsubset! Dmd :Diff != 0 :Tech == "OffRoad"
J.add_pdiff!(Dmd)
DmFrac = J.var("DmFrac", locs; filter=pe, sec) # 
@rsubset! DmFrac :Diff != 0
J.add_pdiff!(DmFrac)

# @. Dmd = Dmd-EE
EE = J.var("EE", locs; filter=pe, sec) # 
@rsubset EE :Diff != 0
# @. Dmd = Dmd-DSMEU
DSMEU = J.var("DSMEU", locs; filter=pe, sec) # 
@rsubset! EE :Diff != 0
# Dmd[enduse,tech,ec,area] = Dmd[enduse,tech,ec,area]+
#                            SqDmd[tech,ec,area]*SqEUTechMap[enduse,tech]
SqDmd = J.var("SqDmd", locs; filter=pe, sec) # 
@rsubset SqDmd :Diff != 0
SqEUTechMap = J.var("SqEUTechMap", locs; filter=pe, sec) # 
@rsubset SqEUTechMap :Diff != 0
#
# @finite_math Dmd[enduse,tech,ec,area] = DER[enduse,tech,ec,area]*
# UMS[enduse,tech,ec,area]*CERSM[enduse,ec,area]*
# CUF[enduse,tech,ec,area]*WCUF[ec,area]*RPEI[enduse,tech,ec,area]/1.0e6*
# (TSLoad[enduse,ec,area]*
# (DDay[enduse,area]/DDayNorm[enduse,area])^DDCoefficient[enduse,ec,area]+
# (1.0-TSLoad[enduse,ec,area]))
DER = J.var("DER", locs; filter=pe, sec) # Problem
J.add_pdiff!(DER)
UMS = J.var("UMS", locs; filter=pe, sec)
CERSM = J.var("CERSM", locs; filter=pe, sec)
CUF = J.var("CUF", locs; filter=pe, sec)
WCUF = J.var("WCUF", locs; filter=pe, sec)
RPEI = J.var("RPEI", locs; filter=pe, sec)
TSLoad = J.var("TSLoad", locs; filter=pe, sec)
DDay = J.var("DDay", locs; filter=pe, sec)
DDayNorm = J.var("DDayNorm", locs; filter=pe, sec)
DDCoefficient = J.var("DDCoefficient", locs; filter=pe, sec)
EEImpact = J.var("EEImpact", locs; filter=pe, sec)
EESat = J.var("EESat", locs; filter=pe, sec)
xDmd = J.var("xDmd", locs; filter=pe, sec)

DERRRExo = J.var("DERRRExo", locs; filter=pe, sec)
PERRRExo = J.var("PERRRExo", locs; filter=pe, sec)
PERRRExo = J.var("PERRRExo", locs; filter=pe, sec)
i = findall(vars.Variable .∈ Ref(["DERV", "DERAV"]) .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
DERV = J.var("DERV", locs; filter=pe, sec) # Problem
select!(DERV, Not([:Enduse, :Tech]))
J.add_pdiff!(DERV)
@rsubset DERV :Vintage == "Vintage 1"
StockAdjustment = J.var("StockAdjustment", locs; filter=pe, sec) # fine
DERAV = J.var("DERAV", locs; filter=pe, sec)
DERA   = J.var("DERA", locs; filter=pe, sec) # problem
J.add_pdiff!(DERA)

PEEBeforeStd = J.var("PEEBeforeStd", locs; filter=pe, sec) # doesn't match
J.add_pdiff!(PEEBeforeStd)

PEE = J.var("PEE", locs; filter=pe, sec) # doesn't match
J.add_pdiff!(PEE)
@rsubset PEE :PDiff != 0
PEESw = J.var("PEESw", locs; filter=pe, sec) # 1

PEEPrice = J.var("PEEPrice", locs; filter=pe, sec) # issue - fixed
J.add_pdiff!(PEEPrice)
@rsubset PEEPrice :PDiff != 0
PEM = J.var("PEM", locs; filter=pe, sec) # matches
J.add_pdiff!(PEM)
PEMM = J.var("PEMM", locs; filter=pe, sec) # matches
# PEEBeforeStd matches PEM = PEM*PEMM which means in Spruce, Electric and Oil are getting the else values
PEPM = J.var("PEPM", locs; filter=pe, sec) # 1
PFPN = J.var("PFPN", locs; filter=pe, sec) # matches
PFTC = J.var("PFTC", locs; filter=pe, sec) # matches
PCCR = J.var("PCCR", locs; filter=pe, sec) # matches
PCCRB = J.var("PCCR", spr_base, tan_base; filter=pe, sec) # matches
PCTC = J.var("PCTC", locs; filter=pe, sec) # matches


H2PipelineMultiplier = J.var("H2PipelineMultiplier", locs; filter=pe, sec) # matches
@rsubset! H2PipelineMultiplier :Diff != 0

ngp_initial = copy(pe)
push!(ngp_initial, :Year => "1985")

# PEM[enduse,ec,area] = maximum(PEEA[enduse,tech,ec,area,InitialYear]*
#   PEMX[enduse,tech,ec,area] for tech in Techs)

PEEA = J.var("PEEA", locs; filter=ngp_initial, sec) # issue - fixed
J.add_pdiff!(PEEA)
PEMX = J.var("PEMX", locs; filter=pe, sec) # 

# @finite_math PEEA[enduse,tech,ec,area,InitialYear] = 
#   FPCI[enduse,tech,ec,area]*DSt[enduse,ec,area,InitialYear]/
#   PER[enduse,tech,ec,area,InitialYear]

FPCI = J.var("FPCI", locs; filter=ngp_initial, sec) # unknown; temp var
Dst = J.var("Dst", locs; filter=ngp_initial, sec) # matches
PER = J.var("PER", locs; filter=ngp_initial, sec) # an issue
J.add_pdiff(PER)

# @finite_math FPCI[enduse,tech,ec,area] = 
#   sum(PCLV[age,ecc,area,InitialYear] for age in Ages)*
#   FPC[enduse,tech,ec,area]/TFPC[enduse,ec,area]
PCLV = J.var("PCLV", locs; filter=ngp_initial, sec) # fine when summed
J.add_pdiff(PCLV)
sum(PCLV.Diff)
FPC = J.var("FPC", locs; filter=ngp_initial, sec) # temp var
TFPC = J.var("TFPC", locs; filter=ngp_initial, sec) # temp var, total of FPC

# @finite_math FPC[enduse,tech,ec,area] = 
#   PER[enduse,tech,ec,area,InitialYear]/DSt[enduse,ec,area,InitialYear]*
#   PDif[enduse,tech,ec,area]
PER = J.var("PER", locs; filter=ngp_initial, sec) # issue - fixed
J.add_pdiff!(PER)
PDif = J.var("PDif", locs; filter=ngp_initial, sec) # Matches

# PER[enduse,tech,ec,area,InitialYear] = 
#   DER[enduse,tech,ec,area,InitialYear]*DEEA[enduse,tech,ec,area,InitialYear]
DER = J.var("DER", locs; filter=ngp_initial, sec) # fixed
J.add_pdiff!(DER)
DEEA = J.var("DEEA", locs; filter=ngp_initial, sec) # issue - fixed
J.add_pdiff(DEEA)

DEE = J.var("DEE", locs; filter=ngp_initial, sec) # issue - fixed

# @finite_math DEE[enduse,tech,ec,area,InitialYear] = 
#   DEM[enduse,tech,ec,area]*DEMM[enduse,tech,ec,area,InitialYear]*
#   (1/(1+(ECFP[enduse,tech,ec,area,InitialYear]/
#   Inflation[area,InitialYear]*DEPM[enduse,tech,ec,area,InitialYear]/
#   DFPN[enduse,tech,ec,area])^DFTC[enduse,tech,ec,area,InitialYear]))

# DEE[enduse,tech,ec,area,InitialYear] = 
#   max(DEE[enduse,tech,ec,area,InitialYear],
#     DEStd[enduse,tech,ec,area,InitialYear],
#     DEStdP[enduse,tech,ec,area,InitialYear])

DEStd = J.var("DEStd", locs; filter=ngp_initial, sec) # match
DEStdP = J.var("DEStdP", locs; filter=ngp_initial, sec) # match

DEM = J.var("DEM", locs; filter=ngp_initial, sec) # match
DEMM = J.var("DEMM", locs; filter=ngp_initial, sec) # missing in Spruce? it's zero
ECFP = J.var("ECFP", locs; filter=ngp_initial, sec) # doesn't match - fixed
J.add_pdiff!(ECFP)
Inflation = J.var("Inflation", locs; filter=ngp_initial, sec) # match
DEPM = J.var("DEPM", locs; filter=ngp_initial, sec) # match
DFPN = J.var("DFPN", locs; filter=ngp_initial, sec) # doesn't match - fixed
J.add_pdiff!(DFPN)
DFTC = J.var("DFTC", locs; filter=ngp_initial, sec) # doesn't match - fixed
J.add_pdiff!(DFTC)

# @finite_math ECFP[enduse,tech,ec,area,curtime] = 
#   sum((FPEC[fuel,ec,area,curtime]+FPCFSNet[fuel,ec,area,curtime])*
#       DmdFuelTechPrior[enduse,fuel,tech,ec,area] for fuel in Fuels)/
#   sum(DmdFuelTechPrior[enduse,fuel,tech,ec,area] for fuel in Fuels)+
#   PCostTech[tech,ec,area,curtime]
# Note that curtime in this instance would be InitialYear
FPEC = J.var("FPEC", locs; filter=ngp_initial, sec) # v off
FPCFSNet = J.var("FPCFSNet", locs; filter=ngp_initial, sec) # zeros
DmdFuelTechPrior = J.var("DmdFuelTechPrior", locs; filter=ngp_initial, sec) # temp var
PCostTech = J.var("PCostTech", locs; filter=ngp_initial, sec) # zeros

# FPEC[fuel,ec,area,curtime] = FPF[fuel,es,area,curtime]
FPF = J.var("FPF", locs; filter=ngp_initial, sec) # v off
push!(ngp_initial, :ES => ["Commercial", "Gas"])
pop!(ngp_initial, :ES)
@rsubset FPF :Diff != 0

pop!(ngp_initial, :ES)
ECESMap = J.var("ECESMap", locs; filter=ngp_initial, sec) # v off

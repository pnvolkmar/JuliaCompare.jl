using Revise
import JuliaCompare as J
import PromulaDBA as P
import SmallModel as M 
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

################################################################################
# Inputs for file locations ####################################################
################################################################################
DATA_FOLDER1 = raw"\\Pink\c\2020CanadaPine\2020Model"
DATA_FOLDER2 = raw"\\Pink\c\2020CanadaRedwood\2020Model"
################################################################################

vars = J.list_vars(DATA_FOLDER1, db_files);
vars_j = J.list_vars(joinpath(DATA_FOLDER2,"database.hdf5"))
loc1 = J.Loc_p(vars, DATA_FOLDER1, "Pine");
loc2 = J.Loc_j(vars_j, joinpath(DATA_FOLDER2, "database.hdf5"), "Redwood");
red_base = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaRedwood\\2020Model\\Base\\database.hdf5", "Redwood");
pn_base = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaPine\\2020Model\\Base", "Pine");
red_ogref = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaRedwood\\2020Model\\OGRef\\database.hdf5", "Redwood");
pn_ogref = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaPine\\2020Model\\OGRef", "Pine");

################################################################################
# Initialize filters and sec for analysis ######################################
################################################################################
sec = 'I'
dimension_filters = Dict{Symbol,Any}()
push!(dimension_filters, :Area => Canada, :Year => string.(1986:2050))
ECC = M.ReadDisk(loc2.HDF5_path, "E2020DB/ECCKey")
push!(dimension_filters, :ECC => ECC[occursin.("Foreign", ECC) .== false])

lng = copy(dimension_filters)
push!(lng, :EC => "LNGProduction",  :ECC => "LNGProduction")
push!(lng, :Poll => "CO2")
push!(lng, :Area => ["BC"])
push!(lng,:Year => string.(2024:2026))
push!(lng, :Fuel => "NaturalGas",  :FuelEP => "NaturalGas")
push!(lng, :Enduse => "Heat")
push!(lng, :Tech => "Gas")

################################################################################
# Analysis of variables: Overview ##############################################
################################################################################

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters, sec) # 
@rsubset! TotPol abs(:Diff) >= 1
J.add_pdiff!(TotPol)
J.plot_diff(TotPol, dim="ECC"; title = "Differences in TotPol")

LNG = @rsubset TotPol :ECC == "LNGProduction"
J.plot_diff(LNG, dim="Area") # BC
J.plot_diff(LNG, dim="Poll") # CO2
@rsubset LNG :Year == 2025

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters=lng, sec) # 
J.add_pdiff!(TotPol)
J.plot_diff(TotPol, dim="Poll")
J.plot_diff(TotPol, dim="Area")

    # TotPol[ecc,poll,area] = EnPol[ecc,poll,area]+NcPol[ecc,poll,area]+
    #   MEPol[ecc,poll,area]+VnPol[ecc,poll,area]+FlPol[ecc,poll,area]+
    #   FuPol[ecc,poll,area]+ORMEPol[ecc,poll,area]+SqPol[ecc,poll,area]
EnPol = J.diff_fast("EnPol", loc1, loc2; dimension_filters=lng, sec) # Problem
NcPol = J.diff_fast("NcPol", loc1, loc2; dimension_filters=lng, sec) # 
MEPol = J.diff_fast("MEPol", loc1, loc2; dimension_filters=lng, sec) # 
VnPol = J.diff_fast("VnPol", loc1, loc2; dimension_filters=lng, sec) # 
FlPol = J.diff_fast("FlPol", loc1, loc2; dimension_filters=lng, sec) # 
FuPol = J.diff_fast("FuPol", loc1, loc2; dimension_filters=lng, sec) # 
ORMEPol = J.diff_fast("ORMEPol", loc1, loc2; dimension_filters=lng, sec) # 
SqPol = J.diff_fast("SqPol", loc1, loc2; dimension_filters=lng, sec) # 

EuFPol = J.diff_fast("EuFPol", loc1, loc2; dimension_filters=lng, sec) # Problem
@rsubset! EuFPol abs(:Diff) >= .1
J.plot_diff(EuFPol, dim = "FuelEP")

Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters=lng, sec = sec) # this is fine
J.add_pdiff!(Polute)
@rsubset! Polute abs(:PDiff) != 0
J.plot_diff(Polute, dim = "Enduse")

i = findall(vars.Variable .== "POCA" .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
# Polute[enduse,fuelep,ec,poll,area] = EuDem[enduse,fuelep,ec,area]*
#     POCA[enduse,fuelep,ec,poll,area]+xPolute[enduse,fuelep,ec,poll,area]

POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters=lng, sec) # 
EuDem = J.diff_fast("EuDem", loc1, loc2; dimension_filters=lng, sec) # Problem
@rsubset! EuDem :Diff != 0
J.add_pdiff!(EuDem)

xPolute = J.diff_fast("xPolute", loc1, loc2; dimension_filters=lng, sec) # fine
J.add_pdiff!(xPolute)
@rsubset! xPolute abs(:PDiff) != 0

# EuDemF[enduse,fuel,ecc,area] = sum(Dmd[enduse,tech,ec,area]*
# DmFrac[enduse,fuel,tech,ec,area] for tech in Techs)
Dmd = J.diff_fast("Dmd", loc1, loc2; dimension_filters=lng, sec) # problem
J.add_pdiff!(Dmd)
@rsubset! Dmd :Diff != 0
DmFrac = J.diff_fast("DmFrac", loc1, loc2; dimension_filters=lng, sec) # fine
@rsubset! DmFrac :Diff != 0
J.add_pdiff!(DmFrac)


# @rsubset! DmFrac :Tanoak != 0
# xDmFrac = J.diff_fast("xDmFrac", loc1, loc2; dimension_filters=zoom, sec) # 
# @rsubset! xDmFrac :Tanoak != 0
# DmFracMin = J.diff_fast("DmFracMin", loc1, loc2; dimension_filters=lng, sec) # 
# @rsubset! DmFracMin :Tanoak != 0
# @rsubset DmFracMin :Area == "QC" :Fuel == "Biodiesel"
# @rsubset DmFracMin :Area == "BC" :Fuel == "Biodiesel"
# DmFracMax = J.diff_fast("DmFracMax", loc1, loc2; dimension_filters=zoom, sec) # 
# @rsubset DmFracMax :Tanoak != 0 :Diff != 0

# @. Dmd = Dmd-EE
EE = J.diff_fast("EE", loc1, loc2; dimension_filters=lng, sec) # 
@rsubset EE :Diff != 0
# @. Dmd = Dmd-DSMEU
DSMEU = J.diff_fast("DSMEU", loc1, loc2; dimension_filters=lng, sec) # not going to worry about now
@rsubset  DSMEU :Diff != 0
# Dmd[enduse,tech,ec,area] = Dmd[enduse,tech,ec,area]+
#                            SqDmd[tech,ec,area]*SqEUTechMap[enduse,tech]
SqDmd = J.diff_fast("SqDmd", loc1, loc2; dimension_filters=lng, sec) #  not going to worry aobut now
@rsubset! SqDmd :Diff != 0
J.add_pdiff!(SqDmd)
SqEUTechMap = J.diff_fast("SqEUTechMap", loc1, loc2; dimension_filters=lng, sec) # 
#
# @finite_math Dmd[enduse,tech,ec,area] = DER[enduse,tech,ec,area]*
# UMS[enduse,tech,ec,area]*CERSM[enduse,ec,area]*
# CUF[enduse,tech,ec,area]*WCUF[ec,area]*RPEI[enduse,tech,ec,area]/1.0e6*
# (TSLoad[enduse,ec,area]*
# (DDay[enduse,area]/DDayNorm[enduse,area])^DDCoefficient[enduse,ec,area]+
# (1.0-TSLoad[enduse,ec,area]))
DER = J.diff_fast("DER", loc1, loc2; dimension_filters=lng, sec) # Problem
J.add_pdiff!(DER)
UMS = J.diff_fast("UMS", loc1, loc2; dimension_filters=lng, sec)
CERSM = J.diff_fast("CERSM", loc1, loc2; dimension_filters=lng, sec)
CUF = J.diff_fast("CUF", loc1, loc2; dimension_filters=lng, sec)
WCUF = J.diff_fast("WCUF", loc1, loc2; dimension_filters=lng, sec)
RPEI = J.diff_fast("RPEI", loc1, loc2; dimension_filters=lng, sec)
TSLoad = J.diff_fast("TSLoad", loc1, loc2; dimension_filters=lng, sec)
DDay = J.diff_fast("DDay", loc1, loc2; dimension_filters=lng, sec)
DDayNorm = J.diff_fast("DDayNorm", loc1, loc2; dimension_filters=lng, sec)
DDCoefficient = J.diff_fast("DDCoefficient", loc1, loc2; dimension_filters=lng, sec)
EEImpact = J.diff_fast("EEImpact", loc1, loc2; dimension_filters=lng, sec)
EESat = J.diff_fast("EESat", loc1, loc2; dimension_filters=lng, sec)
xDmd = J.diff_fast("xDmd", loc1, loc2; dimension_filters=lng, sec)
xDmdTrend = J.diff_fast("xDmdTrend", loc1, loc2; dimension_filters=lng, sec)

DERRRExo = J.diff_fast("DERRRExo", loc1, loc2; dimension_filters=lng, sec) # not going to worry about rn
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters=lng, sec)
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters=lng, sec)
i = findall(vars.Variable .∈ Ref(["DERV", "DERAV"]) .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
DERV = J.diff_fast("DERV", loc1, loc2; dimension_filters=lng, sec) # Problem
select!(DERV, Not([:Enduse, :Tech]))
J.add_pdiff!(DERV)
@rsubset DERV :Vintage == "Vintage 1"
StockAdjustment = J.diff_fast("StockAdjustment", loc1, loc2; dimension_filters=lng, sec) # fine
DERAV = J.diff_fast("DERAV", loc1, loc2; dimension_filters=lng, sec)
DERA   = J.diff_fast("DERA", loc1, loc2; dimension_filters=lng, sec) # problem
J.add_pdiff!(DERA)
@rsubset DERA isnan(:Tanoak) :Year == 2020
#   @. DERA = DERAPC+DERAP+DERAD+DERARC
DERAPC = J.diff_fast("DERAPC", loc1, loc2; dimension_filters=lng, sec) #  problem
J.add_pdiff!(DERAPC)
DERAP = J.diff_fast("DERAP", loc1, loc2; dimension_filters=lng, sec) # fine
DERAD = J.diff_fast("DERAD", loc1, loc2; dimension_filters=lng, sec) # problem
J.add_pdiff!(DERAD)
DERARC = J.diff_fast("DERARC", loc1, loc2; dimension_filters=lng, sec) # fine

RetroSwExo = J.diff_fast("RetroSwExo", loc1, loc2; dimension_filters=lng, sec) # -99

#   @. @finite_math DERAPC = (PERAPC+PERADSt)/DEE
PERAPC = J.diff_fast("PERAPC", loc1, loc2; dimension_filters=lng, sec) # problem
J.add_pdiff!(PERAPC)
PERADSt = J.diff_fast("PERADSt", loc1, loc2; dimension_filters=lng, sec) # fine
@rsubset! PERADSt :Diff != 0
DEE = J.diff_fast("DEE", loc1, loc2; dimension_filters=lng, sec) # fine
J.add_pdiff!(DEE)

# process additions from production capacity (PERAPC)
# @finite_math PERAPC[enduse,tech,ec,area] = EUPCAPC[enduse,tech,New,ec,area]*
#     DSt[enduse,ec,area]/PEE[enduse,tech,ec,area]
push!(lng, :Age => "New")
EUPCAPC = J.diff_fast("EUPCAPC", loc1, loc2; dimension_filters=lng, sec) # big problem
J.add_pdiff!(EUPCAPC)
DSt = J.diff_fast("DSt", loc1, loc2; dimension_filters=lng, sec) # 
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters=lng, sec) # 
J.add_pdiff!(PEE)
# EUPCAPC[enduse,tech,New,ec,area] = PCA[New,ecc,area]*MMSF[enduse,tech,ec,area]
PCA = J.diff_fast("PCA", loc1, loc2; dimension_filters=lng, sec) # 
age = pop!(lng, :Age)
PCA_mid = J.diff_fast("PCA", loc1, loc2; dimension_filters=lng, sec) # 
# PCA[New,ecc,area] = (PCPrior[ecc,area]*(ECGR[ecc,area]+GRPGR[area]))+PCR[Old,ecc,area]
push!(lng, :Age => "Old")
PC = J.diff_fast("PC", loc1, loc2; dimension_filters=lng, sec) # problem
J.add_pdiff!(PC)
ECGR = J.diff_fast("ECGR", loc1, loc2; dimension_filters=lng, sec) # problem
J.add_pdiff!(ECGR)
GRPGR = J.diff_fast("GRPGR", loc1, loc2; dimension_filters=lng, sec) # fine
PCR = J.diff_fast("PCR", loc1, loc2; dimension_filters=lng, sec) # problem
J.add_pdiff!(PCR)
#     @finite_math ECGR[ecc,area] = (PC[ecc,area]-PCPrior[ecc,area])/PCPrior[ecc,area]-GRPGR[area]
# PCMin = J.diff_fast("PCMin", loc1, loc2; dimension_filters=lng, sec) # temp var
xPC = J.diff_fast("xPC", loc1, loc2; dimension_filters=lng, sec) # problem
# @finite_math xPC[ecc,area] = Driver[ecc,area]/xECUF[ecc,area]   
Driver = J.diff_fast("Driver", loc1, loc2; dimension_filters=lng, sec) # 
J.add_pdiff!(Driver)
sort!(Driver, :PDiff)
xECUF = J.diff_fast("xECUF", loc1, loc2; dimension_filters=lng, sec) # 
# Driver[ecc,area] = Driver[ecc,area]*DriverMultiplier[ecc,area]
DriverMultiplier = J.diff_fast("DriverMultiplier", loc1, loc2; dimension_filters=lng, sec) # fine
xDriver = J.diff_fast("xDriver", loc1, loc2; dimension_filters=lng, sec) # issue
DrSwitch = J.diff_fast("DrSwitch", loc1, loc2; dimension_filters=lng, sec) # 21
xGO = J.diff_fast("xGO", loc1, loc2; dimension_filters=lng, sec) # issue
vGO = J.diff_fast("vGO", loc1, loc2; dimension_filters=lng, sec) # doesn't prvide data for CN
    # xGO[ecc,area,year]  =sum(xGOECC[ecc,areatom,year]*
    #   MapAreaTOM[area,areatom] for areatom in areatoms)

# xGOECC = J.diff_fast("xGOECC", loc1, loc2; dimension_filters=lng, sec) # temp var
MapECCfromTOM = J.diff_fast("KInput/MapECCfromTOM", loc1, loc2; dimension_filters=lng, sec='K') # 
@rsubset MapECCfromTOM :Pine != 0 || :Redwood != 0
push!(lng, :ECCfromTOM => ["Fertilizer"])
GY = J.diff_fast("KOutput/GY", pn_base, red_base; dimension_filters=lng, sec='K') # 
@rsubset GY :AreaTOM == "ON"
tech = pop!(lng, :Tech)
MMSF = J.diff_fast("MMSF", loc1, loc2; dimension_filters=lng, sec) # 
J.add_pdiff!(MMSF)
@rsubset! MMSF :Spruce != 0 || :Tanoak != 0
xProcSw = J.diff_fast("xProcSw", loc1, loc2; dimension_filters=lng, sec) # 
@rsubset xProcSw :PI == "MShare"
xMMSF = J.diff_fast("xMMSF", loc1, loc2; dimension_filters=lng, sec) # 

# MMSF[enduse,tech,ec,area] = MMSF[enduse,tech,ec,area]+
# max(MMSF[enduse,tech,ec,area]-MMSFB[enduse,tech,ec,area],0)*
# ETSwitch[tech,area] 
ETSwitch = J.diff_fast("ETSwitch", loc1, loc2; dimension_filters=lng, sec) # 

#   for area in Areas,enduse in Enduses,ec in ECs,tech in Techs
#     if MCFU[enduse,tech,ec,area]  > 0.0 && 
#        MCFU0[enduse,tech,ec,area] > 0.0 &&
#        PEE[enduse,tech,ec,area]   > 0.0 && 
#        PEE0[enduse,tech,ec,area]  > 0.0 &&
#        SPC[ec,area]               > 0.0 && 
#        SPC0[ec,area]              > 0.0 && 
#        SPop[ec,area]              > 0.0 && 
#        SPop[ec,area]              > 0.0 && 
#        MMSM0[enduse,tech,ec,area] > -150.0 

#       MAW[enduse,tech,ec,area] = exp(MMSM0[enduse,tech,ec,area]+
#         log(MSMM[enduse,tech,ec,area])+
#         MMSMI[enduse,tech,ec,area]*(SPC[ec,area]/SPop[ec,area])/
#                                     (SPC0[ec,area]/SPop0[ec,area])+
#         MVF[enduse,tech,ec,area]*
#         log((MCFU[enduse,tech,ec,area]/Inflation[area]/PEE[enduse,tech,ec,area])/
#         (MCFU0[enduse,tech,ec,area]/Inflation0[area]/PEE0[enduse,tech,ec,area])))*
#         MSLimit[enduse,tech,ec,area]
#     else
#       MAW[enduse,tech,ec,area] = 0.0
#     end
#   end

#   for area in Areas,enduse in Enduses,ec in ECs,tech in Techs
#     TMAW[enduse,ec,area] = sum(MAW[enduse,tech,ec,area] for tech in Techs)
#     @finite_math MMSF[enduse,tech,ec,area] = MAW[enduse,tech,ec,area]/TMAW[enduse,ec,area]
#   end
pop!(lng, :Tech)
push!(lng, :Year => string.(2022:2023))
push!(lng, :Tech => ["OffRoad"], :Enduse => ["OffRoad"])
pop!(lng, :Tech)
MMSM0 = J.diff_fast("MMSM0", loc1, loc2; dimension_filters=lng, sec) # doesn't match
MMSM0 = J.diff_fast("MMSM0", spr_base, tan_base; dimension_filters=lng, sec) # doesn't match
@rsubset MMSM0 :Spruce > -170 || :Tanoak > -170
MSMM = J.diff_fast("MSMM", loc1, loc2; dimension_filters=lng, sec) # matches
MVF = J.diff_fast("MVF", loc1, loc2; dimension_filters=lng, sec) # matches
MCFU = J.diff_fast("MCFU", loc1, loc2; dimension_filters=lng, sec) # close
MMSMI = J.diff_fast("MMSMI", loc1, loc2; dimension_filters=lng, sec) # 0
AMSF = J.diff_fast("AMSF", loc1, loc2; dimension_filters=lng, sec) # doesn't match
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters=lng, sec) # doesn't match
@rsubset PEE isnan(:Tanoak)
PEESw = J.diff_fast("PEESw", loc1, loc2; dimension_filters=lng, sec) # 1

            # @finite_math MAW[enduse,tech,ec,area] = exp(MMSMI[enduse,tech,ec,area]*
            #   (SPC[ec,area]/SPop[ec,area])/(SPC0[ec,area]/SPop0[ec,area])+
            #   MVF[enduse,tech,ec,area]*log((MCFU[enduse,tech,ec,area]/Inflation[area]/
            #   PEE[enduse,tech,ec,area])/(MCFU0[enduse,tech,ec,area]/Inflation0[area]/PEE0[enduse,tech,ec,area])))

            MMSMI = J.diff_fast("MMSMI", loc1, loc2; dimension_filters=lng, sec) # 0
MVF = J.diff_fast("MVF", loc1, loc2; dimension_filters=lng, sec) # matches
MCFU = J.diff_fast("MCFU", loc1, loc2; dimension_filters=lng, sec) # close
AMSF = J.diff_fast("AMSF", loc1, loc2; dimension_filters=lng, sec) # doesn't match
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters=lng, sec) # doesn't match



Enduse = pop!(lng, :Enduse)
push!(lng, :Enduse => Enduse)
PEEBeforeStd = J.diff_fast("PEEBeforeStd", loc1, loc2; dimension_filters=lng, sec) # doesn't match
PEEBeforeStd = J.diff_fast("PEEBeforeStd", loc1, loc2; dimension_filters, sec) # doesn't match
J.add_pdiff!(PEEBeforeStd)
@rsubset! PEEBeforeStd :Year >= 2020 :Year <= 2030
sort!(PEEBeforeStd, :PDiff)
lng = copy(dimension_filters)
push!(lng, :Year => "2020")
push!(lng, :Tech => ["Gas", "Oil", "Coal"])
push!(lng, :Enduse => "Steam")
P_issues = @rsubset PEEBeforeStd isnan(:Tanoak) :Year == 2020 :Enduse == "Steam" 
unique(P_issues.Tech)
@rsubset! PEEBeforeStd :Spruce != 0 || :Diff != 0
J.add_pdiff!(PEEBeforeStd)
@rsubset! PEEBeforeStd :PDiff == 100
# ==> implies issue is before std with PEEPrice

#   PEEPrice[enduse,tech,ec,area] = PEM[enduse,ec,area]*
#     PEMM[enduse,tech,ec,area]*
#     (1/(1+(MCFU[enduse,tech,ec,area]/Inflation[area]*
#     PEPM[enduse,tech,ec,area]/PFPN[enduse,tech,ec,area])^PFTC[enduse,tech,ec,area]*
#     (PCCR[enduse,tech,ec,area]/PCCRB[enduse,tech,ec,area])))
PEEPrice = J.diff_fast("PEEPrice", loc1, loc2; dimension_filters=lng, sec) # issue
J.add_pdiff!(PEEPrice)
PEM = J.diff_fast("PEM", loc1, loc2; dimension_filters=lng, sec) # off
J.add_pdiff!(PEM)
PEMM = J.diff_fast("PEMM", loc1, loc2; dimension_filters=lng, sec) # matches
# PEEBeforeStd matches PEM = PEM*PEMM which means in Spruce, Electric and Oil are getting the else values
PEPM = J.diff_fast("PEPM", loc1, loc2; dimension_filters=lng, sec) # 1
PFPN = J.diff_fast("PFPN", loc1, loc2; dimension_filters=lng, sec) # matches
PFTC = J.diff_fast("PFTC", loc1, loc2; dimension_filters=lng, sec) # matches
PCCR = J.diff_fast("PCCR", loc1, loc2; dimension_filters=lng, sec) # matches
PCCRB = J.diff_fast("PCCR", spr_base, tan_base; dimension_filters=lng, sec) # matches
PCTC = J.diff_fast("PCTC", loc1, loc2; dimension_filters=lng, sec) # matches

df = leftjoin(PEEBeforeStd,PCTC, on=[:Enduse,:Tech,:EC,:Area,:Year],  renamecols = "_PEEB" => "_PCTC")
@rsubset df :Spruce_PCTC == 0
PFPN = J.diff_fast("PFPN", loc1, loc2; dimension_filters=lng, sec) # matchesl
df = leftjoin(df, PFPN, on = [:Enduse, :Tech, :EC, :Area], renamecols = "" => "_PFPN")
@rsubset! df :Enduse != "lngRoad"

# Sometime PEE is getting halved and sometimes not. 
# For EC(Food) Enduse(OthSub), Tech(Oil) get halved but Tech(lngRoad) does not.
# Let's findout why
lng = copy(lng)
PEEBeforeStd_f = J.diff_fast("PEEBeforeStd", loc1, loc2; dimension_filters=lng, sec) # doesn't match
PEEBeforeStd # Doesn't match between Oil and lngRoad in Spruce
push!(lng, :EC => "Food", :Enduse => "OthSub", :Tech => ["Oil", "lngRoad"])
PEM = J.diff_fast("PEM", loc1, loc2; dimension_filters=lng, sec) # doesn't have tech
PEMM = J.diff_fast("PEMM", loc1, loc2; dimension_filters=lng, sec) # matches
PCCR = J.diff_fast("PCCR", loc1, loc2; dimension_filters=lng, sec) # matches
PCCRB = J.diff_fast("PCCR", spr_base, tan_base; dimension_filters=lng, sec) # matches
PFTC = J.diff_fast("PFTC", loc1, loc2; dimension_filters=lng, sec) # matches

PEPM = J.diff_fast("PEPM", loc1, loc2; dimension_filters=lng, sec) # 1
PFPN = J.diff_fast("PFPN", loc1, loc2; dimension_filters=lng, sec) # matches

MCFU = J.diff_fast("MCFU", loc1, loc2; dimension_filters=lng, sec) # 1
Inflation = J.diff_fast("Inflation", loc1, loc2; dimension_filters=lng, sec) # matches

InitialDemandYear = J.diff_fast("InitialDemandYear", loc1, loc2; dimension_filters=lng, sec) # matches
lng_initial = copy(lng)
push!(lng_initial, :Year => "2025")

# PEM[enduse,ec,area] = maximum(PEEA[enduse,tech,ec,area,InitialYear]*
#   PEMX[enduse,tech,ec,area] for tech in Techs)
PEEA = J.diff_fast("PEEA", loc1, loc2; dimension_filters=lng_initial, sec) # matches
J.add_pdiff!(PEEA)
pop!(lng_initial, :Tech)
PEMX = J.diff_fast("PEMX", loc1, loc2; dimension_filters=lng_initial, sec) # matches
@rsubset PEMX :Diff != 0

PEEB = J.diff_fast("PEE", pn_base, red_base; dimension_filters=lng, sec)
J.add_pdiff!(PEEB)
PEEAB = J.diff_fast("PEEA", pn_base, red_base; dimension_filters=lng, sec)
J.add_pdiff!(PEEAB)

# @finite_math PEEA[enduse,tech,ec,area] =
#   sum(EUPC[enduse,tech,age,ec,area] for age in Ages)*
#   DSt[enduse,ec,area]/(PER[enduse,tech,ec,area])
EUPC = J.diff_fast("EUPC", pn_base, red_base; dimension_filters=lng, sec)
J.add_pdiff!(EUPC)
DSt = J.diff_fast("DSt", pn_base, red_base; dimension_filters=lng, sec)
J.add_pdiff!(DSt)
PER = J.diff_fast("PER", pn_base, red_base; dimension_filters=lng, sec)
J.add_pdiff!(PER)
PERRRExo = J.diff_fast("PERRRExo", pn_base, red_base; dimension_filters=lng, sec)
J.add_pdiff!(PERRRExo)
RetroSwExo = J.diff_fast("RetroSwExo", pn_base, red_base; dimension_filters=lng, sec) -99
StockAdjustment = J.diff_fast("StockAdjustment", pn_base, red_base; dimension_filters=lng, sec) -99
StockAdjustment_pine = J.var("IInput/StockAdjustment", loc1)
J.subset_dataframe!(StockAdjustment_pine, lng; drop_filtered_cols=false)
StockAdjustment_red = J.var("IInput/StockAdjustment", loc2)
J.subset_dataframe!(StockAdjustment_red, lng; drop_filtered_cols=false)

PERA = J.diff_fast("PERA", pn_base, red_base; dimension_filters=lng, sec)

# @. PERA = PERAPC+PERADSt+PERAP+PERARC
PERAPC = J.diff_fast("PERAPC", pn_base, red_base; dimension_filters=lng, sec) # big problem
PERADSt = J.diff_fast("PERADSt", pn_base, red_base; dimension_filters=lng, sec) # fine
PERAP = J.diff_fast("PERAP", pn_base, red_base; dimension_filters=lng, sec) # problem
PERARC = J.diff_fast("PERARC", pn_base, red_base; dimension_filters=lng, sec) # fine

# New = Select(Age,"New")
# @finite_math PERAPC[enduse,tech,ec,area] = EUPCAPC[enduse,tech,New,ec,area]*
# DSt[enduse,ec,area]/PEE[enduse,tech,ec,area]
push!(lng, :Age => "New")
EUPCAPC = J.diff_fast("EUPCAPC", pn_base, red_base; dimension_filters=lng, sec) # fine
DSt = J.diff_fast("DSt", pn_base, red_base; dimension_filters=lng, sec) # fine
PEE = J.diff_fast("PEE", pn_base, red_base; dimension_filters=lng, sec) # fine
J.add_pdiff!(PEE)

# @. @finite_math PERAP = PERRP*PEEAPrior/PEE

# @finite_math PEEA[enduse,tech,ec,area,InitialYear] = 
# FPCI[enduse,tech,ec,area]*DSt[enduse,ec,area,InitialYear]/
# PER[enduse,tech,ec,area,InitialYear]

xDmd = J.diff_fast("xDmd", loc1, loc2; dimension_filters=lng, sec)

i = findall(vars.Variable .== "DmdFuel" .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Input"

DmdFuel = J.diff_fast("DmdFuel", loc1, loc2; dimension_filters=lng, sec) # IInput2
vDmd = J.diff_fast("vDmd", loc1, loc2; dimension_filters=lng, sec) #
@rsubset vDmd :Diff != 0

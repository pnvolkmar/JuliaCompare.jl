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
sec = 'I'
dimension_filters = Dict{Symbol,Any}()
push!(dimension_filters, :Area => Canada, :Year => string.(1986:2050))
ECC = M.ReadDisk(loc2.HDF5_path, "E2020DB/ECCKey")
push!(dimension_filters, :ECC => ECC[occursin.("Foreign", ECC) .== false])

omm = copy(dimension_filters)
push!(omm, :ECC => "OtherMetalMining", :EC => "OtherMetalMining")
push!(omm, :Poll => "CO2", :Year => string.(2020:2050))
push!(omm, :Area => ["NU"])
push!(omm, :FuelEP => ["Diesel","Gasoline", "Kerosene"], :Fuel => ["Diesel","Gasoline", "Kerosene"])
push!(omm, :Enduse => ["OffRoad", "Heat"])
push!(omm, :Tech => ["OffRoad", "Oil"])

################################################################################
# Analysis of variables: Overview ##############################################
################################################################################

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters, sec) # 
@rsubset! TotPol abs(:Diff) >= 1
J.add_pdiff!(TotPol)
J.plot_diff(TotPol, dim="ECC"; title = "Differences in TotPol")

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters=omm, sec) # 
@rsubset! TotPol abs(:Diff) >= 1
J.add_pdiff!(TotPol)
J.plot_diff(TotPol, dim="Area"; title = "Differences in OtherMetalMining TotPol")

    # TotPol[ecc,poll,area] = EnPol[ecc,poll,area]+NcPol[ecc,poll,area]+
    #   MEPol[ecc,poll,area]+VnPol[ecc,poll,area]+FlPol[ecc,poll,area]+
    #   FuPol[ecc,poll,area]+ORMEPol[ecc,poll,area]+SqPol[ecc,poll,area]
EnPol = J.diff_fast("EnPol", loc1, loc2; dimension_filters=omm, sec) # Problem
NcPol = J.diff_fast("NcPol", loc1, loc2; dimension_filters=omm, sec) # 
MEPol = J.diff_fast("MEPol", loc1, loc2; dimension_filters=omm, sec) # 
VnPol = J.diff_fast("VnPol", loc1, loc2; dimension_filters=omm, sec) # 
FlPol = J.diff_fast("FlPol", loc1, loc2; dimension_filters=omm, sec) # 
FuPol = J.diff_fast("FuPol", loc1, loc2; dimension_filters=omm, sec) # 
ORMEPol = J.diff_fast("ORMEPol", loc1, loc2; dimension_filters=omm, sec) # 
SqPol = J.diff_fast("SqPol", loc1, loc2; dimension_filters=omm, sec) # 

EuFPol = J.diff_fast("EuFPol", loc1, loc2; dimension_filters=omm, sec) # Problem
@rsubset! EuFPol abs(:Diff) >= 1
J.plot_diff(EuFPol, dim = "FuelEP"; title = "EuFPol OtherMetalMining CO2, BC,ON,QC, by FuelEP")

push!(omm, :Area => "ON", :Year => string.(2023:2029), :Enduse => "ommRoad")
push!(omm, :Tech => "ommRoad")
push!(omm, :Enduse => "ommRoad")

Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters=omm, sec = sec) # this is fine
J.add_pdiff!(Polute)
@rsubset! Polute abs(:PDiff) != 0
J.plot_diff(Polute, dim = "Enduse")

i = findall(vars.Variable .== "POCA" .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
# Polute[enduse,fuelep,ec,poll,area] = EuDem[enduse,fuelep,ec,area]*
#     POCA[enduse,fuelep,ec,poll,area]+xPolute[enduse,fuelep,ec,poll,area]

POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters=omm, sec) # 
EuDem = J.diff_fast("EuDem", loc1, loc2; dimension_filters=omm, sec) # Problem
@rsubset! EuDem :Diff != 0
J.add_pdiff!(EuDem)
@rsubset EuDem :Area == "BC" :Year == 2050
xPolute = J.diff_fast("xPolute", loc1, loc2; dimension_filters=omm, sec) # 
J.add_pdiff!(EuDem)
@rsubset! EuDem abs(:PDiff) != 0

EuDemF = J.diff_fast("EuDemF", loc1, loc2; dimension_filters=omm, sec) # not on db

# EuDemF[enduse,fuel,ecc,area] = sum(Dmd[enduse,tech,ec,area]*
# DmFrac[enduse,fuel,tech,ec,area] for tech in Techs)
Dmd = J.diff_fast("Dmd", loc1, loc2; dimension_filters=omm, sec) # problem
@rsubset! Dmd :Diff != 0 :Tech == "ommRoad"
J.add_pdiff!(Dmd)
@rsubset Dmd :PDiff != 0
DmFrac = J.diff_fast("DmFrac", loc1, loc2; dimension_filters=omm, sec) # fine
@rsubset! DmFrac :Diff != 0
J.add_pdiff!(DmFrac)
zoom = copy(omm)
push!(zoom, :Area => "QC", :Year => "2050")
pop!(zoom, :Fuel)
DmFrac = J.diff_fast("DmFrac", loc1, loc2; dimension_filters=zoom, sec) # 
DmFrac = J.diff_fast("DmFrac", spr_base, tan_base; dimension_filters=omm, sec) # 
@rsubset DmFrac :Area == "BC" :Fuel == "Diesel"

@rsubset! DmFrac :Tanoak != 0
xDmFrac = J.diff_fast("xDmFrac", loc1, loc2; dimension_filters=zoom, sec) # 
@rsubset! xDmFrac :Tanoak != 0
DmFracMin = J.diff_fast("DmFracMin", loc1, loc2; dimension_filters=omm, sec) # 
@rsubset! DmFracMin :Tanoak != 0
@rsubset DmFracMin :Area == "QC" :Fuel == "Biodiesel"
@rsubset DmFracMin :Area == "BC" :Fuel == "Biodiesel"
DmFracMax = J.diff_fast("DmFracMax", loc1, loc2; dimension_filters=zoom, sec) # 
@rsubset DmFracMax :Tanoak != 0 :Diff != 0

# @. Dmd = Dmd-EE
EE = J.diff_fast("EE", loc1, loc2; dimension_filters=omm, sec) # 
@rsubset
# @. Dmd = Dmd-DSMEU
DSMEU = J.diff_fast("DSMEU", loc1, loc2; dimension_filters=omm, sec) # 
# Dmd[enduse,tech,ec,area] = Dmd[enduse,tech,ec,area]+
#                            SqDmd[tech,ec,area]*SqEUTechMap[enduse,tech]
SqDmd = J.diff_fast("SqDmd", loc1, loc2; dimension_filters=omm, sec) # 
SqEUTechMap = J.diff_fast("SqEUTechMap", loc1, loc2; dimension_filters=omm, sec) # 
#
# @finite_math Dmd[enduse,tech,ec,area] = DER[enduse,tech,ec,area]*
# UMS[enduse,tech,ec,area]*CERSM[enduse,ec,area]*
# CUF[enduse,tech,ec,area]*WCUF[ec,area]*RPEI[enduse,tech,ec,area]/1.0e6*
# (TSLoad[enduse,ec,area]*
# (DDay[enduse,area]/DDayNorm[enduse,area])^DDCoefficient[enduse,ec,area]+
# (1.0-TSLoad[enduse,ec,area]))
DER = J.diff_fast("DER", loc1, loc2; dimension_filters=omm, sec) # Problem
J.add_pdiff!(DER)
UMS = J.diff_fast("UMS", loc1, loc2; dimension_filters=omm, sec)
CERSM = J.diff_fast("CERSM", loc1, loc2; dimension_filters=omm, sec)
CUF = J.diff_fast("CUF", loc1, loc2; dimension_filters=omm, sec)
WCUF = J.diff_fast("WCUF", loc1, loc2; dimension_filters=omm, sec)
RPEI = J.diff_fast("RPEI", loc1, loc2; dimension_filters=omm, sec)
TSLoad = J.diff_fast("TSLoad", loc1, loc2; dimension_filters=omm, sec)
DDay = J.diff_fast("DDay", loc1, loc2; dimension_filters=omm, sec)
DDayNorm = J.diff_fast("DDayNorm", loc1, loc2; dimension_filters=omm, sec)
DDCoefficient = J.diff_fast("DDCoefficient", loc1, loc2; dimension_filters=omm, sec)
EEImpact = J.diff_fast("EEImpact", loc1, loc2; dimension_filters=omm, sec)
EESat = J.diff_fast("EESat", loc1, loc2; dimension_filters=omm, sec)
xDmd = J.diff_fast("xDmd", loc1, loc2; dimension_filters=omm, sec)
xDmdTrend = J.diff_fast("xDmdTrend", loc1, loc2; dimension_filters=omm, sec)

DERRRExo = J.diff_fast("DERRRExo", loc1, loc2; dimension_filters=omm, sec)
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters=omm, sec)
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters=omm, sec)
i = findall(vars.Variable .∈ Ref(["DERV", "DERAV"]) .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
DERV = J.diff_fast("DERV", loc1, loc2; dimension_filters=omm, sec) # Problem
select!(DERV, Not([:Enduse, :Tech]))
J.add_pdiff!(DERV)
@rsubset DERV :Vintage == "Vintage 1"
StockAdjustment = J.diff_fast("StockAdjustment", loc1, loc2; dimension_filters=omm, sec) # fine
DERAV = J.diff_fast("DERAV", loc1, loc2; dimension_filters=omm, sec)
DERA   = J.diff_fast("DERA", loc1, loc2; dimension_filters=omm, sec) # problem
J.add_pdiff(DERA)
@rsubset DERA isnan(:Tanoak) :Year == 2020
#   @. DERA = DERAPC+DERAP+DERAD+DERARC
DERAPC = J.diff_fast("DERAPC", loc1, loc2; dimension_filters=omm, sec) #  problem
J.add_pdiff(DERAPC)
DERAP = J.diff_fast("DERAP", loc1, loc2; dimension_filters=omm, sec) # missing from Tanoak
DERAD = J.diff_fast("DERAD", loc1, loc2; dimension_filters=omm, sec) # problem
J.add_pdiff(DERAD)
DERARC = J.diff_fast("DERARC", loc1, loc2; dimension_filters=omm, sec) # fine

RetroSwExo = J.diff_fast("RetroSwExo", loc1, loc2; dimension_filters=omm, sec)

#   @. @finite_math DERAPC = (PERAPC+PERADSt)/DEE
PERAPC = J.diff_fast("PERAPC", loc1, loc2; dimension_filters=omm, sec) # big problem: Missing in Tanoak
J.add_pdiff!(PERAPC)
PERADSt = J.diff_fast("PERADSt", loc1, loc2; dimension_filters=omm, sec) # fine
DEE = J.diff_fast("DEE", loc1, loc2; dimension_filters=omm, sec) # fine
J.add_pdiff!(DEE)

# process additions from production capacity (PERAPC)
# @finite_math PERAPC[enduse,tech,ec,area] = EUPCAPC[enduse,tech,New,ec,area]*
#     DSt[enduse,ec,area]/PEE[enduse,tech,ec,area]
push!(omm, :Age => "New")
EUPCAPC = J.diff_fast("EUPCAPC", loc1, loc2; dimension_filters=omm, sec) # big problem
J.add_pdiff(EUPCAPC)
DSt = J.diff_fast("DSt", loc1, loc2; dimension_filters=omm, sec) # 
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters=omm, sec) # 
J.add_pdiff(PEE)
# EUPCAPC[enduse,tech,New,ec,area] = PCA[New,ecc,area]*MMSF[enduse,tech,ec,area]
PCA = J.diff_fast("PCA", loc1, loc2; dimension_filters=omm, sec) # 
age = pop!(omm, :Age)
PCA_mid = J.diff_fast("PCA", loc1, loc2; dimension_filters=omm, sec) # 
# PCA[New,ecc,area] = (PCPrior[ecc,area]*(ECGR[ecc,area]+GRPGR[area]))+PCR[Old,ecc,area]
push!(omm, :Age => "Old")
PC = J.diff_fast("PC", loc1, loc2; dimension_filters=omm, sec) # about 1/4 of the problem
J.add_pdiff!(PC)
ECGR = J.diff_fast("ECGR", loc1, loc2; dimension_filters=omm, sec) # a lot of problem
J.add_pdiff(ECGR)
GRPGR = J.diff_fast("GRPGR", loc1, loc2; dimension_filters=omm, sec) # fine
PCR = J.diff_fast("PCR", loc1, loc2; dimension_filters=omm, sec) # fine
J.add_pdiff(PCR)
#     @finite_math ECGR[ecc,area] = (PC[ecc,area]-PCPrior[ecc,area])/PCPrior[ecc,area]-GRPGR[area]
# PCMin = J.diff_fast("PCMin", loc1, loc2; dimension_filters=omm, sec) # temp var
xPC = J.diff_fast("xPC", loc1, loc2; dimension_filters=omm, sec) # problem
# @finite_math xPC[ecc,area] = Driver[ecc,area]/xECUF[ecc,area]   
Driver = J.diff_fast("Driver", loc1, loc2; dimension_filters=omm, sec) # 
J.add_pdiff!(Driver)
sort!(Driver, :PDiff)
xECUF = J.diff_fast("xECUF", loc1, loc2; dimension_filters=omm, sec) # 
# Driver[ecc,area] = Driver[ecc,area]*DriverMultiplier[ecc,area]
DriverMultiplier = J.diff_fast("DriverMultiplier", loc1, loc2; dimension_filters=omm, sec) # fine
xDriver = J.diff_fast("xDriver", loc1, loc2; dimension_filters=omm, sec) # issue
DrSwitch = J.diff_fast("DrSwitch", loc1, loc2; dimension_filters=omm, sec) # 21
xGO = J.diff_fast("xGO", loc1, loc2; dimension_filters=omm, sec) # issue
vGO = J.diff_fast("vGO", loc1, loc2; dimension_filters=omm, sec) # doesn't prvide data for CN
    # xGO[ecc,area,year]  =sum(xGOECC[ecc,areatom,year]*
    #   MapAreaTOM[area,areatom] for areatom in areatoms)

# xGOECC = J.diff_fast("xGOECC", loc1, loc2; dimension_filters=omm, sec) # temp var
MapECCfromTOM = J.diff_fast("KInput/MapECCfromTOM", loc1, loc2; dimension_filters=omm, sec='K') # 
@rsubset MapECCfromTOM :Spruce != 0 || :Tanoak != 0
push!(omm, :ECCfromTOM => ["CopperMining","GoldOreMining","OtherMetalMining"])
GY = J.diff_fast("KOutput/GY", loc1, loc2; dimension_filters=omm, sec='K') # 
@rsubset GY :AreaTOM == "ON"
tech = pop!(omm, :Tech)
MMSF = J.diff_fast("MMSF", loc1, loc2; dimension_filters=omm, sec) # 
J.add_pdiff!(MMSF)
@rsubset! MMSF :Spruce != 0 || :Tanoak != 0
xProcSw = J.diff_fast("xProcSw", loc1, loc2; dimension_filters=omm, sec) # 
@rsubset xProcSw :PI == "MShare"
xMMSF = J.diff_fast("xMMSF", loc1, loc2; dimension_filters=omm, sec) # 

# MMSF[enduse,tech,ec,area] = MMSF[enduse,tech,ec,area]+
# max(MMSF[enduse,tech,ec,area]-MMSFB[enduse,tech,ec,area],0)*
# ETSwitch[tech,area] 
ETSwitch = J.diff_fast("ETSwitch", loc1, loc2; dimension_filters=omm, sec) # 

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
pop!(omm, :Tech)
push!(omm, :Year => string.(2022:2023))
push!(omm, :Tech => ["OffRoad"], :Enduse => ["OffRoad"])
pop!(omm, :Tech)
MMSM0 = J.diff_fast("MMSM0", loc1, loc2; dimension_filters=omm, sec) # doesn't match
MMSM0 = J.diff_fast("MMSM0", spr_base, tan_base; dimension_filters=omm, sec) # doesn't match
@rsubset MMSM0 :Spruce > -170 || :Tanoak > -170
MSMM = J.diff_fast("MSMM", loc1, loc2; dimension_filters=omm, sec) # matches
MVF = J.diff_fast("MVF", loc1, loc2; dimension_filters=omm, sec) # matches
MCFU = J.diff_fast("MCFU", loc1, loc2; dimension_filters=omm, sec) # close
MMSMI = J.diff_fast("MMSMI", loc1, loc2; dimension_filters=omm, sec) # 0
AMSF = J.diff_fast("AMSF", loc1, loc2; dimension_filters=omm, sec) # doesn't match
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters=omm, sec) # doesn't match
@rsubset PEE isnan(:Tanoak)
PEESw = J.diff_fast("PEESw", loc1, loc2; dimension_filters=omm, sec) # 1

            # @finite_math MAW[enduse,tech,ec,area] = exp(MMSMI[enduse,tech,ec,area]*
            #   (SPC[ec,area]/SPop[ec,area])/(SPC0[ec,area]/SPop0[ec,area])+
            #   MVF[enduse,tech,ec,area]*log((MCFU[enduse,tech,ec,area]/Inflation[area]/
            #   PEE[enduse,tech,ec,area])/(MCFU0[enduse,tech,ec,area]/Inflation0[area]/PEE0[enduse,tech,ec,area])))

            MMSMI = J.diff_fast("MMSMI", loc1, loc2; dimension_filters=omm, sec) # 0
MVF = J.diff_fast("MVF", loc1, loc2; dimension_filters=omm, sec) # matches
MCFU = J.diff_fast("MCFU", loc1, loc2; dimension_filters=omm, sec) # close
AMSF = J.diff_fast("AMSF", loc1, loc2; dimension_filters=omm, sec) # doesn't match
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters=omm, sec) # doesn't match



Enduse = pop!(omm, :Enduse)
push!(omm, :Enduse => Enduse)
PEEBeforeStd = J.diff_fast("PEEBeforeStd", loc1, loc2; dimension_filters=omm, sec) # doesn't match
PEEBeforeStd = J.diff_fast("PEEBeforeStd", loc1, loc2; dimension_filters, sec) # doesn't match
J.add_pdiff!(PEEBeforeStd)
@rsubset! PEEBeforeStd :Year >= 2020 :Year <= 2030
sort!(PEEBeforeStd, :PDiff)
nas = copy(dimension_filters)
push!(nas, :Year => "2020")
push!(nas, :Tech => ["Gas", "Oil", "Coal"])
push!(nas, :Enduse => "Steam")
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
PEEPrice = J.diff_fast("PEEPrice", loc1, loc2; dimension_filters=nas, sec) # issue
J.add_pdiff!(PEEPrice)
PEM = J.diff_fast("PEM", loc1, loc2; dimension_filters=nas, sec) # matches
PEMM = J.diff_fast("PEMM", loc1, loc2; dimension_filters=nas, sec) # matches
# PEEBeforeStd matches PEM = PEM*PEMM which means in Spruce, Electric and Oil are getting the else values
PEPM = J.diff_fast("PEPM", loc1, loc2; dimension_filters=nas, sec) # 1
PFPN = J.diff_fast("PFPN", loc1, loc2; dimension_filters=nas, sec) # matches
PFTC = J.diff_fast("PFTC", loc1, loc2; dimension_filters=nas, sec) # matches
PCCR = J.diff_fast("PCCR", loc1, loc2; dimension_filters=nas, sec) # matches
PCCRB = J.diff_fast("PCCR", spr_base, tan_base; dimension_filters=nas, sec) # matches
PCTC = J.diff_fast("PCTC", loc1, loc2; dimension_filters=nas, sec) # matches

df = leftjoin(PEEBeforeStd,PCTC, on=[:Enduse,:Tech,:EC,:Area,:Year],  renamecols = "_PEEB" => "_PCTC")
@rsubset df :Spruce_PCTC == 0
PFPN = J.diff_fast("PFPN", loc1, loc2; dimension_filters=nas, sec) # matchesl
df = leftjoin(df, PFPN, on = [:Enduse, :Tech, :EC, :Area], renamecols = "" => "_PFPN")
@rsubset! df :Enduse != "ommRoad"

# Sometime PEE is getting halved and sometimes not. 
# For EC(Food) Enduse(OthSub), Tech(Oil) get halved but Tech(ommRoad) does not.
# Let's findout why
food = copy(omm)
PEEBeforeStd_f = J.diff_fast("PEEBeforeStd", loc1, loc2; dimension_filters=food, sec) # doesn't match
PEEBeforeStd # Doesn't match between Oil and ommRoad in Spruce
push!(food, :EC => "Food", :Enduse => "OthSub", :Tech => ["Oil", "ommRoad"])
PEM = J.diff_fast("PEM", loc1, loc2; dimension_filters=food, sec) # doesn't have tech
PEMM = J.diff_fast("PEMM", loc1, loc2; dimension_filters=food, sec) # matches
PCCR = J.diff_fast("PCCR", loc1, loc2; dimension_filters=food, sec) # matches
PCCRB = J.diff_fast("PCCR", spr_base, tan_base; dimension_filters=food, sec) # matches
PFTC = J.diff_fast("PFTC", loc1, loc2; dimension_filters=food, sec) # matches

PEPM = J.diff_fast("PEPM", loc1, loc2; dimension_filters=food, sec) # 1
PFPN = J.diff_fast("PFPN", loc1, loc2; dimension_filters=food, sec) # matches

MCFU = J.diff_fast("MCFU", loc1, loc2; dimension_filters=food, sec) # 1
Inflation = J.diff_fast("Inflation", loc1, loc2; dimension_filters=food, sec) # matches



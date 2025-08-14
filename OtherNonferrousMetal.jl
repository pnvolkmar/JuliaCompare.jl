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
loc1 = J.Loc_p(vars, DATA_FOLDER1, "Pine");
loc2 = J.Loc_j(vars_j, joinpath(DATA_FOLDER2, "database.hdf5"), "Redwood");
tan_base = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaTanoak\\2020Model\\Base\\database.hdf5", "Tan_Base");
spr_base = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaSpruce\\2020Model\\Base", "Spr_Base");
tan_ogref = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaTanoak\\2020Model\\OGRef\\database.hdf5", "Tan_OGRef");
spr_ogref = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaSpruce\\2020Model\\OGRef", "Spr_OGRef");
locs = [locs]
################################################################################
# Initialize filters and sec for analysis ######################################
################################################################################
sec = 'I'
filter = Dict{Symbol,Any}()
push!(filter, :Area => Canada, :Year => string.(1986:2050))
ECC = M.ReadDisk(loc2.HDF5_path, "E2020DB/ECCKey")
push!(filter, :ECC => ECC[occursin.("Foreign", ECC) .== false])

omm = copy(filter)
push!(omm, :ECC => "OtherMetalMining", :EC => "OtherMetalMining")
push!(omm, :Poll => "CO2", :Year => string.(2023:2024))
push!(omm, :Area => ["NU"])
push!(omm, :FuelEP => ["Diesel","Gasoline", "Kerosene"], :Fuel => ["Diesel","Gasoline", "Kerosene"])
push!(omm, :Enduse => ["OffRoad", "Heat"])
push!(omm, :Tech => ["OffRoad", "Oil"])

################################################################################
# Analysis of variables: Overview ##############################################
################################################################################

TotPol = J.var("TotPol", locs; filter=omm, sec, diff = true, pdiff = true) # 
J.plot_sets(TotPol; col = :Pine_minus_Redwood, dim="Area", title = "Differences in OtherMetalMining TotPol")

    # TotPol[ecc,poll,area] = EnPol[ecc,poll,area]+NcPol[ecc,poll,area]+
    #   MEPol[ecc,poll,area]+VnPol[ecc,poll,area]+FlPol[ecc,poll,area]+
    #   FuPol[ecc,poll,area]+ORMEPol[ecc,poll,area]+SqPol[ecc,poll,area]
EnPol = J.var("EnPol", locs; filter=omm, sec) # Problem
NcPol = J.var("NcPol", locs; filter=omm, sec) # 
MEPol = J.var("MEPol", locs; filter=omm, sec) # 
VnPol = J.var("VnPol", locs; filter=omm, sec) # 
FlPol = J.var("FlPol", locs; filter=omm, sec) # 
FuPol = J.var("FuPol", locs; filter=omm, sec) # 
ORMEPol = J.var("ORMEPol", locs; filter=omm, sec) # 
SqPol = J.var("SqPol", locs; filter=omm, sec) # 

EuFPol = J.var("EuFPol", locs; filter=omm, sec) # Problem
@rsubset! EuFPol :Pine != :Redwood
J.plot_diff(EuFPol, dim = "FuelEP"; title = "EuFPol OtherMetalMining CO2, BC,ON,QC, by FuelEP")

Polute = J.var("Polute", locs; filter=omm, sec = sec, diff = true, pdiff = true) # this is fine
@rsubset! Polute abs(:Pine_pdiff_Redwood) != 0
J.plot_diff(Polute, dim = "Enduse")

i = findall(vars.Variable .== "POCA" .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
# Polute[enduse,fuelep,ec,poll,area] = EuDem[enduse,fuelep,ec,area]*
#     POCA[enduse,fuelep,ec,poll,area]+xPolute[enduse,fuelep,ec,poll,area]

POCA = J.var("POCA", locs; filter=omm, sec) # 
EuDem = J.var("EuDem", locs; filter=omm, sec) # Problem
@rsubset! EuDem :Diff != 0
J.add_pdiff!(EuDem)
@rsubset EuDem :Area == "BC" :Year == 2050
xPolute = J.var("xPolute", locs; filter=omm, sec) # 
J.add_pdiff!(EuDem)
@rsubset! EuDem abs(:PDiff) != 0

EuDemF = J.var("EuDemF", locs; filter=omm, sec) # not on db

# EuDemF[enduse,fuel,ecc,area] = sum(Dmd[enduse,tech,ec,area]*
# DmFrac[enduse,fuel,tech,ec,area] for tech in Techs)
Dmd = J.var("Dmd", locs; filter=omm, sec) # problem
@rsubset! Dmd :Diff != 0 :Tech == "ommRoad"
J.add_pdiff!(Dmd)
@rsubset Dmd :PDiff != 0
DmFrac = J.var("DmFrac", locs; filter=omm, sec) # fine
@rsubset! DmFrac :Diff != 0
J.add_pdiff!(DmFrac)
zoom = copy(omm)
push!(zoom, :Area => "QC", :Year => "2050")
pop!(zoom, :Fuel)
DmFrac = J.var("DmFrac", locs; filter=zoom, sec) # 
DmFrac = J.var("DmFrac", spr_base, tan_base; filter=omm, sec) # 
@rsubset DmFrac :Area == "BC" :Fuel == "Diesel"

@rsubset! DmFrac :Tanoak != 0
xDmFrac = J.var("xDmFrac", locs; filter=zoom, sec) # 
@rsubset! xDmFrac :Tanoak != 0
DmFracMin = J.var("DmFracMin", locs; filter=omm, sec) # 
@rsubset! DmFracMin :Tanoak != 0
@rsubset DmFracMin :Area == "QC" :Fuel == "Biodiesel"
@rsubset DmFracMin :Area == "BC" :Fuel == "Biodiesel"
DmFracMax = J.var("DmFracMax", locs; filter=zoom, sec) # 
@rsubset DmFracMax :Tanoak != 0 :Diff != 0

# @. Dmd = Dmd-EE
EE = J.var("EE", locs; filter=omm, sec) # 
@rsubset
# @. Dmd = Dmd-DSMEU
DSMEU = J.var("DSMEU", locs; filter=omm, sec) # 
# Dmd[enduse,tech,ec,area] = Dmd[enduse,tech,ec,area]+
#                            SqDmd[tech,ec,area]*SqEUTechMap[enduse,tech]
SqDmd = J.var("SqDmd", locs; filter=omm, sec) # 
SqEUTechMap = J.var("SqEUTechMap", locs; filter=omm, sec) # 
#
# @finite_math Dmd[enduse,tech,ec,area] = DER[enduse,tech,ec,area]*
# UMS[enduse,tech,ec,area]*CERSM[enduse,ec,area]*
# CUF[enduse,tech,ec,area]*WCUF[ec,area]*RPEI[enduse,tech,ec,area]/1.0e6*
# (TSLoad[enduse,ec,area]*
# (DDay[enduse,area]/DDayNorm[enduse,area])^DDCoefficient[enduse,ec,area]+
# (1.0-TSLoad[enduse,ec,area]))
DER = J.var("DER", locs; filter=omm, sec, pdiff = true) # Problem
@rsubset DER :Pine_pdiff_Redwood != 0
UMS = J.var("UMS", locs; filter=omm, sec)
CERSM = J.var("CERSM", locs; filter=omm, sec)
CUF = J.var("CUF", locs; filter=omm, sec)
WCUF = J.var("WCUF", locs; filter=omm, sec)
RPEI = J.var("RPEI", locs; filter=omm, sec)
TSLoad = J.var("TSLoad", locs; filter=omm, sec)
DDay = J.var("DDay", locs; filter=omm, sec)
DDayNorm = J.var("DDayNorm", locs; filter=omm, sec)
DDCoefficient = J.var("DDCoefficient", locs; filter=omm, sec)
EEImpact = J.var("EEImpact", locs; filter=omm, sec)
EESat = J.var("EESat", locs; filter=omm, sec)
xDmd = J.var("xDmd", locs; filter=omm, sec)
xDmdTrend = J.var("xDmdTrend", locs; filter=omm, sec)

DERRRExo = J.var("DERRRExo", locs; filter=omm, sec)
PERRRExo = J.var("PERRRExo", locs; filter=omm, sec)
PERRRExo = J.var("PERRRExo", locs; filter=omm, sec)
i = findall(vars.Variable .∈ Ref(["DERV", "DERAV"]) .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
DERV = J.var("DERV", locs; filter=omm, sec) # Problem
select!(DERV, Not([:Enduse, :Tech]))
J.add_pdiff!(DERV)
@rsubset DERV :Vintage == "Vintage 1"
StockAdjustment = J.var("StockAdjustment", locs; filter=omm, sec) # fine
DERAV = J.var("DERAV", locs; filter=omm, sec)
DERA   = J.var("DERA", locs; filter=omm, sec) # problem
J.add_pdiff(DERA)
@rsubset DERA isnan(:Tanoak) :Year == 2020
#   @. DERA = DERAPC+DERAP+DERAD+DERARC
DERAPC = J.var("DERAPC", locs; filter=omm, sec) #  problem
J.add_pdiff(DERAPC)
DERAP = J.var("DERAP", locs; filter=omm, sec) # missing from Tanoak
DERAD = J.var("DERAD", locs; filter=omm, sec) # problem
J.add_pdiff(DERAD)
DERARC = J.var("DERARC", locs; filter=omm, sec) # fine

RetroSwExo = J.var("RetroSwExo", locs; filter=omm, sec)

#   @. @finite_math DERAPC = (PERAPC+PERADSt)/DEE
PERAPC = J.var("PERAPC", locs; filter=omm, sec) # big problem: Missing in Tanoak
J.add_pdiff!(PERAPC)
PERADSt = J.var("PERADSt", locs; filter=omm, sec) # fine
DEE = J.var("DEE", locs; filter=omm, sec) # fine
J.add_pdiff!(DEE)

# process additions from production capacity (PERAPC)
# @finite_math PERAPC[enduse,tech,ec,area] = EUPCAPC[enduse,tech,New,ec,area]*
#     DSt[enduse,ec,area]/PEE[enduse,tech,ec,area]
push!(omm, :Age => "New")
EUPCAPC = J.var("EUPCAPC", locs; filter=omm, sec) # big problem
J.add_pdiff(EUPCAPC)
DSt = J.var("DSt", locs; filter=omm, sec) # 
PEE = J.var("PEE", locs; filter=omm, sec) # 
J.add_pdiff(PEE)
# EUPCAPC[enduse,tech,New,ec,area] = PCA[New,ecc,area]*MMSF[enduse,tech,ec,area]
PCA = J.var("PCA", locs; filter=omm, sec) # 
age = pop!(omm, :Age)
PCA_mid = J.var("PCA", locs; filter=omm, sec) # 
# PCA[New,ecc,area] = (PCPrior[ecc,area]*(ECGR[ecc,area]+GRPGR[area]))+PCR[Old,ecc,area]
push!(omm, :Age => "Old")
PC = J.var("PC", locs; filter=omm, sec) # about 1/4 of the problem
J.add_pdiff!(PC)
ECGR = J.var("ECGR", locs; filter=omm, sec) # a lot of problem
J.add_pdiff(ECGR)
GRPGR = J.var("GRPGR", locs; filter=omm, sec) # fine
PCR = J.var("PCR", locs; filter=omm, sec) # fine
J.add_pdiff(PCR)
#     @finite_math ECGR[ecc,area] = (PC[ecc,area]-PCPrior[ecc,area])/PCPrior[ecc,area]-GRPGR[area]
# PCMin = J.var("PCMin", locs; filter=omm, sec) # temp var
xPC = J.var("xPC", locs; filter=omm, sec) # problem
# @finite_math xPC[ecc,area] = Driver[ecc,area]/xECUF[ecc,area]   
Driver = J.var("Driver", locs; filter=omm, sec) # 
J.add_pdiff!(Driver)
sort!(Driver, :PDiff)
xECUF = J.var("xECUF", locs; filter=omm, sec) # 
# Driver[ecc,area] = Driver[ecc,area]*DriverMultiplier[ecc,area]
DriverMultiplier = J.var("DriverMultiplier", locs; filter=omm, sec) # fine
xDriver = J.var("xDriver", locs; filter=omm, sec) # issue
DrSwitch = J.var("DrSwitch", locs; filter=omm, sec) # 21
xGO = J.var("xGO", locs; filter=omm, sec) # issue
vGO = J.var("vGO", locs; filter=omm, sec) # doesn't prvide data for CN
    # xGO[ecc,area,year]  =sum(xGOECC[ecc,areatom,year]*
    #   MapAreaTOM[area,areatom] for areatom in areatoms)

# xGOECC = J.var("xGOECC", locs; filter=omm, sec) # temp var
MapECCfromTOM = J.var("KInput/MapECCfromTOM", locs; filter=omm, sec='K') # 
@rsubset MapECCfromTOM :Spruce != 0 || :Tanoak != 0
push!(omm, :ECCfromTOM => ["CopperMining","GoldOreMining","OtherMetalMining"])
GY = J.var("KOutput/GY", locs; filter=omm, sec='K') # 
@rsubset GY :AreaTOM == "ON"
tech = pop!(omm, :Tech)
MMSF = J.var("MMSF", locs; filter=omm, sec) # 
J.add_pdiff!(MMSF)
@rsubset! MMSF :Spruce != 0 || :Tanoak != 0
xProcSw = J.var("xProcSw", locs; filter=omm, sec) # 
@rsubset xProcSw :PI == "MShare"
xMMSF = J.var("xMMSF", locs; filter=omm, sec) # 

# MMSF[enduse,tech,ec,area] = MMSF[enduse,tech,ec,area]+
# max(MMSF[enduse,tech,ec,area]-MMSFB[enduse,tech,ec,area],0)*
# ETSwitch[tech,area] 
ETSwitch = J.var("ETSwitch", locs; filter=omm, sec) # 

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
MMSM0 = J.var("MMSM0", locs; filter=omm, sec, pdiff = true) 
@rsubset MMSM0 :Pine_pdiff_Redwood != 0# doesn't match
MMSM0 = J.var("MMSM0", spr_base, tan_base; filter=omm, sec) # doesn't match
@rsubset MMSM0 :Spruce > -170 || :Tanoak > -170
MSMM = J.var("MSMM", locs; filter=omm, sec) # matches
MVF = J.var("MVF", locs; filter=omm, sec) # matches
MCFU = J.var("MCFU", locs; filter=omm, sec) # close
omm_init = copy(omm)
push!(omm_init, :Year => "1985")
MCFU0 = J.var("MCFU", locs; filter=omm_init, sec) # close

MMSMI = J.var("MMSMI", locs; filter=omm, sec) # 0
AMSF = J.var("AMSF", locs; filter=omm, sec) # doesn't match
PEE = J.var("PEE", locs; filter=omm, sec) # doesn't match
PEE0 = J.var("PEE", locs; filter=omm_init, sec) # doesn't match
@rsubset PEE isnan(:Tanoak)
MMSM0 = J.var("MMSM0", locs; filter=omm, sec) # doesn't match
PEESw = J.var("PEESw", locs; filter=omm, sec) # 1

            # @finite_math MAW[enduse,tech,ec,area] = exp(MMSMI[enduse,tech,ec,area]*
            #   (SPC[ec,area]/SPop[ec,area])/(SPC0[ec,area]/SPop0[ec,area])+
            #   MVF[enduse,tech,ec,area]*log((MCFU[enduse,tech,ec,area]/Inflation[area]/
            #   PEE[enduse,tech,ec,area])/(MCFU0[enduse,tech,ec,area]/Inflation0[area]/PEE0[enduse,tech,ec,area])))
push!(omm, :Tech => "Electric")
MMSMI = J.var("MMSMI", locs; filter=omm, sec) # 0
MVF = J.var("MVF", locs; filter=omm, sec) # matches
MCFU = J.var("MCFU", locs; filter=omm, sec) # close
AMSF = J.var("AMSF", locs; filter=omm, sec) # doesn't match
PEE = J.var("PEE", locs; filter=omm, sec) # doesn't match



Enduse = pop!(omm, :Enduse)
push!(omm, :Enduse => Enduse)
PEEBeforeStd = J.var("PEEBeforeStd", locs; filter=omm, sec) # doesn't match
PEEBeforeStd = J.var("PEEBeforeStd", locs; filter, sec) # doesn't match
J.add_pdiff!(PEEBeforeStd)
@rsubset! PEEBeforeStd :Year >= 2020 :Year <= 2030
sort!(PEEBeforeStd, :PDiff)
nas = copy(filter)
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
PEEPrice = J.var("PEEPrice", locs; filter=nas, sec) # issue
J.add_pdiff!(PEEPrice)
PEM = J.var("PEM", locs; filter=nas, sec) # matches
PEMM = J.var("PEMM", locs; filter=nas, sec) # matches
# PEEBeforeStd matches PEM = PEM*PEMM which means in Spruce, Electric and Oil are getting the else values
PEPM = J.var("PEPM", locs; filter=nas, sec) # 1
PFPN = J.var("PFPN", locs; filter=nas, sec) # matches
PFTC = J.var("PFTC", locs; filter=nas, sec) # matches
PCCR = J.var("PCCR", locs; filter=nas, sec) # matches
PCCRB = J.var("PCCR", spr_base, tan_base; filter=nas, sec) # matches
PCTC = J.var("PCTC", locs; filter=nas, sec) # matches

df = leftjoin(PEEBeforeStd,PCTC, on=[:Enduse,:Tech,:EC,:Area,:Year],  renamecols = "_PEEB" => "_PCTC")
@rsubset df :Spruce_PCTC == 0
PFPN = J.var("PFPN", locs; filter=nas, sec) # matchesl
df = leftjoin(df, PFPN, on = [:Enduse, :Tech, :EC, :Area], renamecols = "" => "_PFPN")
@rsubset! df :Enduse != "ommRoad"

# Sometime PEE is getting halved and sometimes not. 
# For EC(Food) Enduse(OthSub), Tech(Oil) get halved but Tech(ommRoad) does not.
# Let's findout why
food = copy(omm)
PEEBeforeStd_f = J.var("PEEBeforeStd", locs; filter=food, sec) # doesn't match
PEEBeforeStd # Doesn't match between Oil and ommRoad in Spruce
push!(food, :EC => "Food", :Enduse => "OthSub", :Tech => ["Oil", "ommRoad"])
PEM = J.var("PEM", locs; filter=food, sec) # doesn't have tech
PEMM = J.var("PEMM", locs; filter=food, sec) # matches
PCCR = J.var("PCCR", locs; filter=food, sec) # matches
PCCRB = J.var("PCCR", spr_base, tan_base; filter=food, sec) # matches
PFTC = J.var("PFTC", locs; filter=food, sec) # matches

PEPM = J.var("PEPM", locs; filter=food, sec) # 1
PFPN = J.var("PFPN", locs; filter=food, sec) # matches

MCFU = J.var("MCFU", locs; filter=food, sec) # 1
Inflation = J.var("Inflation", locs; filter=food, sec) # matches



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

################################################################################
# Initialize filters and sec for analysis ######################################
################################################################################
sec = 'T'
dimension_filters = Dict{Symbol,Any}()
push!(dimension_filters, :Area => Canada, :Year => string.(1986:2050))
ECC = M.ReadDisk(loc2.HDF5_path, "E2020DB/ECCKey")
push!(dimension_filters, :ECC => ECC[occursin.("Foreign", ECC) .== false])

frt = copy(dimension_filters)
push!(frt, :EC => "Freight",  :ECC => "Freight")
push!(frt, :Poll => "CO2")

push!(frt, :Area => ["NT"])
push!(frt, :AreaTOM => ["NT"])
push!(frt, :Fuel => "Diesel", :FuelEP => "Diesel")


################################################################################
# Analysis of variables: Overview ##############################################
################################################################################

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters, sec) # 
@rsubset! TotPol abs(:Diff) >= 1
J.add_pdiff!(TotPol)
J.plot_diff(TotPol, dim="ECC"; title = "Differences in TotPol")

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters=frt, sec) # 
J.add_pdiff!(TotPol)
J.plot_diff(TotPol, dim="Poll")
J.plot_diff(TotPol, dim="Area")

    # TotPol[ecc,poll,area] = EnPol[ecc,poll,area]+NcPol[ecc,poll,area]+
    #   MEPol[ecc,poll,area]+VnPol[ecc,poll,area]+FlPol[ecc,poll,area]+
    #   FuPol[ecc,poll,area]+ORMEPol[ecc,poll,area]+SqPol[ecc,poll,area]
EnPol = J.diff_fast("EnPol", loc1, loc2; dimension_filters=frt, sec) # Problem
NcPol = J.diff_fast("NcPol", loc1, loc2; dimension_filters=frt, sec) # 
MEPol = J.diff_fast("MEPol", loc1, loc2; dimension_filters=frt, sec) # 
VnPol = J.diff_fast("VnPol", loc1, loc2; dimension_filters=frt, sec) # 
FlPol = J.diff_fast("FlPol", loc1, loc2; dimension_filters=frt, sec) # 
FuPol = J.diff_fast("FuPol", loc1, loc2; dimension_filters=frt, sec) # 
ORMEPol = J.diff_fast("ORMEPol", loc1, loc2; dimension_filters=frt, sec) # 
SqPol = J.diff_fast("SqPol", loc1, loc2; dimension_filters=frt, sec) # 

EuFPol = J.diff_fast("EuFPol", loc1, loc2; dimension_filters=frt, sec) # Problem
@rsubset! EuFPol abs(:Diff) >= 1
J.plot_diff(EuFPol, dim = "FuelEP")

Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters=frt, sec = sec) # this is fine
J.add_pdiff!(Polute)
@rsubset! Polute abs(:PDiff) != 0
J.plot_diff(Polute, dim = "Tech")

i = findall(vars.Variable .== "POCA" .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
# Polute[enduse,fuelep,ec,poll,area] = EuDem[enduse,fuelep,ec,area]*
#     POCA[enduse,fuelep,ec,poll,area]+xPolute[enduse,fuelep,ec,poll,area]

POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters=frt, sec) # 
@rsubset! POCA :Diff != 0
EuDem = J.diff_fast("EuDem", loc1, loc2; dimension_filters=frt, sec) # Problem
@rsubset! EuDem abs(:Diff) >= 0.0001
J.add_pdiff!(EuDem)

xPolute = J.diff_fast("xPolute", loc1, loc2; dimension_filters=frt, sec) # doesn't exist in T
J.add_pdiff!(xPolute)
@rsubset! xPolute abs(:PDiff) != 0

# EuDemF[enduse,fuel,ecc,area] = sum(Dmd[enduse,tech,ec,area]*
# DmFrac[enduse,fuel,tech,ec,area] for tech in Techs)
Dmd = J.diff_fast("Dmd", loc1, loc2; dimension_filters=frt, sec) # problem
J.add_pdiff!(Dmd)
@rsubset! Dmd :PDiff != 0
xDmd = J.diff_fast("xDmd", loc1, loc2; dimension_filters=frt, sec) # problem
J.add_pdiff!(xDmd)
@rsubset! xDmd :Redwood != 0
DmFrac = J.diff_fast("DmFrac", loc1, loc2; dimension_filters=frt, sec) # fine
@rsubset! DmFrac :Diff != 0
J.add_pdiff!(DmFrac)


# @rsubset! DmFrac :Tanoak != 0
# xDmFrac = J.diff_fast("xDmFrac", loc1, loc2; dimension_filters=zoom, sec) # 
# @rsubset! xDmFrac :Tanoak != 0
# DmFracMin = J.diff_fast("DmFracMin", loc1, loc2; dimension_filters=frt, sec) # 
# @rsubset! DmFracMin :Tanoak != 0
# @rsubset DmFracMin :Area == "QC" :Fuel == "Biodiesel"
# @rsubset DmFracMin :Area == "BC" :Fuel == "Biodiesel"
# DmFracMax = J.diff_fast("DmFracMax", loc1, loc2; dimension_filters=zoom, sec) # 
# @rsubset DmFracMax :Tanoak != 0 :Diff != 0

# @. Dmd = Dmd-EE
EE = J.diff_fast("EE", loc1, loc2; dimension_filters=frt, sec) # 
@rsubset EE :Diff != 0
# @. Dmd = Dmd-DSMEU
DSMEU = J.diff_fast("DSMEU", loc1, loc2; dimension_filters=frt, sec) # not going to worry about now
@rsubset  DSMEU :Diff != 0
# Dmd[enduse,tech,ec,area] = Dmd[enduse,tech,ec,area]+
#                            SqDmd[tech,ec,area]*SqEUTechMap[enduse,tech]
SqDmd = J.diff_fast("SqDmd", loc1, loc2; dimension_filters=frt, sec) #  not going to worry aobut now
@rsubset! SqDmd :Diff != 0
J.add_pdiff!(SqDmd)
SqEUTechMap = J.diff_fast("SqEUTechMap", loc1, loc2; dimension_filters=frt, sec) # 
#
# @finite_math Dmd[enduse,tech,ec,area] = DER[enduse,tech,ec,area]*
# UMS[enduse,tech,ec,area]*CERSM[enduse,ec,area]*
# CUF[enduse,tech,ec,area]*WCUF[ec,area]*RPEI[enduse,tech,ec,area]/1.0e6*
# (TSLoad[enduse,ec,area]*
# (DDay[enduse,area]/DDayNorm[enduse,area])^DDCoefficient[enduse,ec,area]+
# (1.0-TSLoad[enduse,ec,area]))
DER = J.diff_fast("DER", loc1, loc2; dimension_filters=frt, sec) # Problem
J.add_pdiff!(DER)
@rsubset! DER :PDiff != 0
UMS = J.diff_fast("UMS", loc1, loc2; dimension_filters=frt, sec)
CERSM = J.diff_fast("CERSM", loc1, loc2; dimension_filters=frt, sec)
CUF = J.diff_fast("CUF", loc1, loc2; dimension_filters=frt, sec)
WCUF = J.diff_fast("WCUF", loc1, loc2; dimension_filters=frt, sec)
RPEI = J.diff_fast("RPEI", loc1, loc2; dimension_filters=frt, sec)
TSLoad = J.diff_fast("TSLoad", loc1, loc2; dimension_filters=frt, sec)
DDay = J.diff_fast("DDay", loc1, loc2; dimension_filters=frt, sec)
DDayNorm = J.diff_fast("DDayNorm", loc1, loc2; dimension_filters=frt, sec)
DDCoefficient = J.diff_fast("DDCoefficient", loc1, loc2; dimension_filters=frt, sec)
EEImpact = J.diff_fast("EEImpact", loc1, loc2; dimension_filters=frt, sec)
EESat = J.diff_fast("EESat", loc1, loc2; dimension_filters=frt, sec)
xDmd = J.diff_fast("xDmd", loc1, loc2; dimension_filters=frt, sec)
xDmdTrend = J.diff_fast("xDmdTrend", loc1, loc2; dimension_filters=frt, sec)

DERRRExo = J.diff_fast("DERRRExo", loc1, loc2; dimension_filters=frt, sec) # not going to worry about rn
J.add_pdiff!(DERRRExo)
@rsubset! DERRRExo :PDiff != 0
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters=frt, sec)
J.add_pdiff!(PERRRExo)
@rsubset! PERRRExo :PDiff != 0
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters=frt, sec)
J.add_pdiff!(PERRRExo)
@rsubset! PERRRExo :PDiff != 0
i = findall(vars.Variable .∈ Ref(["DERV", "DERAV"]) .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
DERV = J.diff_fast("DERV", loc1, loc2; dimension_filters=frt, sec) # Problem
select!(DERV, Not([:Enduse, :Tech]))
J.add_pdiff!(DERV)
@rsubset DERV :Vintage == "Vintage 1"
StockAdjustment = J.diff_fast("StockAdjustment", loc1, loc2; dimension_filters=frt, sec) # fine
DERAV = J.diff_fast("DERAV", loc1, loc2; dimension_filters=frt, sec)
DERA   = J.diff_fast("DERA", loc1, loc2; dimension_filters, sec) # problem
J.add_pdiff!(DERA)
issue_DERA = @rsubset DERA isnan(:Pine) 
unique(issue_DERA.EC)
unique(issue_DERA.Area)
@by issue_DERA [:EC,:Area] :first_yr = minimum(:Year)
#   @. DERA = DERAPC+DERAP+DERAD+DERARC
push!(frt, :EC => ["Freight","AirFreight","ForeignFreight"])
DERAPC = J.diff_fast("DERAPC", loc1, loc2; dimension_filters=frt, sec) #  problem
J.add_pdiff!(DERAPC)
DERAP = J.diff_fast("DERAP", loc1, loc2; dimension_filters=frt, sec) # problem
DERAD = J.diff_fast("DERAD", loc1, loc2; dimension_filters=frt, sec) # problem
J.add_pdiff!(DERAD)
DERARC = J.diff_fast("DERARC", loc1, loc2; dimension_filters=frt, sec) # fine

RetroSwExo = J.diff_fast("RetroSwExo", loc1, loc2; dimension_filters=frt, sec) # -99

#   @. @finite_math DERAPC = (PERAPC+PERADSt)/DEE
PERAPC = J.diff_fast("PERAPC", loc1, loc2; dimension_filters=frt, sec) # big problem: Missing in Tanoak
J.add_pdiff!(PERAPC)
PERADSt = J.diff_fast("PERADSt", loc1, loc2; dimension_filters=frt, sec) # fine
@rsubset! PERADSt :Diff != 0
DEE = J.diff_fast("DEE", loc1, loc2; dimension_filters=frt, sec) # fine
J.add_pdiff!(DEE)

# process additions from production capacity (PERAPC)
# @finite_math PERAPC[enduse,tech,ec,area] = EUPCAPC[enduse,tech,New,ec,area]*
#     DSt[enduse,ec,area]/PEE[enduse,tech,ec,area]
push!(frt, :Age => "New")
EUPCAPC = J.diff_fast("EUPCAPC", loc1, loc2; dimension_filters=frt, sec) # big problem
J.add_pdiff!(EUPCAPC)
@rsubset EUPCAPC isnan(:PDiff)
DSt = J.diff_fast("DSt", loc1, loc2; dimension_filters=frt, sec) # 
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters=frt, sec) # 
J.add_pdiff!(PEE)
# EUPCAPC[enduse,tech,New,ec,area] = PCA[New,ecc,area]*MMSF[enduse,tech,ec,area]
PCA = J.diff_fast("PCA", loc1, loc2; dimension_filters=frt, sec) # 
@rsubset PCA isnan(:Pine)
MMSF = J.diff_fast("MMSF", loc1, loc2; dimension_filters=frt, sec) # 
@rsubset MMSF isnan(:Pine)
age = pop!(frt, :Age)
PCA_mid = J.diff_fast("PCA", loc1, loc2; dimension_filters=frt, sec) # 
# PCA[New,ecc,area] = (PCPrior[ecc,area]*(ECGR[ecc,area]+GRPGR[area]))+PCR[Old,ecc,area]
push!(frt, :Age => "Old")
PC = J.diff_fast("PC", loc1, loc2; dimension_filters=frt, sec) # problem
J.add_pdiff!(PC)
ECGR = J.diff_fast("ECGR", loc1, loc2; dimension_filters=frt, sec) # problem
J.add_pdiff!(ECGR)
GRPGR = J.diff_fast("GRPGR", loc1, loc2; dimension_filters=frt, sec) # fine
PCR = J.diff_fast("PCR", loc1, loc2; dimension_filters=frt, sec) # problem
J.add_pdiff!(PCR)
#     @finite_math ECGR[ecc,area] = (PC[ecc,area]-PCPrior[ecc,area])/PCPrior[ecc,area]-GRPGR[area]
# PCMin = J.diff_fast("PCMin", loc1, loc2; dimension_filters=frt, sec) # temp var
xPC = J.diff_fast("xPC", loc1, loc2; dimension_filters=frt, sec) # problem
# @finite_math xPC[ecc,area] = Driver[ecc,area]/xECUF[ecc,area]   
Driver = J.diff_fast("Driver", loc1, loc2; dimension_filters=frt, sec) # 
J.add_pdiff!(Driver)
sort!(Driver, :PDiff)
xECUF = J.diff_fast("xECUF", loc1, loc2; dimension_filters=frt, sec) # 
# Driver[ecc,area] = Driver[ecc,area]*DriverMultiplier[ecc,area]
DriverMultiplier = J.diff_fast("DriverMultiplier", loc1, loc2; dimension_filters=frt, sec) # fine
xDriver = J.diff_fast("xDriver", loc1, loc2; dimension_filters=frt, sec) # issue
DrSwitch = J.diff_fast("DrSwitch", loc1, loc2; dimension_filters=frt, sec) # 21
xGO = J.diff_fast("xGO", loc1, loc2; dimension_filters=frt, sec) # issue
vGO = J.diff_fast("vGO", loc1, loc2; dimension_filters=frt, sec) # doesn't prvide data for CN
    # xGO[ecc,area,year]  =sum(xGOECC[ecc,areatom,year]*
    #   MapAreaTOM[area,areatom] for areatom in areatoms)

# xGOECC = J.diff_fast("xGOECC", loc1, loc2; dimension_filters=frt, sec) # temp var
MapECCfromTOM = J.diff_fast("KInput/MapECCfromTOM", loc1, loc2; dimension_filters=frt, sec='K') # 
@rsubset MapECCfromTOM :Pine != 0 || :Redwood != 0
push!(frt, :ECCfromTOM => ["Fertilizer"])
GY = J.diff_fast("KOutput/GY", loc1, loc2; dimension_filters=frt, sec='K') # 
@rsubset GY :AreaTOM == "ON"
tech = pop!(frt, :Tech)
MMSF = J.diff_fast("MMSF", loc1, loc2; dimension_filters=frt, sec) # 
J.add_pdiff!(MMSF)
@rsubset! MMSF :Spruce != 0 || :Tanoak != 0
xProcSw = J.diff_fast("xProcSw", loc1, loc2; dimension_filters=frt, sec) # 
@rsubset xProcSw :PI == "MShare" :Year ∈ [2007,2008] # 1
xMMSF = J.diff_fast("xMMSF", loc1, loc2; dimension_filters=frt, sec) # 
@rsubset xMMSF isnan(:Pine)
# MMSF[enduse,tech,ec,area] = MMSF[enduse,tech,ec,area]+
# max(MMSF[enduse,tech,ec,area]-MMSFB[enduse,tech,ec,area],0)*
# ETSwitch[tech,area] 
ETSwitch = J.diff_fast("ETSwitch", loc1, loc2; dimension_filters=frt, sec) # 

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
pop!(frt, :Tech)
push!(frt, :Year => string.(2008))
push!(frt, :Tech => ["OffRoad"], :Enduse => ["OffRoad"])
pop!(frt, :Tech)
MMSM0 = J.diff_fast("MMSM0", loc1, loc2; dimension_filters=frt, sec) # doesn't match
MMSM0 = J.diff_fast("MMSM0", spr_base, tan_base; dimension_filters=frt, sec) # doesn't match
@rsubset MMSM0 isnan(:Redwood)
MSMM = J.diff_fast("MSMM", loc1, loc2; dimension_filters=frt, sec) # matches
sum(isnan.(MSMM.Pine))
MVF = J.diff_fast("MVF", loc1, loc2; dimension_filters=frt, sec) # matches
sum(isnan.(MVF.Pine))

MCFU = J.diff_fast("MCFU", loc1, loc2; dimension_filters=frt, sec) # 1 problem
@rsubset! MCFU isnan(:Pine)
MMSMI = J.diff_fast("MMSMI", loc1, loc2; dimension_filters=frt, sec) # 0
@rsubset! MMSMI isnan(:Pine)
AMSF = J.diff_fast("AMSF", loc1, loc2; dimension_filters=frt, sec) # 57 problems
@rsubset! AMSF isnan(:Pine)
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters=frt, sec) # 0
sum(isnan.(PEE.Pine))
@rsubset PEE isnan(:Tanoak)
PEESw = J.diff_fast("PEESw", loc1, loc2; dimension_filters=frt, sec) # 1
sum(isnan.(MVF.PEESw))

            # @finite_math MAW[enduse,tech,ec,area] = exp(MMSMI[enduse,tech,ec,area]*
            #   (SPC[ec,area]/SPop[ec,area])/(SPC0[ec,area]/SPop0[ec,area])+
            #   MVF[enduse,tech,ec,area]*log((MCFU[enduse,tech,ec,area]/Inflation[area]/
            #   PEE[enduse,tech,ec,area])/(MCFU0[enduse,tech,ec,area]/Inflation0[area]/PEE0[enduse,tech,ec,area])))

            MMSMI = J.diff_fast("MMSMI", loc1, loc2; dimension_filters=frt, sec) # 0
MVF = J.diff_fast("MVF", loc1, loc2; dimension_filters=frt, sec) # matches
MCFU = J.diff_fast("MCFU", loc1, loc2; dimension_filters=frt, sec) # NaN
@rsubset MCFU isnan(:Pine)
MCFUPolicy = J.diff_fast("MCFUPolicy", loc1, loc2; dimension_filters=frt, sec) # 0
AMSF = J.diff_fast("AMSF", loc1, loc2; dimension_filters=frt, sec) # doesn't match
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters=frt, sec) # doesn't match

# MCFU=DCCR*DCC+DOMC+ECFP/DEE+IdrtCost*Inflation
DCCR = J.diff_fast("DCCR", loc1, loc2; dimension_filters=frt, sec) #
DCCRPolicy = J.diff_fast("DCCRPolicy", loc1, loc2; dimension_filters=frt, sec) # 0
@rsubset DCCR isnan(:Pine)
DCC = J.diff_fast("DCC", loc1, loc2; dimension_filters=frt, sec) # issue
xDCC = J.diff_fast("xDCC", loc1, loc2; dimension_filters=frt, sec) # issue
xDCCPolicy = J.diff_fast("xDCCPolicy", loc1, loc2; dimension_filters=frt, sec) # issue
@rsubset DCC isnan(:Pine)
ECFP = J.diff_fast("ECFP", loc1, loc2; dimension_filters=frt, sec) #
@rsubset ECFP isnan(:Pine)
DEE = J.diff_fast("DEE", loc1, loc2; dimension_filters=frt, sec) #
@rsubset DEE isnan(:Pine)
IdrtCost = J.diff_fast("IdrtCost", loc1, loc2; dimension_filters=frt, sec) #
@rsubset IdrtCost isnan(:Pine)
Inflation = J.diff_fast("Inflation", loc1, loc2; dimension_filters=frt, sec) #
@rsubset Inflation isnan(:Pine)

push!(frt, :Tech => "TrainElectric")
# DCC=DCCBefore*(1-DEEPolicyMSF)+(DCCPolicy-DCCSubsidy*Inflation)*DEEPolicyMSF
DCCFullCost = J.diff_fast("DCCFullCost", loc1, loc2; dimension_filters=frt, sec) # fine
DCCBefore = J.diff_fast("DCCBefore", loc1, loc2; dimension_filters=frt, sec) # fine
@rsubset DCCBefore isnan(:Pine)
DEEPolicyMSF = J.diff_fast("DEEPolicyMSF", loc1, loc2; dimension_filters=frt, sec) # fine
@rsubset DEEPolicyMSF isnan(:Pine)
DCCPolicy = J.diff_fast("DCCPolicy", loc1, loc2; dimension_filters=frt, sec) # fine
xDCCPolicy = J.diff_fast("xDCCPolicy", loc1, loc2; dimension_filters=frt, sec) # fine
@rsubset DCCPolicy isnan(:Pine)
DCCSubsidy = J.diff_fast("DCCSubsidy", loc1, loc2; dimension_filters=frt, sec) # fine
@rsubset DCCSubsidy isnan(:Pine)

DEESw = J.diff_fast("DEESw", loc1, loc2; dimension_filters=frt, sec) # fine
DCCPrice = J.diff_fast("DCCPrice", loc1, loc2; dimension_filters=frt, sec) # fine
DCCBeforeStd = J.diff_fast("DCCBeforeStd", loc1, loc2; dimension_filters=frt, sec) # fine
# DCCPrice=DCCN*DCMM*(DEM*DEMM/DEEPrice-1)**(1/DCTC)*
#         (1+STX)*(1-DGF)*Inflation    
DCCN = J.diff_fast("DCCN", loc1, loc2; dimension_filters=frt, sec) # fine
DCMM = J.diff_fast("DCMM", loc1, loc2; dimension_filters=frt, sec) # NaN
DEM = J.diff_fast("DEM", loc1, loc2; dimension_filters=frt, sec) # fine
DEMM = J.diff_fast("DEMM", loc1, loc2; dimension_filters=frt, sec) # fine
DEEPrice = J.diff_fast("DEEPrice", loc1, loc2; dimension_filters=frt, sec) # fine
DCTC = J.diff_fast("DCTC", loc1, loc2; dimension_filters=frt, sec) # fine
STX = J.diff_fast("STX", loc1, loc2; dimension_filters=frt, sec) # fine
DGF = J.diff_fast("DGF", loc1, loc2; dimension_filters=frt, sec) # fine

xDCMM = J.diff_fast("xDCMM", loc1, loc2; dimension_filters=frt, sec) # 1
xDCC = J.diff_fast("xDCC", loc1, loc2; dimension_filters=frt, sec) # 1
# DCMM=xDCMM*DCMMPrior/xDCMMPrior
push!(frt, :Year => "2007")
xDCMMPrior = J.diff_fast("xDCMM", loc1, loc2; dimension_filters=frt, sec) # 1
DCMMPrior = J.diff_fast("DCMM", loc1, loc2; dimension_filters=frt, sec) # 1
push!(frt, :Year => "2008")

Enduse = pop!(frt, :Enduse)
push!(frt, :Enduse => Enduse)
PEEBeforeStd = J.diff_fast("PEEBeforeStd", loc1, loc2; dimension_filters=frt, sec) # doesn't match
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
@rsubset! df :Enduse != "frtRoad"

# Sometime PEE is getting halved and sometimes not. 
# For EC(Food) Enduse(OthSub), Tech(Oil) get halved but Tech(frtRoad) does not.
# Let's findout why
food = copy(frt)
PEEBeforeStd_f = J.diff_fast("PEEBeforeStd", loc1, loc2; dimension_filters=food, sec) # doesn't match
PEEBeforeStd # Doesn't match between Oil and frtRoad in Spruce
push!(food, :EC => "Food", :Enduse => "OthSub", :Tech => ["Oil", "OffRoad"])
PEM = J.diff_fast("PEM", loc1, loc2; dimension_filters=food, sec) # doesn't have tech
PEMM = J.diff_fast("PEMM", loc1, loc2; dimension_filters=food, sec) # matches
PCCR = J.diff_fast("PCCR", loc1, loc2; dimension_filters=food, sec) # matches
PCCRB = J.diff_fast("PCCR", spr_base, tan_base; dimension_filters=food, sec) # matches
PFTC = J.diff_fast("PFTC", loc1, loc2; dimension_filters=food, sec) # matches

PEPM = J.diff_fast("PEPM", loc1, loc2; dimension_filters=food, sec) # 1
PFPN = J.diff_fast("PFPN", loc1, loc2; dimension_filters=food, sec) # matches

MCFU = J.diff_fast("MCFU", loc1, loc2; dimension_filters=food, sec) # 1
Inflation = J.diff_fast("Inflation", loc1, loc2; dimension_filters=food, sec) # matches



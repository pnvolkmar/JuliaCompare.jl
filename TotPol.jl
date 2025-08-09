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
loc1 = J.Loc_p(vars, raw"\\Pink\c\2020CanadaPineAqua\2020Model\Ref25", "Pine");
loc2 = J.Loc_j(vars_j, joinpath(raw"\\Pink\c\2020CanadaRedwoodAqua\2020Model\Ref25", "database.hdf5"), "Redwood");
red_base = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaRedwood\\2020Model\\Base\\database.hdf5", "Redwood");
pine_base = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaPine\\2020Model\\Base", "Pine");
red_ogref = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaRedwood\\2020Model\\OGRef\\database.hdf5", "Redwood");
pine_ogref = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaPine\\2020Model\\OGRef", "Pine");
red_aqua = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaRedwoodAqua\\2020Model\\Ref25\\database.hdf5", "RedAqua");
pine_aqua = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaPineAqua\\2020Model\\Ref25", "PineAqua");

################################################################################
# Initialize filters and sec for analysis ######################################
################################################################################
sec = 'T'
dimension_filters = Dict{Symbol,Any}()
push!(dimension_filters, :Area => Canada, :Year => string.(1986:2050))
ECC = M.ReadDisk(loc2.HDF5_path, "E2020DB/ECCKey")
push!(dimension_filters, :ECC => ECC[occursin.("Foreign", ECC) .== false])

pass = copy(dimension_filters)
push!(pass, :Year => string.([1986; 2023:2030]))
push!(pass, :EC => "Passenger",  :ECC => "Passenger")
push!(pass, :Poll => "CO2")
push!(pass, :Area => ["ON"]) # Non exhaustive
push!(pass, :AreaTOM => ["NS", "BC","ON"])
push!(pass, :Fuel => ["Biodiesel","Diesel","Gasoline"], :FuelEP => ["Biodiesel","Diesel","Gasoline"])
push!(pass, :Enduse => "Carriage")
push!(pass, :Tech => ["BusGasoline","BusElectric"])
################################################################################
# Analysis of variables: Overview ##############################################
################################################################################

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters, sec) # 
TotPol_aqua = J.diff_fast("TotPol", pine_aqua, red_aqua; dimension_filters, sec) # 
TotPol = leftjoin(TotPol,TotPol_aqua, on = [:ECC,:Poll,:Area,:Year], renamecols = "" => "Aqua")
@rsubset! TotPol abs(:Diff) >= 1
@rsubset! TotPol :Diff != :DiffAqua
@rsubset TotPol :ECC == "Freight" :Year == 2050 :Area  == "ON"
@rsubset TotPol :ECC == "Passenger" :Year == 2050 :Area  == "ON"
J.add_pdiff!(TotPol)
J.plot_diff(TotPol, dim="ECC"; title = "Differences in TotPol")

xMMSF = J.diff_fast("xMMSF", loc1, loc2; dimension_filters, sec) # 
xMMSF_aqua = J.diff_fast("xMMSF", pine_aqua, red_aqua; dimension_filters, sec) # 
xMMSF = leftjoin(xMMSF,xMMSF_aqua, on = [:Enduse,:Tech,:EC,:Area,:Year], renamecols = "" => "Aqua")
@rsubset xMMSF :EC == "Passenger" :Year == 2050 :Area  == "ON"
J.plot_diff(xMMSF, dim="ECC"; title = "Differences in xMMSF")

xMMSF.PineDiff = xMMSF.Pine - xMMSF.PineAquaAqua
@rsubset xMMSF abs(:PineDiff) > 0
xMMSF = J.diff_fast("xMMSF", loc1, pine_aqua; dimension_filters, sec) # 

MMSM0 = J.diff_fast("MMSM0", loc1, loc2; dimension_filters, sec) # 
MMSM0_aqua = J.diff_fast("MMSM0", pine_aqua, red_aqua; dimension_filters, sec) # 
MMSM0 = leftjoin(MMSM0,MMSM0_aqua, on = [:Enduse,:Tech,:EC,:Area,:Year], renamecols = "" => "Aqua")
@rsubset MMSM0 :EC == "Passenger" :Year == 2050 :Area  == "ON"
J.plot_diff(MMSM0, dim="ECC"; title = "Differences in xMMSF")

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters=pass, sec) # 
J.add_pdiff!(TotPol)
@rsubset! TotPol :Year > 2022
J.plot_diff(TotPol, dim="Poll")
J.plot_diff(TotPol, dim="Area")

    # TotPol[ecc,poll,area] = EnPol[ecc,poll,area]+NcPol[ecc,poll,area]+
    #   MEPol[ecc,poll,area]+VnPol[ecc,poll,area]+FlPol[ecc,poll,area]+
    #   FuPol[ecc,poll,area]+ORMEPol[ecc,poll,area]+SqPol[ecc,poll,area]
EnPol = J.diff_fast("EnPol", loc1, loc2; dimension_filters=pass, sec) # Problem
NcPol = J.diff_fast("NcPol", loc1, loc2; dimension_filters=pass, sec) # Problem
MEPol = J.diff_fast("MEPol", loc1, loc2; dimension_filters=pass, sec) # Problem
VnPol = J.diff_fast("VnPol", loc1, loc2; dimension_filters=pass, sec) # 
FlPol = J.diff_fast("FlPol", loc1, loc2; dimension_filters=pass, sec) # 
FuPol = J.diff_fast("FuPol", loc1, loc2; dimension_filters=pass, sec) # 
ORMEPol = J.diff_fast("ORMEPol", loc1, loc2; dimension_filters=pass, sec) # 
SqPol = J.diff_fast("SqPol", loc1, loc2; dimension_filters=pass, sec) # Problem

EuFPol = J.diff_fast("EuFPol", loc1, loc2; dimension_filters=pass, sec) # Problem
@rsubset! EuFPol abs(:Diff) >= 1
J.add_pdiff!(EuFPol)
@rsubset! EuFPol :Year == 2050
J.plot_diff(EuFPol, dim = "FuelEP")

Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters=pass, sec = sec) # prob
J.add_pdiff!(Polute)
@rsubset! Polute abs(:PDiff) != 0
J.plot_diff(Polute, dim = "Tech")

i = findall(vars.Variable .== "POCA" .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
# Polute[enduse,fuelep,ec,poll,area] = EuDem[enduse,fuelep,ec,area]*
#     POCA[enduse,fuelep,ec,poll,area]+xPolute[enduse,fuelep,ec,poll,area]

POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters=pass, sec) # 
@rsubset! POCA :Diff != 0
EuDem = J.diff_fast("EuDem", loc1, loc2; dimension_filters=pass, sec) # Problem
@rsubset! EuDem :Diff != 0
J.add_pdiff!(EuDem)

xPolute = J.diff_fast("xPolute", loc1, loc2; dimension_filters=pass, sec) # DNE in transportation
J.add_pdiff!(xPolute)
@rsubset! xPolute abs(:PDiff) != 0

# EuDemF[enduse,fuel,ecc,area] = sum(Dmd[enduse,tech,ec,area]*
# DmFrac[enduse,fuel,tech,ec,area] for tech in Techs)
Dmd = J.diff_fast("Dmd", loc1, loc2; dimension_filters=pass, sec) # problem
J.add_pdiff!(Dmd)
@rsubset! Dmd :Diff != 0
DmFrac = J.diff_fast("DmFrac", loc1, loc2; dimension_filters=pass, sec) # fine
@rsubset! DmFrac :Diff != 0
J.add_pdiff!(DmFrac)


# @rsubset! DmFrac :Tanoak != 0
# xDmFrac = J.diff_fast("xDmFrac", loc1, loc2; dimension_filters=zoom, sec) # 
# @rsubset! xDmFrac :Tanoak != 0
# DmFracMin = J.diff_fast("DmFracMin", loc1, loc2; dimension_filters=pass, sec) # 
# @rsubset! DmFracMin :Tanoak != 0
# @rsubset DmFracMin :Area == "QC" :Fuel == "Biodiesel"
# @rsubset DmFracMin :Area == "BC" :Fuel == "Biodiesel"
# DmFracMax = J.diff_fast("DmFracMax", loc1, loc2; dimension_filters=zoom, sec) # 
# @rsubset DmFracMax :Tanoak != 0 :Diff != 0

# @. Dmd = Dmd-EE
EE = J.diff_fast("EE", loc1, loc2; dimension_filters=pass, sec) # 
@rsubset EE :Diff != 0
# @. Dmd = Dmd-DSMEU
DSMEU = J.diff_fast("DSMEU", loc1, loc2; dimension_filters=pass, sec) # not going to worry about now
@rsubset  DSMEU :Diff != 0
# Dmd[enduse,tech,ec,area] = Dmd[enduse,tech,ec,area]+
#                            SqDmd[tech,ec,area]*SqEUTechMap[enduse,tech]
SqDmd = J.diff_fast("SqDmd", loc1, loc2; dimension_filters=pass, sec) #  not going to worry aobut now
@rsubset! SqDmd :Diff != 0
J.add_pdiff!(SqDmd)
SqEUTechMap = J.diff_fast("SqEUTechMap", loc1, loc2; dimension_filters=pass, sec) # 
#
# @finite_math Dmd[enduse,tech,ec,area] = DER[enduse,tech,ec,area]*
# UMS[enduse,tech,ec,area]*CERSM[enduse,ec,area]*
# CUF[enduse,tech,ec,area]*WCUF[ec,area]*RPEI[enduse,tech,ec,area]/1.0e6*
# (TSLoad[enduse,ec,area]*
# (DDay[enduse,area]/DDayNorm[enduse,area])^DDCoefficient[enduse,ec,area]+
# (1.0-TSLoad[enduse,ec,area]))
DER = J.diff_fast("DER", loc1, loc2; dimension_filters=pass, sec) # Problem
J.add_pdiff!(DER)
UMS = J.diff_fast("UMS", loc1, loc2; dimension_filters=pass, sec)
CERSM = J.diff_fast("CERSM", loc1, loc2; dimension_filters=pass, sec)
CUF = J.diff_fast("CUF", loc1, loc2; dimension_filters=pass, sec)
WCUF = J.diff_fast("WCUF", loc1, loc2; dimension_filters=pass, sec)
RPEI = J.diff_fast("RPEI", loc1, loc2; dimension_filters=pass, sec)
TSLoad = J.diff_fast("TSLoad", loc1, loc2; dimension_filters=pass, sec)
DDay = J.diff_fast("DDay", loc1, loc2; dimension_filters=pass, sec)
DDayNorm = J.diff_fast("DDayNorm", loc1, loc2; dimension_filters=pass, sec)
DDCoefficient = J.diff_fast("DDCoefficient", loc1, loc2; dimension_filters=pass, sec)
EEImpact = J.diff_fast("EEImpact", loc1, loc2; dimension_filters=pass, sec)
EESat = J.diff_fast("EESat", loc1, loc2; dimension_filters=pass, sec)
xDmd = J.diff_fast("xDmd", loc1, loc2; dimension_filters=pass, sec)
xDmdTrend = J.diff_fast("xDmdTrend", loc1, loc2; dimension_filters=pass, sec)

DERRRExo = J.diff_fast("DERRRExo", loc1, loc2; dimension_filters=pass, sec) # not going to worry about rn
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters=pass, sec)
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters=pass, sec)
i = findall(vars.Variable .∈ Ref(["DERV", "DERAV"]) .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
DERV = J.diff_fast("DERV", loc1, loc2; dimension_filters=pass, sec) # Problem
select!(DERV, Not([:Enduse, :Tech]))
J.add_pdiff!(DERV)
@rsubset DERV :Vintage == "Vintage 1"
StockAdjustment = J.diff_fast("StockAdjustment", loc1, loc2; dimension_filters=pass, sec) # fine
DERAV = J.diff_fast("DERAV", loc1, loc2; dimension_filters=pass, sec)
DERA   = J.diff_fast("DERA", loc1, loc2; dimension_filters=pass, sec) # problem
J.add_pdiff!(DERA)
@rsubset DERA isnan(:Tanoak) :Year == 2020
#   @. DERA = DERAPC+DERAP+DERAD+DERARC
DERAPC = J.diff_fast("DERAPC", loc1, loc2; dimension_filters=pass, sec) #  problem
J.add_pdiff!(DERAPC)
DERAP = J.diff_fast("DERAP", loc1, loc2; dimension_filters=pass, sec) # Fine
DERAD = J.diff_fast("DERAD", loc1, loc2; dimension_filters=pass, sec) # fine
J.add_pdiff!(DERAD)
DERARC = J.diff_fast("DERARC", loc1, loc2; dimension_filters=pass, sec) # prob

RetroSwExo = J.diff_fast("RetroSwExo", loc1, loc2; dimension_filters=pass, sec) # -99

#   @. @finite_math DERAPC = (PERAPC+PERADSt)/DEE
PERAPC = J.diff_fast("PERAPC", loc1, loc2; dimension_filters=pass, sec) # big problem: Missing in Tanoak
J.add_pdiff!(PERAPC)
PERADSt = J.diff_fast("PERADSt", loc1, loc2; dimension_filters=pass, sec) # fine
@rsubset! PERADSt :Diff != 0
DEE = J.diff_fast("DEE", loc1, loc2; dimension_filters=pass, sec) # fine
J.add_pdiff!(DEE)

# process additions from production capacity (PERAPC)
# @finite_math PERAPC[enduse,tech,ec,area] = EUPCAPC[enduse,tech,New,ec,area]*
#     DSt[enduse,ec,area]/PEE[enduse,tech,ec,area]
push!(pass, :Age => "New")
EUPCAPC = J.diff_fast("EUPCAPC", loc1, loc2; dimension_filters=pass, sec) # big problem
J.add_pdiff!(EUPCAPC)
DSt = J.diff_fast("DSt", loc1, loc2; dimension_filters=pass, sec) # 
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters=pass, sec) # 
J.add_pdiff!(PEE)
# EUPCAPC[enduse,tech,New,ec,area] = PCA[New,ecc,area]*MMSF[enduse,tech,ec,area]
PCA = J.diff_fast("PCA", loc1, loc2; dimension_filters=pass, sec) # 
age = pop!(pass, :Age)
PCA_mid = J.diff_fast("PCA", loc1, loc2; dimension_filters=pass, sec) # 
# PCA[New,ecc,area] = (PCPrior[ecc,area]*(ECGR[ecc,area]+GRPGR[area]))+PCR[Old,ecc,area]
push!(pass, :Age => "Old")
PC = J.diff_fast("PC", loc1, loc2; dimension_filters=pass, sec) # problem
J.add_pdiff!(PC)
ECGR = J.diff_fast("ECGR", loc1, loc2; dimension_filters=pass, sec) # problem
J.add_pdiff!(ECGR)
GRPGR = J.diff_fast("GRPGR", loc1, loc2; dimension_filters=pass, sec) # fine
PCR = J.diff_fast("PCR", loc1, loc2; dimension_filters=pass, sec) # problem
J.add_pdiff!(PCR)
#     @finite_math ECGR[ecc,area] = (PC[ecc,area]-PCPrior[ecc,area])/PCPrior[ecc,area]-GRPGR[area]
# PCMin = J.diff_fast("PCMin", loc1, loc2; dimension_filters=pass, sec) # temp var
xPC = J.diff_fast("xPC", loc1, loc2; dimension_filters=pass, sec) # problem
# @finite_math xPC[ecc,area] = Driver[ecc,area]/xECUF[ecc,area]   
Driver = J.diff_fast("Driver", loc1, loc2; dimension_filters=pass, sec) # 
J.add_pdiff!(Driver)
sort!(Driver, :PDiff)
xECUF = J.diff_fast("xECUF", loc1, loc2; dimension_filters=pass, sec) # 
# Driver[ecc,area] = Driver[ecc,area]*DriverMultiplier[ecc,area]
DriverMultiplier = J.diff_fast("DriverMultiplier", loc1, loc2; dimension_filters=pass, sec) # fine
xDriver = J.diff_fast("xDriver", loc1, loc2; dimension_filters=pass, sec) # issue
DrSwitch = J.diff_fast("DrSwitch", loc1, loc2; dimension_filters=pass, sec) # 21
xGO = J.diff_fast("xGO", loc1, loc2; dimension_filters=pass, sec) # issue
vGO = J.diff_fast("vGO", loc1, loc2; dimension_filters=pass, sec) # doesn't prvide data for CN
    # xGO[ecc,area,year]  =sum(xGOECC[ecc,areatom,year]*
    #   MapAreaTOM[area,areatom] for areatom in areatoms)

# xGOECC = J.diff_fast("xGOECC", loc1, loc2; dimension_filters=pass, sec) # temp var
MapECCfromTOM = J.diff_fast("KInput/MapECCfromTOM", loc1, loc2; dimension_filters=pass, sec='K') # 
@rsubset MapECCfromTOM :Pine != 0 || :Redwood != 0
push!(pass, :ECCfromTOM => ["Fertilizer"])
GY = J.diff_fast("KOutput/GY", pn_base, red_base; dimension_filters=pass, sec='K') # 
@rsubset GY :AreaTOM == "ON"
tech = pop!(pass, :Tech)
MMSF = J.diff_fast("MMSF", loc1, loc2; dimension_filters=pass, sec) # 
J.add_pdiff!(MMSF)
@rsubset! MMSF :Spruce != 0 || :Tanoak != 0
xProcSw = J.diff_fast("xProcSw", loc1, loc2; dimension_filters=pass, sec) # 
@rsubset xProcSw :PI == "MShare"
xDmd = J.diff_fast("xDmd", loc1, loc2; dimension_filters=pass, sec) # 
xMMSF = J.diff_fast("xMMSF", loc1, loc2; dimension_filters=pass, sec) # 
xMMSFB = J.diff_fast("xMMSF", pine_base, red_base; dimension_filters=pass, sec) # 
J.add_pdiff!(xMMSFB)
xDmdB = J.diff_fast("xDmd", pine_base, red_base; dimension_filters=pass, sec) # 
J.add_pdiff!(xDmdB)

DmdFuel = J.diff_fast("DmdFuel", pine_base, red_base; dimension_filters=pass, sec) # IInput2
@by DmdFuel by = [:Tech, :EC, :Year] :Pine = sum(:Pine) :Redwood = sum(:Redwood)
vDmd = J.diff_fast("vDmd", pine_base, red_base; dimension_filters=pass, sec) #
@rsubset vDmd :Year == 2015
vTrDmd = J.diff_fast("vTrDmd", pine_base, red_base; dimension_filters=pass, sec) #
@rsubset vTrDmd :Year == 2015
push!(pass, :TechTrans => pass[:Tech])
push!(pass, :ECTrans => pass[:EC])

# MMSF[enduse,tech,ec,area] = MMSF[enduse,tech,ec,area]+
# max(MMSF[enduse,tech,ec,area]-MMSFB[enduse,tech,ec,area],0)*
# ETSwitch[tech,area] 
ETSwitch = J.diff_fast("ETSwitch", loc1, loc2; dimension_filters=pass, sec) # 

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
pop!(pass, :Tech)
push!(pass, :Year => string.(2022:2023))
push!(pass, :Tech => ["OffRoad"], :Enduse => ["OffRoad"])
pop!(pass, :Tech)
CFraction = J.diff_fast("CFraction", loc1, loc2; dimension_filters=pass, sec) # doesn't match
MMSM0 = J.diff_fast("MMSM0", loc1, loc2; dimension_filters=pass, sec) # doesn't match
MMSM0 = J.diff_fast("MMSM0", spr_base, tan_base; dimension_filters=pass, sec) # doesn't match
@rsubset MMSM0 :Spruce > -170 || :Tanoak > -170
MSMM = J.diff_fast("MSMM", loc1, loc2; dimension_filters=pass, sec) # matches
MVF = J.diff_fast("MVF", loc1, loc2; dimension_filters=pass, sec) # matches
MCFU = J.diff_fast("MCFU", loc1, loc2; dimension_filters=pass, sec) # close
MMSMI = J.diff_fast("MMSMI", loc1, loc2; dimension_filters=pass, sec) # 0
AMSF = J.diff_fast("AMSF", loc1, loc2; dimension_filters=pass, sec) # doesn't match
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters=pass, sec) # doesn't match
J.add_pdiff!(PEE)
PEESw = J.diff_fast("PEESw", loc1, loc2; dimension_filters=pass, sec) # 1

            # @finite_math MAW[enduse,tech,ec,area] = exp(MMSMI[enduse,tech,ec,area]*
            #   (SPC[ec,area]/SPop[ec,area])/(SPC0[ec,area]/SPop0[ec,area])+
            #   MVF[enduse,tech,ec,area]*log((MCFU[enduse,tech,ec,area]/Inflation[area]/
            #   PEE[enduse,tech,ec,area])/(MCFU0[enduse,tech,ec,area]/Inflation0[area]/PEE0[enduse,tech,ec,area])))

            MMSMI = J.diff_fast("MMSMI", loc1, loc2; dimension_filters=pass, sec) # 0
MVF = J.diff_fast("MVF", loc1, loc2; dimension_filters=pass, sec) # matches
MCFU = J.diff_fast("MCFU", loc1, loc2; dimension_filters=pass, sec) # close
AMSF = J.diff_fast("AMSF", loc1, loc2; dimension_filters=pass, sec) # doesn't match
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters=pass, sec) # doesn't match

# MAW=exp(ln(MSMM)+MVF*ln((MCFU/xInflation/PEE)/(MCFU0/Inflation0/PEE0)))
# TAW=sum(Tech)(MAW(Tech))
# MU=MMSFx/MAW
# MMSM0(TE,Y)=LN(MU(TE)/max(Tech)(MU(Tech)))    
YMMSM = J.diff_fast("YMMSM", loc1, loc2; dimension_filters=pass, sec) # doesn't match
EUPC = J.diff_fast("EUPC", loc1, loc2; dimension_filters=pass, sec) # doesn't match
J.add_pdiff!(EUPC)

# EUPC=EUPCPrior+DT*(EUPCA-EUPCR)
# EUPCAdj=EUPC*StockAdjustment
# EUPC=EUPC+EUPCAdj
push!(pass, :Year => string.(2023:2024))
StockAdjustment = J.diff_fast("StockAdjustment", loc1, loc2; dimension_filters=pass, sec) # doesn't match
DT = J.diff_fast("DT", loc1, loc2; dimension_filters=pass, sec) # doesn't match
EUPCA = J.diff_fast("EUPCA", loc1, loc2; dimension_filters=pass, sec) # doesn't match
EUPCR = J.diff_fast("EUPCR", loc1, loc2; dimension_filters=pass, sec) # doesn't match
J.add_pdiff!(EUPC)
J.add_pdiff!(StockAdjustment)
J.add_pdiff!(EUPCA)
J.add_pdiff!(EUPCR)

# EUPCA=EUPCAPC+EUPCAC
# EUPCR=EUPCRPC+EUPCRC
EUPCAPC = J.diff_fast("EUPCAPC", loc1, loc2; dimension_filters=pass, sec) # doesn't match
EUPCAC = J.diff_fast("EUPCAC", loc1, loc2; dimension_filters=pass, sec) # doesn't match
EUPCRPC = J.diff_fast("EUPCRPC", loc1, loc2; dimension_filters=pass, sec) # doesn't match
EUPCRC = J.diff_fast("EUPCRC", loc1, loc2; dimension_filters=pass, sec) # doesn't match
J.add_pdiff!(EUPCR)


Enduse = pop!(pass, :Enduse)
push!(pass, :Enduse => Enduse)
PEEBeforeStd = J.diff_fast("PEEBeforeStd", loc1, loc2; dimension_filters=pass, sec) # doesn't match
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
food = copy(pass)
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



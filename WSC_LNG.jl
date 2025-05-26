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
SCENARIO1 = "Ref24"
SCENARIO2 = "Ref24"
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
tan_ogref = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaTanoak\\2020Model\\OGRef\\database.hdf5", "Tan_OGRef");

################################################################################
# Initialize variables for analysis ############################################
################################################################################
sec = 'S'
dimension_filters = Dict{Symbol,Any}()

push!(dimension_filters, :Area => ["WSC"])
push!(dimension_filters, :Process => ["LNGProduction"])
push!(dimension_filters, :ECC => ["LNGProduction"])
push!(dimension_filters, :EC => ["LNGProduction"])
push!(dimension_filters, :Year => string.(2020:2021))
push!(dimension_filters, :OGUnit => ["WSC_LNG_0001"])
################################################################################
# Analysis of variables ########################################################
################################################################################
GAProd = J.diff_fast("GAProd", loc1, loc2; dimension_filters, sec) # bad
GProd = J.diff_fast("GProd", loc1, loc2; dimension_filters, sec) # bad
@rsubset! GAProd isnan(:Tanoak)
@by GAProd [:Process, :Area] :Count = length(:Diff)
xGAProd = J.diff_fast("xGAProd", loc1, loc2; dimension_filters, sec) # good 
GOMult = J.diff_fast("GOMult", loc1, loc2; dimension_filters, sec) # good
Pd = J.diff_fast("Pd", loc1, loc2; dimension_filters, sec) # 0 which isn't right
LNGAProd = J.diff_fast("LNGAProd", loc1, loc2; dimension_filters, sec) # 0 which isn't right
@rsubset! Pd isnan(:Tanoak)
xPd = J.diff_fast("xPd", loc1, loc2; dimension_filters, sec) # good
#       Pd[ogunit] = min(RsDevPrior[ogunit]*PdRate[ogunit],PdMax[ogunit])
RsDev = J.diff_fast("RsDev", loc1, loc2; dimension_filters, sec) # good for 2020 and we're using prior bad 2021
PdRate = J.diff_fast("PdRate", loc1, loc2; dimension_filters, sec) # bad
PdMax = J.diff_fast("PdMax", loc1, loc2; dimension_filters, sec) # good
# PdRate[ogunit] = xPdRate[ogunit]*PdRateM[ogunit]
xPdRate = J.diff_fast("xPdRate", loc1, loc2; dimension_filters, sec) # good
PdRateM = J.diff_fast("PdRateM", loc1, loc2; dimension_filters, sec) # good
PdRateMRef = J.diff_fast("PdRateM", loc1, tan_ogref; dimension_filters, sec) # good
PdSw = J.diff_fast("PdSw", loc1, loc2; dimension_filters, sec) # good
OAPrEOR = J.diff_fast("OAPrEOR", loc1, loc2; dimension_filters, sec) # 0 as expected

# GAProd[lngproduction,wsc] = Exports[ngfuelep,us]-ExportsPipeline[us]-LNGProdMin[us]
push!(dimension_filters, :FuelEP => "NaturalGas")
push!(dimension_filters, :Nation => "US")
Exports = J.diff_fast("Exports", loc1, loc2; dimension_filters, sec) # bad
ExportsPipeline = J.diff_fast("ExportsPipeline", loc1, loc2; dimension_filters, sec) # bad
LNGProdMin = J.diff_fast("LNGProdMin", loc1, loc2; dimension_filters, sec) # good
# Exports[ngfuelep,nation] = ExportsMin[ngfuelep,nation]+LNGProdMin[nation]+max(GMarket[nation]-
#   GDemand[nation]+Imports[ngfuelep,nation]-ExportsMin[ngfuelep,nation]-
#   LNGProdMin[nation]+SupplyAdjustments[ngfuelep,nation],0)
ExportsMin = J.diff_fast("ExportsMin", loc1, loc2; dimension_filters, sec) # good
GMarket = J.diff_fast("GMarket", loc1, loc2; dimension_filters, sec) # bad
GDemand = J.diff_fast("GDemand", loc1, loc2; dimension_filters, sec) # good
Imports = J.diff_fast("Imports", loc1, loc2; dimension_filters, sec) # bad
ExportsMin = J.diff_fast("ExportsMin", loc1, loc2; dimension_filters, sec) # good
SupplyAdjustments = J.diff_fast("SupplyAdjustments", loc1, loc2; dimension_filters, sec) # bad
#  GMarket[nation] = (TProd[nation]-GRaw[nation])*GMMult[nation]+(VnGProd[nation]+
#    FuGProd[nation]+FlGProd[nation])
pop!(dimension_filters, :Area)
GRaw = J.diff_fast("GRaw", loc1, loc2; dimension_filters, sec) # v small but good
GARaw = J.diff_fast("GARaw", loc1, loc2; dimension_filters, sec) # v small but good
GMMult = J.diff_fast("GMMult", loc1, loc2; dimension_filters, sec) # good
VnGProd = J.diff_fast("VnGProd", loc1, loc2; dimension_filters, sec) # good
FuGProd = J.diff_fast("FuGProd", loc1, loc2; dimension_filters, sec) # good
FlGProd = J.diff_fast("FlGProd", loc1, loc2; dimension_filters, sec) # good 
# TProd[nation] = sum(GProd[process,nation]*GasProductionMap[process] for process in Processes)
GProd = J.diff_fast("GProd", loc1, loc2; dimension_filters, sec) # good 
pop!(dimension_filters, :Process)
GasProductionMap = J.diff_fast("GasProductionMap", loc1, loc2; dimension_filters, sec) # good 
push!(dimension_filters, :Process => ["ConventionalGasProduction", "UnconventionalGasProduction", "AssociatedGasProduction"])
GProd = J.diff_fast("GProd", loc1, loc2; dimension_filters, sec) # good 
push!(dimension_filters, :Process => ["ConventionalGasProduction"])
pop!(dimension_filters, :Area)
GAProd = J.diff_fast("GAProd", loc1, loc2; dimension_filters, sec) # good 
@rsubset! GAProd isnan(:Tanoak)
pop!(dimension_filters, :Nation)
# @finite_math PdRateM[ogunit] = max(min((1+PdVar[ogunit])/(PdVar[ogunit]+
#   (PdROI[ogunit]/PdROIRef[ogunit])^PdVF[ogunit]),
#   PdMaxM[ogunit]),PdMinM[ogunit])
# PdRate[ogunit] = xPdRate[ogunit]*PdRateM[ogunit]
# Pd[ogunit] = min(RsDevPrior[ogunit]*PdRate[ogunit],PdMax[ogunit])
PdVar = J.diff_fast("PdVar", loc1, loc2; dimension_filters, sec) # good
PdROI = J.diff_fast("PdROI", loc1, loc2; dimension_filters, sec) # bad
PdROIRef = J.diff_fast("PdROI", loc1, tan_ogref; dimension_filters, sec) # bad
PdVF = J.diff_fast("PdVF", loc1, loc2; dimension_filters, sec) # discrepancies
PdMaxM = J.diff_fast("PdMaxM", loc1, loc2; dimension_filters, sec) # good
PdMinM = J.diff_fast("PdMinM", loc1, loc2; dimension_filters, sec) # good
#   @finite_math PdROI[ogunit] = PdNtInc[ogunit] /(SusCap[ogunit])
PdNtInc = J.diff_fast("PdNtInc", loc1, loc2; dimension_filters, sec) # bad
SusCap = J.diff_fast("SusCap", loc1, loc2; dimension_filters, sec) # bad
# PdNtInc[ogunit] = OGRev[ogunit]-PdExp[ogunit] -PdITax[ogunit]
OGRev = J.diff_fast("OGRev", loc1, loc2; dimension_filters, sec) # good
PdExp = J.diff_fast("PdExp", loc1, loc2; dimension_filters, sec) # bad
PdITax = J.diff_fast("PdITax", loc1, loc2; dimension_filters, sec) # bad
#   PdExp[ogunit] = OGOpCosts[ogunit]+OpWrkCap[ogunit]+OGAbCosts[ogunit]+RyLev[ogunit]+SusDep[ogunit]
OGOpCosts = J.diff_fast("OGOpCosts", loc1, loc2; dimension_filters, sec) # bad in 2022
OpWrkCap = J.diff_fast("OpWrkCap", loc1, loc2; dimension_filters, sec) # good
OGAbCosts = J.diff_fast("OGAbCosts", loc1, loc2; dimension_filters, sec) # bad
RyLev = J.diff_fast("RyLev", loc1, loc2; dimension_filters, sec) # bad
SusDep = J.diff_fast("SusDep", loc1, loc2; dimension_filters, sec) # bad
#   OGAbCosts[ogunit] = DevCap[ogunit]*OGAbCFr[ogunit]
DevCap = J.diff_fast("DevCap", loc1, loc2; dimension_filters, sec) # bad
OGAbCFr = J.diff_fast("OGAbCFr", loc1, loc2; dimension_filters, sec) # good
#  @finite_math DevCap[ogunit] = xDevCap[ogunit]*DevDM[ogunit]/DevDMRef[ogunit]*DevLCM[ogunit]*
# DevIM[ogunit]*CCMDem[ogunit]*InflationOGUnit[ogunit]
xDevCap = J.diff_fast("xDevCap", loc1, loc2; dimension_filters, sec) # good
DevDM = J.diff_fast("DevDM", loc1, loc2; dimension_filters, sec) # bad in 2022
DevDMRef = J.diff_fast("DevDM", loc1, tan_ogref; dimension_filters, sec) # good
DevLCM = J.diff_fast("DevLCM", loc1, loc2; dimension_filters, sec) # good
DevIM = J.diff_fast("DevIM", loc1, loc2; dimension_filters, sec) # good
CCMDem = J.diff_fast("CCMDem", loc1, loc2; dimension_filters, sec) # bad
InflationOGUnit = J.diff_fast("InflationOGUnit", loc1, loc2; dimension_filters, sec) # good
# if (CCMDemSw[ogunit] == 1.0) && (OGNation[ogunit] == "CN")
#   CCMDem[ogunit] = DemCCMultPrior[ecc,area]
# else
#   CCMDem[ogunit] = 1.0
# end
CCMDemSw = J.diff_fast("CCMDemSw", loc1, loc2; dimension_filters, sec) # 1
DemCCMult = J.diff_fast("DemCCMult", loc1, loc2; dimension_filters, sec) # bad
# from IDemand: 
#         @finite_math DemCCMult[ecc,area] = DemCC[ecc,area]/DemCCRef[ecc,area]
DemCC = J.diff_fast("DemCC", loc1, loc2; dimension_filters, sec) # bad
DemCCRef = J.diff_fast("DemCC", loc1, loc2; dimension_filters, sec) # bad which takes PCCRef which is PCC
# DemCC[ecc,area] = sum(DCC[enduse,tech,ec,area]/Inflation[area]*
#   PERRef[enduse,tech,ec,area] for tech in Techs, enduse in Enduses)
sec = 'I'
DCC = J.diff_fast("DCC", loc1, loc2; dimension_filters, sec) # good
@rsubset DCC isnan(:Tanoak)
Inflation = J.diff_fast("Inflation", loc1, loc2; dimension_filters, sec) # good
PERRef = J.diff_fast("PER", loc1, loc2; dimension_filters, sec) # good
@rsubset PERRef isnan(:Tanoak)
#         DemCC[ecc,area] = DemCC[ecc,area]+
#           sum(PCC[enduse,tech,ec,area]/Inflation[area]*
#           EUPCTemp[enduse,tech,ec,area] for tech in Techs, enduse in enduses)
PCC = J.diff_fast("PCC", loc1, loc2; dimension_filters, sec) # good
@rsubset PCC isnan(:Tanoak) # bad for enduse heat and techs LPG in NB and Gas in NS all the way back to 1986
PEESw = J.diff_fast("PEESw", loc1, loc2; dimension_filters, sec) # 1
push!(dimension_filters, :Enduse => "Heat")
push!(dimension_filters, :Tech => ["LPG", "Gas"])
PCCPrice = J.diff_fast("PCCPrice", loc1, loc2; dimension_filters, sec) # bad

# @finite_math PCCPrice[enduse,tech,ec,area] = PCCN[enduse,tech,ec,area]*
#   PCCMM[enduse,tech,ec,area]*Inflation[area]*
#   (1+STX[area])*(PEM[enduse,ec,area]*PEMM[enduse,tech,ec,area]/
#   PEEPrice[enduse,tech,ec,area]-1)^(1/PCTC[enduse,tech,ec,area])

PCCMM = J.diff_fast("PCCMM", loc1, loc2; dimension_filters, sec) # bad
Inflation = J.diff_fast("Inflation", loc1, loc2; dimension_filters, sec) # bad
STX = J.diff_fast("STX", loc1, loc2; dimension_filters, sec) # bad
PEM = J.diff_fast("PEM", loc1, loc2; dimension_filters, sec) # bad
PEMM = J.diff_fast("PEMM", loc1, loc2; dimension_filters, sec) # bad
PEEPrice = J.diff_fast("PEEPrice", loc1, loc2; dimension_filters, sec) # bad
PCCN = J.diff_fast("PCCN", loc1, loc2; dimension_filters, sec) # zeros
PCTC = J.diff_fast("PCTC", loc1, loc2; dimension_filters, sec) # zeros
@rsubset PCTC :Tanoak == 0


InitialDemandYear = J.diff_fast("InitialDemandYear", loc1, loc2; dimension_filters, sec) # 1985

push!(dimension_filters, :Year => "1985")
Driver = J.diff_fast("Driver", loc1, loc2; dimension_filters, sec) #  < 0.0001, meaning average coef as used

# Driver subset issues were likely causing the issues in PCCN

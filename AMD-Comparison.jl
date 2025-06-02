using Revise
import JuliaCompare as J
import PromulaDBA as P
import SmallModel as M 
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta, Plots

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
# Electric Utility ###############################################################
################################################################################
sec = 'E'
UtilityGen = Dict{Symbol,Any}()
PolConv = J.var("PolConv", loc1)
PollKey = J.var("PollKey", loc1)

push!(UtilityGen, :Area => Canada)
push!(UtilityGen, :ECC => "UtilityGen")
push!(UtilityGen, :Poll => PollKey.PollKey[PolConv.PolConv.>0])

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters = UtilityGen, sec)
J.add_pdiff!(TotPol)
@rsubset! TotPol :Diff != 0

leftjoin!(TotPol, PolConv, on = :Poll)
TotPol.Spruce = TotPol.Spruce .* TotPol.PolConv
TotPol.Tanoak = TotPol.Tanoak .* TotPol.PolConv
TotPol = @by TotPol [:Area,:Year] :Spruce = sum(:Spruce) :Tanoak = sum(:Tanoak)
TotPol.Diff = TotPol.Spruce - TotPol.Tanoak
J.plot_diff(TotPol; dim = "Area", title = "Utility Gen TotPol Differences (GHG Tonnes) by Area")

################################################################################
# Transportation ###############################################################
################################################################################
T = 'T'
Transportation = Dict{Symbol,Any}()
T_ECCs = M.ReadDisk(loc2.HDF5_path, "TInput/ECKey")
push!(Transportation, :Area => Canada)
push!(Transportation, :ECC => T_ECCs[occursin.("Foreign", T_ECCs) .== false])

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters = Transportation, sec = T)
J.add_pdiff!(TotPol)
@rsubset! TotPol abs(:PDiff) >= 0.01
leftjoin!(TotPol, PolConv, on = :Poll)
TotPol.Spruce = TotPol.Spruce .* TotPol.PolConv
TotPol.Tanoak = TotPol.Tanoak .* TotPol.PolConv
TotPol = @by TotPol [:ECC, :Area,:Year] :Spruce = sum(:Spruce) :Tanoak = sum(:Tanoak)
TotPol.Diff = TotPol.Spruce - TotPol.Tanoak

J.plot_diff(TotPol; dim = "ECC", title = "Transportation TotPol Differences (GHG Tonnes) by ECC")

################################################################################
# Buildings ###############################################################
################################################################################
Buildings = Dict{Symbol,Any}()
eccs = [M.ReadDisk(loc2.HDF5_path, "RInput/ECKey"); M.ReadDisk(loc2.HDF5_path, "CInput/ECKey")]
push!(Buildings, :Area => Canada)
push!(Buildings, :ECC => eccs)

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters = Buildings, sec = sec)
J.add_pdiff!(TotPol)
@rsubset! TotPol abs(:PDiff) >= 0.01
leftjoin!(TotPol, PolConv, on = :Poll)
TotPol.Spruce = TotPol.Spruce .* TotPol.PolConv
TotPol.Tanoak = TotPol.Tanoak .* TotPol.PolConv
TotPol = @by TotPol [:ECC, :Area,:Year] :Spruce = sum(:Spruce) :Tanoak = sum(:Tanoak)
TotPol.Diff = TotPol.Spruce - TotPol.Tanoak

J.plot_diff(TotPol; dim = "ECC", title = "Buildings TotPol Differences (GHG Tonnes) by ECC")

################################################################################
# Industry #####################################################################
################################################################################
Industry = Dict{Symbol,Any}()
eccs = M.ReadDisk(loc2.HDF5_path, "IInput/ECKey")
push!(Industry, :Area => Canada)
push!(Industry, :ECC => eccs)

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters = Industry, sec = sec)
J.add_pdiff!(TotPol)
@rsubset! TotPol abs(:PDiff) >= 0.01
leftjoin!(TotPol, PolConv, on = :Poll)
TotPol.Spruce = TotPol.Spruce .* TotPol.PolConv
TotPol.Tanoak = TotPol.Tanoak .* TotPol.PolConv
TotPol = @by TotPol [:ECC, :Area,:Year] :Spruce = sum(:Spruce) :Tanoak = sum(:Tanoak)
TotPol.Diff = TotPol.Spruce - TotPol.Tanoak
@rsubset! TotPol :Year > 2020
J.plot_diff(TotPol; dim = "ECC", title = "Industry TotPol Differences (GHG Tonnes) by ECC")

################################################################################
# Waste #####################################################################
################################################################################
Waste = Dict{Symbol,Any}()
eccs = J.var("ECCKey", loc1)
push!(Waste, :Area => Canada)
push!(Waste, :ECC => ["SolidWaste", "Wastewater", "Incineration"])

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters = Waste, sec = sec)
J.add_pdiff!(TotPol)
@rsubset! TotPol abs(:PDiff) >= 0.01
leftjoin!(TotPol, PolConv, on = :Poll)
TotPol.Spruce = TotPol.Spruce .* TotPol.PolConv
TotPol.Tanoak = TotPol.Tanoak .* TotPol.PolConv
TotPol = @by TotPol [:ECC, :Area,:Year] :Spruce = sum(:Spruce) :Tanoak = sum(:Tanoak)
TotPol.Diff = TotPol.Spruce - TotPol.Tanoak
@rsubset! TotPol :Year > 2020
J.plot_diff(TotPol; dim = "Area", title = "TotPol Differences (GHG Tonnes) for Waste ECCs by Area")

################################################################################
# Transportation Detail ########################################################
################################################################################

Dmd = J.diff_fast("Dmd", loc1, loc2; dimension_filters = Transportation, sec = T)
LDV_Gasoline = @rsubset Dmd :Tech == "LDVGasoline"
LDV_Gasoline = @by LDV_Gasoline :Year :Spruce = sum(:Spruce) :Tanoak = sum(:Tanoak) :Diff = sum(abs.(:Diff))

LDV_Diesel = @rsubset Dmd :Tech == "LDVDiesel"
LDV_Diesel = @by LDV_Diesel :Year :Spruce = sum(:Spruce) :Tanoak = sum(:Tanoak) :Diff = sum(abs.(:Diff))

p1 = plot(LDV_Gasoline.Year, [LDV_Gasoline.Spruce LDV_Gasoline.Tanoak], 
          labels=["Spruce" "Tanoak"], title="LDVGasoline Demand")
p2 = plot(LDV_Diesel.Year, [LDV_Diesel.Spruce LDV_Diesel.Tanoak], 
          labels=["Spruce" "Tanoak"], title="LDVDiesel Demand")

plot(p1, p2, layout=(2,1))


################################################################################
# Energy Demands ###############################################################
################################################################################
Energy = Dict{Symbol,Any}()
ECCs = M.ReadDisk(loc2.HDF5_path, "E2020DB/ECCKey")
push!(Energy, :Area => Canada)
push!(Energy, :ECC => ECCs[occursin.("Foreign", ECCs) .== false])
Aviation = copy(Energy)
push!(Aviation, :Fuel => "AviationGasoline", :FuelEP => "AviationGasoline")

Geothermal = copy(Energy)
push!(Geothermal, :Fuel => "Geothermal", :FuelEP => "Geothermal")

Electric = copy(Energy)
push!(Electric, :Fuel => "Electric", :FuelEP => "Electric")

EuFPol_av = J.diff_fast("EuFPol", loc1, loc2; dimension_filters = Aviation, sec = T)
@rsubset! EuFPol_av :Diff != 0
J.plot_diff(EuFPol_av; dim = "ECC", num = 2, title = "GHG Emissions differences for Aviation Gasoline by ECC")
J.plot_diff(EuFPol_av; dim = "Poll")
@rsubset EuFPol_av :Year == 2030 :Poll == "CO2" abs(:Diff) > 1
av_g = @by EuFPol_av :Year :Spruce = sum(:Spruce) :Tanoak = sum(:Tanoak) :Diff = sum(abs.(:Diff))
J.add_pdiff!(av_g)

Dmd_el = J.diff_fast("Dmd", loc1, loc2; dimension_filters = Electric, sec = 'I')
@rsubset! Dmd_el :Tech == "Electric" abs(:Diff) > 1e-3
J.plot_diff(Dmd_el; dim = "EC", title = "TotPol Differences (GHG Tonnes) for Waste ECCs by Area")
@by Dmd_el :Year :Spruce = sum(:Spruce)

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

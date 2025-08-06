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
SCENARIO2 = ""
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
# Electric Utility #############################################################
################################################################################
sec = 'S'
EFlows = Dict{Symbol,Any}()

push!(EFlows, :Area => Canada)
push!(EFlows, :Fuel => "NaturalGas", :FuelEP => "NaturalGas")
push!(EFlows, :Year => ["2022", "2023","2024"])

Exports = J.diff_fast("Exports", loc1, loc2; dimension_filters = EFlows, sec)
@rsubset! Exports :Diff != 0 :Nation == "CN"
J.plot_diff(Exports; dim = "FuelEP", title = "Exports Differences by Fuel")

Imports = J.diff_fast("Imports", loc1, loc2; dimension_filters = EFlows, sec)
@rsubset! Imports :Diff != 0 :Nation == "CN"
J.plot_diff(Imports; dim = "FuelEP", title = "Imports Differences by Fuel")
sort(Imports, :Diff)

# Imports=ImportsMin+xmax(GDemand-GMarket-ImportsMin+ExportsMin+LNGProdMin-SupplyAdjustments,0)
Imports = J.diff_fast("Imports", loc1, loc2; dimension_filters = EFlows, sec)
xImports = J.diff_fast("xImports", loc1, loc2; dimension_filters = EFlows, sec) # good
ImportsMin = J.diff_fast("ImportsMin", loc1, loc2; dimension_filters = EFlows, sec) # match
ExportsMin = J.diff_fast("ExportsMin", loc1, loc2; dimension_filters = EFlows, sec) # match
LNGProdMin = J.diff_fast("LNGProdMin", loc1, loc2; dimension_filters = EFlows, sec) # match
SupplyAdjustments = J.diff_fast("SupplyAdjustments", loc1, loc2; dimension_filters = EFlows, sec) # don't match throughout
xSupplyAdjustments = J.diff_fast("xSupplyAdjustments", loc1, loc2; dimension_filters = EFlows, sec) # don't match throughout
GDemand = J.diff_fast("GDemand", loc1, loc2; dimension_filters = EFlows, sec) # missing?
GMarket = J.diff_fast("GMarket", loc1, loc2; dimension_filters = EFlows, sec) # matches
test = M.ReadDisk(loc2.HDF5_path, "SOutput/GDemand") # missing
sum(test)
test = M.ReadDisk(loc2.HDF5_path, "SpOutput/GADemand") # missing
sum(test)
ImportsMin.Tanoak[5]+max(GDemand.Tanoak[5]-GMarket.Tanoak[5]-ImportsMin.Tanoak[5]+ExportsMin.Tanoak[5]+LNGProdMin.Tanoak[5]-SupplyAdjustments.Tanoak[5],0)
Imports.Tanoak[5] # 2943.1 which is not right matches the value above. Assumes GDemand is 0. 
ImportsMin.Tanoak[6]+max(GDemand.Tanoak[6]-GMarket.Tanoak[6]-ImportsMin.Tanoak[6]+ExportsMin.Tanoak[6]+LNGProdMin.Tanoak[6]-SupplyAdjustments.Tanoak[6],0)
Imports.Tanoak[6] # 1372.4 which is not right matches the value above. Assumes GDemand is 0. 
# fuel = Select(Fuel,"NaturalGas")
#       GADemand[area] = sum(TotDemand[fuel,ecc,area] for ecc in ECCs)
TotDemand = J.diff_fast("TotDemand", loc1, loc2; dimension_filters = EFlows, sec) # good
@rsubset! TotDemand abs(:Diff) >= 0.001
J.plot_diff(TotDemand; dim = "Area", title = "Total Demand Differences by ECC")

# GMarket # NaNs throughout
# SupplyAdjustments # NaNs in 2022

    # GMarket[nation] = (TProd[nation]-GRaw[nation])*GMMult[nation]+
    #                   (VnGProd[nation]+FuGProd[nation]+FlGProd[nation])
      TProd = zeros(Float32,length(Nation))

  for nation in Nations
    TProd[nation] = sum(GProd[process,nation]*GasProductionMap[process] for process in Processes)
  end
  GProd = J.diff_fast("GProd", loc1, loc2; dimension_filters = EFlows, sec)
  @rsubset GProd isnan(:Tanoak) # Conventional Gas Production
  PCCN = J.diff_fast("PCCN", loc1, loc2; dimension_filters = EFlows, sec='I')
  @rsubset PCCN isnan(:Tanoak)
  dimension_filters = EFlows

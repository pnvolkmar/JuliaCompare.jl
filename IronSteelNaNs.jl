using Revise
import JuliaCompare as J
import PromulaDBA as P
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

BASE_FOLDER = raw"\\Silver\c\2020CanadaWalnut"
BASE_FOLDER2 = raw"\\Silver\c\24.11.15 2020CanadaWalnut"
SCENARIO1 = "StartBase"
SCENARIO2 = "StartBase"

CODE_FOLDER = joinpath(BASE_FOLDER, "Engine")
DATA_FOLDER1 = joinpath(BASE_FOLDER, "2020Model", SCENARIO1)
DATA_FOLDER2 = joinpath(BASE_FOLDER2, "2020Model", SCENARIO2)
E2020_Folder = joinpath(BASE_FOLDER, "2020Model")

vars = J.list_vars(CODE_FOLDER, DATA_FOLDER1, db_files);
loc1 = J.Loc_p(vars, DATA_FOLDER1, "Updated");
loc2 = J.Loc_p(vars, DATA_FOLDER2, "Official");

nm = "IOutput/DERAPC"
fname = joinpath(DATA_FOLDER1, "IOutput.dba")
vname = string("DERAPC")
out = P.dataframe(fname, vname; skip_zeros=false)

df = J.diff("PInv", loc1, loc2)
@rsubset! df :Diff != 0
J.plot_sets(df; dim="ECC", num=3)
@by(df, :ECC, :Count = sum(isreal.(:Diff)), :NACount = sum(isnan.(:Diff)))



@rsubset df3 :Year == 2050 :Area == "AB" :ECC == "Passenger"

# Row │ Fuel        ECC        Poll    Area    Year   Promula   Julia     Diff    
# ────┼───────────────────────────────────────────────────────────────────────────
#   1 │ Lubricants  Passenger  CO2     AB       2050   3801.08       0.0  3801.08

df = J.var("TOutput/FsPol", loc3) # No difference in lubricants co2
df2 = J.var("TOutput/FsPol", loc2)
@rsubset! df :Year > 2020 :Year < 2051 :EC ∈ ["Passenger", "Freight"] :Fuel == "Lubricants" :Poll == "CO2"
@rsubset! df2 :Year > 2020 :Year < 2051 :EC ∈ ["Passenger", "Freight"] :Fuel == "Lubricants" :Poll == "CO2"
df3 = J.diff(df, df2; name1="Promula", name2="Julia")
@rsubset df3 :Year == 2050 :Area == "AB" :EC == "Passenger"
# julia> @rsubset df3 :Year == 2050 :Area == "AB" :EC == "Passenger" :Promula != 0
# 0×9 DataFrame
#  Row │ Fuel    Tech    EC      Poll    Area    Year   Promula   Julia     Diff    
#      │ String  String  String  String  String  Int64  Float64?  Float64?  Float64 
# ─────┴────────────────────────────────────────────────────────────────────────────
# Nothing to see here!
#
J.plot_diff(df3; dim="Poll", num=3)
J.plot_diff(df3; dim="Area", num=3)
J.plot_diff(df3; dim="H2Tech", num=3)

totpol_j = J.var("TotPol", loc1)
totpol_p = J.var("SOutput/TotPol", loc2)
totpol = J.diff(totpol_j, totpol_p; name1="Julia", name2="Promula")
@rsubset! totpol :Diff != 0 :Year == 2050
@by(totpol, by = ["ECC"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))

polconv = J.var("PolConv", loc1)
ECCKey = J.var("ECCKey", loc1)

EnPol = J.diff("SOutput/EnPol", loc1, loc2; name1="Promula", name2="Julia")
@rsubset! EnPol :Year == 2050 :ECC == "Steam"
@by(EnPol, by = ["ECC"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))

MEPol = J.diff("SOutput/MEPol", loc1, loc2; name1="Promula", name2="Julia")
@rsubset! MEPol :Year == 2050 :ECC == "Fertilizer" :Area ∈ Canada
@by(MEPol, by = ["ECC"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
J.plot_diff(MEPol; dim="Area", num=3) # AB, ON, SK
J.plot_diff(MEPol; dim="Poll", num=1) # N2O

MEReduce = J.diff("SOutput/MEReduce", loc1, loc2; name1="Promula", name2="Julia")
@rsubset! MEReduce :Year == 2050 :ECC == "Fertilizer" :Area ∈ Canada
@by(MEReduce, by = ["ECC"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
J.plot_diff(MEReduce; dim="Area", num=3) # AB, ON, SK
J.plot_diff(MEReduce; dim="Poll", num=1) # N2O

MERP = J.diff("MEOutput/MERP", loc1, loc2; name1="Promula", name2="Julia")
@rsubset! MERP :ECC == "Fertilizer" :Area ∈ Canada :Year >= 2020
@by(MERP, by = ["ECC"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
J.plot_diff(MERP; dim="Area", num=10) # Most Areas
J.plot_diff(MERP; dim="Poll", num=5) # N2O and HFC
@by(MERP, by = ["Year"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
@rsubset MERP :Area == "AB" :Poll == "N2O"

MEC0 = J.diff("MEInput/MEC0", loc1, loc2; name1="Promula", name2="Julia")
@rsubset! MEC0 :ECC == "Fertilizer" :Area ∈ Canada :Year >= 2020
@by(MEC0, by = ["ECC"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
J.plot_diff(MEC0; dim="Area", num=10) # Most Areas
J.plot_diff(MEC0; dim="Poll", num=5) # N2O and HFC
@by(MEC0, by = ["Year"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
@rsubset MEC0 :Area == "AB" :Poll == "N2O"

MEPrice = J.diff("MEOutput/MEPrice", loc1, loc2; name1="Promula", name2="Julia")
@rsubset! MEPrice :ECC == "Fertilizer" :Area ∈ Canada :Year >= 2020
@by(MEPrice, by = ["ECC"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
J.plot_diff(MEPrice; dim="Area", num=10) # Most Areas
J.plot_diff(MEPrice; dim="Poll", num=5) # N2O and HFC
@by(MEPrice, by = ["Year"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
@rsubset MEPrice :Area == "AB" :Poll ∈ ["N2O", "HFC"]

MEPriceSw = J.diff("MEInput/MEPriceSw", loc1, loc2; name1="Promula", name2="Julia")
@rsubset! MEPriceSw :ECC == "Fertilizer" :Area ∈ Canada :Year >= 2020
@by(MEPriceSw, by = ["ECC"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
J.plot_diff(MEPriceSw; dim="Area", num=10) # Most Areas
J.plot_diff(MEPriceSw; dim="Poll", num=5) # N2O and HFC
@by(MEPriceSw, by = ["Year"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
@rsubset MEPriceSw :Area == "AB" :Poll ∈ ["N2O", "HFC"]

PolConv = J.diff("SInput/PolConv", loc1, loc2; name1="Promula", name2="Julia")
sum(abs.(PolConv.Diff))
PCost = J.diff("SOutput/PCost", loc1, loc2; name1="Promula", name2="Julia")
@rsubset! PCost :ECC == "Fertilizer" :Area ∈ Canada :Year >= 2020
@by(PCost, by = ["ECC"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
J.plot_diff(PCost; dim="Area", num=10) # Most Areas
J.plot_diff(PCost; dim="Poll", num=5) # N2O and HFC
@by(PCost, by = ["Year"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
@rsubset PCost :Area == "AB" :Poll ∈ ["N2O", "HFC"]

ECoverage = J.diff("SInput/ECoverage", loc1, loc2; name1="Promula", name2="Julia")
@rsubset! ECoverage :ECC == "Fertilizer" :Area ∈ Canada :Year >= 2020
@by(ECoverage, by = ["ECC"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
J.plot_diff(ECoverage; dim="Area", num=10) # Most Areas
J.plot_diff(ECoverage; dim="Poll", num=5) # N2O and HFC
@by(ECoverage, by = ["Year"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
@rsubset ECoverage :Area == "AB" :Poll ∈ ["N2O", "HFC"] :PCov == "Process"

PCostExo = J.diff("SInput/PCostExo", loc1, loc2; name1="Promula", name2="Julia")
@rsubset! PCostExo :ECC == "Fertilizer" :Area ∈ Canada :Year >= 2020
sum(abs.(PCostExo.Diff))
@by(PCostExo, by = ["ECC"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
J.plot_diff(PCostExo; dim="Area", num=10) # Most Areas
J.plot_diff(PCostExo; dim="Poll", num=5) # N2O and HFC
@by(PCostExo, by = ["Year"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
@rsubset PCostExo :Area == "AB" :Poll ∈ ["N2O", "HFC"]

FPCFSObligated = J.diff("SOutput/FPCFSObligated", loc1, loc2; name1="Promula", name2="Julia")
@rsubset! FPCFSObligated :ECC == "Fertilizer" :Area ∈ Canada :Year >= 2020
sum(abs.(FPCFSObligated.Diff))
@by(FPCFSObligated, by = ["ECC"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
J.plot_diff(FPCFSObligated; dim="Area", num=10) # Most Areas
J.plot_diff(FPCFSObligated; dim="Poll", num=5) # N2O and HFC
@by(FPCFSObligated, by = ["Year"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
@rsubset FPCFSObligated :Area == "AB" :Poll ∈ ["N2O", "HFC"]

MEIRP = J.diff("MEOutput/MEIRP", loc1, loc2; name1="Promula", name2="Julia")
@rsubset! MEIRP :ECC == "Fertilizer" :Area == "AB" :Year == 2050
@by(MEIRP, by = ["ECC"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
J.plot_diff(MEIRP; dim="Area", num=10) # Most Areas
J.plot_diff(MEIRP; dim="Poll", num=5) # N2O and HFC
@by(MEIRP, by = ["Year"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))

MERPStd = J.diff("MERPStd", loc1, loc2; name1="Promula", name2="Julia")
@rsubset! MERPStd :ECC == "Fertilizer" :Area ∈ Canada :Year >= 2020
@by(MERPStd, by = ["ECC"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
J.plot_diff(MERPStd; dim="Area", num=10) # Most Areas
J.plot_diff(MERPStd; dim="Poll", num=5) # N2O and HFC
@by(MERPStd, by = ["Year"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))

RPolicy = J.diff("SOutput/RPolicy", loc1, loc2; name1="Promula", name2="Julia")
@rsubset! RPolicy :ECC == "Fertilizer" :Area ∈ Canada :Year >= 2020
@by(RPolicy, by = ["ECC"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(:Diff))
J.plot_diff(RPolicy; dim="Area", num=10) # Most Areas
J.plot_diff(RPolicy; dim="Poll", num=5) # N2O and HFC
@by(RPolicy, by = ["Year"], :Julia = sum(:Julia), :Promula = sum(:Promula), :Diff = sum(abs(:Diff)))

MEIRP
PCC = J.var("IOutput/PCC", loc1)
df = @rsubset PCC isnan(:PCC)
@by(df, by = ["Tech", "EC", "Area"], :PCC = sum(:PCC))

PCCPrice = J.var("IOutput/PCCPrice", loc1)
df = @rsubset PCCPrice isnan(:PCCPrice)
@by(df, by = ["Tech", "EC", "Area"], :PCCPrice = sum(:PCCPrice))

#PCCPoll, PCTC, and PEEPrice don't have NaNs
PCTC = J.var("IOutput/PCTC", loc1)
df = @rsubset PCTC isnan(:PCTC)
PEEPrice = J.var("IOutput/PEEPrice", loc1)
df = @rsubset PEEPrice isnan(:PEEPrice)

J.filter_ns(PCTC)
J.filter_ne(PCTC)

J.var("IECFPFuel", loc1)
ECFP = J.diff("IOutput/ECFPFuel", loc1, loc2; name1="Ref24", name2="Ref24A")
@rsubset ECFP :Fuel ∈ ["Hydrogen", "NaturalGas"] :EC == "Petroleum" :Year == 2023 :Area == "AB"

J.var("FPCFS", loc1)
FPCFS = J.diff("IOutput/FPCFS", loc1, loc2; name1="Ref24", name2="Ref24A")
@rsubset FPCFS :Fuel ∈ ["Hydrogen", "NaturalGas"] :EC == "Petroleum" :Year == 2023 :Area == "AB"

J.var("DmFrac", loc1)
DmFrac = J.diff("IOutput/DmFrac", loc1, loc2; name1="Ref24", name2="Ref24A")
@rsubset DmFrac :Fuel ∈ ["Hydrogen", "NaturalGas"] :EC == "Petroleum" :Year == 2023 :Area == "AB"

J.var("DmFracMSF", loc1)
DmFracMSF = J.diff("IOutput/DmFracMSF", loc1, loc2; name1="Ref24", name2="Ref24A")
@rsubset DmFracMSF :Fuel ∈ ["Hydrogen", "NaturalGas"] :EC == "Petroleum" :Year == 2023 :Area == "AB"

undmd = J.var("undmd", loc1)
uncogen = J.var("uncogen", loc1)
unsector = J.var("unsector", loc1)
unarea = J.var("unarea", loc1)
unxsw = J.var("unxsw", loc1)

undmd = J.join_vars(undmd, uncogen, unsector, unarea)
undmd = J.join_vars(undmd, unxsw)

xundmd = J.var("xundmd", loc1)
df = J.join_vars(undmd, xundmd)
vundmd = J.var("vundmd", loc1)
df = J.join_vars(df, vundmd)
@rsubset df :UnArea == "NL" :UnSector == "HeavyOilMining" :Year == 2022 :UnCogen == 1 :FuelEP == "NaturalGasRaw"
@rsubset df :UnArea == "NL" :UnSector == "FrontierOilMining" :Year == 2022 :UnCogen == 1 :FuelEP == "NaturalGasRaw"
@rsubset df :UnArea == "ON" :UnSector == "Offices" :Year == 2022 :UnCogen == 1 :FuelEP == "RNG"
@rsubset @by(undmd, by = :Year, :UnXSw = sum(:UnXSw)) :Year > 2019
@rsubset df :UnArea == "NL" :Year == 2022 :UnCogen == 0 :FuelEP == "NaturalGasRaw"

@transform! df begin
  :xdiff = :xUnDmd - :UnDmd
  :vdiff = :vUnDmd - :UnDmd
end

df4 = @rsubset df :Year == 2022 :vdiff != 0
unique(df4.UnArea)

cgdmd = J.var("ioutput/cgdmd", loc1)
xcgdmd = J.var("iinput/xcgdmd", loc1)
J.var(vcgdmd, loc1)
df2 = J.join_vars(cgdmd, xcgdmd)
@rsubset df2 :Area == "NL" :EC == "HeavyOilMining" :Year == 2022
@rsubset df2 :Area == "NL" :EC == "FrontierOilMining" :Year == 2022

xeudmd = J.var("xeudmd", loc1)
veudmd = J.var("veudmd", loc1)
eudemand = J.var("eudemand", loc1)
@rsubset! eudemand :ECC == "UtilityGen"
select!(eudemand, Not(:ECC))
xeudemand = J.var("xeudemand", loc1)
@rsubset! xeudemand :ECC == "UtilityGen"
select!(xeudemand, Not(:ECC))
xeudmd.xEUDmd .*= 1054.61
veudmd.vEUDmd .*= 1054.61
eudemand.EuDemand .*= 1054.61
xeudemand.xEuDemand .*= 1054.61
df3 = J.join_vars(veudmd, xeudmd, xeudemand, eudemand)
@rsubset df3 :Year == 2022 :Area ∈ ["SK", "AB", "ON"] :Fuel == "NaturalGas"
@rsubset df3 :Year == 2022 :Area ∈ ["NU", "YT", "NL"] :Fuel == "Diesel"

xeudmd = J.var("xeudmd", loc2)
veudmd = J.var("veudmd", loc2)
eudemand = J.var("eudemand", loc2)
@rsubset! eudemand :ECC == "UtilityGen"
select!(eudemand, Not(:ECC))
xeudemand = J.var("xeudemand", loc2)
@rsubset! xeudemand :ECC == "UtilityGen"
select!(xeudemand, Not(:ECC))
xeudmd.xEUDmd .*= 1054.61
veudmd.vEUDmd .*= 1054.61
eudemand.EuDemand .*= 1054.61
xeudemand.xEuDemand .*= 1054.61
df3 = J.join_vars(veudmd, xeudmd, xeudemand, eudemand)
@rsubset df3 :Year == 2022 :Area ∈ ["SK", "AB", "ON"] :Fuel == "NaturalGas"
@rsubset df3 :Year == 2022 :Area ∈ ["NU", "YT", "NL"] :Fuel == "Diesel"


pol_w = J.var("TotPol", loc1)
pol_wa = J.var("TotPol", loc2)

pol_d = J.diff("TotPol", loc1, loc2; name1="Walnut", name2="WAspen")
J.plot_diff(pol_d; dim="ECC", num=10, title="Difference in TotPol for all Areas")
pol_cn = @rsubset pol_d :Area in Canada

J.plot_diff(pol_cn; dim="ECC", num=10, title="Difference in Canadian TotPol")
@by(pol_cn, by = "Year", :Walnut = sum(:Walnut), :WAspen = sum(:WAspen), :Diff = sum(:Diff))
@by(pol_d, by = "Year", :Walnut = sum(:Walnut), :WAspen = sum(:WAspen), :Diff = sum(:Diff))

totdemand = J.var("TotDemand", loc1)
@rsubset! totdemand :Area == "NU" :ECC == "OtherMetalMining"
J.plot_sets(totdemand, dim="Fuel")

eudemand = J.var("EuDemand", loc1)
@rsubset! eudemand :Area == "NU" :ECC == "OtherMetalMining"
J.plot_sets(eudemand, dim="Fuel")

cgdemand = J.var("cgdemand", loc1)
@rsubset! cgdemand :Area == "NU" :ECC == "OtherMetalMining"
J.plot_sets(cgdemand, dim="Fuel")

fsdemand = J.var("fsdemand", loc1)
@rsubset! fsdemand :Area == "NU" :ECC == "OtherMetalMining"
J.plot_sets(fsdemand, dim="Fuel")

cgec = J.var("cgec", loc1)
@rsubset! cgec :Area == "NU" :ECC == "OtherMetalMining"
eeconv = P.data(joinpath(DATA_FOLDER1, "SInput.dba"), "EEConv")
cgec.CgEC .= cgec.CgEC .* eeconv ./ 1e6
J.plot_sets(cgec, dim="ECC")

psoecc = J.var("PSoECC", loc1)
@rsubset! psoecc :Area == "NU" :ECC == "OtherMetalMining"
psoecc.PSoECC .= psoecc.PSoECC .* eeconv ./ 1e6
J.plot_sets(psoecc, dim="ECC")



totdemand.Year .+= 1984
df = @rsubset totdemand :Fuel == "RNG" :TotDemand != 0
df = @rsubset df :Year >= 2020 :Year <= 2030
J.plot_sets(df, dim="ECC")

DmFrac_Ash = J.var("IOutput/DMFrac", loc1)
@rsubset! DmFrac_Ash :Area .== "BC" :Year > 2023 - 1985 + 1 :Year < 2027 - 1985 + 1
@rsubset! DmFrac_Ash :Fuel .== "RNG"
@rsubset DmFrac_Ash :DmFrac != 0

ECCs = J.var("ECCDS", loc1)
Year = J.var("YrKey", loc1)

leftjoin!(DmFrac_Ash, Year, on="Year")
select!(DmFrac_Ash, Not("Year"))
dims = names(DmFrac_Ash)[1:end-1]
rename!(DmFrac_Ash, [dims; "Year"])

DmFrac_Ash.Year = "Y" .* DmFrac_Ash.Year
df = unstack(DmFrac_Ash, :Year, :DmFrac, combine=sum)
@rsubset df :Y2026 > 0

@rsubset! xGoalPol :xGoalPol != 0
@rsubset! xGoalPol :Market == "Market(161)"

plot(xGoalPol.Year, xGoalPol.xGoalPol)

xAdjust = J.var("xRPPAdjustArea", loc1)
vAdjust = J.var("vRPPAdjustmentsArea", loc1)
vStock = J.var("vRPPStockChangeArea", loc1)
vTransfer = J.var("vRPPTransfersArea", loc1)
df = J.var("xRPPProdArea", loc1)

rename!(vAdjust, Dict("vArea" => "Area"))
unique(vAdjust.Area)
unique(xAdjust.Area)
@rsubset! xAdjust :Area in canada
df = J.diff(xAdjust, vAdjust)
@rsubset df :Diff != 0
@rsubset df :old != 0

J.plot_diff(@rsubset df :Fuel == "Diesel"; dim="Area")


J.plot_sets(@rsubset xAdjust :Area == "AB"; dim="Fuel")
J.plot_sets(@rsubset vAdjust :vArea == "AB"; dim="Fuel")
J.plot_sets(vStock; dim="Fuel")
J.plot_sets(vTransfer; dim="Fuel")


totpol = J.var("TotPol", loc1)
@subset! totpol :Poll .== "CO2"
totpol0 = deepcopy(totpol)
totpol0.TotPol .= 0

df = J.diff(totpol, totpol0; name1="ref", name2="zero");
# diffs = J.comparedata(vars1, vars2, fnames)

anmap = J.var("anmap", loc1)
@subset! anmap :ANMap .== 1
select!(anmap, Not(:ANMap))

@subset! df :Year .>= 2020
leftjoin!(df, anmap, on=:Area)
@subset! df :Nation .== "CN"

classes = J.var("eccclmap", loc1)
@subset! classes :ECCCLMap .== 1
select!(classes, Not(:ECCCLMap))
df = leftjoin(df, classes, on=:ECC)
@subset! df :Class .== "Com"

J.plot_diff(@subset df :Area .== "AB"; dim="ECC")

EuDem = J.var("COutput/EuDem", loc1)
@rsubset! EuDem :Area in ["ON", "AB"]
@rsubset! EuDem :EC in ["Offices", "Health", "Education"]
@rsubset! EuDem :Year >= 2020
J.plot_sets(EuDem; dim="EC")
J.plot_sets(EuDem; dim="FuelEP")
J.plot_sets(EuDem; dim="Enduse")

EuDem = J.var("IOutput/EuDem", loc1)
@rsubset! EuDem :EC == "PulpPaperMills"
leftjoin!(EuDem, anmap, on=:Area)
@rsubset! EuDem :Nation == "CN"
select!(EuDem, Not(:Nation))

J.plot_sets(EuDem; dim="Area")
dropmissing!(EuDem)
dropmissing!(EuDem, :Year)

unega = J.var("ungc", loc1)
uncogen = J.var("uncogen", loc1)
unplant = J.var("unplant", loc1)

df = rightjoin(uncogen, unega, on=:Unit)
@rsubset! df :UnCogen == 1
df = rightjoin(unplant, df, on=:Unit)

# @rsubset! df :UnCogen == 1 :UnPlant == "BiomassCCS"

unsector = J.var("unsector", loc1)
df = rightjoin(unsector, df, on=:Unit)

@rsubset! df :Year > 2020
@rsubset! df :UnSector == "PulpPaperMills"
J.plot_sets(df; dim="UnPlant")

# 25.	Air Passenger - emissions higher in NZ than Ref23
df = J.diff("TotPol", Ref23, loc1; name1="Ref23", name2="NZ_CT")

# df = @by(df, [:ECC, :Area, :Year], :Ref23 = sum(:Ref23), :NZ_CT = sum(:NZ_CT))
@rsubset! df :Area in canada :Poll == "CO2"
df = @by(df, [:ECC, :Year], :Ref23 = sum(:Ref23), :NZ_CT = sum(:NZ_CT))
@rsubset! df :Year == 2050
@. df.Diff = df.Ref23 - df.NZ_CT
@rsubset df :Diff < 0
@rsubset df :ECC == "AirPassenger"

ref23 = J.var("TotPol", Ref23)
@rsubset! ref23 :Area in canada :Poll == "CO2" :Year == 2050
ref23 = @by(ref23, [:ECC], :TotPol = sum(:TotPol))
@rsubset ref23 :ECC == "AirPassenger"

ref23a = J.var("TotPol", loc2)
@rsubset! ref23a :Area in canada :Poll == "CO2" :Year == 2050
ref23a = @by(ref23a, [:ECC], :TotPol = sum(:TotPol))
@rsubset ref23a :ECC == "AirPassenger"

nz = J.var("totpol", loc1)
@rsubset! nz :Area in canada :Poll == "CO2" :Year == 2050
# nz = @by(nz, [:ECC], :TotPol = sum(:TotPol))
J.plot_sets(@rsubset nz :ECC == "AirPassenger"; dim="Area")
@rsubset ref23a :ECC == "AirPassenger"
@rsubset ref23 :ECC == "AirPassenger"

totdemand = J.var("totdemand", loc1)
df = @rsubset totdemand :Year >= 2020 :ECC == "AirPassenger" :Area in canada
J.plot_sets(df; dim="Fuel", num=10)

totdemand = J.var("totdemand", loc2)
df = @rsubset totdemand :Year >= 2020 :ECC == "AirPassenger" :Area in canada
J.plot_sets(df; dim="Fuel", num=10)

v = J.var("voaprod", loc1)
x = J.var("xoaprod", loc1)
o = J.var("oaprod", loc1)

dims = names(x)[1:end-1]
rename!(v, [dims; "OAProd"])
rename!(x, [dims; "OAProd"])
d = J.diff(v, x)
@rsubset d :Diff != 0 :Area in canada

unretire = J.var("UnRetire", Ref23)
uncogen = J.var("uncogen", Ref23)
unsector = J.var("unsector", Ref23)
unplant = J.var("unplant", Ref23)
unonline = J.var("unonline", Ref23)


@rsubset! unretire :Unit != ""
leftjoin!(unretire, uncogen, on=:Unit)
leftjoin!(unretire, unsector, on=:Unit)
leftjoin!(unretire, unplant, on=:Unit)
leftjoin!(unretire, unonline, on=:Unit)

@rsubset! unretire :UnCogen == 1 :UnSector == "PulpPaperMills"
@rsubset! unretire :UnPlant in ["BiomassCCS"; "Biomass"]
@rsubset! unretire :Year == 2030
@rsubset! unretire :UnRetire > 2023
@rsubset unretire :UnRetire != 2200

J.var("CD", Ref23)



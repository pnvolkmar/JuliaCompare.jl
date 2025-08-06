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
SCENARIO1 = ""
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
spr_base = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaSpruce\\2020Model\\Base", "Spr_Base");
loc2 = J.Loc_j(vars_j, HDF5_path, "Tanoak");
tan_base = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaTanoak\\2020Model\\Base\\database.hdf5", "Tan_Base");

################################################################################
# Initialize variables for analysis ############################################
################################################################################
sec = 'T'
dimension_filters = Dict{Symbol,Any}()
push!(dimension_filters, :Area => Canada)
push!(dimension_filters, :ECC => ["CommercialOffRoad","Passenger","ForeignPassenger"])
push!(dimension_filters, :EC => ["CommercialOffRoad","Passenger","ForeignPassenger"])
offroad = copy(dimension_filters)
push!(offroad, :ECC => ["CommercialOffRoad"])
push!(offroad, :EC => ["CommercialOffRoad"])
push!(offroad, :Area => "ON")

OMExp = J.diff_fast("OMExp", loc1, loc2; dimension_filters, sec)
@rsubset! OMExp :Diff != 0
J.plot_diff(OMExp; dim="Area", num=10, title="OMExp diffs by Area")
J.plot_diff(OMExp; dim="ECC", num=10, title="OMExp diffs by ECC")

DOMExp = J.diff_fast("DOMExp", loc1, loc2; dimension_filters=offroad, sec)
@rsubset! DOMExp :Diff != 0
J.plot_diff(DOMExp; dim="Area", num=10, title="OMExp diffs by Area")

CgOMExp = J.diff_fast("CgOMExp", loc1, loc2; dimension_filters, sec)
@rsubset! CgOMExp :Diff != 0
J.plot_diff(CgOMExp; dim="Area", num=10, title="OMExp diffs by Area")

      # DOMExp[ecc,area] = sum(DCCFullCost[enduse,tech,ec,area]*DOCF[enduse,tech,ec,area]*
      #   PER[enduse,tech,ec,area]*Inflation[area] for tech in Techs,enduse in Enduses)/1000000
DCCFullCost = J.diff_fast("DCCFullCost", loc1, loc2; dimension_filters=offroad, sec)
DOCF = J.diff_fast("DOCF", loc1, loc2; dimension_filters=offroad, sec)
PER = J.diff_fast("PER", loc1, loc2; dimension_filters=offroad, sec)
Inflation = J.diff_fast("Inflation", loc1, loc2; dimension_filters=offroad, sec)
@rsubset! DCCFullCost abs(:Diff) >= 1e-3
@rsubset! DOCF :Diff != 0
@rsubset! PER :Diff != 0
@rsubset! Inflation :Diff != 0

unique(DCCFullCost.Tech)
unique(DOCF.Tech)
# @finite_math DCCFullCost[enduse,tech,ec,area] = (DCC[enduse,tech,ec,area]+
#   DCCSubsidy[enduse,tech,ec,area]*Inflation[area]*
#   DEEPolicyMSF[enduse,tech,ec,area])*((1+STXB[area])/(1+STX[area]))/
#   (1-DGF[enduse,tech,ec,area])

      
sec = 'C'
push!(dimension_filters, :ECC => "Wholesale")
push!(dimension_filters, :EC => "Wholesale")
push!(dimension_filters, :Year => "2021")

Inflation = J.diff_fast("Inflation", loc1, loc2; dimension_filters, sec) # fine
@rsubset Inflation abs(:Diff) > 1e-4
CgDmd = J.diff_fast("CgDmd", loc1, loc2; dimension_filters, sec) # fine
@rsubset CgDmd abs(:Diff) > 1e-4
CgDC = J.diff_fast("CgDC", loc1, loc2; dimension_filters, sec) # fine
@rsubset CgDC abs(:Diff) > 1e-4
CgOF = J.diff_fast("CgOF", loc1, loc2; dimension_filters, sec) # fine
@rsubset CgOF abs(:Diff) > 1e-4
CgCC = J.diff_fast("CgCC", loc1, loc2; dimension_filters, sec) # fine
@rsubset CgCC abs(:Diff) > 1e-4
DCCFullCost = J.diff_fast("DCCFullCost", loc1, loc2; dimension_filters, sec) # fine
@rsubset DCCFullCost abs(:Diff) > 1e-4
DOCF = J.diff_fast("DOCF", loc1, loc2; dimension_filters, sec) # fine
@rsubset! DOCF abs(:Diff) > 1e-4
PER = J.diff_fast("PER", loc1, loc2; dimension_filters, sec) # fine
@rsubset! PER abs(:Diff) > 1e-4

@rsubset! EuFPol :ECC == "CommercialOffRoad" :Poll == "CO2" :Year >= 2020
EuFPol_g = @by EuFPol :Year begin
  :Spruce = sum(:Spruce)
  :Tanoak = sum(:Tanoak)
  :Diff = sum(abs.(:Diff))
end 
EuFPol_g.PDiff = EuFPol_g.Diff ./ EuFPol_g.Spruce

# Almost all Gasoline and Diesel
J.plot_diff(EuFPol; dim="FuelEP", num=10, title="EuFPol diffs by FuelEP") 
push!(dimension_filters, :FuelEP => ["Gasoline", "Diesel"])
push!(dimension_filters, :ECC => "CommercialOffRoad")
push!(dimension_filters, :EC => "CommercialOffRoad")

# CO2 only
@rsubset! EuFPol :Spruce != 0 && :Tanoak != 0
J.plot_diff(EuFPol; dim="Poll", num=10, title="EuFPol diffs by Poll") 
push!(dimension_filters, :Poll => "CO2")
# ON, AB, BC, QC are the biggest contributors, though many are present
J.plot_diff(EuFPol; dim="Area", num=10, title="EuFPol diffs by Area")
push!(dimension_filters, :Area => "AB")


Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters, sec)
J.plot_diff(Polute; dim="Tech", num=10, title="Polute diffs by Tech")
push!(dimension_filters, :Tech => ["OffRoad"])
@rsubset! Polute :Diff != 0 :Year > 2020
J.add_pdiff!(Polute)
@by(Polute, :Year, :Diff = sum(:Diff))
push!(dimension_filters, :Year => ["2022", "2023", "2024", "2029", "2030"])

POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters, sec)
i = findall(vars.Variable .== "POCA" .&& first.(vars.Database) .== 'T')
i = i[1]
vars[i,:Database] = "TOutput2"
POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters, sec)
@rsubset POCA :Diff != 0

DmdFEPTech = J.diff_fast("DmdFEPTech", loc1, loc2; dimension_filters, sec)
@rsubset DmdFEPTech !isapprox(:Diff, 0)

DmFrac = J.diff_fast("DmFrac", loc1, loc2; dimension_filters, sec)
@rsubset! DmFrac :Diff != 0
J.add_pdiff!(DmFrac)
@rsubset! DmFrac :PDiff != 0
@rsubset DmFrac :Fuel == "Ethanol"
Dmd = J.diff_fast("Dmd", loc1, loc2; dimension_filters, sec)
push!(dimension_filters, :Fuel => "Ethanol", :FuelEP => "Ethanol")
# DmFrac's Ethanol is causing the issue

DmFracMSM0 = J.diff_fast("DmFracMSM0", loc1, loc2; dimension_filters, sec)
DmFracMSF = J.diff_fast("DmFracMSF", loc1, loc2; dimension_filters, sec)
J.add_pdiff(DmFracMSF)
DmFracMarginal = J.diff_fast("DmFracMarginal", loc1, loc2; dimension_filters, sec)
J.add_pdiff(DmFracMarginal)
DmFracMin = J.diff_fast("DmFracMin", loc1, loc2; dimension_filters, sec)
J.add_pdiff(DmFracMin)
DmFracMax = J.diff_fast("DmFracMax", loc1, loc2; dimension_filters, sec)
J.add_pdiff(DmFracMax)
DmFracMax_b = J.diff_fast("DmFracMax", spr_base, tan_base; dimension_filters, sec)
J.add_pdiff(DmFracMax)

tech_check = copy(dimension_filters)
push!(tech_check, :Year => string.(2022:2030))
DmFracMax
sets = M.ReadSets(loc2.HDF5_path,"TOutput/DmFrac")
sets.Tech
DmFracMax = J.diff_fast("DmFracMax", loc1, loc2; dimension_filters=tech_check, sec)
J.add_pdiff(DmFracMax)

EEImpact = J.diff_fast("EEImpact", loc1, loc2; dimension_filters, sec)
EESat = J.diff_fast("EESat", loc1, loc2; dimension_filters, sec)
xDmd = J.diff_fast("xDmd", loc1, loc2; dimension_filters, sec)
DSMEU = J.diff_fast("DSMEU", loc1, loc2; dimension_filters, sec)
SqDmd = J.diff_fast("SqDmd", loc1, loc2; dimension_filters, sec)
SqEUTechMap = J.diff_fast("SqEUTechMap", loc1, loc2; dimension_filters, sec)
RPEI = J.diff_fast("RPEI", loc1, loc2; dimension_filters, sec)
CERSM = J.diff_fast("CERSM", loc1, loc2; dimension_filters, sec)
TSLoad = J.diff_fast("TSLoad", loc1, loc2; dimension_filters, sec)
UMS = J.diff_fast("UMS", loc1, loc2; dimension_filters, sec)
DDay = J.diff_fast("DDay", loc1, loc2; dimension_filters, sec)
CUF = J.diff_fast("CUF", loc1, loc2; dimension_filters, sec)
WCUF = J.diff_fast("WCUF", loc1, loc2; dimension_filters, sec)
DDCoefficient = J.diff_fast("DDCoefficient", loc1, loc2; dimension_filters, sec)
DDayNorm = J.diff_fast("DDayNorm", loc1, loc2; dimension_filters, sec)


DER = J.diff_fast("DER", loc1, loc2; dimension_filters, sec)
DERRRExo = J.diff_fast("DERRRExo", loc1, loc2; dimension_filters, sec)
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters, sec)
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters, sec)
DERV = J.diff_fast("DERV", loc1, loc2; dimension_filters, sec)
i = findall(vars.Variable .∈ Ref(["DERV", "DERAV"]) .&& first.(vars.Database) .== 'T')
i = i[1]
vars[i,:Database] = "TOutput2"
vars[i,:]
DERV = J.diff_fast("DERV", loc1, loc2; dimension_filters, sec)
StockAdjustment = J.diff_fast("StockAdjustment", loc1, loc2; dimension_filters, sec)
DERAV = J.diff_fast("DERAV", loc1, loc2; dimension_filters, sec)
dimension_filters[:Year] = string.(2020:2040)
DERA   = J.diff_fast("DERA", loc1, loc2; dimension_filters, sec = 'T')
DERAPC = J.diff_fast("DERAPC", loc1, loc2; dimension_filters, sec = 'T')
DERAP = J.diff_fast("DERAP", loc1, loc2; dimension_filters, sec = 'T')
DERAD = J.diff_fast("DERAD", loc1, loc2; dimension_filters, sec = 'T')
DERARC = J.diff_fast("DERARC", loc1, loc2; dimension_filters, sec = 'T')

@rsubset DERA   :Year == 2020
@rsubset DERAPC :Year == 2020
@rsubset DERAP  :Year == 2020
@rsubset DERAD  :Year == 2020
@rsubset DERARC :Year == 2020

# DERARC is the largest contributor to the difference in DERA in 2020

@rsubset DERA   :Year == 2039
@rsubset DERAPC :Year == 2039
@rsubset DERAP  :Year == 2039
@rsubset DERAD  :Year == 2039
@rsubset DERARC :Year == 2039

# DERAP is the largest contributor to the difference in DERA in 2039

df2020 = copy(dimension_filters)
df2020[:Year] = "2020"
df2020[:Year] = ["2019","2020"]

DERRRC = J.diff_fast("DERRRC", loc1, loc2; dimension_filters = df2020, T)
CMSF = J.diff_fast("CMSF", loc1, loc2; dimension_filters = df2020, T)
DEE = J.diff_fast("DEE", loc1, loc2; dimension_filters = df2020, T)
DEEA = J.diff_fast("DEEA", loc1, loc2; dimension_filters = df2020, T)
CnvrtEU = J.diff_fast("CnvrtEU", loc1, loc2; dimension_filters = df2020, T)
xCMSF = J.diff_fast("xCMSF", loc1, loc2; dimension_filters = df2020, T)
DEEA = J.diff_fast("DEEA", loc1, loc2; dimension_filters = df2020, T)
df2020[:Tech] = ["LDVGasoline","LDVDiesel"]
MMSF = J.diff_fast("MMSF", loc1, loc2; dimension_filters = df2020, T)
MMSM0 = J.diff_fast("MMSM0", loc1, loc2; dimension_filters, T)

dfDiesel = copy(dimension_filters)
dfDiesel[:Tech] = "LDVDiesel"
MMSM0 = J.diff_fast("MMSM0", loc1, loc2; dimension_filters, T)
# MMSF is off

CMSM0 = J.diff_fast("CMSM0", loc1, loc2; dimension_filters = df2020, T) # bad
@rsubset CMSM0 abs(:Diff) > 1e-5
CVF = J.diff_fast("CVF", loc1, loc2; dimension_filters = df2020, T) # good
MCFU = J.diff_fast("MCFU", loc1, loc2; dimension_filters = df2020, T) # good
DCCR = J.diff_fast("DCCR", loc1, loc2; dimension_filters = df2020, T) # good
FDCC = J.diff_fast("FDCC", loc1, loc2; dimension_filters = df2020, T) # good
FDCCU = J.diff_fast("FDCCU", loc1, loc2; dimension_filters = df2020, T) # good
MCFU0 = J.diff_fast("MCFU", loc1, loc2; dimension_filters = dfFirst, T) # good
CMSMI = J.diff_fast("CMSMI", loc1, loc2; dimension_filters = df2020, T) # good 
SPC = J.diff_fast("SPC", loc1, loc2; dimension_filters = df2020, T) # didn't work
df2020[:ECC] = "Passenger"
SPC = J.diff_fast("PC", loc1, loc2; dimension_filters = df2020, T) # didn't work
SPC0 = J.diff_fast("PC", loc1, loc2; dimension_filters = dfFirst, sec='M') # didn't work 
SPop = J.diff_fast("Pop", loc1, loc2; dimension_filters = df2020, sec='M') # didn't work 
SPop0 = J.diff_fast("Pop", loc1, loc2; dimension_filters = dfFirst, sec='M') # didn't work 
dfFirst = copy(df2020)
dfFirst[:Year] = "1986"

xMMSF = J.diff_fast("xMMSF", loc1, loc2; dimension_filters = df2020, sec) # didn't work
CMMSF = J.diff_fast("CMM", loc1, loc2; dimension_filters = df2020, sec) # didn't work

MMSM0 = J.diff_fast("MMSM0", loc1, loc2; dimension_filters = df2020, sec) # bad
@rsubset MMSM0 abs(:Diff) > 1e-5

dimension_filters[:Tech] = "BusGasoline"
MMSM0 = J.diff_fast("MMSM0", loc1, loc2; dimension_filters, sec = s) # bad
MMSM0_base = J.diff_fast("MMSM0", loc_s_base, loc_t_base; dimension_filters, sec = s) # bad
@rsubset MMSM0 abs(:Diff) > 1e-5

loc_s_base = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaSpruce\\2020Model\\Base", "SpruceBase")
loc_t_base = J.Loc_j(vars_j, HDF5_path, "TanoakBase")

MMSM0_base = J.diff_fast("MMSM0", loc_s_base, loc_t_base; dimension_filters = df2020, sec = s) # bad
push!(df2020, :Area => Canada)
push!(df2020, :Tech => ["BusGasoline", "LDVGasoline"])
@rsubset MMSM0_base abs(:Diff) > 1e-5
CUF_base = J.diff_fast("CUF", loc_s_base, loc_t_base; dimension_filters = df2020, sec = s) # bad
@rsubset CUF_base abs(:Diff) > 1e-5

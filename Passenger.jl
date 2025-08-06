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

################################################################################
# Initialize variables for analysis ############################################
################################################################################
sec = 'T'
dimension_filters = Dict{Symbol,Any}()
ECC = M.ReadDisk(loc2.HDF5_path, "E2020DB/ECCKey")
push!(dimension_filters, :ECC => ECC[ECC .∉ Ref(["ForeignPassenger", "ForeignFreight"])])
push!(dimension_filters, :Area => Canada)

EuFPol = J.diff_fast("EuFPol", loc1, loc2; dimension_filters, sec)
J.plot_diff(EuFPol; dim="ECC", num=10, title="EuFPol diffs by ECC")
@rsubset! EuFPol :Diff != 0
J.plot_diff(EuFPol, dim = "ECC")

@rsubset! EuFPol :ECC == "CommercialOffRoad" :Poll == "CO2" :Year >= 2020
EuFPol_g = @by EuFPol :Year begin
  :Spruce = sum(:Spruce)
  :Tanoak = sum(:Tanoak)
  :Diff = sum(abs.(:Diff))
end 
EuFPol_g.PDiff = EuFPol_g.Diff ./ EuFPol_g.Spruce

# Almost all Gasoline and Diesel
J.plot_diff(EuFPol; dim="FuelEP", num=10, title="EuFPol diffs by FuelEP") 
push!(dimension_filters, :FuelEP => "Gasoline")
push!(dimension_filters, :ECC => "Passenger")
push!(dimension_filters, :EC => "Passenger")

# CO2 mainly with trace COX
@rsubset! EuFPol :Spruce != 0 && :Tanoak != 0
J.plot_diff(EuFPol; dim="Poll", num=10, title="EuFPol diffs by Poll") 
push!(dimension_filters, :Poll => "CO2")
# ON, AB, BC are the biggest contributors, though many are present
J.plot_diff(EuFPol; dim="Area", num=10, title="EuFPol diffs by Area")
push!(dimension_filters, :Area => "ON")


Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters)
s = 'T'
Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters, s)
J.plot_diff(Polute; dim="Tech", num=10, title="Polute diffs by Tech")
push!(dimension_filters, :Tech => ["LDVGasoline"])
push!(dimension_filters, :Year => "2030")

POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters, s)
i = findall(vars.Variable .== "POCA" .&& first.(vars.Database) .== 'T')
i = i[1]
vars[i,:Database] = "TOutput2"
push!(dimension_filters, :EC => "Passenger")
POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters, s)
@rsubset POCA :Diff != 0

DmdFEPTech = J.diff_fast("DmdFEPTech", loc1, loc2; dimension_filters, s)
@rsubset DmdFEPTech !isapprox(:Diff, 0)

DmdFuelTech = J.diff_fast("DmdFuelTech", loc1, loc2; dimension_filters, s)
dimension_filters[:FuelEP] = "Gasoline"
dimension_filters[:Tech] = "LDVGasoline"
push!(dimension_filters, :Fuel => "Gasoline")
DmdFuelTech = J.diff_fast("DmdFuelTech", loc1, loc2; dimension_filters, s)

DmFrac = J.diff_fast("DmFrac", loc1, loc2; dimension_filters, s)
Dmd = J.diff_fast("Dmd", loc1, loc2; dimension_filters, s)

EEImpact = J.diff_fast("EEImpact", loc1, loc2; dimension_filters, s)
EESat = J.diff_fast("EESat", loc1, loc2; dimension_filters, s)
xDmd = J.diff_fast("xDmd", loc1, loc2; dimension_filters, s)
DSMEU = J.diff_fast("DSMEU", loc1, loc2; dimension_filters, s)
SqDmd = J.diff_fast("SqDmd", loc1, loc2; dimension_filters, s)
SqEUTechMap = J.diff_fast("SqEUTechMap", loc1, loc2; dimension_filters, s)
RPEI = J.diff_fast("RPEI", loc1, loc2; dimension_filters, s)
CERSM = J.diff_fast("CERSM", loc1, loc2; dimension_filters, s)
TSLoad = J.diff_fast("TSLoad", loc1, loc2; dimension_filters, s)
UMS = J.diff_fast("UMS", loc1, loc2; dimension_filters, s)
DDay = J.diff_fast("DDay", loc1, loc2; dimension_filters, s)
CUF = J.diff_fast("CUF", loc1, loc2; dimension_filters, s)
WCUF = J.diff_fast("WCUF", loc1, loc2; dimension_filters, s)
DDCoefficient = J.diff_fast("DDCoefficient", loc1, loc2; dimension_filters, s)
DDayNorm = J.diff_fast("DDayNorm", loc1, loc2; dimension_filters, s)


DER = J.diff_fast("DER", loc1, loc2; dimension_filters, s)
DERRRExo = J.diff_fast("DERRRExo", loc1, loc2; dimension_filters, s)
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters, s)
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters, s)
DERV = J.diff_fast("DERV", loc1, loc2; dimension_filters, s)
i = findall(vars.Variable .∈ Ref(["DERV", "DERAV"]) .&& first.(vars.Database) .== 'T')
i = i[1]
vars[i,:Database] = "TOutput2"
vars[i,:]
DERV = J.diff_fast("DERV", loc1, loc2; dimension_filters, s)
StockAdjustment = J.diff_fast("StockAdjustment", loc1, loc2; dimension_filters, s)
DERAV = J.diff_fast("DERAV", loc1, loc2; dimension_filters, s)
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

xMMSF = J.diff_fast("xMMSF", loc1, loc2; dimension_filters = df2020, s) # didn't work
CMMSF = J.diff_fast("CMM", loc1, loc2; dimension_filters = df2020, s) # didn't work

MMSM0 = J.diff_fast("MMSM0", loc1, loc2; dimension_filters = df2020, s) # bad
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

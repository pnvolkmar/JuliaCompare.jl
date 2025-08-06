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
sec = 'T'
dimension_filters = Dict{Symbol,Any}()
push!(dimension_filters, :Area => Canada, :Year => string.(1986:2050))
push!(dimension_filters, :Poll => "SOX")

onf = copy(dimension_filters)
push!(onf, :ECC => ["Passenger","Freight"], :EC => ["Passenger","Freight"])
push!(onf, :Area => "ON", :Year => string.(2022:2023))

################################################################################
# Analysis of variables: Overview ########################################################
################################################################################
EC = J.var("ECKey", loc1, sec)
EC = J.var_id(3906, loc1)
EuFPol = J.diff_fast("EuFPol", loc1, loc2; dimension_filters=onf, sec) # 
J.add_pdiff!(EuFPol)
@rsubset! EuFPol :PDiff != 0 :ECC != "ForeignPassenger" :ECC != "ForeignFreight"
J.plot_diff(EuFPol, dim="ECC")
J.plot_diff(EuFPol, dim="Area")
TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters=onf, sec) # 
@rsubset! TotPol :Diff != 0 :ECC != "ForeignPassenger" :ECC != "ForeignFreight"
J.plot_diff(TotPol, dim="ECC")
J.plot_diff(TotPol, dim="Area")

    # TotPol[ecc,poll,area] = EnPol[ecc,poll,area]+NcPol[ecc,poll,area]+
    #   MEPol[ecc,poll,area]+VnPol[ecc,poll,area]+FlPol[ecc,poll,area]+
    #   FuPol[ecc,poll,area]+ORMEPol[ecc,poll,area]+SqPol[ecc,poll,area]
EnPol = J.diff_fast("EnPol", loc1, loc2; dimension_filters=onf, sec) # Problem
NcPol = J.diff_fast("NcPol", loc1, loc2; dimension_filters=onf, sec) # 
MEPol = J.diff_fast("MEPol", loc1, loc2; dimension_filters, sec) # 
@rsubset! MEPol abs(:Diff) >= 0.001
J.plot_diff(MEPol, dim ="ECC"; title = "MEPol differences after fix")
VnPol = J.diff_fast("VnPol", loc1, loc2; dimension_filters=onf, sec) # 
FlPol = J.diff_fast("FlPol", loc1, loc2; dimension_filters=onf, sec) # 
FuPol = J.diff_fast("FuPol", loc1, loc2; dimension_filters=onf, sec) # 
ORMEPol = J.diff_fast("ORMEPol", loc1, loc2; dimension_filters=onf, sec) # 
SqPol = J.diff_fast("SqPol", loc1, loc2; dimension_filters=onf, sec) # 

EuFPol = J.diff_fast("EuFPol", loc1, loc2; dimension_filters=onf, sec) # Problem
@rsubset! EuFPol abs(:Diff) >= 1
J.add_pdiff!(EuFPol)
J.plot_diff(EuFPol, dim = "FuelEP"; title = "EuFPol OtherMetalMining CO2, BC,ON,QC, by FuelEP")

push!(onf, :FuelEP => unique(EuFPol.FuelEP))

Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters=onf, sec = sec) # problem
J.add_pdiff!(Polute)
@rsubset! Polute abs(:PDiff) != 0
J.plot_diff(Polute, dim = "Enduse")

i = findall(vars.Variable .== "POCA" .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
# Polute[enduse,fuelep,ec,poll,area] = EuDem[enduse,fuelep,ec,area]*
#     POCA[enduse,fuelep,ec,poll,area]+xPolute[enduse,fuelep,ec,poll,area]

POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters=onf, sec) # Problem
J.add_pdiff!(POCA)
@rsubset! POCA :PDiff != 0
EuDem = J.diff_fast("EuDem", loc1, loc2; dimension_filters=onf, sec) # 
@rsubset! EuDem :Diff != 0
J.add_pdiff!(EuDem)
@rsubset EuDem :Area == "BC" :Year == 2050
xPolute = J.diff_fast("xPolute", loc1, loc2; dimension_filters=onf, sec) # 
J.add_pdiff!(EuDem)
@rsubset! EuDem abs(:PDiff) != 0

PolSw = J.diff_fast("PolSw", tan_base, loc2; dimension_filters=onf, sec) # matches, 2.0 everywhere

#   POCA[enduse,fuelep,tech,ec,poll,area] = POCX[enduse,fuelep,tech,ec,poll,area]*
#     RM[tech,ec,poll,area]

POCX = J.diff_fast("POCX", loc1, loc2; dimension_filters=onf, sec) # matches, everywhere
@rsubset! POCX abs(:Diff) > 1e-3
@rsubset! POCX :Year ∈ [2022:2023]
RM = J.diff_fast("RM", loc1, loc2; dimension_filters=onf, sec) # issue

#   @. RM = (1.0-RP)*(1.0-VR)*xRM

xRM = J.diff_fast("xRM", loc1, loc2; dimension_filters=onf, sec) # issue
@rsubset xRM abs(:Diff) > 1e-1
RP = J.diff_fast("RP", loc1, loc2; dimension_filters=onf, sec) # match
@rsubset RP :Tanoak == 1

i = findall(vars.Variable .== "VR" .&& first.(vars.Database) .== 'I')
i = i[1]
vars[i,:Database] = "IOutput2"

VR = J.diff_fast("VR", loc1, loc2; dimension_filters=onf, sec) # match
@rsubset VR abs(:Diff) > 1e-1



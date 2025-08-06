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
red_base = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaRedwood\\2020Model\\Base\\database.hdf5", "Redwood");
pn_base = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaPine\\2020Model\\Base", "Pine");
red_ogref = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaRedwood\\2020Model\\OGRef\\database.hdf5", "Redwood");
pn_ogref = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaPine\\2020Model\\OGRef", "Pine");

################################################################################
# Initialize filters and sec for analysis ######################################
################################################################################
sec = 'I'
dimension_filters = Dict{Symbol,Any}()
push!(dimension_filters, :Area => Canada, :Year => string.(1986:2050))
ECC = M.ReadDisk(loc2.HDF5_path, "E2020DB/ECCKey")
push!(dimension_filters, :ECC => ECC[occursin.("Foreign", ECC) .== false])

cnst = copy(dimension_filters)
push!(cnst, :EC => "Construction",  :ECC => "Construction")
push!(cnst, :Poll => ["CO2","PMT","PM10"])
push!(cnst, :Area => ["AB"])
push!(cnst, :AreaTOM => ["AB"])
push!(cnst, :Fuel => "NaturalGas", :FuelEP => "NaturalGas")
push!(cnst, :Enduse => "Heat")
push!(cnst, :Tech => ["Electric","Gas","LPG"])
################################################################################
# Analysis of variables: Overview ##############################################
################################################################################

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters, sec) # 
@rsubset! TotPol abs(:Diff) >= 1
J.add_pdiff!(TotPol)
J.plot_diff(TotPol, dim="ECC"; title = "Differences in TotPol")

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters=cnst, sec) # 
J.add_pdiff!(TotPol)
J.plot_diff(TotPol, dim="Poll")
J.plot_diff(TotPol, dim="Area")

    # TotPol[ecc,poll,area] = EnPol[ecc,poll,area]+NcPol[ecc,poll,area]+
    #   MEPol[ecc,poll,area]+VnPol[ecc,poll,area]+FlPol[ecc,poll,area]+
    #   FuPol[ecc,poll,area]+ORMEPol[ecc,poll,area]+SqPol[ecc,poll,area]
EnPol = J.diff_fast("EnPol", loc1, loc2; dimension_filters=cnst, sec) # 10 percent off 
J.add_pdiff!(EnPol)
NcPol = J.diff_fast("NcPol", loc1, loc2; dimension_filters=cnst, sec) # Problem
J.add_pdiff!(NcPol)
MEPol = J.diff_fast("MEPol", loc1, loc2; dimension_filters=cnst, sec) # Problem
VnPol = J.diff_fast("VnPol", loc1, loc2; dimension_filters=cnst, sec) # 
FlPol = J.diff_fast("FlPol", loc1, loc2; dimension_filters=cnst, sec) # 
FuPol = J.diff_fast("FuPol", loc1, loc2; dimension_filters=cnst, sec) # 
ORMEPol = J.diff_fast("ORMEPol", loc1, loc2; dimension_filters=cnst, sec) # 
SqPol = J.diff_fast("SqPol", loc1, loc2; dimension_filters=cnst, sec) # Problem

MEReduce = J.diff_fast("MEReduce", loc1, loc2; dimension_filters=cnst, sec) # Fine

J.add_pdiff!(MEPol)
# MEPol[ecc,poll,area] = (MEPOCA[ecc,poll,area]*MEDriver[ecc,area]*MEPolSwitch[ecc,poll,area])+
#   (xMEPol[ecc,poll,area]*(1-MEPolSwitch[ecc,poll,area]))
MEPOCA = J.diff_fast("MEPOCA", loc1, loc2; dimension_filters=cnst, sec) # 
J.add_pdiff!(MEPOCA) # fine
MEDriver = J.diff_fast("MEDriver", loc1, loc2; dimension_filters=cnst, sec) # matches
Driver = J.diff_fast("Driver", loc1, loc2; dimension_filters=cnst, sec) # matches
xDriver = J.diff_fast("xDriver", loc1, loc2; dimension_filters=cnst, sec) # matches
DriverMultiplier = J.diff_fast("DriverMultiplier", loc1, loc2; dimension_filters=cnst, sec) # matches
GRPAdj = J.diff_fast("GRPAdj", loc1, loc2; dimension_filters=cnst, sec) # matches
MEPolSwitch = J.diff_fast("MEPolSwitch", loc1, loc2; dimension_filters=cnst, sec) # matches, 1
xMEPol = J.diff_fast("xMEPol", loc1, loc2; dimension_filters=cnst, sec) # matches

# MEPOCA[ecc,poll,area] = MEPOCX[ecc,poll,area]*MERM[ecc,poll,area]*MEPOCM[ecc,poll,area]
MEPOCX = J.diff_fast("MEPOCX", loc1, loc2; dimension_filters=cnst, sec = 'M') # issue
J.add_pdiff!(MEPOCX) 
@rsubset! MEPOCX :Diff != 0

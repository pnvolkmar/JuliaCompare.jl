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
sec = 'I'
dimension_filters = Dict{Symbol,Any}()
push!(dimension_filters, :Area => Canada, :Year => string.(1986:2050))
geo = copy(dimension_filters)
push!(geo, :Fuel => "Geothermal")
push!(geo, :Tech => "Geothermal")
pop!(geo, :Tech)
pop!(geo, :Poll)
push!(geo, :Year => "2000")
push!(geo, :Area => "ON", :ECC => "Food", :EC => "Food")
################################################################################
# Analysis of variables: Overview ########################################################
################################################################################
EuDemand = J.diff_fast("EuDemand", loc1, loc2; dimension_filters=geo, sec) # 
@rsubset! EuDemand :Diff != 0 :ECC != "ForeignPassenger" :ECC != "ForeignFreight"
J.plot_diff(EuDemand, dim="ECC")

FPF = J.diff_fast("FPF", loc1, loc2; dimension_filters=geo, sec) # 
J.add_pdiff!(FPF)
@rsubset! FPF abs(:Diff) > 1 abs(:PDiff) > 1
sort!(FPF, :PDiff)

EGFA = J.diff_fast("EGFA", loc1, loc2; dimension_filters=geo, sec) # 
@rsubset! EGFA :Diff != 0

Dmd = J.diff_fast("Dmd", loc1, loc2; dimension_filters=geo, sec) # 
unique(Dmd.EC)

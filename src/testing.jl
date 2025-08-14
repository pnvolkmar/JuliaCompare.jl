
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

cfilename = "TInput"
CODE_FOLDER = raw"\\Pink\c\2020CanadaPine\Engine"
DATA_FOLDER = raw"\\Pink\c\2020CanadaPine\2020Model"
J.list_var(cfilename, CODE_FOLDER, DATA_FOLDER)


vars = J.list_vars(DATA_FOLDER1, db_files);
vars_j = J.list_vars(joinpath(DATA_FOLDER2,"database.hdf5"))
loc1 = J.Loc_p(vars, raw"\\Pink\c\2020CanadaPineAqua\2020Model\Ref25", "Pine");
loc2 = J.Loc_j(vars_j, joinpath(raw"\\Pink\c\2020CanadaRedwoodAqua\2020Model\Ref25", "database.hdf5"), "Redwood");

sec = 'T'
filter = Dict{Symbol,Any}()
push!(filter, :Area => Canada, :Year => string.(1986:2050))


vname2, dbname2 = J.lookup_database("ECKey", loc2; sec = 'T')
arr2, sets2 = J.arr_set(vname2, dbname2, loc2)
vname1, dbname1 = J.lookup_database("RPPMap", loc1; sec = 'T')
J.arr_set(vname1, dbname1, loc1)
vattr = M.data_attrs(loc2.HDF5_path, string(dbname2,"/", vname2))

J.diff_fast("TotPol", loc1, loc2; filter, sec)
J.diff_fast("ECKey", loc1, loc2; filter, sec)
df = J.var("ECKey", loc1; filter, sec)
df = J.var("TotPol", [loc1, loc2]; filter, sec, diff = :from_first, pdiff = :from_first)
J.plot_sets(df; col = :Redwood_minus_Pine, dim = "Area", title = "Peter's Chart", units = "Pounds")
df
J.plot_lines(df, [loc1, loc2]; title = "Peter", units = "test")

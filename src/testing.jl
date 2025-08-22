
using Revise
import JuliaCompare as J
import PromulaDBA as P
import SmallModel as M 
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

################################################################################
# Inputs for file locations ####################################################
################################################################################
loc1 = J.loc(raw"\\Pink\c\2020CanadaPine\2020Model\Ref25", "Pine");
loc2 = J.loc(raw"\\Pink\c\2020CanadaRedwood\2020Model\Ref25", "Redwood");

vnames = ["MEInput/FuA0","MEInput/FuB0"]
variables, summary = J.compare_vars(vnames, [loc1, loc2]; diff = true, pdiff = true)
test = J.var("RInput/DmFracMin", [loc1, loc2]; diff = true, pdiff = true)
vnames = ["RInput/DmFracMin","CInput/DmFracMin"]
variables, summary = J.compare_vars(vnames, [loc1, loc2]; diff = true, pdiff = true)
locations = [redwood, pine]
redwood = J.loc("C:/2020CanadaRedwood/2020Model/Ref25", "Redwood")
pine = J.loc("C:/2020CanadaPine/2020Model/Ref25", "Pine")
DmFrac = variables["RInput/DmFracMin"]
J.plot_sets(DmFrac; col = "Pine_minus_Redwood", dim = "EC")

vnames = ["ROutput/DERRRExo","COutput/DERRRExo", "EGInput/xUnFlFr"]
J.lookup_database("COutput/DERRRExo", loc1)
J.lookup_database("ROutput/DERRRExo", loc2)
locations = [loc1,loc2]
variables, summary = J.compare_vars(vnames,locations;diff=true,pdiff=true)
summary
df = variables["ROutput/DERRRExo"]

red_unit = M.ReadDisk(loc2.HDF5_path,"EGInput/Unit")
red_code = M.ReadDisk(loc2.HDF5_path,"EGInput/UnCode")
red_unit[red_unit.!="Null"]
red_code[red_code.!="Null"]
pine_code = P.data(joinpath(loc1.DATA_FOLDER,"EGInput.dba"), "UnCode")
pine_unit = P.data(joinpath(loc1.DATA_FOLDER,"EGInput.dba"), "UnitKey")
red_unit == red_code
setdiff(red_code, red_unit)
setdiff(red_code, pine_code)

red_unit = M.ReadDisk(loc2.HDF5_path,"MInput/OGUnit")
red_code = M.ReadDisk(loc2.HDF5_path,"EGInput/UnCode")
red_unit[red_unit.!="Null"]
red_code[red_code.!="Null"]
pine_code = P.data(joinpath(loc1.DATA_FOLDER,"2020DB.dba"), "OGCode")
J.lookup_database("OGCode", loc1)
pine_unit = P.data(joinpath(loc1.DATA_FOLDER,"EGInput.dba"), "UnitKey")
red_unit == red_code
setdiff(red_code, red_unit)
setdiff(red_code, pine_code)


describe(df)
@rsubset df :Redwood_minus_Pine != 0
sort!(df, :Redwood_minus_Pine)
J.plot_sets(df, col = :Redwood_minus_Pine, dim = "EC")

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

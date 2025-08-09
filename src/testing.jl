
test with ECKey and TotPol
vname, dbname = J.lookup_database("RPPMap", loc2; sec = 'T')
arr, set = J.arr_set(vname, dbname, loc2)
SmallModel.data_attrs(loc2.HDF5_path, string(dbname,"/", vname))



function arr_set(vname::String, dbname::String, loc::Loc_j)
  set = SmallModel.ReadSets(loc.HDF5_path, string(dbname,"/", vname))
  arr = ReadDisk(loc.HDF5_path, string(dbname,"/", vname))
  if :Unit ∈ keys(sets)
    UnCode = ReadDisk(loc.HDF5_path, "EGInput/UnCode")
    set.Unit[:] = UnCode[:]
  end
  return(arr, set)
end

function arr_set(vname::String, dbname::String, loc::Loc_p)
  arr = P.data(joinpath(loc.DATA_FOLDER,string(dbname,".dba")), vname)
  idx = findfirst(loc.vars.Variable .== vname .&& loc.vars.Database .== dbname)
  dim_names = split(loc.vars.Dimensions[idx], ",")
  dim_pairs = loc.vars.DPairs[idx]
  set = [P.data(key,value) for (key, value) in dim_pairs]
  names(set) = dim_names
  return(arr, set)
end


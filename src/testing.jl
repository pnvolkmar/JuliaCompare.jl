
test with ECKey and TotPol
vname2, dbname2 = J.lookup_database("RPPMap", loc2; sec = 'T')
arr2, sets2 = J.arr_set(vname2, dbname2, loc2)
vname1, dbname1 = J.lookup_database("RPPMap", loc1; sec = 'T')
J.arr_set(vname1, dbname1, loc1)
SmallModel.data_attrs(loc2.HDF5_path, string(dbname,"/", vname))

J.diff_fast("TotPol", loc1, loc2; dimension_filters, sec)

using Revise
import JuliaCompare as J
import PromulaDBA as P
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

BASE_FOLDER = raw"\\Pink\c\2020CanadaSpruce"
BASE_FOLDER2 = raw"\\Pink\c\2020CanadaTanoak"
SCENARIO1 = "CalibCom"
SCENARIO2 = "CalibCom"

CODE_FOLDER = joinpath(BASE_FOLDER, "Engine")
DATA_FOLDER1 = joinpath(BASE_FOLDER, "2020Model", SCENARIO1)
DATA_FOLDER2 = joinpath(BASE_FOLDER2, "2020Model")
E2020_Folder = joinpath(BASE_FOLDER, "2020Model")
HDF5_path = joinpath(DATA_FOLDER2, "database.hdf5")

vars = J.list_vars(CODE_FOLDER, DATA_FOLDER1, db_files);
vars_j = J.list_vars(HDF5_path)
loc1 = J.Loc_p(vars, DATA_FOLDER1, "Spruce");
loc2 = J.Loc_j(vars_j, HDF5_path, "Tanoak");
loc3 = J.Loc_j(vars_j, joinpath(BASE_FOLDER2, "2020Model","CalibPrice3","database.hdf5"), "CalibPrice3");

df = J.var("CInput/xDmd", loc1) # No difference in lubricants co2
df2 = J.var("CInput/xDmd", loc2)
df3 = J.diff(df, df2; name1="Spruce", name2="Tanoak")
@rsubset df3 abs(:Diff) > 1e-5

for n in names(df3)[1:(end-3)]
  J.hist_diff(df3; dim = n)
end

for n in names(df3)[1:(end-4)]
  J.plot_diff(df3; dim=n, num=3, title="Dmd difference by $n")
end

df4 = @rsubset df3 :Area == "MX"
J.plot_diff(df4; dim="EC", num = 5, title = "MX's Dmd diff by EC")
df5 = @rsubset df3 :EC == "Wholesale"
J.plot_diff(df5; dim="Area", num = 5, title = "Wholesale's Dmd diff by Area")


df4 = @rsubset df3 :Area == "ENC"
J.plot_diff(df4; dim="EC", num = 5, title = "ENC's Dmd diff by EC")
df5 = @rsubset df3 :EC == "Health"
J.plot_diff(df5; dim="Area", num = 5, title = "Health's Dmd diff by Area")


df5 = @rsubset df3 :EC == "OtherCommercial"
J.plot_diff(df5; dim="Area", num = 5, title = "OtherCommercial's Dmd diff by Area")

xDEE = J.diff("CInput/xDEE", loc1, loc2)
J.filter_ns(xDEE) # -99
DEEBeforeStd = J.diff("COutput/DEEBeforeStd", loc1, loc2)
J.filter_ns(DEEBeforeStd) # big differences
DEStd = J.diff("CInput/DEStd", loc1, loc2)
J.filter_ns(DEStd)
DEStdP = J.diff("CInput/DEStdP", loc1, loc2)
J.filter_ns(DEStdP)
DEE = J.diff("COutput/DEE", loc1, loc2)
J.filter_ns(DEE)

DEM = J.diff("CInput/DEM", loc1, loc2)
J.filter_ns(DEM)
DEMM = J.diff("CCalDB/DEMM", loc1, loc2)
J.filter_ns(DEMM) # big differences

xDEMM = J.diff("CInput/xDEMM", loc1, loc2)
J.filter_ns(xDEMM) # big differences

DEEThermalMax = J.diff("CInput/DEEThermalMax", loc1, loc2)
J.filter_ns(DEEThermalMax) # big differences

ECFP = J.diff("COutput/ECFP", loc1, loc2)
J.filter_ns(ECFP) # same

DCTC = J.diff("COutput/DCTC", loc1, loc2)
J.filter_ns(DCTC) # same

DEESw = J.diff("CInput/DEESw", loc1, loc2)
J.filter_ns(DEESw) # same

EUPC = J.diff("COutput/EUPC", loc1, loc2)
J.filter_85(EUPC) # differ on LPG

PCERF = J.diff("SInput/PCERF", loc1, loc2)
J.filter_85(PCERF) # agree on LPG

ECFP = J.diff("COutput/ECFP", loc1, loc2)
J.filter_85(ECFP) # agree

PCLV = J.diff("MOutput/PCLV", loc1, loc2)
J.filter_85(PCLV) # agree

PER = J.diff("COutput/PER", loc1, loc2)
J.filter_85(PER) # agree

DSt = J.diff("COutput/DSt", loc1, loc2)
J.filter_85(DSt) # agree

PDif = J.diff("CInput/PDif", loc1, loc2)
J.filter_85(PDif) # agree

using SmallModel
import SmallModel as M 
M.ReadDisk(HDF5_path, "CInput/CurTime")

xDEE = J.diff("CInput/xDEE", loc1, loc2)
J.filter_bm(xDEE) # -99

DEE = J.diff("COutput/DEE", loc1, loc2)
J.filter_bm(DEE) # different

  # CurTime = 16 

  DEStd = J.diff("CInput/DEStd", loc1, loc2)
  J.filter_bm(DEStd) # same

  DEStdP = J.diff("CInput/DEStdP", loc1, loc2)
  J.filter_bm(DEStdP) # same

  DEM = J.diff("CInput/DEM", loc1, loc2)
  J.filter_bm(DEM) # same

  DEPM = J.diff("CInput/DEPM", loc1, loc2)
  J.filter_bm(DEPM) # same
  
  DFTC = J.diff("COutput/DFTC", loc1, loc2)
  J.filter_bm(DFTC) # same
  
  Inflation = J.diff("MOutput/Inflation", loc1, loc2)
  J.filter_bm(Inflation) # same
  
DEMM = J.diff("CCalDB/DEMM", loc1, loc2)
J.filter_bm(DEMM) # different because DEE is zero and it's used to limit DEMM

  xDEMM = J.diff("CInput/xDEMM", loc1, loc2)
  J.filter_bm(xDEMM) # same

DFPN = J.diff("COutput/DFPN", loc1, loc2)
J.filter_bm(DFPN) # different

  DCCRN = J.diff("COutput/DCCR", loc1, loc2)
  J.filter_bm(DCCRN) # basically the same

  DCTC = J.diff("COutput/DCTC", loc1, loc2)
  J.filter_bm(DCTC) # basically the same

DCCN = J.diff("COutput/DCCN", loc1, loc2)
J.filter_bm(DCCN) # v different

DCCN = J.diff("COutput/DCCN", loc3, loc2)
J.filter_bm(DCCN) # v different

  xDCC = J.diff("CInput/xDCC", loc1, loc2)
  J.filter_bm(xDCC) # same curtime

  J.filter_bm(DEM) # good
  J.filter_bm(DCTC) # good

  DEEStd = J.diff("COutput/DEEStd", loc1, loc2)
  J.filter_bm(DEEStd) # unknown

# DEEStd = J.diff("COutput/DEEStd", loc1, loc2)
# J.filter_bm(DEEStd) # unknown
  J.filter_bm(xDCC) # good
  J.filter_bm(DEStd) # good
  J.filter_bm(DEStdP) # good
  




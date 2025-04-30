using Revise
import JuliaCompare as J
import PromulaDBA as P
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

BASE_FOLDER = raw"\\Pink\c\2020CanadaSpruce"
BASE_FOLDER2 = raw"\\Pink\c\2020CanadaTanoak"
SCENARIO1 = "Calib2"
SCENARIO2 = "Calib2"

CODE_FOLDER = joinpath(BASE_FOLDER, "Engine")
DATA_FOLDER1 = joinpath(BASE_FOLDER, "2020Model", SCENARIO1)
DATA_FOLDER2 = joinpath(BASE_FOLDER2, "2020Model", SCENARIO2)
E2020_Folder = joinpath(BASE_FOLDER, "2020Model")
HDF5_path = joinpath(DATA_FOLDER2, "database.hdf5")

vars = J.list_vars(CODE_FOLDER, DATA_FOLDER1, db_files);
vars_j = J.list_vars(HDF5_path)
loc1 = J.Loc_p(vars, DATA_FOLDER1, "Spruce");
loc2 = J.Loc_j(vars_j, HDF5_path, "Tanoak");
loc3 = J.Loc_j(vars_j, joinpath(BASE_FOLDER2, "2020Model","CalibPrice3","database.hdf5"), "CalibPrice3");

BaseAdj = J.diff("SCalDB/BaseAdj", loc1, loc2); J.f_on(BaseAdj)
xPkLoad = J.diff("SInput/xPkLoad", loc1, loc2); J.f_on(xPkLoad)
xPkSavECC = J.diff("SInput/xPkSavECC", loc1, loc2); J.f_on(xPkSavECC)
PkLoad = J.diff("SOutput/PkLoad", loc1, loc2); J.f_on(PkLoad)
TDEF = J.diff("SInput/TDEF", loc1, loc2); J.f_on(TDEF)
CgNetPk = J.diff("SOutput/CgNetPk", loc1, loc2); J.f_on(CgNetPk)
HoursPerMonth = J.diff("SInput/HoursPerMonth", loc1, loc2); J.f_on(HoursPerMonth)
CgLDCECC = J.diff("SOutput/CgLDCECC", loc1, loc2); J.f_on(CgLDCECC)
CgLDCSoldECC = J.diff("SOutput/CgLDCSoldECC", loc1, loc2); J.f_on(CgLDCSoldECC)
J.f_on(xPkSavECC)
sum(xPkSavECC.Spruce)
sum(xPkSavECC.Tanoak)
SLDC = J.diff("SOutput/SLDC", loc1, loc2); J.f_on(SLDC)
xSalSw = J.diff("SInput/xSalSw", loc1, loc2); J.f_on(xSalSw)
ECCCLMap_t = J.var("E2020DB/ECCCLMap", loc2);
ECCCLMap_s = J.var("2020DB/ECCCLMap", loc1);
ECCCLMap = J.diff(ECCCLMap_s, ECCCLMap_t); J.f_on(ECCCLMap)
CLSF = J.diff("SCalDB/CLSF", loc1, loc2); J.f_on(CLSF)
xCLSF = J.diff("SInput/xCLSF", loc1, loc2); J.f_on(xCLSF)
LDCECC = J.diff("SOutput/LDCECC", loc1, loc2); J.f_on(LDCECC)
SaEC = J.diff("SOutput/SaEC", loc1, loc2); J.f_on(SaEC)
xSaEC = J.diff("SInput/xSaEC", loc1, loc2); J.f_on(xSaEC)
CgLDCECC = J.diff("SOutput/CgLDCECC", loc1, loc2); J.f_on(CgLDCECC)
LDCEUECC = J.diff("SOutput/LDCEUECC", loc1, loc2); J.f_on(LDCEUECC)
ESales = J.diff("IOutput/ESales", loc1, loc2); J.f_on(ESales)
LSF = J.diff("ICalDB/LSF", loc1, loc2); J.f_on(LSF)
HPKM = J.diff("SCalDB/HPKM", loc1, loc2); J.f_on(HPKM)
TSLoad = J.diff("IInput/TSLoad", loc1, loc2); J.f_on(TSLoad)
xPkSav = J.diff("IInput/xPkSav", loc1, loc2); J.f_on(xPkSav)

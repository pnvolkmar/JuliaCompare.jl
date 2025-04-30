using Revise
import JuliaCompare as J
import PromulaDBA as P
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

BASE_FOLDER = raw"\\Pink\c\2020CanadaSpruce"
BASE_FOLDER2 = raw"\\Pink\c\25.04.14_2020CanadaTanoak"
SCENARIO1 = "Base"
SCENARIO2 = "Base"

CODE_FOLDER = joinpath(BASE_FOLDER, "Engine")
DATA_FOLDER1 = joinpath(BASE_FOLDER, "2020Model")
DATA_FOLDER2 = joinpath(BASE_FOLDER2, "2020Model", SCENARIO2)
E2020_Folder = joinpath(BASE_FOLDER, "2020Model")
HDF5_path = joinpath(DATA_FOLDER2, "database.hdf5")

vars = J.list_vars(CODE_FOLDER, DATA_FOLDER1, db_files);
vars_j = J.list_vars(HDF5_path)
loc1 = J.Loc_p(vars, DATA_FOLDER1, "Spruce");
loc2 = J.Loc_j(vars_j, HDF5_path, "Tanoak");
loc3 = J.Loc_j(vars_j, joinpath(BASE_FOLDER2, "2020Model","CalibPrice3","database.hdf5"), "CalibPrice3");
S_Ref24 = J.Loc_p(vars, joinpath(E2020_Folder,"OGRef"), "S_Ref24")
T_Ref24 = J.Loc_j(vars_j, joinpath(raw"\\Pink\c\2020CanadaTanoak", "2020Model","Ref24","database.hdf5"), "T_Ref24");
UnOR = J.var("EGInput/UnORTM", loc1)
UnOR_j = J.var("EGInput/UnORTM", loc2)
# Extract the number from each TimeP value and add parentheses around it
UnOR_j.TimeP = map(tp -> "TimeP($(match(r"\d+", tp).match))", UnOR_j.TimeP)
UnOR = J.diff(UnOR,UnOR_j; name1 = "Spruce", name2 = "Tanoak")


UnArea = J.var("UnArea", loc1)
UnPlant = J.var("UnPlant", loc1)
UnOR = leftjoin(UnOR, UnArea, on = :Unit, makeunique=true)
UnOR = leftjoin(UnOR, UnPlant, on = :Unit, makeunique=true)
# Reorder columns to move UnArea to be the second column
new_order = [:Unit, :UnArea, :UnPlant, :TimeP, :Month, :Year, :Spruce, :Tanoak, :Diff]
UnOR = UnOR[:, new_order]
@rsubset UnOR abs(:Diff) > 1e-5 :Year == 2025 :TimeP == "TimeP(1)" :Month == "Summer" :Unit != "Null" :UnArea ∈ Canada
@rsubset UnOR :Year == 2023 :TimeP == "TimeP(1)" :Month == "Summer" :Unit != "Null" :UnArea == "NB" :Diff > 1e-5


df = J.diff("SOutput/DemRq", loc1,loc2) # 
@rsubset df :Area == "AB" :ECC == "ConventionalGasProduction" :Year == 2025 abs(:Diff) > 1e-5

DemRq = J.diff("SOutput/DemRq", loc1,loc2) # 
J.f_og(DemRq)

DmdRq = J.diff("IOutput/DmdRq", loc1,loc2) # 
J.f_og(DmdRq)

ECUF = J.diff("MOutput/ECUF", loc1,loc2) # 
J.f_og(ECUF)

RPEI = J.diff("IOutput/RPEI", loc1,loc2) # 
J.f_og(RPEI)

DEE = J.diff("IOutput/DEE", loc1,loc2) # 
J.f_og(DEE)
DEEPoll = J.diff("IOutput/DEEPoll", loc1,loc2) # 
J.f_og(DEEPoll)
DEESw = J.diff("IInput/DEESw", loc1,loc2) # 
J.f_og(DEESw)

PEE = J.diff("IOutput/PEE", loc1,loc2) # 
J.f_og(PEE)

  PEMM = J.diff("ICalDB/PEMM", loc1,loc2) # 
  J.f_og(PEMM)
  PEMMM = J.diff("IOutput/PEMMM", loc1,loc2) # 
  J.f_og(PEMMM)
  PER = J.diff("IOutput/PER", loc1,loc2) # 
  J.f_og(PER)
  CHRM = J.diff("IInput/CHRM", loc1,loc2) # 
  J.f_og(CHRM)
  CHR = J.diff("ICalDB/CHR", loc1,loc2) # 
  J.f_og(CHR)

DSt = J.diff("IOutput/DSt", loc1,loc2) #  
J.f_og(DSt)

CERSM = J.diff("ICalDB/CERSM", loc1,loc2) # 
J.f_og(CERSM)

AMSF = J.diff("IOutput/AMSF", loc1,loc2) # 
J.f_og(AMSF)

DCCR = J.diff("IOutput/DCCR", loc1,loc2) # 
J.f_og(DCCR)

FPCFSTech = J.diff("IOutput/FPCFSTech", loc1,loc2) # 
J.f_og(FPCFSTech)

DEM = J.diff("IInput/DEM", loc1,loc2) # 
J.f_og(DEM)

DFTC = J.diff("IOutput/DFTC", loc1,loc2) # 
J.f_og(DFTC)

DCTC = J.diff("IOutput/DCTC", loc1,loc2) # 
J.f_og(DCTC)

PCostTech = J.diff("IOutput/PCostTech", loc1,loc2) # Has an issue with Tech Storage!!
J.f_og(PCostTech)

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
  
  PHEG0 = J.diff("EGOutput/PHEG0", loc1, loc2)
  @rsubset PHEG0 :Node == "QC" :Year == 2023
  
  PHPDP0 = J.diff("EGOutput/PHPDP0", loc1, loc2)
  @rsubset PHPDP0 :Node == "QC" :Year == 2023

  J.var("EGOutput/PHPDP0", T_Ref24)
  PHPDP1 = J.diff("EGOutput/PHPDP1", S_Ref24, T_Ref24)
  @rsubset PHPDP1 :Node == "QC" :Year == 2023
  
  sum(PHPDP1.S_Ref24)
  sum(PHPDP1.T_Ref24)

  BaseAdj = J.diff("SCalDB/BaseAdj", S_Ref24, T_Ref24); J.f_on(BaseAdj)
  xPkLoad = J.diff("SInput/xPkLoad", S_Ref24, T_Ref24); J.f_on(xPkLoad)
  xPkSavECC = J.diff("SInput/xPkSavECC", S_Ref24, T_Ref24); J.f_on(xPkSavECC)
  PkLoad = J.diff("SOutput/PkLoad", S_Ref24, T_Ref24); J.f_on(PkLoad)
  TDEF = J.diff("SInput/TDEF", S_Ref24, T_Ref24); J.f_on(TDEF)
  CgNetPk = J.diff("SOutput/CgNetPk", S_Ref24, T_Ref24); J.f_on(CgNetPk)
  HoursPerMonth = J.diff("SInput/HoursPerMonth", S_Ref24, T_Ref24); J.f_on(HoursPerMonth)
  CgLDCECC = J.diff("SOutput/CgLDCECC", S_Ref24, T_Ref24); J.f_on(CgLDCECC)
  CgLDCSoldECC = J.diff("SOutput/CgLDCSoldECC", S_Ref24, T_Ref24); J.f_on(CgLDCSoldECC)
  J.f_on(xPkSavECC)
  sum(xPkSavECC.S_Ref24)
  sum(xPkSavECC.T_Ref24)

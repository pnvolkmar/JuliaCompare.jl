using Revise
import JuliaCompare as J
import PromulaDBA as P
import SmallModel as M
import JuliaCompare: db_files, Canada
using CSV, DataFrames, DataFramesMeta

################################################################################
# Inputs for file locations ####################################################
################################################################################
BASE_FOLDER = raw"\\Silver\c\2020CanadaSpruce"
BASE_FOLDER2 = raw"\\Silver\c\2020CanadaTanoak"
SCENARIO1 = "Ref24"
SCENARIO2 = "Ref24"
################################################################################

CODE_FOLDER = joinpath(BASE_FOLDER, "Engine")
DATA_FOLDER1 = joinpath(BASE_FOLDER, "2020Model")
DATA_FOLDER2 = joinpath(BASE_FOLDER2, "2020Model")

DATA_FOLDER1 = SCENARIO1 == "" ? DATA_FOLDER1 : joinpath(DATA_FOLDER1,SCENARIO1)
DATA_FOLDER2 = SCENARIO2 == "" ? DATA_FOLDER2 : joinpath(DATA_FOLDER2,SCENARIO1)
E2020_Folder = joinpath(BASE_FOLDER, "2020Model")
HDF5_path = joinpath(DATA_FOLDER2, "database.hdf5")

vars = J.list_vars(CODE_FOLDER, DATA_FOLDER1, db_files);
vars_j = J.list_vars(HDF5_path)
loc1 = J.Loc_p(vars, DATA_FOLDER1, "Spruce");
loc2 = J.Loc_j(vars_j, HDF5_path, "Tanoak");

################################################################################
# Initialize variables for analysis ############################################
################################################################################
sec = 'E'
dimension_filters = Dict{Symbol,Any}()

CalibLTime = J.diff_fast("CalibLTime", loc1, loc2; dimension_filters, sec)
maximum(CalibLTime.Spruce)
maximum(CalibLTime.Tanoak)
EGCalSw = J.diff_fast("EGCalSw", loc1, loc2; dimension_filters, sec)
@rsubset EGCalSw :Spruce == 1
push!(dimension_filters, :Area => ["ON"])
push!(dimension_filters, :Node => ["ON"])
push!(dimension_filters, :ECC => "UtilityGen")
push!(dimension_filters, :Fuel => "NaturalGas")
push!(dimension_filters, :FuelEP => "NaturalGas")
push!(dimension_filters, :Year => string.(2021:2025))
push!(dimension_filters, :Poll => "CO2")
push!(dimension_filters, :Unit => "ON_New_021")

zoom = copy(dimension_filters)
push!(dimension_filters, :Year => "2024", :Month => "Summer", :TimeP => "TimeP6")

################################################################################
# Analysis of variables ########################################################
################################################################################
TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters, sec)
J.plot_diff(TotPol; dim="Poll", num=10, title="TotPol diffs by Poll")
J.plot_diff(TotPol; dim="Area", num=10, title="TotPol diffs by Area")


EuFPol = J.diff_fast("EuFPol", loc1, loc2; dimension_filters, sec)
@rsubset! EuFPol abs(:Diff) > 1e-46 :Year >= 2020
J.plot_diff(EuFPol; dim="Area", num=2, title="EuFPol diffs by Area")
J.plot_diff(EuFPol; dim="FuelEP", num=2, title="EuFPol diffs by FuelEP")
@by EuFPol [:Area, :Year] :Diff = sum(:Diff)
J.add_pdiff!(EuFPol)

i = findall(vars.Variable .== "UnPolGross")
i = i[1]
vars[i,:Database] = "EGOutput3"
pop!(dimension_filters, :Unit)
UnPolGross = J.diff_fast("UnPolGross", loc1, loc2; dimension_filters, sec)
@rsubset! UnPolGross :Unit != "Null"
UnArea = M.ReadDisk(loc2.HDF5_path, "EGInput/UnArea")
UnCode = M.ReadDisk(loc2.HDF5_path, "EGInput/UnCode")
UnArea = DataFrame(Unit = UnCode, UnArea = UnArea)
@rsubset! UnArea :Unit != "Null"
UnPlant = M.ReadDisk(loc2.HDF5_path, "EGInput/UnPlant")
UnPlant = DataFrame(Unit = UnCode, UnPlant = UnPlant)
@rsubset! UnPlant :Unit != "Null"

sort(UnPolGross, :Diff)

UnEG = J.diff_fast("UnEG", loc1, loc2; dimension_filters, sec) # bad
@rsubset! UnEG :Unit != "Null"
leftjoin!(UnEG, UnArea, on = :Unit)
leftjoin!(UnEG, UnPlant, on = :Unit)
@rsubset! UnEG :UnArea == "ON" abs(:Diff) >= 1e-5
sort!(UnEG, :Diff)
J.add_pdiff!(UnEG)

UnGC = J.diff_fast("UnGC", loc1, loc2; dimension_filters, sec) # bad
@rsubset! UnGC :Unit != "Null"
leftjoin!(UnGC, UnArea, on = :Unit)
leftjoin!(UnGC, UnPlant, on = :Unit)
@rsubset! UnGC :UnArea == "ON" abs(:Diff) >= 1e-5
sort!(UnGC, :Diff)
J.add_pdiff!(UnGC)

UnEAV = J.diff_fast("UnEAV", loc1, loc2; dimension_filters, sec) # bad
@rsubset! UnEAV :Unit != "Null"
leftjoin!(UnEAV, UnArea, on = :Unit)
leftjoin!(UnEAV, UnPlant, on = :Unit)
@rsubset! UnEAV :UnArea == "ON" :Spruce != 0 && :Tanoak != 0
sort!(UnEAV, :Diff)
J.add_pdiff!(UnEAV)

function get_me(Str)
  df = J.diff_fast(Str, loc1, loc2; dimension_filters, sec) # bad
  @rsubset! df :Unit != "Null"
  leftjoin!(df, UnArea, on = :Unit)
  leftjoin!(df, UnPlant, on = :Unit)
  @rsubset! df :UnArea == "ON" abs(:Diff) >= 1e-5
  sort!(df, :Diff)
  J.add_pdiff!(df)
  return df
end


#  UnEAV[unit,timep,month] = UnGC[unit]*(1-UnOUREG[unit])*UnEAF[unit,month]*(1-UnOOR[unit])*HoursPerMonth[month]/1000
UnOOR = get_me("UnOOR") # bad
@by(UnOOR, [:UnPlant, :Year], :Diff = sum(:Diff)/length(:Diff))
UnOR = get_me("UnOR") # good
HDADPwithStorage = J.diff_fast("HDADPwithStorage", loc1, loc2; dimension_filters, sec)
UnEG_slice
@by(UnEG_slice, :Year, :Spruce = sum(:Spruce), :Tanoak = sum(:Tanoak), :Diff = sum(abs.(:Diff)))

pop!(dimension_filters, :Node)
LLVC = J.diff_fast("LLVC", loc1, loc2; dimension_filters, sec)
@rsubset! LLVC abs(:Diff) >= 1e-6
LLMax = J.diff_fast("LLMax", loc1, loc2; dimension_filters, sec)
@rsubset! LLMax abs(:Diff) >= 1e-6



# we're seeing the small numbers issue with UnPolGross
# UnPolGross[unit,fuelep,poll] = UnDmd[unit,fuelep]*UnPOCA[unit,fuelep,poll]
UnDmd = J.diff_fast("UnDmd", loc1, loc2; dimension_filters, sec) # bad
UnPOCA = J.diff_fast("UnPOCA", loc1, loc2; dimension_filters, sec) # good
#    UnDmd[unit,fuelep] = max(UnEGGross[unit]*UnHRt[unit]/1e6*UnFlFr[unit,fuelep],0)
UnEGGross = J.diff_fast("UnEGGross", loc1, loc2; dimension_filters, sec) # bad
UnHRt = J.diff_fast("UnHRt", loc1, loc2; dimension_filters, sec) # good
UnFlFr = J.diff_fast("UnFlFr", loc1, loc2; dimension_filters, sec) # good
#     @finite_math UnEGGross[unit] = UnEGA[unit]/(1-UnOUREG[unit])
UnEGA = J.diff_fast("UnEGA", loc1, loc2; dimension_filters, sec) # bad
UnOUREG = J.diff_fast("UnOUREG", loc1, loc2; dimension_filters, sec) # good
UUnEGA = J.diff_fast("UUnEGA", loc1, loc2; dimension_filters, sec) # bad
# UUnEGA[unit] = sum(UnEG[unit,timep,month]+UnStorCurtailed[unit,timep,month]
#                    for month in Months, timep in TimePs)
@rsubset UnEG abs(:Diff) >= 1e-5
UnStorCurtailed = J.diff_fast("UnStorCurtailed", loc1, loc2; dimension_filters, sec) # good
@rsubset UnStorCurtailed :Diff != 0

HDPDP = J.diff_fast("HDPDP", loc1, loc2; dimension_filters, sec)
@rsubset! HDPDP :Spruce != 0 && :Tanoak != 0
J.plot_diff(HDPDP; dim="TimeP", num=10, title="HDPDP diffs by TimeP")
J.add_pdiff!(HDPDP)

PkLoad = J.diff_fast("PkLoad", loc1, loc2; dimension_filters, sec = 'E')
J.add_pdiff!(PkLoad)
J.plot_diff(PkLoad; dim="Area", num=10, title="PkLoad diffs by Area")
PkLoad_yr = @by PkLoad [:Area, :Year] begin
  :Diff = sum(:Diff)
end
@rsubset PkLoad_yr abs(:Diff) > 1 :Year > 2020
# let's focus on ON
push!(dimension_filters, :Node => "ON", :Area => "ON", :Year => string.(2020:2027))
# SLDC[hour,day,month,area] = sum(LDCECC[ecc,hour,day,month,area] for ecc in ECCs)
# /TDEF[electric,area]
SLDC = J.diff_fast("SLDC", loc1, loc2; dimension_filters, sec = 'E')
J.add_pdiff(SLDC)
pop!(dimension_filters, :ECC)
LDCECC = J.diff_fast("LDCECC", loc1, loc2; dimension_filters, sec = 'E')
J.add_pdiff!(LDCECC)
@rsubset! LDCECC abs(:PDiff) > .01
J.plot_diff(LDCECC; dim = "ECC", num = 10, title = "LDCECC diffs by ECC")
TDEF = J.diff_fast("TDEF", loc1, loc2; dimension_filters, sec = 'E')
@rsubset TDEF :Diff != 0
J.plot_diff(LDCECC; dim="ECC", num=10, title="LDCECC diffs by ECC")
push!(dimension_filters, :ECC => "SingleFamilyDetached")
push!(dimension_filters, :EC => "SingleFamilyDetached")

# LDCECC  = SaEC[ecc,area]*CLSF[class,hour,average,month,area]/8760*1000
SaEC = J.diff_fast("SaEC", loc1, loc2; dimension_filters, sec = 'E')

# ElecDmd[ecc,area] = sum(ESales[enduse,ec,area] for enduse in Enduses)
# SaEC[ecc,area] = ElecDmd[ecc,area]-CgEC[ecc,area]+PSoECC[ecc,area]
ESales = J.diff_fast("ESales", loc1, loc2; dimension_filters, sec = 'R')
CgEC = J.diff_fast("CgEC", loc1, loc2; dimension_filters, sec = 'R')
PSoECC = J.diff_fast("PSoECC", loc1, loc2; dimension_filters, sec = 'R')

#     CgEC[ecc,area] = sum(CgGen[fuel,ecc,area] for fuel in Fuels)
CgGen = J.diff_fast("CgGen", loc1, loc2; dimension_filters, sec = 'I')
@rsubset CgGen abs(:Diff) > 0.01


Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters, sec = 'R')
@rsubset! Polute abs(:Diff) >= 1e-1
@rsubset Polute :Year == 2027
sec = 'T'
Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters, sec)
J.plot_diff(Polute; dim="Tech", num=10, title="Polute diffs by Tech")
push!(dimension_filters, :Tech => ["LDVDiesel","LDVGasoline"])
push!(dimension_filters, :Year => "2030")

POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters, sec)
i = findall(vars.Variable .== "POCA" .&& first.(vars.Database) .== 'T')
i = i[1]
vars[i,:Database] = "TOutput2"
push!(dimension_filters, :EC => "Passenger")
POCA = J.diff_fast("POCA", loc1, loc2; dimension_filters, sec)
@rsubset POCA :Diff != 0

DmdFEPTech = J.diff_fast("DmdFEPTech", loc1, loc2; dimension_filters, sec)
@rsubset DmdFEPTech isapprox(:Diff, 0)

DmdFuelTech = J.diff_fast("DmdFuelTech", loc1, loc2; dimension_filters, sec)
dimension_filters[:FuelEP] = "Gasoline"
dimension_filters[:Tech] = "LDVGasoline"
push!(dimension_filters, :Fuel => "Gasoline")
DmdFuelTech = J.diff_fast("DmdFuelTech", loc1, loc2; dimension_filters, sec)
sec = 'R'
push!(dimension_filters, :Enduse => "Heat", :Fuel => "Electric")
DmFrac = J.diff_fast("DmFrac", loc1, loc2; dimension_filters, sec)
@rsubset DmFrac abs(:Diff) > 1e-45

DmFracMSM0 = J.diff_fast("DmFracMSM0", loc1, loc2; dimension_filters, sec)
DmFrac_T = select(DmFrac, [:Enduse,:Fuel,:Tech,:EC,:Area,:Year,:Tanoak])
DmFracMSM0_T = select(DmFracMSM0, [:Enduse,:Fuel,:Tech,:EC,:Area,:Year,:Tanoak])
rename!(DmFrac_T, :Tanoak => :DmFrac)
rename!(DmFracMSM0_T, :Tanoak => :DmFracMSM0)
DmFrac_T.DmFracMSM0 = DmFracMSM0_T.DmFracMSM0
@rsubset DmFrac_T :DmFrac > 1e-46 :Year == 2020
outer_join
@rsubset! DmFrac abs(:Diff) >= 1e-45

Dmd = J.diff_fast("Dmd", loc1, loc2; dimension_filters, sec)
@rsubset! Dmd abs(:Diff) > 0.001
push!(dimension_filters, :Tech => "Gas")
# Modify your DataFrame directly
J.add_pdiff!(Dmd)

push!(dimension_filters, :Year => "2027")
DER = J.diff_fast("DER", loc1, loc2; dimension_filters, sec)
@rsubset! DER abs(:Diff) > 0.001
# Modify your DataFrame directly
J.add_pdiff!(DER)

EEImpact = J.diff_fast("EEImpact", loc1, loc2; dimension_filters, sec)
EESat = J.diff_fast("EESat", loc1, loc2; dimension_filters, sec)
xDmd = J.diff_fast("xDmd", loc1, loc2; dimension_filters, sec)
DSMEU = J.diff_fast("DSMEU", loc1, loc2; dimension_filters, sec)
SqDmd = J.diff_fast("SqDmd", loc1, loc2; dimension_filters, sec)
SqEUTechMap = J.diff_fast("SqEUTechMap", loc1, loc2; dimension_filters, sec)
RPEI = J.diff_fast("RPEI", loc1, loc2; dimension_filters, sec)
CERSM = J.diff_fast("CERSM", loc1, loc2; dimension_filters, sec)
TSLoad = J.diff_fast("TSLoad", loc1, loc2; dimension_filters, sec)
vars[2625,:Database] = "ROutput2"
UMS = J.diff_fast("UMS", loc1, loc2; dimension_filters, sec)
DDay = J.diff_fast("DDay", loc1, loc2; dimension_filters, sec)
CUF = J.diff_fast("CUF", loc1, loc2; dimension_filters, sec)
WCUF = J.diff_fast("WCUF", loc1, loc2; dimension_filters, sec)
DDCoefficient = J.diff_fast("DDCoefficient", loc1, loc2; dimension_filters, sec)
DDayNorm = J.diff_fast("DDayNorm", loc1, loc2; dimension_filters, sec)

sec = 'R'
DER = J.diff_fast("DER", loc1, loc2; dimension_filters, sec)
@rsubset DER abs(:Diff) >= 1
J.add_pdiff!(DER)
DERRRExo = J.diff_fast("DERRRExo", loc1, loc2; dimension_filters, sec)
J.add_pdiff!(DERRRExo)
@rsubset DERRRExo abs(:Diff) >= 1
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters, sec)
J.add_pdiff!(PERRRExo)
@rsubset PERRRExo abs(:Diff) >= 1
RetroSwExo = J.diff_fast("RetroSwExo", loc1, loc2; dimension_filters, sec)
#     DER[enduse,tech,ec,area] = sum(DERV[enduse,tech,ec,area,vintage] for vintage in Vintages)


DERV = J.diff_fast("DERV", loc1, loc2; dimension_filters, sec)
i = findall(vars.Variable .∈ Ref(["DERV", "DERAV"]) .&& first.(vars.Database) .== 'R')
i = i[1]
vars[i,:Database] = "ROutput2"
vars[i,:]
DERV = J.diff_fast("DERV", loc1, loc2; dimension_filters, sec)

    # DERV[enduse,tech,ec,area,vintage] =
    #   DERV[enduse,tech,ec,area,vintage]*(1+StockAdjustment[enduse,tech,ec,area])

StockAdjustment = J.diff_fast("StockAdjustment", loc1, loc2; dimension_filters, sec)

@rsubset! DERV abs(:Diff) >= 1e-5
J.add_pdiff!(DERV)
@by DERV :Vintage begin
  :Diff = sum(:Diff)
  :PDiff = sum(:PDiff)
end
@rsubset DERV :Vintage == "Vintage 1"

    # DERV[enduse,tech,ec,area,vintage] = DERV[enduse,tech,ec,area,vintage]+DT*
    #   (DERAV[enduse,tech,ec,area,vintage]-DERRV[enduse,tech,ec,area,vintage])

DERAV = J.diff_fast("DERAV", loc1, loc2; dimension_filters, sec)
dimension_filters[:Year] = string.(2020:2040)
DERA   = J.diff_fast("DERA", loc1, loc2; dimension_filters, sec)
@rsubset! DERA abs(:Diff) >= 1e-5
J.add_pdiff!(DERA)
sum(DERA.PDiff) # -6.18

# DERA=DERAPC+DERAP+DERAD+DERARC

DERAPC = J.diff_fast("DERAPC", loc1, loc2; dimension_filters, sec)
@rsubset! DERAPC abs(:Diff) >= 1e-5
J.add_pdiff!(DERAPC)
sum(DERAPC.PDiff) # -7.14

DERAP = J.diff_fast("DERAP", loc1, loc2; dimension_filters, sec)
@rsubset! DERAP abs(:Diff) >= 1e-5
J.add_pdiff!(DERAP)
sum(DERAP.PDiff) # 7.08

DERAD = J.diff_fast("DERAD", loc1, loc2; dimension_filters, sec)
@rsubset! DERAD abs(:Diff) >= 1e-5
J.add_pdiff!(DERAD)
sum(DERAD.PDiff) # 0

DERARC = J.diff_fast("DERARC", loc1, loc2; dimension_filters, sec)
@rsubset! DERARC abs(:Diff) >= 1e-5
J.add_pdiff!(DERARC)
sum(DERARC.PDiff) # .73

# DERAPC is the largest contributor to the difference in DERA in 2027
# @. @finite_math DERAPC = (PERAPC+PERADSt)/DEE
PERAPC = J.diff_fast("PERAPC", loc1, loc2; dimension_filters, sec)
@rsubset! PERAPC abs(:Diff) >= 1e-5
J.add_pdiff!(PERAPC)
sum(PERAPC.PDiff) # -5.59

PERADSt = J.diff_fast("PERADSt", loc1, loc2; dimension_filters, sec)
@rsubset! PERADSt abs(:Diff) >= 1e-5
J.add_pdiff!(PERADSt)
sum(PERADSt.PDiff) # -3.22

DEE = J.diff_fast("DEE", loc1, loc2; dimension_filters, sec)
@rsubset! DEE abs(:Diff) >= 1e-4
J.add_pdiff!(DEE)
sum(DEE.PDiff) # 11.15 but it's mainly an issue for Heat/HW Enduse with DuelHPump Tech

# PERAPC
# New = Select(Age,"New")
# for area in Areas,ec in ECs,enduse in Enduses,tech in Techs
#   @finite_math PERAPC[enduse,tech,ec,area] = EUPCAPC[enduse,tech,New,ec,area]*
#     DSt[enduse,ec,area]/PEE[enduse,tech,ec,area]
# end
push!(dimension_filters, :Age => "New")

# EUPCAPC
EUPCAPC = J.diff_fast("EUPCAPC", loc1, loc2; dimension_filters, sec)
@rsubset! EUPCAPC abs(:Diff) >= 1e-4
J.add_pdiff!(EUPCAPC)
sum(EUPCAPC.PDiff) # -6.03

# DSt
DSt = J.diff_fast("DSt", loc1, loc2; dimension_filters, sec)
@rsubset! DSt abs(:Diff) >= 1e-4
J.add_pdiff!(DSt)
sum(DSt.PDiff) # 0

# PEE
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters, sec)
@rsubset! PEE abs(:Diff) >= 1e-7
J.add_pdiff!(PEE)
sum(PEE.PDiff) #

# EUPCAPC
  # EUPCAPC[enduse,tech,New,ec,area] = PCA[New,ecc,area]*MMSF[enduse,tech,ec,area]
# PCA
PCA = J.diff_fast("PCA", loc1, loc2; dimension_filters, sec)
@rsubset! PCA abs(:Diff) >= 1e-7
J.add_pdiff!(PCA)
sum(PCA.PDiff) #

# MMSF
MMSF = J.diff_fast("MMSF", loc1, loc2; dimension_filters, sec)
@rsubset! MMSF abs(:Diff) != 0
J.add_pdiff!(MMSF)
sum(MMSF.PDiff) # 6.03

# @finite_math MMSF[enduse,tech,ec,area] = MAW[enduse,tech,ec,area]/TMAW[enduse,ec,area]

  # MAW[enduse,tech,ec,area] = exp(MMSM0[enduse,tech,ec,area]+
  #   log(MSMM[enduse,tech,ec,area])+
  #   MMSMI[enduse,tech,ec,area]*(SPC[ec,area]/SPop[ec,area])/
  #                               (SPC0[ec,area]/SPop0[ec,area])+
  #   MVF[enduse,tech,ec,area]*
  #   log((MCFU[enduse,tech,ec,area]/Inflation[area]/PEE[enduse,tech,ec,area])/
  #   (MCFU0[enduse,tech,ec,area]/Inflation0[area]/PEE0[enduse,tech,ec,area])))
MMSM0 = J.diff_fast("MMSM0", loc1, loc2; dimension_filters, sec)
@rsubset! MMSM0 abs(:Diff) != 0
J.add_pdiff!(MMSM0)
@rsubset! MMSM0 :Spruce > -170 || :Tanoak > -170
sum(MMSM0.PDiff) # 6.03

MCFU = J.diff_fast("MCFU", loc1, loc2; dimension_filters, sec)
@rsubset! MCFU abs(:Diff) != 0
J.add_pdiff!(MCFU)
@rsubset! MCFU :Spruce > -170 || :Tanoak > -170
sum(MCFU.PDiff) # 83

MSMM = J.diff_fast("MSMM", loc1, loc2; dimension_filters, sec)
@rsubset! MSMM abs(:Diff) != 0

    # @finite_math MCFU[enduse,tech,ec,area] =
    #   DCCR[enduse,tech,ec,area]*DCC[enduse,tech,ec,area]+DOMC[enduse,tech,ec,area]+
    #   ECFP[enduse,tech,ec,area]/
    #   DEE[enduse,tech,ec,area]+
    #   IdrtCost[enduse,tech,ec,area]*Inflation[area]

# ECFP
ECFP = J.diff_fast("ECFP", loc1, loc2; dimension_filters, sec)
@rsubset! ECFP abs(:Diff) != 0
J.add_pdiff!(ECFP)
@rsubset! ECFP :Spruce > -170 || :Tanoak > -170
sum(ECFP.PDiff) # 8.6

# DOMC = DOCF*DCCFullCost
DCCFullCost = J.diff_fast("DCCFullCost", loc1, loc2; dimension_filters, sec)
@rsubset! DCCFullCost abs(:Diff) != 0
J.add_pdiff!(DCCFullCost)
@rsubset! DCCFullCost :Spruce > -170 || :Tanoak > -170
sum(DCCFullCost.PDiff) # 8.6

DOCF = J.diff_fast("DOCF", loc1, loc2; dimension_filters, sec)
@rsubset! DOCF abs(:Diff) != 0
J.add_pdiff!(DOCF)
@rsubset! DOCF :Spruce > -170 || :Tanoak > -170
sum(DOCF.PDiff) # 0

DCC = J.diff_fast("DCC", loc1, loc2; dimension_filters, sec)
@rsubset! DCC abs(:Diff) != 0
J.add_pdiff!(DCC)
@rsubset! DCC :Spruce > -170 || :Tanoak > -170
sum(DCC.PDiff) # 8.6

DCCR = J.diff_fast("DCCR", loc1, loc2; dimension_filters, sec)
@rsubset! DCCR abs(:Diff) != 0
J.add_pdiff!(DCCR)
@rsubset! DCCR :Spruce > -170 || :Tanoak > -170
sum(DCCR.PDiff) # 0


# Dict{Symbol, Any} with 6 entries:
#   :Area => "ON"
#   :Node => "ON"
#   :ECC  => "SingleFamilyDetached"
#   :Year => "2027"
#   :Age  => "New"
#   :EC   => "SingleFamilyDetached"
small = copy(dimension_filters)
push!(small, Enduse)


# Another issue is with outage rates in NS Coal Units:
sec = 'E'
ns_issue = Dict{Symbol,Any}()
push!(ns_issue, :Unit => [
  "NS00008004102",
  "NS00008004101",
  "NS00008004001",
  "NS00008003901"
])
push!(ns_issue, :Node => "NS", :Area => "NS")

OOR = J.diff_fast("UnOOR", loc1, loc2; dimension_filters = ns_issue, sec)
@rsubset! OOR abs(:Diff) > .01

UnPlant = M.ReadDisk(loc2.HDF5_path, "EGInput/UnPlant")
UnNode = M.ReadDisk(loc2.HDF5_path, "EGInput/UnNode")
UnCode = M.ReadDisk(loc2.HDF5_path, "EGInput/UnCode")
UnInfo = DataFrame(UnCode = UnCode, UnPlant = UnPlant, UnNode = UnNode)
@rsubset! UnInfo :UnPlant == "BaseHydro" :UnNode == "ON"

on_bh = Dict{Symbol,Any}()
push!(on_bh, :Unit => UnInfo.UnCode, :Year => ["2024", "2025"])
UnGC = J.diff_fast("UnGC", loc1, loc2, dimension_filters = on_bh, sec = 'E')
@rsubset! UnGC !isapprox(:Diff, 0)

push!(on_bh, :Node => "ON", :Area => "ON", :Plant => "BaseHydro")
vGC = J.diff_fast("GC", loc1, loc2, dimension_filters = on_bh, sec = 'E')
@rsubset! vGC !isapprox(:Diff, 0)

# GCap[plant,genco]
GCap = J.diff_fast("GCap", loc1, loc2, dimension_filters = on_bh, sec = 'E')
@rsubset! GCap !isapprox(:Diff, 0)

Capacity = J.diff_fast("Capacity", loc1, loc2, dimension_filters = on_bh, sec = 'E')
@rsubset! Capacity !isapprox(:Diff, 0)

xCapacity = J.diff_fast("xCapacity", loc1, loc2, dimension_filters = on_bh, sec = 'E')
@rsubset! xCapacity !isapprox(:Diff, 0)

push!(on_bh, :Area => ["ON", "NL"])
xRnImports = J.diff_fast("xRnImports", loc1, loc2, dimension_filters = on_bh, sec = 'E')
@rsubset! xRnImports !isapprox(:Diff, 0)



# DERARC is the largest contributor to the difference in DERA in 2027

df2020 = copy(dimension_filters)
df2020[:Year] = "2020"
df2020[:Year] = ["2019","2020"]

DERRRC = J.diff_fast("DERRRC", loc1, loc2; dimension_filters = df2020, sec)
CMSF = J.diff_fast("CMSF", loc1, loc2; dimension_filters = df2020, sec)
DEE = J.diff_fast("DEE", loc1, loc2; dimension_filters = df2020, sec)
DEEA = J.diff_fast("DEEA", loc1, loc2; dimension_filters = df2020, sec)
CnvrtEU = J.diff_fast("CnvrtEU", loc1, loc2; dimension_filters = df2020, sec)
xCMSF = J.diff_fast("xCMSF", loc1, loc2; dimension_filters = df2020, sec)
DEEA = J.diff_fast("DEEA", loc1, loc2; dimension_filters = df2020, sec)
df2020[:Tech] = ["LDVGasoline","LDVDiesel"]
MMSF = J.diff_fast("MMSF", loc1, loc2; dimension_filters = df2020, sec)
MMSM0 = J.diff_fast("MMSM0", loc1, loc2; dimension_filters, sec)

dfDiesel = copy(dimension_filters)
dfDiesel[:Tech] = "LDVDiesel"
MMSM0 = J.diff_fast("MMSM0", loc1, loc2; dimension_filters, sec)
# MMSF is off

CMSM0 = J.diff_fast("CMSM0", loc1, loc2; dimension_filters = df2020, sec) # good
@rsubset! CMSM0 abs(:Diff) > 1e-7
CVF = J.diff_fast("CVF", loc1, loc2; dimension_filters = df2020, sec) # good
@rsubset! CVF abs(:Diff) > 1e-7
MCFU = J.diff_fast("MCFU", loc1, loc2; dimension_filters = df2020, sec) # good
DCCR = J.diff_fast("DCCR", loc1, loc2; dimension_filters = df2020, sec) # good
FDCC = J.diff_fast("FDCC", loc1, loc2; dimension_filters = df2020, sec) # good
FDCCU = J.diff_fast("FDCCU", loc1, loc2; dimension_filters = df2020, sec) # good
@rsubset! FDCCU :Diff != 0
MCFU0 = J.diff_fast("MCFU", loc1, loc2; dimension_filters, sec) # good
CMSMI = J.diff_fast("CMSMI", loc1, loc2; dimension_filters = df2020, sec) # good
SPC = J.diff_fast("SPC", loc1, loc2; dimension_filters = df2020, sec) # didn't work
df2020[:ECC] = "Passenger"
SPC = J.diff_fast("PC", loc1, loc2; dimension_filters = df2020, sec) # didn't work
SPC0 = J.diff_fast("PC", loc1, loc2; dimension_filters = dfFirst, sec='M') # didn't work
SPop = J.diff_fast("Pop", loc1, loc2; dimension_filters = df2020, sec='M') # didn't work
SPop0 = J.diff_fast("Pop", loc1, loc2; dimension_filters = dfFirst, sec='M') # didn't work
dfFirst = copy(df2020)
dfFirst[:Year] = "1986"

xMMSF = J.diff_fast("xMMSF", loc1, loc2; dimension_filters = df2020, sec) # didn't work
CMMSF = J.diff_fast("CMM", loc1, loc2; dimension_filters = df2020, sec) # didn't work

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
tan_CInd = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaTanoak\\2020Model\\CalibInd\\database.hdf5", "Tan_CInd");
spr_CInd = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaSpruce\\2020Model\\CalibInd", "Spr_CInd");

################################################################################
# Initialize filters and sec for analysis ######################################
################################################################################
sec = 'I'
dimension_filters = Dict{Symbol,Any}()
push!(dimension_filters, :Area => Canada, :Year => string.(1986:2050))

off = copy(dimension_filters)
push!(off, :ECC => "OnFarmFuelUse", :EC => "OnFarmFuelUse")
push!(off, :Year => string.(1985:2050))
push!(off, :Area => ["AB"])
# push!(off, :FuelEP => ["Diesel","Biodiesel"], :Fuel => ["Diesel","Biodiesel"])
push!(off, :Tech => ["Electric","Coal"])
push!(off, :Enduse => ["OthNSub","Miscellaneous"])

################################################################################
# Analysis of variables: Overview ##############################################
################################################################################

Dmd = J.diff_fast("Dmd", loc1, loc2; dimension_filters=off, sec) # 
J.plot_diff(Dmd, dim = :Enduse)
@rsubset! Dmd :Diff != 0 :Tech == "OffRoad"
J.add_pdiff!(Dmd)
xDmd = J.diff_fast("xDmd", loc1, loc2; dimension_filters=off, sec) # 
unique(xDmd.Enduse)
vDmd = J.diff_fast("vDmd", loc1, loc2; dimension_filters=off, sec) # 
unique(vDmd.vEnduse)
vEUMap = J.diff_fast("vEUMap", loc1, loc2; dimension_filters=off, sec) # 

J.plot_diff(xDmd, dim = :Enduse)
@rsubset! xDmd :Tech == "Electric" :Spruce != 0
J.add_pdiff!(xDmd)
DmFrac = J.diff_fast("DmFrac", loc1, loc2; dimension_filters=off, sec) # 
@rsubset DmFrac :Year == 2050
@rsubset! DmFrac :Diff != 0
J.add_pdiff!(DmFrac)
zoom = copy(off)
push!(zoom, :Area => "QC", :Year => "2050")
pop!(zoom, :Fuel)
DmFrac = J.diff_fast("DmFrac", loc1, loc2; dimension_filters=zoom, sec) # 
DmFrac = J.diff_fast("DmFrac", spr_base, tan_base; dimension_filters=off, sec) # 
@rsubset DmFrac :Area == "BC" :Fuel == "Diesel"

@rsubset! DmFrac :Tanoak != 0
xDmFrac = J.diff_fast("xDmFrac", loc1, loc2; dimension_filters=zoom, sec) # 
@rsubset! xDmFrac :Tanoak != 0
DmFracMin = J.diff_fast("DmFracMin", loc1, loc2; dimension_filters=off, sec) # 
@rsubset! DmFracMin :Tanoak != 0
@rsubset DmFracMin :Area == "QC" :Fuel == "Biodiesel"
@rsubset DmFracMin :Area == "BC" :Fuel == "Biodiesel"
DmFracMax = J.diff_fast("DmFracMax", loc1, loc2; dimension_filters=zoom, sec) # 
@rsubset DmFracMax :Tanoak != 0 :Diff != 0

# @. Dmd = Dmd-EE
EE = J.diff_fast("EE", loc1, loc2; dimension_filters=off, sec) # 
@rsubset
# @. Dmd = Dmd-DSMEU
DSMEU = J.diff_fast("DSMEU", loc1, loc2; dimension_filters=off, sec) # 
# Dmd[enduse,tech,ec,area] = Dmd[enduse,tech,ec,area]+
#                            SqDmd[tech,ec,area]*SqEUTechMap[enduse,tech]
SqDmd = J.diff_fast("SqDmd", loc1, loc2; dimension_filters=off, sec) # 
SqEUTechMap = J.diff_fast("SqEUTechMap", loc1, loc2; dimension_filters=off, sec) # 
#
# @finite_math Dmd[enduse,tech,ec,area] = DER[enduse,tech,ec,area]*
# UMS[enduse,tech,ec,area]*CERSM[enduse,ec,area]*
# CUF[enduse,tech,ec,area]*WCUF[ec,area]*RPEI[enduse,tech,ec,area]/1.0e6*
# (TSLoad[enduse,ec,area]*
# (DDay[enduse,area]/DDayNorm[enduse,area])^DDCoefficient[enduse,ec,area]+
# (1.0-TSLoad[enduse,ec,area]))
DER = J.diff_fast("DER", loc1, loc2; dimension_filters=off, sec) # Problem
J.add_pdiff!(DER)
UMS = J.diff_fast("UMS", loc1, loc2; dimension_filters=off, sec)
CERSM = J.diff_fast("CERSM", loc1, loc2; dimension_filters=off, sec)
CUF = J.diff_fast("CUF", loc1, loc2; dimension_filters=off, sec)
WCUF = J.diff_fast("WCUF", loc1, loc2; dimension_filters=off, sec)
RPEI = J.diff_fast("RPEI", loc1, loc2; dimension_filters=off, sec)
TSLoad = J.diff_fast("TSLoad", loc1, loc2; dimension_filters=off, sec)
DDay = J.diff_fast("DDay", loc1, loc2; dimension_filters=off, sec)
DDayNorm = J.diff_fast("DDayNorm", loc1, loc2; dimension_filters=off, sec)
DDCoefficient = J.diff_fast("DDCoefficient", loc1, loc2; dimension_filters=off, sec)
EEImpact = J.diff_fast("EEImpact", loc1, loc2; dimension_filters=off, sec)
EESat = J.diff_fast("EESat", loc1, loc2; dimension_filters=off, sec)
xDmd = J.diff_fast("xDmd", loc1, loc2; dimension_filters=off, sec)

DERRRExo = J.diff_fast("DERRRExo", loc1, loc2; dimension_filters=off, sec)
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters=off, sec)
PERRRExo = J.diff_fast("PERRRExo", loc1, loc2; dimension_filters=off, sec)
i = findall(vars.Variable .∈ Ref(["DERV", "DERAV"]) .&& first.(vars.Database) .== sec)
vars[i,:Database] .= sec*"Output2"
DERV = J.diff_fast("DERV", loc1, loc2; dimension_filters=off, sec) # Problem
select!(DERV, Not([:Enduse, :Tech]))
J.add_pdiff!(DERV)
@rsubset DERV :Vintage == "Vintage 1"
StockAdjustment = J.diff_fast("StockAdjustment", loc1, loc2; dimension_filters=off, sec) # fine
DERAV = J.diff_fast("DERAV", loc1, loc2; dimension_filters=off, sec)
DERA   = J.diff_fast("DERA", loc1, loc2; dimension_filters=off, sec) # problem
J.add_pdiff(DERA)
@rsubset DERA isnan(:Tanoak) :Year == 2020
#   @. DERA = DERAPC+DERAP+DERAD+DERARC
DERAPC = J.diff_fast("DERAPC", loc1, loc2; dimension_filters=off, sec) #  problem
J.add_pdiff(DERAPC)
DERAP = J.diff_fast("DERAP", loc1, loc2; dimension_filters=off, sec) # fine
DERAD = J.diff_fast("DERAD", loc1, loc2; dimension_filters=off, sec) # problem
J.add_pdiff(DERAD)
DERARC = J.diff_fast("DERARC", loc1, loc2; dimension_filters=off, sec) # fine

RetroSwExo = J.diff_fast("RetroSwExo", loc1, loc2; dimension_filters=off, sec)

#   @. @finite_math DERAPC = (PERAPC+PERADSt)/DEE
PERAPC = J.diff_fast("PERAPC", loc1, loc2; dimension_filters=off, sec) # big problem
J.add_pdiff!(PERAPC)
PERADSt = J.diff_fast("PERADSt", loc1, loc2; dimension_filters=off, sec) # 
DEE = J.diff_fast("DEE", loc1, loc2; dimension_filters=off, sec) # 
J.add_pdiff!(DEE)

# process additions from production capacity (PERAPC)
# @finite_math PERAPC[enduse,tech,ec,area] = EUPCAPC[enduse,tech,New,ec,area]*
#     DSt[enduse,ec,area]/PEE[enduse,tech,ec,area]
push!(off, :Age => "New")
EUPCAPC = J.diff_fast("EUPCAPC", loc1, loc2; dimension_filters=off, sec) # big problem
J.add_pdiff(EUPCAPC)
DSt = J.diff_fast("DSt", loc1, loc2; dimension_filters=off, sec) # 
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters=off, sec) # 
J.add_pdiff(PEE)
# EUPCAPC[enduse,tech,New,ec,area] = PCA[New,ecc,area]*MMSF[enduse,tech,ec,area]
PCA = J.diff_fast("PCA", loc1, loc2; dimension_filters=off, sec) # 
age = pop!(off, :Age)
PCA_mid = J.diff_fast("PCA", loc1, loc2; dimension_filters=off, sec) # 
# PCA[New,ecc,area] = (PCPrior[ecc,area]*(ECGR[ecc,area]+GRPGR[area]))+PCR[Old,ecc,area]
push!(off, :Age => "Old")
PC = J.diff_fast("PC", loc1, loc2; dimension_filters=off, sec) # about 1/4 of the problem
J.add_pdiff!(PC)
ECGR = J.diff_fast("ECGR", loc1, loc2; dimension_filters=off, sec) # a lot of problem
J.add_pdiff(ECGR)
GRPGR = J.diff_fast("GRPGR", loc1, loc2; dimension_filters=off, sec) # fine
PCR = J.diff_fast("PCR", loc1, loc2; dimension_filters=off, sec) # fine
J.add_pdiff(PCR)
#     @finite_math ECGR[ecc,area] = (PC[ecc,area]-PCPrior[ecc,area])/PCPrior[ecc,area]-GRPGR[area]
# PCMin = J.diff_fast("PCMin", loc1, loc2; dimension_filters=off, sec) # temp var
xPC = J.diff_fast("xPC", loc1, loc2; dimension_filters=off, sec) # problem
# @finite_math xPC[ecc,area] = Driver[ecc,area]/xECUF[ecc,area]   
Driver = J.diff_fast("Driver", loc1, loc2; dimension_filters=off, sec) # 
J.add_pdiff!(Driver)
sort!(Driver, :PDiff)
xECUF = J.diff_fast("xECUF", loc1, loc2; dimension_filters=off, sec) # 
# Driver[ecc,area] = Driver[ecc,area]*DriverMultiplier[ecc,area]
DriverMultiplier = J.diff_fast("DriverMultiplier", loc1, loc2; dimension_filters=off, sec) # fine
xDriver = J.diff_fast("xDriver", loc1, loc2; dimension_filters=off, sec) # issue
DrSwitch = J.diff_fast("DrSwitch", loc1, loc2; dimension_filters=off, sec) # 21
xGO = J.diff_fast("xGO", loc1, loc2; dimension_filters=off, sec) # issue
vGO = J.diff_fast("vGO", loc1, loc2; dimension_filters=off, sec) # doesn't prvide data for CN
    # xGO[ecc,area,year]  =sum(xGOECC[ecc,areatom,year]*
    #   MapAreaTOM[area,areatom] for areatom in areatoms)

# xGOECC = J.diff_fast("xGOECC", loc1, loc2; dimension_filters=off, sec) # temp var
MapECCfromTOM = J.diff_fast("KInput/MapECCfromTOM", loc1, loc2; dimension_filters=off, sec='K') # 
@rsubset MapECCfromTOM :Spruce != 0 || :Tanoak != 0
push!(off, :ECCfromTOM => ["CopperMining","GoldOreMining","OtherMetalMining"])
GY = J.diff_fast("KOutput/GY", loc1, loc2; dimension_filters=off, sec='K') # 
@rsubset GY :AreaTOM == "ON"
MMSF = J.diff_fast("MMSF", loc1, loc2; dimension_filters=off, sec) # problem
MMSFB = J.diff_fast("MMSF", spr_base, tan_base; dimension_filters=off, sec) # problem
J.add_pdiff!(MMSF)
MMSFB = J.diff_fast("MMSF", spr_base, tan_base; dimension_filters=off, sec) # problem
J.add_pdiff!(MMSFB)
# xMMSF is none zero fro Coal in Base and ProcSw is 0 so MMSF should get that value....but it is zero. 
@rsubset! MMSF :Spruce != 0 || :Tanoak != 0
xProcSw = J.diff_fast("xProcSw", loc1, loc2; dimension_filters=off, sec) # 
@rsubset xProcSw :PI == "MShare" # 0 is exogenous
xProcSw = J.diff_fast("xProcSw", spr_base, tan_base; dimension_filters=off, sec) # 
@rsubset xProcSw :PI == "MShare" # 0 is exogenous
xMMSF = J.diff_fast("xMMSF", loc1, loc2; dimension_filters=off, sec) # different
push!(off, :Year => string.(2025:2050))
push!(off, :Tech => ["Electric", "Coal"])
push!(off, :Area => "AB")
push!(off, :Enduse => "OthNSub")
@rsubset xMMSF :Year == 2050 :Diff != 0
xMMSFB = J.diff_fast("xMMSF", spr_base, tan_base; dimension_filters=off, sec) # different
df = @rsubset xMMSFB :Tech == "Coal"
lines(df.Year, df.Tan_Base)
lines(df.Year, df.Spr_Base)
xMMSFC = J.diff_fast("xMMSF", spr_CInd, tan_CInd; dimension_filters=off, sec) # different
df = @rsubset xMMSFC :Tech == "Coal"
# MMSF[enduse,tech,ec,area] = MMSF[enduse,tech,ec,area]+
# max(MMSF[enduse,tech,ec,area]-MMSFB[enduse,tech,ec,area],0)*
# ETSwitch[tech,area] 

# Create the figure and axis
fig = Figure(resolution = (800, 600))
ax = Axis(fig[1, 1], 
    xlabel = "Year", 
    ylabel = "Index Value",
    title = "xMMSF for AB, OthNSub, Coal, OnFarmFuelUse from CalibInd"
)

# Plot the two lines
lines!(ax, df.Year, df.Spr_CInd, label = "Spruce", color = :blue, linewidth = 2)
lines!(ax, df.Year, df.Tan_CInd, label = "Tanoak", color = :red, linewidth = 2)

# Add legend
axislegend(ax, position = :rt)

# Display the figure
fig

ETSwitch = J.diff_fast("ETSwitch", loc1, loc2; dimension_filters=off, sec) # 

#   for area in Areas,enduse in Enduses,ec in ECs,tech in Techs
#     if MCFU[enduse,tech,ec,area]  > 0.0 && 
#        MCFU0[enduse,tech,ec,area] > 0.0 &&
#        PEE[enduse,tech,ec,area]   > 0.0 && 
#        PEE0[enduse,tech,ec,area]  > 0.0 &&
#        SPC[ec,area]               > 0.0 && 
#        SPC0[ec,area]              > 0.0 && 
#        SPop[ec,area]              > 0.0 && 
#        SPop[ec,area]              > 0.0 && 
#        MMSM0[enduse,tech,ec,area] > -150.0 

#       MAW[enduse,tech,ec,area] = exp(MMSM0[enduse,tech,ec,area]+
#         log(MSMM[enduse,tech,ec,area])+
#         MMSMI[enduse,tech,ec,area]*(SPC[ec,area]/SPop[ec,area])/
#                                     (SPC0[ec,area]/SPop0[ec,area])+
#         MVF[enduse,tech,ec,area]*
#         log((MCFU[enduse,tech,ec,area]/Inflation[area]/PEE[enduse,tech,ec,area])/
#         (MCFU0[enduse,tech,ec,area]/Inflation0[area]/PEE0[enduse,tech,ec,area])))*
#         MSLimit[enduse,tech,ec,area]
#     else
#       MAW[enduse,tech,ec,area] = 0.0
#     end
#   end

#   for area in Areas,enduse in Enduses,ec in ECs,tech in Techs
#     TMAW[enduse,ec,area] = sum(MAW[enduse,tech,ec,area] for tech in Techs)
#     @finite_math MMSF[enduse,tech,ec,area] = MAW[enduse,tech,ec,area]/TMAW[enduse,ec,area]
#   end
pop!(off, :Tech)
push!(off, :Year => string.(2024:2025))
push!(off, :Tech => ["Electric", "Oil", "OffRoad"])
MMSM0 = J.diff_fast("MMSM0", loc1, loc2; dimension_filters=off, sec) # doesn't match
@rsubset MMSM0 :Spruce > -170 || :Tanoak > -170
MSMM = J.diff_fast("MSMM", loc1, loc2; dimension_filters=off, sec) # matches
MVF = J.diff_fast("MVF", loc1, loc2; dimension_filters=off, sec) # matches
MCFU = J.diff_fast("MCFU", loc1, loc2; dimension_filters=off, sec) # close
MMSMI = J.diff_fast("MMSMI", loc1, loc2; dimension_filters=off, sec) # 0
AMSF = J.diff_fast("AMSF", loc1, loc2; dimension_filters=off, sec) # doesn't match
PEE = J.diff_fast("PEE", loc1, loc2; dimension_filters=off, sec) # doesn't match
@rsubset PEE isnan(:Tanoak)
PEESw = J.diff_fast("PEESw", loc1, loc2; dimension_filters=off, sec) # 1
Enduse = pop!(off, :Enduse)
push!(off, :Enduse => Enduse)
PEEBeforeStd = J.diff_fast("PEEBeforeStd", loc1, loc2; dimension_filters=off, sec) # doesn't match
PEEBeforeStd = J.diff_fast("PEEBeforeStd", loc1, loc2; dimension_filters, sec) # doesn't match
J.add_pdiff!(PEEBeforeStd)
@rsubset! PEEBeforeStd :Year >= 2020 :Year <= 2030
sort!(PEEBeforeStd, :PDiff)
nas = copy(dimension_filters)
push!(nas, :Year => "2020")
push!(nas, :Tech => ["Gas", "Oil", "Coal"])
push!(nas, :Enduse => "Steam")
P_issues = @rsubset PEEBeforeStd isnan(:Tanoak) :Year == 2020 :Enduse == "Steam" 
unique(P_issues.Tech)
@rsubset! PEEBeforeStd :Spruce != 0 || :Diff != 0
J.add_pdiff!(PEEBeforeStd)
@rsubset! PEEBeforeStd :PDiff == 100
# ==> implies issue is before std with PEEPrice

#   PEEPrice[enduse,tech,ec,area] = PEM[enduse,ec,area]*
#     PEMM[enduse,tech,ec,area]*
#     (1/(1+(MCFU[enduse,tech,ec,area]/Inflation[area]*
#     PEPM[enduse,tech,ec,area]/PFPN[enduse,tech,ec,area])^PFTC[enduse,tech,ec,area]*
#     (PCCR[enduse,tech,ec,area]/PCCRB[enduse,tech,ec,area])))
PEEPrice = J.diff_fast("PEEPrice", loc1, loc2; dimension_filters=nas, sec) # issue
J.add_pdiff!(PEEPrice)
PEM = J.diff_fast("PEM", loc1, loc2; dimension_filters=nas, sec) # matches
PEMM = J.diff_fast("PEMM", loc1, loc2; dimension_filters=nas, sec) # matches
# PEEBeforeStd matches PEM = PEM*PEMM which means in Spruce, Electric and Oil are getting the else values
PEPM = J.diff_fast("PEPM", loc1, loc2; dimension_filters=nas, sec) # 1
PFPN = J.diff_fast("PFPN", loc1, loc2; dimension_filters=nas, sec) # matches
PFTC = J.diff_fast("PFTC", loc1, loc2; dimension_filters=nas, sec) # matches
PCCR = J.diff_fast("PCCR", loc1, loc2; dimension_filters=nas, sec) # matches
PCCRB = J.diff_fast("PCCR", spr_base, tan_base; dimension_filters=nas, sec) # matches
PCTC = J.diff_fast("PCTC", loc1, loc2; dimension_filters=nas, sec) # matches

df = leftjoin(PEEBeforeStd,PCTC, on=[:Enduse,:Tech,:EC,:Area,:Year],  renamecols = "_PEEB" => "_PCTC")
@rsubset df :Spruce_PCTC == 0
PFPN = J.diff_fast("PFPN", loc1, loc2; dimension_filters=nas, sec) # matchesl
df = leftjoin(df, PFPN, on = [:Enduse, :Tech, :EC, :Area], renamecols = "" => "_PFPN")
@rsubset! df :Enduse != "OffRoad"

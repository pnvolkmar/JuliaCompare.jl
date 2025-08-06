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
tan_c3 = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaTanoak\\2020Model\\Calib3\\database.hdf5", "Tanoak");
spr_c3 = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaSpruce\\2020Model\\Calib3", "Spruce");
tan_c4 = J.Loc_j(vars_j, "\\\\Pink\\c\\2020CanadaTanoak\\2020Model\\Calib4\\database.hdf5", "Tanoak");
spr_c4 = J.Loc_p(vars, "\\\\Pink\\c\\2020CanadaSpruce\\2020Model\\Calib4", "Spruce");

################################################################################
# Initialize filters and sec for analysis ######################################
################################################################################
sec = 'I'
dimension_filters = Dict{Symbol,Any}()
push!(dimension_filters, :Area => Canada, :Year => string.(1986:2050))

is = copy(dimension_filters)
push!(is, :ECC => "OtherMetalMining", :EC => "OtherMetalMining")
push!(is, :Poll => "CO2", :Year => string(2050))
push!(is, :Area => ["NU"])
push!(is, :FuelEP => ["Diesel","Gasoline", "Kerosene"], :Fuel => ["Diesel","Gasoline", "Kerosene"])
push!(is, :Enduse => ["OffRoad", "Heat"])
push!(is, :Tech => ["OffRoad", "Oil"])

is = copy(dimension_filters)
push!(is, :Poll => "NOX", :Year => string(2015))
push!(is, :ECC => "IronSteel", :EC => "IronSteel")
push!(is, :Area => ["ON"])
push!(is, :Poll => ["SOX", "NOX", "VOC"], :Year => string.(2013:2050))

################################################################################
# Analysis of variables: Overview ##############################################
################################################################################

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters=is, sec) # 
@by TotPol by = :Year :Spruce = sum(:Spruce) :Tanoak == sum(:Tanoak)
@rsubset! TotPol abs(:Diff) >= 1 :Year == 2050
J.add_pdiff!(TotPol)
J.plot_diff(TotPol, dim="Area"; title = "Differences in OtherMetalMining TotPol")

TotPol = J.diff_fast("TotPol", loc1, loc2; dimension_filters=is, sec) # 
@by TotPol by = :Year :Spruce = sum(:Spruce) :Tanoak == sum(:Tanoak)
@rsubset! TotPol abs(:Diff) >= 1 :Year == 2050
J.add_pdiff!(TotPol)
J.plot_diff(TotPol, dim="Area"; title = "Differences in OtherMetalMining TotPol")
@rsubset TotPol :Area == "ON" :Year == 2050

@by TotPol by = :Area :PDiff = maximum(:PDiff)
    # TotPol[ecc,poll,area] = EnPol[ecc,poll,area]+NcPol[ecc,poll,area]+
    #   MEPol[ecc,poll,area]+VnPol[ecc,poll,area]+FlPol[ecc,poll,area]+
    #   FuPol[ecc,poll,area]+ORMEPol[ecc,poll,area]+SqPol[ecc,poll,area]
EnPol = J.diff_fast("EnPol", loc1, loc2; dimension_filters=is, sec) # small problem
NcPol = J.diff_fast("NcPol", loc1, loc2; dimension_filters=is, sec) # 
MEPol = J.diff_fast("MEPol", loc1, loc2; dimension_filters=is, sec) # problem
VnPol = J.diff_fast("VnPol", loc1, loc2; dimension_filters=is, sec) # 
FlPol = J.diff_fast("FlPol", loc1, loc2; dimension_filters=is, sec) # 
FuPol = J.diff_fast("FuPol", loc1, loc2; dimension_filters=is, sec) # 
ORMEPol = J.diff_fast("ORMEPol", loc1, loc2; dimension_filters=is, sec) # 
SqPol = J.diff_fast("SqPol", loc1, loc2; dimension_filters=is, sec) # 

MEReduce = J.diff_fast("MEReduce", loc1, loc2; dimension_filters=is, sec) # Fine

J.add_pdiff!(MEPol)
# MEPol[ecc,poll,area] = (MEPOCA[ecc,poll,area]*MEDriver[ecc,area]*MEPolSwitch[ecc,poll,area])+
#   (xMEPol[ecc,poll,area]*(1-MEPolSwitch[ecc,poll,area]))
MEPOCA = J.diff_fast("MEPOCA", loc1, loc2; dimension_filters=is, sec) # 
J.add_pdiff!(MEPOCA) # issue
MEDriver = J.diff_fast("MEDriver", loc1, loc2; dimension_filters=is, sec) # matches
Driver = J.diff_fast("Driver", loc1, loc2; dimension_filters=is, sec) # matches
xDriver = J.diff_fast("xDriver", loc1, loc2; dimension_filters=is, sec) # matches
DriverMultiplier = J.diff_fast("DriverMultiplier", loc1, loc2; dimension_filters=is, sec) # matches
GRPAdj = J.diff_fast("GRPAdj", loc1, loc2; dimension_filters=is, sec) # matches
MEPolSwitch = J.diff_fast("MEPolSwitch", loc1, loc2; dimension_filters=is, sec) # matches, 1
xMEPol = J.diff_fast("xMEPol", loc1, loc2; dimension_filters=is, sec) # matches

# MEPOCA[ecc,poll,area] = MEPOCX[ecc,poll,area]*MERM[ecc,poll,area]*MEPOCM[ecc,poll,area]
MEPOCX = J.diff_fast("MEPOCX", loc1, loc2; dimension_filters=is, sec = 'M') # issue
J.add_pdiff!(MEPOCX) 
@rsubset! MEPOCX :Diff != 0
# issue starts in 2019, only in NS, Petroleum, VOC

using Plots

# Create the plot with separate lines for each poll
plot(MEPOCX.Year, MEPOCX.Tanoak, 
     group = MEPOCX.Poll,
     xlabel = "Year",
     ylabel = "Tanoak",
     title = "Tanoak by Year for Different Polls",
     legend = :topleft,
     linewidth = 2)
     
@by MEPOCX :Year :PDiff = sum(:PDiff)

MEPOCXB = J.diff_fast("MEPOCX", spr_base, tan_base; dimension_filters=is, sec = 'M') # issue

rename!(MEPOCXB, Dict("Spr_Base" => "Spruce", "Tan_Base" => "Tanoak"))
J.add_pdiff!(MEPOCXB) 
# issue starts in 2019, only in NS, Petroleum, VOC

using Plots

# Create the plot with separate lines for each poll
plot(MEPOCXB.Year, MEPOCXB.Tanoak, 
     group = MEPOCXB.Poll,
     xlabel = "Year",
     ylabel = "Tanoak",
     title = "Tanoak by Year for Different Polls",
     legend = :topleft,
     linewidth = 2)
     
@by MEPOCXB :Year :PDiff = sum(:PDiff)

# Calib3 is good
MEPOCXC = J.diff_fast("MEPOCX", spr_c4, tan_c4; dimension_filters=is, sec = 'M') # issue
J.add_pdiff!(MEPOCXC) 

zoom = copy(is)
push!(zoom, :Year => "2017")
pop!(zoom, :Poll)
MEPOCXC = J.diff_fast("MEPOCX", spr_c4, tan_c4; dimension_filters=zoom, sec = 'M') # issue
J.add_pdiff!(MEPOCXC)
sort!(MEPOCXC, :Spruce)
@rsubset! MEPOCXC :PDiff != 0
push!(zoom, :Poll => "NOX")
# Just ON

area  = pop!(zoom, :Area)
push!(zoom, :Area => area)
ecc = pop!(zoom, :ECC)
push!(zoom, :ECC => ecc)



# @finite_math if abs(MisPol[ec,poll,area] / sum(xEnFPol[fep,ecc,poll,area,year] + xOREnFPol[fep,ecc,poll,area,year] for fep in FuelEPs)) < 0.0001
  #   MisPol[ec,poll,area] = 0
  
  
  # @finite_math MEPOCX[ecc,poll,area,year] = (xMEPol[ecc,poll,area,year] + MisPol[ec,poll,area]) / MEDriver[ecc,area,year]
  xMEPol = J.diff_fast("xMEPol", spr_c4, tan_c4; dimension_filters=zoom, sec = 'M') # fine``
  MEDriver = J.diff_fast("MEDriver", spr_c4, tan_c4; dimension_filters=zoom, sec = 'M') # fine
  MisPol = J.diff_fast("MisPol", spr_c4, tan_c4; dimension_filters=zoom, sec = 'M') # missing
  
  # MisPol[ec,poll,area] = sum(xEnFPol[fep,ecc,poll,area,year] + xOREnFPol[fep,ecc,poll,area,year] + xXCgFPol[fep,ecc,poll,area] -
  # (sum(Polute[eu,fep,ec,poll,area,year] for eu in Enduses) + CgPolEC[fep,ec,poll,area,year]) for fep in FuelEPs)
  
  xEnFPol = J.diff_fast("xEnFPol", loc1, loc2; dimension_filters=zoom, sec = 'M') # matches
  @by xEnFPol :Year :Spruce = sum(:Spruce) :Tanoak = sum(:Tanoak)
  xOREnFPol = J.diff_fast("xOREnFPol", loc1, loc2; dimension_filters=zoom, sec = 'M') # matches
  @by xOREnFPol :Year :Spruce = sum(:Spruce) :Tanoak = sum(:Tanoak)
  @rsubset! xOREnFPol :Diff != 0
  xXCgFPol = J.diff_fast("xXCgFPol", loc1, loc2; dimension_filters=zoom, sec = 'S') # missing
  xCgFPol = J.diff_fast("xCgFPol", loc1, loc2; dimension_filters=zoom, sec = 'S') # matches
  @rsubset! xCgFPol :Diff != 0
  Polute = J.diff_fast("Polute", loc1, loc2; dimension_filters=zoom, sec = 'I') # almost matches
  @rsubset! Polute :Diff != 0
  J.add_pdiff!(Polute)
  CgPolEC = J.diff_fast("CgPolEC", loc1, loc2; dimension_filters=zoom, sec = 'I') # matches
  @rsubset! CgPolEC :Diff != 0
  
  #   ecc = Select(ECC,EC[ec])
  # xXCgFPol[fuelep,ecc,poll,area] = max((xCgFPol[fuelep,ecc,poll,area,year] - CgPolNPRI[fuelep,ecc,poll,area]), 0)
  CgPolNPRI = J.diff_fast("CgPolNPRI", spr_c4, tan_c4; dimension_filters=zoom, sec = 'S') # matches
  @rsubset! CgPolNPRI :Diff != 0

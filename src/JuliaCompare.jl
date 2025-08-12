module JuliaCompare

import CSV
using DataFrames, Chain, DataFramesMeta, HDF5
import PromulaDBA as P
using Makie, CairoMakie
using Colors, CategoricalArrays
import SmallModel: ReadDisk, ReadSets, data_attrs

export db_files, Canada
include("UnCodeMapping.jl")
include("TidyingArrays.jl")

greet() = print("Hello Randy")

const open_databases = Dict{String, HDF5.File}()
function CloseAllDatabases()
  for (db, file) in open_databases
    close(file)
  end
  empty!(open_databases)
  return length(open_databases)
end

abstract type Location end

struct Loc_p <: Location
  vars::DataFrame
  DATA_FOLDER::String
  name::String
end
struct Loc_j <: Location
  vars::DataFrame
  HDF5_path::String
  name::String
end

const db_files = ["2020DB", "CCalDB", "CInput", "COutput", "ECalDB", "EGCalDB",
  "EGInput", "EGOutput", "EInput", "EOutput", "ICalDB", "IInput", "IOutput",
  "MCalDB", "MEInput", "MEOutput", "Minput", "Moutput", "RCalDB", "RInput",
  "ROutput", "SCalDB", "Sinput", "SOutput", "SpInput", "SpOutput", "TCalDB",
  "TInput", "TOutput", "VBInput", "vData_ElectricUnits"
]

function list_var(cfilename, CODE_FOLDER, DATA_FOLDER, verbose=false)
  if verbose
    print("cfilename is: ", cfilename, "\n")
  end

  code_path = joinpath(CODE_FOLDER, join([cfilename, ".src"]))
  df = CSV.read(code_path, DataFrame, delim="\t", comment="*", header=[:x])
  a = occursin.("Define Variable", df.x)
  e = findall(a)
  e = reshape(e, 2, :)
  e[1,:] .+= 1
  e[2,:] .-= 1
  rows = map(:,e[1,:],e[2,:])
  rows = map(collect,rows)
  rows = vcat(rows...)
  df = df[rows, :]
  df = df[Not(occursin.(r"^\s*$", df.x)), :] # remove blank rows

  # Identify dimensions of variable from definition lines
  dd1 = findfirst.(''', df.x)
  dd1 = ifelse.(isnothing.(dd1), length.(df.x), dd1)
  name_dim = SubString.(df.x, 1, dd1 .- 1)
  name_dim = replace.(name_dim, r"\s" => "")
  name_dim = replace.(name_dim, ")" => "")

  ss1 = findfirst.('(', name_dim)
  ss1 = ifelse.(isnothing.(ss1), length.(name_dim), ss1)

  dims = string.(SubString.(name_dim, ss1 .+ 1))
  vars = string.(SubString.(name_dim, 1, ss1 .- 1))

  # Identify description and units of variable from definition lines
  dd2 = findlast.(''', df.x)
  desc_units = SubString.(df.x, dd1 .+ 1, dd2 .- 1)
  dd3 = findfirst.('(', desc_units)
  dd4 = findlast.(')', desc_units)
  dd4 = ifelse.(isnothing.(dd4), length.(desc_units), dd4)

  for i in 1:length(dd3)
    if isnothing(dd3[i])
      j = length(desc_units[i]) + 2
      dd3[i] = j
      dd4[i] = dd3[i] + 1
    end
  end
  units = SubString.(desc_units, dd3 .+ 1, dd4 .- 1)
  desc = string.(SubString.(desc_units, 1, dd3 .- 2))

  split_dims = split.(dims, ",")

  sets = ["EC", "Enduse", "Tech"]
  e2020db = joinpath(DATA_FOLDER, "2020db.dba")
  eg = joinpath(DATA_FOLDER, "EGInput.dba")

  if contains(lowercase(cfilename), "input") | contains(lowercase(cfilename), "output") | contains(lowercase(cfilename), "caldb")
    letter = SubString(cfilename, 1, 1)
    cfile = joinpath(DATA_FOLDER, join([letter, "input.dba"]))
  else
    cfile = e2020db
  end
  print("cfilename is: ", cfilename, " cfile is: ", cfile, "\n")


  makepairs = function (d, cfile, e2020db)
    if d in sets
      Pair(cfile, join([d, "Key"]))
    elseif d == "Unit"
      Pair(eg, "UnCode")
    else
      Pair(e2020db, join([d, "Key"]))
    end
  end

  dim_pairs = [[makepairs(d, cfile, e2020db) for d in ds] for ds in split_dims]

  DataFrame(Variable=vars, Dimensions=dims, Description=desc, Database=cfilename, DPairs=dim_pairs)
end

function list_vars(DATA_FOLDER, db_files)
  pathv = splitpath(DATA_FOLDER)
  pathidx = only(findall(pathv.=="2020Model"))
  CODE_FOLDER = joinpath(joinpath(pathv[1:pathidx-1]), "Engine")
  dfs = map(x -> list_var(x, CODE_FOLDER, DATA_FOLDER, true), db_files)
  vars = vcat(dfs...)

  vars.RowID = 1:size(vars, 1)
  vars = @chain vars begin
    select!(_, :RowID, Not(:RowID)) # Move ID column to the front
  end
  return (vars)
end

function list_vars(CODE_FOLDER, DATA_FOLDER, db_files)
  dfs = map(x -> list_var(x, CODE_FOLDER, DATA_FOLDER, true), db_files)
  vars = vcat(dfs...)

  vars.RowID = 1:size(vars, 1)
  vars = @chain vars begin
    select!(_, :RowID, Not(:RowID)) # Move ID column to the front
  end
  return (vars)
end

function list_vars(file::String)
  data = Vector{NamedTuple{(:Variable, :Database),Tuple{String,String}}}()
  sizehint!(data, 5000)

  h5open(file, "r") do fid
    for name in keys(fid)
      obj = fid[name]
      if isa(obj, HDF5.Dataset)
        push!(data, (Variable=name, Database="/"))
      elseif isa(obj, HDF5.Group)
        for gname in keys(obj)
          push!(data, (Variable=gname, Database=name))
        end
      end
    end
  end

  return DataFrame(data)
end

function var_id(i, vars, DATA_FOLDER)
  fname = joinpath(DATA_FOLDER, string(vars.Database[i], ".dba"))

  vname = string(vars.Variable[i])
  opts = vars.DPairs[i]
  print("fname is: ", fname)
  print("\n")
  # out = P.dataframe(fname, vname, opts..., skip_zeros=false)
  out = P.dataframe(fname, vname, opts...; skip_zeros=false)
  rename!(out,
    Dict(zip(names(out)[1:end-1], split(vars.Dimensions[i], ',')))
  )
  if contains(vars.Dimensions[i], "Year")
    out.Year = parse.(Int64, out.Year)
    if minimum(out.Year) < 1900
      out.Year .+= 1984
    end
  end
  if contains(vars.Dimensions[i], "OGUnit")
    @rsubset! out :OGUnit != ""
  elseif contains(vars.Dimensions[i], "Unit")
    @rsubset! out :Unit != ""
  end
  return (out)
end

function var_id(i::Int, loc::Loc_p)
  var_id(i, loc.vars, loc.DATA_FOLDER)
end

function get_data(i, vars, DATA_FOLDER)
  fname = joinpath(DATA_FOLDER, string(vars.Database[i], ".dba"))
  vname = string(vars.Variable[i])
  P.data(fname, vname)
end

function find_var(needle, vars; exact::Bool=false)
  if exact
    i = findall(lowercase.(vars.Variable) .== lowercase(needle))
    j = findall(lowercase.(vars.Database) .* "/" .* lowercase.(vars.Variable) .== lowercase(needle))
    i = [i; j]
  else
    i = findall(occursin.(lowercase(needle), lowercase.(vars.Variable)))
  end
  vars[i, :]
end

const Canada = ["ON", "QC", "BC", "AB", "MB",
  "SK", "NB", "NS", "NL",
  "PI", "YT", "NT", "NV"
]

function var(needle, vars, DATA_FOLDER)
  temp2 = find_var(needle, vars; exact=true)
  if nrow(temp2) == 1
    return (var_id(temp2.RowID[1], vars, DATA_FOLDER))
  end
  temp = find_var(needle, vars)
  if nrow(temp) == 1
    return (var_id(temp.RowID[1], vars, DATA_FOLDER))
  elseif nrow(temp) == 0
    println(needle, " not found in vars, perhaps your variable is in the list below\n")
    temp2 = findall(occursin.(lowercase.(vars.Variable), lowercase(needle)))
    println(vars[temp2,:])
    error("Variable not found")
    return (vars[temp2,1:5])
  elseif unique(temp, [:Variable, :Dimensions]) == 1
    print("Combining ", nrow(temp), "variables")
    return (map(x -> var_id(x, vars, DATA_FOLDER), temp.RowID))
  else
    println(needle, " could have several values, use var_id to select one of the below\n")
    return (temp)
  end
end

function var(needle, loc::Loc_p)
  (; vars, DATA_FOLDER) = loc
  return (var(needle, vars, DATA_FOLDER))
end

function var(needle, loc::Loc_j)
  (; HDF5_path, vars) = loc
  if contains(needle, "/")
    df = ReadDisk(DataFrame, loc.HDF5_path, needle)
    if "TimeP" ∈ names(df)
      if '(' ∉ df.TimeP[1]
        df.TimeP = map(tp -> "TimeP($(match(r"\d+", tp).match))", df.TimeP)
      end
    end
    if "TimeA" ∈ names(df)
      if '(' ∉ df.TimeA[1]
        df.TimeA = map(tp -> "TimeA($(match(r"\d+", tp).match))", df.TimeA)
      end
    end
    return (df)
  else
    i = findall(lowercase.(vars.Variable) .== lowercase(needle))
    if length(i) == 1
      database = vars.Database[i][1]
      needle = string(database,"/",needle)
      println(needle)
      var(needle, loc)
    else
      print("Multiple variables exist, please specify database")
      df = vars[i, :]
      print(df)
      error("Please specify a database")
    end
  end
end

function diff(df1, df2; name1="new", name2="old")
  dims = names(df1)[1:end-1]
  dims != names(df2)[1:end-1] && error("Dimensions don't match")
  rename!(df1, [dims; name1])
  rename!(df2, [dims; name2])
  if size(df1) == size(df2)
    df = hcat(df1, df2[!, Symbol(name2)])
    rename!(df, [dims; name1; name2])
    println("Assuming DataFrame rows match")
  else
    df = outerjoin(df1, df2, on=dims, makeunique=true)
    replace!(df[!, Symbol(name1)], missing => 0)
    replace!(df[!, Symbol(name2)], missing => 0)
  end
  if eltype(df[!,name1]) == eltype(df[!,name2]) == String
    df.Diff = df[:, Symbol(name1)] .!= df[:, Symbol(name2)]
  elseif eltype(df[!, name1]) <: Number && eltype(df[!, name2]) <: Number
    df.Diff = df[:, Symbol(name1)] .= df[:, Symbol(name2)]
  else
    error("Not able to take differences of these element types: ",
    eltype(df[!,name1]), " and ",
    eltype(df[!,name2])
    )
  end
  # @transform!(df, :Diff = name1 - name2)
  return (df)
end

function diff(name, loc1, loc2; name1=loc1.name, name2=loc2.name)
  df1 = var(name, loc1)
  df2 = var(name, loc2)
  df = diff(df1, df2; name1, name2)
end

function lookup_database(name, loc; sec::Char="")
  if contains(name, "/")
    name_vec = string.(split(name, "/"))
    return (name_vec[2], name_vec[1])
  else
    vars = loc.vars
    temp2 = find_var(name, vars; exact=true)
    n = nrow(temp2)
    if n == 1
      return (temp2.Variable[1], temp2.Database[1])
    elseif n > 1 && sec == ""
      println(name, " could have several values, select one of the below\n")
      println(temp2)
      error("Variable not found")
    elseif n > 1 && sec != "" 
      i = findall(first.(temp2.Database) .== sec)
      println("Found these options")
      println(temp2)
      println("I think you want this one: ", i)
      if length(i) == 1
        ind = i[1]
        return (string(temp2[ind,:Variable]), temp2[ind,:Database])
      else
        println(name, " could have several values, select one of the below\n")
        println(temp2)
        error("Variable not found")
      end
    else
      temp = find_var(name, vars)
      m = nrow(temp)
      if m == 1
        return (temp.Variable[1], temp.Database[1])
      elseif m == 0
        println(name, " not found in vars, perhaps your variable is in the list below\n")
        temp2 = findall(occursin.(lowercase.(vars.Variable), lowercase(name)))
        println(vars[temp2,[:Variable, :Database]])
        error("Variable not found")
        return (vars[temp2,[:Variable, :Database]])
      end
    end
  end
end

function arr_set(vname::String, dbname::String, loc::Loc_j)
  arr = ReadDisk(loc.HDF5_path, string(dbname,"/", vname))
  vattr = data_attrs(loc.HDF5_path, string(dbname,"/", vname))
  if vattr["type"] == "variable"
    sets = ReadSets(loc.HDF5_path, string(dbname,"/", vname))
    if :Unit ∈ keys(sets)
      UnCode = ReadDisk(loc.HDF5_path, "EGInput/UnCode")
      sets.Unit[:] = UnCode[:]
    end
  else
    if vattr["type"] != "set"
      println(string(dbname,"/", vname), " is of type ", vattr["type"], ". This is unusual.")
    end
    values = copy(arr)
    idx = findfirst("Key", vname)
    idx_test = findfirst("key", vname)
    if isnothing(idx)
      key = [Symbol(vname)]
    else
      key = [Symbol(vname[setdiff(1:length(vname), idx)])]
    end
    sets =  (; zip(key, [values])...)
  end
  if :TimeP ∈ keys(sets)
    if '(' ∉ sets[:TimeP][1]
      push!(sets(map(tp -> "TimeP($(match(r"\d+", tp).match))", sets[:TimeP])))
    end
  end
  if :TimeA ∈ keys(sets)
    if '(' ∉ sets[:TimeA][1]
      push!(sets(map(tp -> "TimeA($(match(r"\d+", tp).match))", sets[:TimeA])))
    end
  end
  return(arr, sets)
end

function arr_set(vname::String, dbname::String, loc::Loc_p)
  arr = P.data(joinpath(loc.DATA_FOLDER,string(dbname,".dba")), vname)
  idx = findfirst(loc.vars.Variable .== vname .&& loc.vars.Database .== dbname)
  keys = Symbol.(split(loc.vars.Dimensions[idx], ","))
  dim_pairs = loc.vars.DPairs[idx]
  values = [P.data(key,value) for (key, value) in dim_pairs]
  sets =  (; zip(keys, values)...)
  return(arr, sets)
end


function add_differences!(df::DataFrame, locs::Vector{<:Location}, diff::Union{Bool, Symbol, Vector{Int}})
  value_cols = [loc.name for loc in locs]
  
  if diff === true || diff === :all
    # Calculate all pairwise differences
    for i in 1:length(value_cols)
      for j in (i+1):length(value_cols)
        col_name = Symbol("$(value_cols[i])_minus_$(value_cols[j])")
        df[!, col_name] = df[!, value_cols[i]] .- df[!, value_cols[j]]
      end
    end
    
  elseif diff === :sequential
    # Calculate sequential differences (2-1, 3-2, etc.)
    for i in 2:length(value_cols)
      col_name = Symbol("$(value_cols[i])_minus_$(value_cols[i-1])")
      df[!, col_name] = df[!, value_cols[i]] .- df[!, value_cols[i-1]]
    end
    
  elseif diff === :from_first
    # Calculate differences from the first location
    for i in 2:length(value_cols)
      col_name = Symbol("$(value_cols[i])_minus_$(value_cols[1])")
      df[!, col_name] = df[!, value_cols[i]] .- df[!, value_cols[1]]
    end
    
  elseif isa(diff, Vector{Int})
    # Calculate specific pairwise differences based on indices
    if length(diff) != 2
      error("Vector for diff must contain exactly 2 indices")
    end
    i, j = diff[1], diff[2]
    if i < 1 || i > length(value_cols) || j < 1 || j > length(value_cols)
      error("Indices in diff vector are out of range")
    end
    col_name = Symbol("$(value_cols[i])_minus_$(value_cols[j])")
    df[!, col_name] = df[!, value_cols[i]] .- df[!, value_cols[j]]
  end
  
  return df
end

function safe_percent_diff(new_val, old_val)
  if old_val == 0 && new_val == 0
    return 0.0  # Both zero: no change = 0%
  elseif old_val == 0
    return Inf  # Division by zero (could also return missing or a large number)
  else
    return ((new_val - old_val) / old_val) * 100
  end
end

function add_percent_differences!(df::DataFrame, locs::Vector{<:Location}, percent_diff::Union{Bool, Symbol, Vector{Int}})
  value_cols = [loc.name for loc in locs]
  value_cols_sym = [Symbol(col) for col in value_cols]
  
  # Check that all expected columns exist
  missing_cols = [col for col in value_cols if !hasproperty(df, Symbol(col))]
  if !isempty(missing_cols)
    error("Missing columns in DataFrame: $(missing_cols)")
  end
  
  if percent_diff === true || percent_diff === :all
    # Calculate all pairwise percent differences
    for i in 1:length(value_cols_sym)
      for j in (i+1):length(value_cols_sym)
        col_name = Symbol("$(value_cols[i])_pdiff_$(value_cols[j])")
        df[!, col_name] = safe_percent_diff.(df[!, value_cols_sym[i]], df[!, value_cols_sym[j]])
      end
    end
    
  elseif percent_diff === :sequential
    # Calculate sequential percent differences (2 vs 1, 3 vs 2, etc.)
    for i in 2:length(value_cols_sym)
      col_name = Symbol("$(value_cols[i])_pdiff_$(value_cols[i-1])")
      df[!, col_name] = safe_percent_diff.(df[!, value_cols_sym[i]], df[!, value_cols_sym[i-1]])
    end
    
  elseif percent_diff === :from_first
    # Calculate percent differences from the first location (baseline)
    for i in 2:length(value_cols_sym)
      col_name = Symbol("$(value_cols[i])_pdiff_$(value_cols[1])")
      df[!, col_name] = safe_percent_diff.(df[!, value_cols_sym[i]], df[!, value_cols_sym[1]])
    end
    
  elseif percent_diff === :relative_to_mean
    # Calculate percent difference relative to the row-wise mean
    row_means = mean.(eachrow(df[!, value_cols_sym]))
    for i in 1:length(value_cols_sym)
      col_name = Symbol("$(value_cols[i])_pct_diff_mean")
      df[!, col_name] = safe_percent_diff.(df[!, value_cols_sym[i]], row_means)
    end
    
  elseif isa(percent_diff, Vector{Int})
    # Calculate specific pairwise percent difference based on indices
    if length(percent_diff) != 2
      error("Vector for percent_diff must contain exactly 2 indices")
    end
    i, j = percent_diff[1], percent_diff[2]
    if i < 1 || i > length(value_cols) || j < 1 || j > length(value_cols)
      error("Indices in percent_diff vector are out of range")
    end
    col_name = Symbol("$(value_cols[i])_pdiff_$(value_cols[j])")
    df[!, col_name] = safe_percent_diff.(df[!, value_cols_sym[i]], df[!, value_cols_sym[j]])
  end
  
  return df
end

function var(vname::String, loc; 
                   filter::Dict{Symbol, <:Any}=Dict{Symbol,Any}(),
                   sec::Char="")
  vname, dbname = lookup_database(vname, loc; sec)
  arr, sets = arr_set(vname, dbname, loc)
  # set == set2 ? sets = set : error("Sets Don't Match")
  if !isempty(filter)
    arr, sets = subset_array(arr, sets, filter)
  end
  df = to_tidy_dataframe(arr, sets)
  return df
end

function var(vname::String, locs::Vector{<:Location}; 
             filter::Dict{Symbol, <:Any}=Dict{Symbol,Any}(),
             sec::Char="",
             diff::Union{Bool, Symbol, Vector{Int}}=false,
             pdiff::Union{Bool, Symbol, Vector{Int}}=false)
  # 
  arrs = Vector{AbstractArray}()  # or AbstractArray[]
  sets = Vector{NamedTuple}()     # or NamedTuple[]
  # 
  loc_names = [loc.name for loc in locs]
  if length(unique(loc_names)) != length(loc_names)
      error("Location names must be unique. Found duplicates: $(loc_names)")
  end  
  #
  for loc in locs
    vname, dbname = lookup_database(vname, loc; sec) 
    arr, set = arr_set(vname, dbname, loc)
    push!(arrs, arr)
    push!(sets, set)
  end
  
  allequal(sets) ? (set = sets[1]) : error("Sets Don't Match")
  
  if !isempty(filter)
    # Apply filter to all arrays and get the filtered set (same for all)
    set_filtered = nothing
    for i in eachindex(arrs)
      arrs[i], current_set = subset_array(arrs[i], set, filter)
      if set_filtered === nothing
        set_filtered = current_set  # Store the first filtered set
      end
    end
    set = set_filtered  # Update set to the filtered version
  end
  
  df = to_tidy_dataframe(arrs, set, loc_names)
  
  if diff !== false
    df = add_differences!(df, locs, diff)
  end
  
  if pdiff !== false
    df = add_percent_differences!(df, locs, pdiff)
  end
  
  return df  # Don't forget to return the result!
end

function diff_fast(vname::String, loc1, loc2; 
                   filter::Dict{Symbol, <:Any}=Dict{Symbol,Any}(),
                   sec::Char="")
  vname1, dbname1 = lookup_database(vname, loc1; sec)
  vname2, dbname2 = lookup_database(vname, loc2; sec)
  arr1, set1 = arr_set(vname1, dbname1, loc1)
  arr2, set2 = arr_set(vname2, dbname2, loc2)
  set1 == set2 ? sets = set1 : error("Sets Don't Match")
  arr1, set1 = subset_array(arr1, sets, filter)
  arr2, set2 = subset_array(arr2, sets, filter)
  df1 = to_tidy_dataframe(arr1, set1)
  df2 = to_tidy_dataframe(arr2, set2)
  df = diff(df1, df2; name1 = loc1.name, name2 = loc2.name)
end

function join_vars(df1::DataFrame, df2::DataFrame)
  data = names(df1)[end]
  dims = intersect(names(df1), names(df2))
  dims = dims[dims.!="Value"]
  df = leftjoin(df1, df2, on=dims, makeunique=true)
  values = [setdiff(names(df2), dims); setdiff(names(df1), dims)]
  for c ∈ eachcol(df[!, Symbol.(values)])
    replace!(c, missing => 0)
  end
  # select!(df, names(df)[names(df).!=data], data)
  return (df)
end

function join_vars(df1::DataFrame, df2::DataFrame, df3::DataFrame)
  df = join_vars(df1, df2)
  df = join_vars(df, df3)
  return (df)
end

function join_vars(df1::DataFrame, df2::DataFrame, df3::DataFrame, df4::DataFrame)
  df = join_vars(df1, df2, df3)
  df = join_vars(df, df4)
  return (df)
end

function hist_diff(df; dim="ECC")
  df2 = @by(df, [Symbol(dim)], :Diff = sum(abs.(:Diff)))
  fig = Figure(title="Histogram Showing Concentration of Differences")
  ax = Axis(fig[1, 1]; xlabel="Total Difference (absolute value)", ylabel="$dim Count")
  hist!(ax, df2.Diff)
  display(fig)
end

function tops(df; dim="ECC", num=10)
  df2 = @by(df, [Symbol(dim)], :Diff = sum(abs.(:Diff)))
  df2 = @orderby(df2, -:Diff)
  lst = first(df2[!, dim], num)
  return (lst)
end

function plot_diff(data; dim="ECC", num=10, title="New Plot")
  dim = dim isa String ? dim : string(dim)
  df = deepcopy(data)
  ss = tops(df; dim, num)
  others = setdiff(df[:, dim], ss)
  df[in.(df[:, dim], Ref(others)), dim] .= "Other"
  df = @by(df, [Symbol(dim), :Year], :Diff = sum(:Diff))
  cats = categorical(df[:, dim])
  colors = distinguishable_colors(length(unique(cats)))
  fig = Figure()
  ax = Axis(fig[1, 1]; title=title)
  barplot!(ax, df.Year, df.Diff, stack=levelcode.(cats), color=colors[levelcode.(cats)])
  labels = levels(cats)
  # return (labels)
  elements = [PolyElement(polycolor=colors[i]) for i in 1:length(labels)]
  Legend(fig[1, 2], elements, labels, dim)
  display(fig)
end

function plot_sets(data::DataFrame; 
  col::Union{String,Symbol} = "", 
  dim::Union{String,Symbol}="ECC", 
  num::Integer=10, 
  title::String="New Plot")
  df = deepcopy(data)
  if col == ""
    l = last(names(df))
  else 
    l = col
  end
  l_symbol = Symbol(l)
  @rsubset! df $l != 0
  df2 = @by(df, [Symbol(dim)], :Value = sum(abs.($l)))
  df2 = @orderby(df2, -:Value)
  ss = first(df2[!, dim], num)
  others = setdiff(df[:, dim], ss)
  df[in.(df[:, dim], Ref(others)), dim] .= "Other"
  df = @by(df, [Symbol(dim), :Year], $l = sum($l))
  cats = categorical(df[:, dim])
  colors = distinguishable_colors(length(unique(cats)))
  fig = Figure()
  ax = Axis(fig[1, 1]; title=title)
  barplot!(ax, df.Year, df[:, l], stack=levelcode.(cats), color=colors[levelcode.(cats)])
  labels = String.(levels(cats))   # return (labels)
  elements = [PolyElement(polycolor=colors[i]) for i in 1:length(labels)]
  Legend(fig[1, 2], elements, labels, dim)
  display(fig)
end

function plot_lines(df, cols::AbstractVector{<:Union{Symbol,String,Location}};
  title = "")
  
  if typeof(cols[1]) <: Location
    cols = [Symbol(loc.name) for loc in cols]
  elseif typeof(cols[1]) == String
    cols = [Symbol(col) for col in cols]
  end
  
  dfg = combine(groupby(df, :Year), cols .=> sum, renamecols=false)
  dfg = stack(dfg, cols)
  sort!(dfg, :Year)
  
  fig = Figure()
  ax = Axis(fig[1, 1], 
           xlabel = "Year",
           ylabel = "Tonnes of Pollutants", 
           title = "TotPol across versions")
  
  # Plot each group separately
  for var in unique(dfg.variable)
    subset_df = filter(row -> row.variable == var, dfg)
    lines!(ax, subset_df.Year, subset_df.value, 
           label = var,
           linewidth = 2,
           alpha = 0.5)
  end
  
  axislegend(ax, position = :lt)  # :lt = left top
  
  return fig
end


function comparedata(data, data_b, fnames)
  df = DataFrame()
  for i in 1:length(data)
    j = data[i]
    j_b = data_b[i]
    name = fnames.var[i]
    print(rpad(i, 6), name)
    if eltype(j) == String
      print(" is a string\n")
    else
      print(" \n")
      Ref = sum(j)
      Other = sum(j_b)
      Diff = sum(abs.(j .- j_b))
      LDiff = maximum(abs.(j .- j_b))
      PDiff = Diff == 0 ? 0 : Diff / maximum(abs.([Ref, Other]))
      push!(df, (Var=name, Ref, Other, Diff, PDiff, LDiff); promote=true)
    end
  end
  df[sortperm(df.PDiff, rev=true), :]
  return (df)
end # compare function

function unzip_dbas(DATA_FOLDER, dbas)
  dir, subdir = splitdir(DATA_FOLDER)
  # zip = joinpath(DATA_FOLDER, "dba.zip")
  zip = joinpath(subdir, "dba.zip")
  dbas = lowercase.(dbas)
  dba = dbas[1]
  n = length(dba)
  suffix = dba[(n-3):end]
  if suffix != ".dba"
    dbas .*= ".dba"
  end
  if sum(dbas .== "2020db.dba") == 0
    push!(dbas, "2020db.dba")
  end
  other_sets = contains.(lowercase.(dbas), "output")
  address = unique(lowercase.(dbas[other_sets]))
  for a in address
    letter = SubString(a, 1, 1)
    cfile = join([letter, "input.dba"])
    push!(dbas, cfile)
  end
  for dba in dbas
    if !isfile(joinpath(dir, subdir, dba))
      cmd = Cmd(`wzunzip -d -o $zip $subdir $dba`, dir=dir)
      run(cmd)
    end
  end
end

function archive_list(Year::Int)
  if Year == 2024
    fileloc = raw"\\Scarlet\g\Archives\2020 Canada\2024 thru August"
  elseif Year == 2023
    fileloc = raw"\\Blue\e\Archives\2020 Canada\2023"
  end
  list = readdir(fileloc)
end

function archive_list()
  archive_list(2024)
end

function archive_search(Year::Int, Needle::String)
  list = archive_list(Year::Int)
  return (list[occursin.(lowercase(Needle), lowercase.(list))])
end

function archive_search(Needle::String)
  archive_search(2024, Needle)
end

"""
    subset_array(array::AbstractArray, sets::NamedTuple, dimension_filters::Dict{Symbol, <:Any})

Filters a multi-dimensional array based on dimension names and values.

# Arguments
- `array::AbstractArray`: The array to filter
- `sets::NamedTuple`: Named tuple containing the dimension values for each axis
- `dimension_filters::Dict{Symbol, <:Any}`: Dictionary mapping dimension names to filter values

# Returns
- `Tuple{Array, NamedTuple}`: A tuple containing the filtered array and the filtered dimension sets
"""
function subset_array(array::AbstractArray, sets::NamedTuple, dimension_filters::Dict{Symbol, <:Any})
    # Get the dimensions of the array
    dims = ndims(array)
    
    # Get the named dimensions - preserve the original order
    dim_names = collect(propertynames(sets))
    
    # Check that dimensions match
    if dims != length(dim_names)
        error("Number of dimensions in array ($(dims)) doesn't match number of sets ($(length(dim_names)))")
        end
        
    # Create selections for each dimension
    indices_by_dim = []
    filtered_sets_dict = Dict{Symbol, Vector}()
        
    for dim_name in dim_names
        dim_values = getproperty(sets, dim_name)
        
        if haskey(dimension_filters, dim_name)
            filter_value = dimension_filters[dim_name]
                
            # Find the indices that match the filter
            if filter_value isa Function
                # Function filter
                indices = findall(filter_value, dim_values)
            elseif filter_value isa AbstractArray
                # Array of values
                indices = findall(x -> x in filter_value, dim_values)
            else
                # Single value (exact match)
                indices = findall(x -> x == filter_value, dim_values)
                        end
            
            if isempty(indices)
                error("No values found for filter $(dim_name) => $(filter_value)")
            end
            
            push!(indices_by_dim, indices)
            filtered_sets_dict[dim_name] = dim_values[indices]
                else
            # No filter for this dimension, select all
            indices = collect(1:length(dim_values))
            push!(indices_by_dim, indices)
            filtered_sets_dict[dim_name] = dim_values
                    end
                end
                
    # Create indexing expressions for each dimension
    indexing = Tuple(indices_by_dim)
    
    # Extract the subarray using the indices
    result = array[indexing...]
    
    # Convert the filtered sets dictionary to a NamedTuple with the SAME order as original sets
    filtered_sets = NamedTuple{Tuple(dim_names)}(Tuple(filtered_sets_dict[name] for name in dim_names))
    
    return result, filtered_sets
end

"""
    subset_dataframe(df::DataFrame, filters::Dict; drop_filtered_cols::Bool=false)

Filter a DataFrame based on conditions specified in a dictionary.

# Arguments
- `df`: The DataFrame to filter
- `filters`: Dictionary where keys are column names and values are the filter conditions
- `drop_filtered_cols`: If true, removes the filtered columns from the result (default: false)

# Filter conditions can be:
- A single value (equality condition)
- A vector of values (in condition)
- A range (for numeric columns)
- A function that takes a value and returns a boolean
- A Regex for string columns (matches if value contains the pattern)

# Returns
- A filtered DataFrame that meets all conditions

# Example
```julia
df = DataFrame(
    Year = [2020, 2021, 2022, 2023],
    Region = ["North", "South", "North", "East"],
    Value = [100, 150, 120, 200]
)

# Filter rows where Year is 2022 or 2023 AND Region is "North"
result = subset_dataframe(df, Dict(
    "Year" => [2022, 2023],  # Multiple values
    "Region" => "North"      # Single value
))

# Using ranges and functions
result = subset_dataframe(df, Dict(
    "Year" => 2021:2023,            # Range
    "Value" => x -> x > 100,        # Function condition
    "Region" => r"^[NS]"            # Regex (starts with N or S)
))

# Drop filtered columns after filtering
result = subset_dataframe(df, Dict("Year" => 2022), drop_filtered_cols=true)
```
"""
function subset_dataframe(df_in::DataFrame, filters::Dict{<:Any,<:Any}; drop_filtered_cols::Bool=false)
  df = copy(df_in)
  # Check if the DataFrame is empty
  isempty(df) && return df
  
  # Keep track of columns we've filtered on
  filtered_cols = Symbol[]
  
  # Apply each filter
  for (col_key, condition) in filters
      # Handle different key types (Symbol, String, etc.)
      col_sym = col_key isa Symbol ? col_key : Symbol(string(col_key))
      # Check if the column exists in the DataFrame
      if !(col_sym in Symbol.(names(df)))
        # Skip this filter if column doesn't exist
        println("Skipping col ", col_sym)
        continue
      end
      println("Using col ", col_sym)
      if col_sym == :Year
        if condition isa AbstractVector
          condition = [parse(Int,c) for c in condition]
        else 
          condition = parse(Int,condition)
        end
      end
      println("Condition is: ", condition)
      # Add to list of filtered columns
      push!(filtered_cols, col_sym)
      
      # Apply the filter based on the condition type
      if condition isa AbstractVector
          # Vector of allowed values
          filter!(row -> row[col_sym] in condition, df)
      elseif condition isa AbstractRange
          # Range of values
          filter!(row -> row[col_sym] in condition, df)
      elseif condition isa Function
          # Function that returns boolean
          filter!(row -> condition(row[col_sym]), df)
      elseif condition isa Regex && eltype(df[!, col_sym]) <: AbstractString
          # Regex pattern for string columns
          filter!(row -> !isnothing(match(condition, row[col_sym])), df)
      else
          # Single value equality check
          filter!(row -> row[col_sym] == condition, df)
      end
      println("df has nrows: ", size(df))
  end
  
  # Drop the filtered columns if requested
  if drop_filtered_cols && !isempty(filtered_cols)
      select!(df, Not(filtered_cols))
  end
  
  return df
end

"""
    subset_dataframe!(df::DataFrame, filters::Dict; drop_filtered_cols::Bool=false)

In-place version of subset_dataframe that modifies the input DataFrame directly.

# Arguments and behavior are the same as subset_dataframe
"""
function subset_dataframe!(df::DataFrame, filters::Dict; drop_filtered_cols::Bool=false)
  # Check if the DataFrame is empty
  isempty(df) && return df
  
  # Keep track of columns we've filtered on
  filtered_cols = Symbol[]
  
  # Apply each filter
  for (col_key, condition) in filters
      # Handle different key types (Symbol, String, etc.)
      col_sym = col_key isa Symbol ? col_key : Symbol(string(col_key))
      # Check if the column exists in the DataFrame
      if !(col_sym in Symbol.(names(df)))
        # Skip this filter if column doesn't exist
        println("Skipping col ", col_sym)
        continue
      end
      println("Using col ", col_sym)
      if col_sym == :Year
        if condition isa AbstractVector
          condition = [parse(Int,c) for c in condition]
        else 
          condition = parse(Int,condition)
        end
      end
      println("Condition is: ", condition)
      # Add to list of filtered columns
      push!(filtered_cols, col_sym)
      
      # Apply the filter based on the condition type
      if condition isa AbstractVector
          # Vector of allowed values
          filter!(row -> row[col_sym] in condition, df)
      elseif condition isa AbstractRange
          # Range of values
          filter!(row -> row[col_sym] in condition, df)
      elseif condition isa Function
          # Function that returns boolean
          filter!(row -> condition(row[col_sym]), df)
      elseif condition isa Regex && eltype(df[!, col_sym]) <: AbstractString
          # Regex pattern for string columns
          filter!(row -> !isnothing(match(condition, row[col_sym])), df)
      else
          # Single value equality check
          filter!(row -> row[col_sym] == condition, df)
      end
      println("df has nrows: ", size(df))
  end
  
  # Drop the filtered columns if requested
  if drop_filtered_cols && !isempty(filtered_cols)
      select!(df, Not(filtered_cols))
  end
  
  return df
end
# Performance optimized version for large DataFrames
"""
    subset_dataframe_optimized(df::DataFrame, filters::Dict; drop_filtered_cols::Bool=false)

Optimized version of subset_dataframe that uses DataFrames.jl's built-in filtering
capabilities more efficiently. This version is faster for large DataFrames.

# Arguments and behavior are the same as subset_dataframe
"""
function subset_dataframe_optimized(df::DataFrame, filters::Dict; drop_filtered_cols::Bool=false)
    # Check if the DataFrame is empty
    isempty(df) && return df
    
    result_df = df
    filtered_cols = String[]
    
    for (col, condition) in filters
        # Check if the column exists in the DataFrame
        col_sym = Symbol(col)
        if !(col in names(result_df) || col_sym in names(result_df))
            continue
        end
        
        push!(filtered_cols, string(col))
        
        # Use DataFrame's built-in filtering with ByRow for better performance
        if condition isa AbstractVector
            result_df = filter(col_sym => ByRow(x -> x in condition), result_df)
        elseif condition isa AbstractRange
            result_df = filter(col_sym => ByRow(x -> x in condition), result_df)
        elseif condition isa Function
            result_df = filter(col_sym => ByRow(condition), result_df)
        elseif condition isa Regex
            result_df = filter(col_sym => ByRow(x -> 
                typeof(x) <: AbstractString && !isnothing(match(condition, x))), result_df)
        else
            result_df = filter(col_sym => ByRow(x -> x == condition), result_df)
        end
    end
    
    # Drop the filtered columns if requested
    if drop_filtered_cols && !isempty(filtered_cols)
        result_df = select(result_df, Not(Symbol.(filtered_cols)))
    end
    
    return result_df
end

"""
    to_tidy_dataframe(array::AbstractArray, sets::NamedTuple)

Converts a multi-dimensional array and its dimension sets into a tidy DataFrame format.

# Arguments
- `array::AbstractArray`: The array to convert
- `sets::NamedTuple`: Named tuple containing the dimension values for each axis

# Returns
- `DataFrame`: A tidy DataFrame with one row per cell in the array

# Example
```julia
filtered_array, filtered_sets = subset_array(HDGCCI_j, sets, dimension_filters)
df = to_tidy_dataframe(filtered_array, filtered_sets)
```
"""

function to_tidy_dataframe(array::AbstractArray, sets::NamedTuple)
  # Get the dimension names
  dim_names = collect(propertynames(sets))
  # Check that dimensions match
  if ndims(array) != length(dim_names)
    error("Number of dimensions in array ($(ndims(array))) doesn't match number of sets ($(length(dim_names)))")
  end
  
  # Create arrays for each dimension value
  values_by_dim = []
  for dim_name in dim_names
    push!(values_by_dim, getproperty(sets, dim_name))
                    end
  
  # Create the result DataFrame
  df = DataFrame()
  
  # Use Base.product to iterate through all combinations of indices
  indices_iter = Iterators.product([1:length(vals) for vals in values_by_dim]...)
  
  # Initialize arrays to hold the data
  n_elements = length(array)
  dim_columns = [Vector{eltype(values_by_dim[i])}(undef, n_elements) for i in 1:length(dim_names)]
  value_column = Vector{eltype(array)}(undef, n_elements)
  
  # Fill in the data
  for (i, idx) in enumerate(indices_iter)
    # Get the actual value
    value_column[i] = array[idx...]
    
    # Get the dimension values
    for (d, dim_idx) in enumerate(idx)
      dim_columns[d][i] = values_by_dim[d][dim_idx]
    end
  end
  
  # Add the columns to the DataFrame
  for (i, name) in enumerate(dim_names)
    df[!, Symbol(name)] = dim_columns[i]
  end
  
  # Add the value column
  df[!, :Value] = value_column
  
  if :Year ∈ dim_names
    df.Year = parse.(Int64, df.Year)
  end
  
  return df
end

function to_tidy_dataframe(arrays::Vector{<:AbstractArray}, sets::NamedTuple, value_names::Vector{String})
  # Check that we have at least one array
  if isempty(arrays)
    error("At least one array must be provided")
  end
  
  # Check that all arrays have the same dimensions
  ref_size = size(arrays[1])
  ref_ndims = ndims(arrays[1])
  
  for (i, arr) in enumerate(arrays)
    if size(arr) != ref_size
      error("Array $i has size $(size(arr)), but expected $(ref_size)")
    end
    if ndims(arr) != ref_ndims
      error("Array $i has $(ndims(arr)) dimensions, but expected $(ref_ndims)")
    end
  end
  
  # Check that number of value names matches number of arrays
  if length(value_names) != length(arrays)
    error("Number of value names ($(length(value_names))) doesn't match number of arrays ($(length(arrays)))")
  end
  
  # Get the dimension names
  dim_names = collect(propertynames(sets))
  
  # Check that dimensions match the sets
  if ref_ndims != length(dim_names)
    error("Number of dimensions in arrays ($(ref_ndims)) doesn't match number of sets ($(length(dim_names)))")
  end
  
  # Create arrays for each dimension value
  values_by_dim = []
  for dim_name in dim_names
    push!(values_by_dim, getproperty(sets, dim_name))
  end
  
  # Create the result DataFrame
  df = DataFrame()
  
  # Use Base.product to iterate through all combinations of indices
  indices_iter = Iterators.product([1:length(vals) for vals in values_by_dim]...)
  
  # Initialize arrays to hold the data
  n_elements = length(arrays[1])
  dim_columns = [Vector{eltype(values_by_dim[i])}(undef, n_elements) for i in 1:length(dim_names)]
  
  # Initialize value columns for each array
  value_columns = [Vector{eltype(arr)}(undef, n_elements) for arr in arrays]
  
  # Fill in the data
  for (i, idx) in enumerate(indices_iter)
    # Get the dimension values (same for all arrays)
    for (d, dim_idx) in enumerate(idx)
      dim_columns[d][i] = values_by_dim[d][dim_idx]
    end
    
    # Get the values from each array
    for (a, arr) in enumerate(arrays)
      value_columns[a][i] = arr[idx...]
    end
  end
  
  # Add the dimension columns to the DataFrame
  for (i, name) in enumerate(dim_names)
    df[!, Symbol(name)] = dim_columns[i]
  end
  
  # Add the value columns with custom names
  for (i, value_name) in enumerate(value_names)
    df[!, Symbol(value_name)] = value_columns[i]
  end
  
  # Parse Year column if it exists
  if :Year ∈ dim_names
    df.Year = parse.(Int64, df.Year)
  end
  
  return df
end

# Convenience method with default value names
function to_tidy_dataframe(arrays::Vector{<:AbstractArray}, sets::NamedTuple)
  value_names = ["Value_$i" for i in 1:length(arrays)]
  return to_tidy_dataframe(arrays, sets, value_names)
end
"""
    ReadDiskRaw(db::String, name::String; dimension_filters::Dict = Dict())

Reads the dataset named `name` from the HDF5 file specified by `db` and returns
the raw array with filters applied at read time.

# Arguments
- `db::String`: Path to the HDF5 file
- `name::String`: Name of the dataset to read
- `dimension_filters::Dict=Dict()`: A dictionary mapping dimension names to filter criteria

# Returns
- `Tuple{Array, Dict}`: A tuple containing the filtered array and a dictionary 
  mapping dimension names to their filtered values

# Example
```julia
# Read only data for specific years
data, dim_info = ReadDiskRaw("data.h5", "my_dataset", 
                             dimension_filters = Dict("Year" => [2020, 2030]))
```
"""
function ReadDiskRaw(db::String, name::String; dimension_filters::Dict = Dict())
  h5open(db, "r") do f
  if !haskey(f, name)
  throw(HDF5DataSetNotFoundException(db, name))
  end
    dataset = f[name]
    attr = Dict(attrs(dataset))
    
    # Get the dimensions
    if !haskey(attr, "dims")
        error("Dataset $name does not have dimension attributes")
    end
    
    dims = attr["dims"]
    
    # Process dimension data and filters (same as in the DataFrame version)
    dim_values = Dict{String, Vector}()
    dim_indices = Dict{String, Vector{Int}}()
    
    for (i, dim) in enumerate(dims)
        # [Same dimension processing code as in the DataFrame version]
        dim_path = "$(dirname(name))/$dim"
        if haskey(f, dim_path)
            values = read(f[dim_path])
            # If first value is empty, use indices
            if length(values) > 0 && first(values) == ""
                values = collect(1:length(values))
            end
        else
            # If dimension data not available, use indices
            values = collect(1:size(dataset, i))
        end
        
        dim_values[dim] = values
        
        # Apply filter for this dimension if provided
        selected_indices = collect(1:length(values))
        
        if haskey(dimension_filters, dim)
            filter_criterion = dimension_filters[dim]
            
            if filter_criterion isa Function
                # Apply function to values
                selected_indices = findall(i -> filter_criterion(values[i]), selected_indices)
            elseif filter_criterion isa AbstractRange
                # Find indices where values are in the range
                if eltype(values) <: Number && eltype(filter_criterion) <: Number
                    selected_indices = findall(i -> values[i] in filter_criterion, selected_indices)
                else
                    # For string values, need to match exactly
                    selected_indices = findall(i -> string(values[i]) in string.(filter_criterion), selected_indices)
                end
            elseif filter_criterion isa AbstractArray
                # Find indices where values match any in the array
                selected_indices = findall(i -> values[i] in filter_criterion, selected_indices)
            else
                # Single value - find exact matches
                selected_indices = findall(i -> values[i] == filter_criterion, selected_indices)
            end
        end
        
        dim_indices[dim] = selected_indices
    end
    
    # Create the selection for HDF5 hyperslab
    selection = Vector{Any}(undef, length(dims))
    for (i, dim) in enumerate(dims)
        selection[i] = dim_indices[dim]
    end
    
    # Read only the selected data
    filtered_array = dataset[selection...]
    
    # Create a dictionary with the filtered dimension values
    filtered_dim_values = Dict()
    for (dim, indices) in dim_indices
        filtered_dim_values[dim] = dim_values[dim][indices]
    end
    
    # Check if the array has indefinite values
    throw_error_if_arr_contains_indefinite_values(filtered_array, name)
    
    return (filtered_array, filtered_dim_values)
  end
end

"""
    add_pdiff(df::DataFrame)

Add a percent difference (PDiff) column to a DataFrame with Spruce and Tanoak columns.
Sets PDiff to zero when Diff is zero.

# Arguments
- `df`: DataFrame with columns Spruce, Tanoak, and optionally Diff

# Returns
- The modified DataFrame with a new PDiff column added

# Example
```julia
# Assumes Dmd has Spruce, Tanoak, and Diff columns
add_pdiff(Dmd)
```
"""
function add_pdiff(df)
    # Check if required columns exist
    required_cols = [:Spruce, :Tanoak]
    
    for col in required_cols
        if !(col in propertynames(df))
            error("Column '$col' not found in DataFrame")
        end
    end
    
    # Create a copy
    result = copy(df)
    
    # Check if Diff column exists
    has_diff = :Diff in propertynames(df)
    
    # Calculate PDiff
    if has_diff
        # Calculate raw percent differences
        raw_pdiff = (df.Spruce .- df.Tanoak) ./ df.Tanoak .* 100.0
        
        # Round to 2 decimal places
        raw_pdiff = round.(raw_pdiff, digits=2)
        
        # Apply the rule: if Diff is zero, PDiff is zero
        result.PDiff = [isapprox(diff, 0, atol=1e-10) ? 0.0 : pdiff 
                      for (diff, pdiff) in zip(df.Diff, raw_pdiff)]
    else
        # Just calculate percent difference without the zero check
        result.PDiff = round.((df.Spruce .- df.Tanoak) ./ df.Tanoak .* 100.0, digits=2)
    end
    
    return result
end

"""
    add_pdiff!(df::DataFrame)

In-place version that adds a PDiff column directly to the input DataFrame.
"""
function add_pdiff!(df)
    # Check if required columns exist
    col_diff = findfirst(names(df).=="Diff")
    required_cols = Symbol.(names(df)[(col_diff-2):(col_diff-1)])

    for col in required_cols
        if !(col in propertynames(df))
            error("Column '$col' not found in DataFrame")
        end
    end
    
    # Check if Diff column exists
    has_diff = :Diff in propertynames(df)
    
    # Calculate PDiff
    if has_diff
        # Calculate raw percent differences
        raw_pdiff = (df[:,required_cols[1]] .- df[:,required_cols[2]]) ./ df[:,required_cols[2]] .* 100.0
        
        # Round to 2 decimal places
        raw_pdiff = round.(raw_pdiff, digits=2)
        
        # Apply the rule: if Diff is zero, PDiff is zero
        df.PDiff = [isapprox(diff, 0, atol=1e-10) ? 0.0 : pdiff 
                  for (diff, pdiff) in zip(df.Diff, raw_pdiff)]
    else
        # Just calculate percent difference without the zero check
        df.PDiff = round.((df[:,required_cols[1]] .- df[:,required_cols[2]]) ./ df[:,required_cols[2]] .* 100.0, digits=2)
    end
    
    return df
end

# Direct access versions using column properties
"""
    direct_add_pdiff!(df)

A direct implementation that accesses columns directly by property.
May help avoid column lookup issues.
"""
function direct_add_pdiff!(df)
    # Calculate percent differences
    pdiff_values = Vector{Float64}(undef, size(df, 1))
    
    for i in 1:size(df, 1)
        spruce_val = df[i, :Spruce]
        tanoak_val = df[i, :Tanoak]
        
        # Calculate raw percent difference
        raw_pdiff = (spruce_val - tanoak_val) / tanoak_val * 100.0
        
        # Round to 2 decimal places
        raw_pdiff = round(raw_pdiff, digits=2)
        
        # Check if Diff is zero (if Diff column exists)
        if :Diff in propertynames(df)
            diff_val = df[i, :Diff]
            pdiff_values[i] = isapprox(diff_val, 0, atol=1e-10) ? 0.0 : raw_pdiff
        else
            pdiff_values[i] = raw_pdiff
        end
    end
    
    # Add the column
    df.PDiff = pdiff_values
    
    return df
end

end # module JuliaCompare

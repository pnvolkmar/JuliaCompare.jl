module JuliaCompare

import CSV
using DataFrames, Chain, DataFramesMeta, HDF5
import PromulaDBA as P
using Makie, CairoMakie
using Colors, CategoricalArrays
import SmallModel: ReadDisk

include("UnCodeMapping.jl")

greet() = print("Hello Randy")

function f_on(df_in)
  df = copy(df_in)
  if "Area" ∈ names(df)
    @rsubset! df :Area ∈ ["ON"]
  end
  if "Fuel" ∈ names(df)
    @rsubset! df :Fuel == "Ethanol"
  end
  if "FuelEP" ∈ names(df)
    @rsubset! df :FuelEP == "Ethanol"
  end
  if "EC" ∈ names(df)
    @rsubset! df :EC == "ResidentialOffRoad"
  end
  if "ECC" ∈ names(df)
    @rsubset! df :ECC == "ResidentialOffRoad"
  end
  if "Poll" ∈ names(df)
    @rsubset! df :Poll ∈ ["CO2", "COX"]
  end
  if "Year" ∈ names(df)
    @rsubset! df :Year ∈ [2020, 2021]
  end
  if "Diff" ∈ names(df)
    @rsubset! df abs(:Diff) > 1e-7
  end
  return (df)
end

function f_al(df_in)
  df = copy(df_in)
  if "Area" ∈ names(df)
    @rsubset! df :Area ∈ ["QC"]
  end
  if "Year" ∈ names(df)
    @rsubset! df :Year ∈ [1990]
  end
  if "EC" ∈ names(df)
    @rsubset! df :EC ∈ ["Aluminum"]
  end
  if "ECC" ∈ names(df)
    @rsubset! df :ECC ∈ ["Aluminum"]
  end
  return (df)
end

function f_fp(df_in)
  df = copy(df_in)
  if "Area" ∈ names(df)
    @rsubset! df :Area ∈ ["ON"]
  end
  if "FuelEP" ∈ names(df)
    @rsubset! df :FuelEP == "JetFuel"
  end
  if "Fuel" ∈ names(df)
    @rsubset! df :Fuel == "JetFuel"
  end
  if "ECC" ∈ names(df)
    @rsubset! df :ECC == "ForeignPassenger"
  end
  if "EC" ∈ names(df)
    @rsubset! df :EC == "ForeignPassenger"
  end
  if "Poll" ∈ names(df)
    @rsubset! df :Poll == "CO2"
  end
  if "Year" ∈ names(df)
    @rsubset! df :Year > 2007 :Year < 2030 # ∈ [2007,2008]
  end
  if "Diff" ∈ names(df)
    @rsubset! df abs(:Diff) >= 1e-7
  end
  return (df)
end

function filter_bm(df_in)
  df = copy(df_in)
  if "Area" ∈ names(df)
    @rsubset! df :Area ∈ ["ON"]
  end
  if "Enduse" ∈ names(df)
    @rsubset! df :Enduse ∈ ["HW","OthSub"]
  end
  if "EC" ∈ names(df)
    @rsubset! df :EC == "Offices"
  end
  if "ECC" ∈ names(df)
    @rsubset! df :ECC == "Offices"
  end
  if "Tech" ∈ names(df)
    @rsubset! df :Tech ∈ ["Biomass"]
  end
  if "Year" ∈ names(df)
    @rsubset! df :Year ∈ [1985,2000]
  end
  return (df)
end

function f_og(df_in)
  df = copy(df_in)
  if "Area" ∈ names(df)
    @rsubset! df :Area ∈ ["AB"]
  end
  # if "Enduse" ∈ names(df)
  #   @rsubset! df :Enduse ∈ ["HW","OthSub"]
  # end
  if "EC" ∈ names(df)
    @rsubset! df :EC == "ConventionalGasProduction"
  end
  if "ECC" ∈ names(df)
    @rsubset! df :ECC == "ConventionalGasProduction"
  end
  # if "Tech" ∈ names(df)
  #   @rsubset! df :Tech ∈ ["Biomass"]
  # end
  if "Year" ∈ names(df)
    @rsubset! df :Year ∈ [2000]
  end
  if "Diff" ∈ names(df)
    @rsubset! df abs(:Diff) >= 1e-7
  end
  return (df)
end


struct Loc_p
  vars::DataFrame
  DATA_FOLDER::String
  name::String
end
struct Loc_j
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

function list_vars(CODE_FOLDER, DATA_FOLDER, db_files)
  dfs = map(x -> list_var(x, CODE_FOLDER, DATA_FOLDER, true), db_files)
  vars = vcat(dfs...)

  vars.RowID = 1:size(vars, 1)
  vars = @chain vars begin
    select!(_, :RowID, Not(:RowID)) # Move ID column to the front
  end
  return (vars)
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
  df = outerjoin(df1, df2, on=dims, makeunique=true)
  replace!(df[!, Symbol(name1)], missing => 0)
  replace!(df[!, Symbol(name2)], missing => 0)
  df.Diff = df[:, Symbol(name1)] .- df[:, Symbol(name2)]
  # @transform!(df, :Diff = name1 - name2)
  return (df)
end

function diff(name, loc1, loc2; name1=loc1.name, name2=loc2.name)
  df1 = var(name, loc1)
  df2 = var(name, loc2)
  df = diff(df1, df2; name1, name2)
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

function plot_sets(data; dim="ECC", num=10, title="New Plot")
  df = deepcopy(data)
  l = last(names(df))
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
  labels = levels(cats)
  # return (labels)
  elements = [PolyElement(polycolor=colors[i]) for i in 1:length(labels)]
  Legend(fig[1, 2], elements, labels, dim)
  display(fig)
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
    ReadDisk(::Type{DataFrame}, db::String, name::String; 
             skip_zeros = false,
             dimension_filters::Dict = Dict())

Reads the dataset named `name` from the HDF5 file specified by `db` and returns a DataFrame.
Applies dimension-based filters to minimize memory usage.

# Arguments
- `::Type{DataFrame}`: The type parameter specifying a DataFrame return
- `db::String`: Path to the HDF5 file
- `name::String`: Name of the dataset to read
- `skip_zeros::Bool=false`: Whether to skip rows with zero values in the final result
- `dimension_filters::Dict=Dict()`: A dictionary mapping dimension names to filter criteria

# Returns
- `DataFrame`: The dataset converted to a DataFrame with only the filtered data
"""
function ReadDisk(::Type{DataFrame}, db::String, name::String; 
                  skip_zeros = false,
                  dimension_filters::Dict = Dict())
    h5open(db, "r") do f
        if !haskey(f, name)
            throw(HDF5DataSetNotFoundException(db, name))
        end
        
        dataset = f[name]
        attr = Dict(attrs(dataset))
        
        if get(attr, "type", "") != "variable"
            error("`ReadDisk(DataFrame, db, \"$name\")` not supported for sets")
        end
        
        # Get the dimensions
        dims = get(attr, "dims", [])
        if isempty(dims)
            error("Dataset $name does not have dimension attributes")
        end
        
        # For each dimension, get its values and apply filters
        filtered_dim_values = Dict{Symbol, Vector}()
        
        for (i, dim) in enumerate(dims)
            dim_path = "$(dirname(name))/$dim"
            
            # Get dimension values
            if haskey(f, dim_path)
                values = read(f[dim_path])
                if length(values) > 0 && first(values) == ""
                    values = collect(1:length(values))
                end
            else
                values = collect(1:size(dataset, i))
            end
            
            # Filter values if needed
            if haskey(dimension_filters, dim)
                filter_criterion = dimension_filters[dim]
                
                filtered_indices = Int[]
                
                if filter_criterion isa Function
                    # Apply function to values
                    for (idx, val) in enumerate(values)
                        if filter_criterion(val)
                            push!(filtered_indices, idx)
                        end
                    end
                elseif filter_criterion isa AbstractRange
                    # Check if values are in range
                    for (idx, val) in enumerate(values)
                        if val in filter_criterion
                            push!(filtered_indices, idx)
                        end
                    end
                elseif filter_criterion isa AbstractArray
                    # Check if values match any in array
                    for (idx, val) in enumerate(values)
                        if val in filter_criterion
                            push!(filtered_indices, idx)
                        end
                    end
                else
                    # Single value match
                    for (idx, val) in enumerate(values)
                        if val == filter_criterion
                            push!(filtered_indices, idx)
                        end
                    end
                end
                
                # Store filtered values
                filtered_dim_values[Symbol(dim)] = values[filtered_indices]
            else
                # No filter, use all values
                filtered_dim_values[Symbol(dim)] = values
            end
        end
        
        # Create DataFrame with all combinations of filtered dimensions
        df = allcombinations(DataFrame, [(dim => vals) for (dim, vals) in filtered_dim_values]...)
        
        # Now read the values for each row in the DataFrame
        n_rows = nrow(df)
        values = Vector{eltype(dataset)}(undef, n_rows)
        
        for i in 1:n_rows
            # Create the selection indices for this row
            selection = Any[]
            
            for dim in dims
                dim_sym = Symbol(dim)
                # Find the index of this dimension's value in the original dataset
                dim_val = df[i, dim_sym]
                
                dim_path = "$(dirname(name))/$dim"
                if haskey(f, dim_path)
                    orig_values = read(f[dim_path])
                    if length(orig_values) > 0 && first(orig_values) == ""
                        # Using indices
                        idx = dim_val
                    else
                        # Using values
                        idx = findfirst(==(dim_val), orig_values)
                    end
                else
                    # Using indices
                    idx = dim_val
                end
                
                push!(selection, idx)
            end
            
            # Read the value for this combination
            values[i] = dataset[selection...]
        end
        
        # Add values to DataFrame
        df[!, :Value] = values
        
        # Convert Year column to Int if it exists
        if "Year" in names(df)
            df.Year = parse.(Int, string.(df.Year))
        end
        
        # Add metadata
        metadata!(df, "variable", basename(name); style = :note)
        metadata!(df, "group", dirname(name); style = :note)
        metadata!(df, "name", name; style = :note)
        metadata!(df, "units", get(attr, "units", ""); style = :note)
        
        # Apply zero value filtering if requested
        if skip_zeros
            subset!(df, :Value => ByRow(x -> !isapprox(x, 0.0)))
        end
        
        return df
    end
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

end # module JuliaCompare

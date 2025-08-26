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
include("CleanSets.jl")
include("Diagnostics.jl")

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
  "KInput", "KOutput",
  "MCalDB", "MEInput", "MEOutput", "MInput", "MOutput", "RCalDB", "RInput",
  "ROutput", "SCalDB", "SInput", "SOutput", "SpInput", "SpOutput", "TCalDB",
  "TInput", "TOutput", "VBInput", "vData_ElectricUnits"
]

function loc(DATA_FOLDER::String, name::String)
  if endswith(lowercase(DATA_FOLDER), "hdf5")
    HDF5_path = DATA_FOLDER
  else
    HDF5_path = joinpath(DATA_FOLDER,"database.hdf5")
  end
  if isfile(HDF5_path)
    @info "Creating a Julia Location"
    vars = list_vars(HDF5_path)
    return(Loc_j(vars, HDF5_path, name))
  elseif isfile(joinpath(DATA_FOLDER,"2020db.dba"))
    @info "Creating a Promula Location"
    vars = list_vars(DATA_FOLDER, db_files)
    return(Loc_p(vars,DATA_FOLDER, name))
  else
    @error "Neither a database.hdf5 nor a 2020db.dba was found at $DATA_FOLDER"
  end
end

function fltr()
  fltr = Dict{Symbol,Any}()
  push!(fltr, :Area => Canada, :Year => string.(1986:2050))
  return fltr
end

function list_var(cfilename, CODE_FOLDER, DATA_FOLDER, verbose=false)
  if verbose
    print("cfilename is: ", cfilename, "\n")
  end

  code_path = joinpath(CODE_FOLDER, join([cfilename, ".src"]))
  df = CSV.read(code_path, DataFrame, delim="\t", comment="*", header=[:x])
  
  # Find all "Open Input" lines to track database changes
  open_input_lines = findall(occursin.(r"Open.*\.dba", df.x))
  # println(open_input_lines)
  # Extract database names from "Open Input" lines
  db_names = String[]
  for line_idx in open_input_lines
    line = df.x[line_idx]
    # Extract database name between quotes
   match_result = match(r"Open\s+\w+\d*\s+\"([^\"]+)\"", line)
  #  println(match_result)
    if match_result !== nothing
      db_name = replace(match_result.captures[1], ".dba" => "")
      push!(db_names, db_name)
    end
  end
  
  # Find variable definition blocks
  a = occursin.("Define Variable", df.x)
  define_var_lines = findall(a)
  end_define_lines = findall(occursin.("End Define Variable", df.x))
  setdiff!(define_var_lines, end_define_lines)
  
  # Pair up Define Variable with End Define Variable
  var_blocks = []
  for i in 1:length(define_var_lines)
    start_line = define_var_lines[i] + 1
    end_line = end_define_lines[i] - 1
    
    # Find which database this block belongs to
    current_db = cfilename  # default
    for (j, open_line) in enumerate(open_input_lines)
      if open_line < define_var_lines[i] && j <= length(db_names)
        current_db = db_names[j]
      end
    end
    
    push!(var_blocks, (start_line:end_line, current_db))
  end
  
  # Extract all variable definition rows with their corresponding databases
  all_rows = Int[]
  db_assignments = String[]
  
  for (rows, db_name) in var_blocks
    for row in rows
      push!(all_rows, row)
      push!(db_assignments, db_name)
    end
  end
  
  # Filter dataframe to only variable definition rows
  df_vars = df[all_rows, :]
  df_vars = df_vars[Not(occursin.(r"^\s*$", df_vars.x)), :] # remove blank rows

  # Filter corresponding database assignments
  non_empty_mask = .!occursin.(r"^\s*$", df[all_rows, :x])
  db_assignments = db_assignments[non_empty_mask]
  
  # Rest of your existing parsing logic
  dd1 = findfirst.(''', df_vars.x)
  dd1 = ifelse.(isnothing.(dd1), length.(df_vars.x), dd1)
  name_dim = SubString.(df_vars.x, 1, dd1 .- 1)
  name_dim = replace.(name_dim, r"\s" => "")
  name_dim = replace.(name_dim, ")" => "")

  ss1 = findfirst.('(', name_dim)
  ss1 = ifelse.(isnothing.(ss1), length.(name_dim), ss1)

  dims = string.(SubString.(name_dim, ss1 .+ 1))
  vars = string.(SubString.(name_dim, 1, ss1 .- 1))

  # Identify description and units of variable from definition lines
  dd2 = findlast.(''', df_vars.x)
  desc_units = SubString.(df_vars.x, dd1 .+ 1, dd2 .- 1)
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
    kinput = joinpath(DATA_FOLDER,"KInput.dba")
    if d in sets
      Pair(cfile, join([d, "Key"]))
    elseif d == "Unit"
      Pair(eg, "UnCode")
    elseif d == "OGUnit"
      Pair(e2020db, "OGCode")
    elseif contains(d, "TOM")
      Pair(kinput, join([d, "Key"]))
    else
      Pair(e2020db, join([d, "Key"]))
    end
  end

  dim_pairs = [[makepairs(d, cfile, e2020db) for d in ds] for ds in split_dims]

  DataFrame(Variable=vars, Dimensions=dims, Description=desc, Database=db_assignments, DPairs=dim_pairs)
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
  "PE", "YT", "NT", "NU"
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
    @warn needle, " not found in vars, perhaps your variable is in the list below\n"
    temp2 = findall(occursin.(lowercase.(vars.Variable), lowercase(needle)))
    println(vars[temp2,:])
    error("Variable not found")
    return (vars[temp2,1:5])
  elseif unique(temp, [:Variable, :Dimensions]) == 1
    print("Combining ", nrow(temp), "variables")
    return (map(x -> var_id(x, vars, DATA_FOLDER), temp.RowID))
  else
    @error needle, " could have several values, use var_id to select one of the below\n"
    return (temp)
  end
end

function var_legacy(needle, loc::Loc_p)
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
      # println(needle)
      var(needle, loc)
    else
      @warn "Multiple variables exist, please specify database"
      df = vars[i, :]
      print(df)
      @error "Please specify a database"
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

function lookup_database(name, loc; sec::Char=' ', verbose = false)
  verbose && println("Looking up: '$name'")
  
  if contains(name, "/")
    verbose && println("'$name' contains a slash")
    name_parts = string.(split(name, "/"))
    variable_name = name_parts[2]
    database_name = name_parts[1]
    println("  -> Variable: '$variable_name', Database: '$database_name'")
    return (variable_name, database_name)
  else
    verbose && println("'$name' does not contain a slash")
    vars = loc.vars
    temp2 = find_var(name, vars; exact=true)
    n = nrow(temp2)
    
    if n == 1
      result = (temp2.Variable[1], temp2.Database[1])
      verbose && println("  -> Found exact match: $result")
      return result
      
    elseif n > 1 && sec == ' '  # Changed from "" to ' '
      @warn "'$name' could have several values, select one of the below:"
      println(temp2)
      @error "Variable not found - multiple matches without sector specification"
      
    elseif n > 1 && sec != ' '  # Changed from "" to ' '
      i = findall(first.(temp2.Database) .== sec)
      verbose && println("Found these options:")
      verbose && println(temp2)
      verbose && println("Looking for sector '$sec', found indices: $i")
      
      if length(i) == 1
        ind = i[1]
        result = (string(temp2[ind, :Variable]), temp2[ind, :Database])
        verbose && println("  -> Selected: $result")
        return result
      else
        @warn "'$name' could have several values, select one of the below:"
        println(temp2)
        @error "Variable not found - multiple or no matches for sector '$sec'"
      end
      
    else
      # Try fuzzy search
      temp = find_var(name, vars)
      m = nrow(temp)
      
      if m == 1
        result = (temp.Variable[1], temp.Database[1])
        @warn "  -> Found fuzzy match: $result"
        return result
        
      elseif m == 0
        @warn "'$name' not found in vars, perhaps your variable is in the list below:"
        temp2_indices = findall(occursin.(lowercase(name), lowercase.(vars.Variable)))
        println(vars[temp2_indices, [:Variable, :Database]])
        @error "Variable not found - no matches"
        
      else
        @warn "'$name' has multiple fuzzy matches:"
        println(temp)
        @error "Variable not found - multiple fuzzy matches"
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
            UnCode[UnCode .== "Null"] .= ""
            sets = merge(sets, (Unit = UnCode,))
        end
        
        # Clean up sets and compress array
        cleaned_sets, index_mapping = clean_duplicate_dimensions(sets)
        arr = compress_array(arr, index_mapping)
        
        # Convert all set components to categorical
        categorical_sets = NamedTuple{keys(cleaned_sets)}(
            Tuple(categorical(getproperty(cleaned_sets, k)) for k in keys(cleaned_sets))
        )
        sets = categorical_sets
    else
        if vattr["type"] != "set"
            @warn string(dbname,"/", vname), " is of type ", vattr["type"], ". This is unusual."
        end
        values = copy(arr)
        idx = findfirst("Key", vname)
        if isnothing(idx)
            key = [Symbol(vname)]
        else
            key = [Symbol(vname[setdiff(1:length(vname), idx)])]
        end
        
        # For set types, also clean duplicates
        unique_values = unique(values)
        sets = NamedTuple{Tuple(key)}(Tuple([categorical(unique_values)]))
        # Update arr to only contain unique values
        arr = unique_values
    end
    
    # Handle TimeP and TimeA transformations while keeping categorical
    if :TimeP ∈ keys(sets)
        timeP_vals = getproperty(sets, :TimeP)
        if length(timeP_vals) > 0 && '(' ∉ timeP_vals[1]
            new_timeP = categorical([string("TimeP(", match(r"\d+", tp).match, ")") for tp in timeP_vals])
            sets = merge(sets, (TimeP = new_timeP,))
        end
    end
    
    if :TimeA ∈ keys(sets)
        timeA_vals = getproperty(sets, :TimeA)
        if length(timeA_vals) > 0 && '(' ∉ timeA_vals[1]
            new_timeA = categorical([string("TimeA(", match(r"\d+", ta).match, ")") for ta in timeA_vals])
            sets = merge(sets, (TimeA = new_timeA,))
        end
    end
    
    return(arr, sets)
end

function arr_set(vname::String, dbname::String, loc::Loc_p)
    arr = P.data(joinpath(loc.DATA_FOLDER, string(dbname, ".dba")), vname)
    idx = findfirst(loc.vars.Variable .== vname .&& loc.vars.Database .== dbname)
    keys = Symbol.(split(loc.vars.Dimensions[idx], ","))
    dim_pairs = loc.vars.DPairs[idx]
    
    # Get the dimension values
    values = [P.data(key, value) for (key, value) in dim_pairs]
    
    # Clean up dimensions and compress array
    cleaned_values, index_mapping = clean_duplicate_dimensions_from_vectors(values)
    arr = compress_array(arr, index_mapping)
    
    # Convert to categorical
    categorical_values = [categorical(vals) for vals in cleaned_values]
    sets = NamedTuple{Tuple(keys)}(Tuple(categorical_values))
    
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
                   fltr::Dict{Symbol, <:Any}=Dict{Symbol,Any}(),
                   sec::Char=' ')
  vname, dbname = lookup_database(vname, loc; sec)
  arr, sets = arr_set(vname, dbname, loc)
  # set == set2 ? sets = set : error("Sets Don't Match")
  if !isempty(fltr)
    arr, sets = subset_array(arr, sets, fltr)
  end
  df = to_tidy_dataframe(arr, sets; Value = vname)
  return df
end

function var(vname::String, locs::Vector{<:Location}; 
             fltr::Dict{Symbol, <:Any}=Dict{Symbol,Any}(),
             sec::Char=' ',
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
    vname_processed, dbname = lookup_database(vname, loc; sec) 
    arr, set = arr_set(vname_processed, dbname, loc)
    push!(arrs, arr)
    push!(sets, set)
  end
  
  if allequal(sets) 
    (set = sets[1]) 
  
    if !isempty(fltr)
      # Apply filter to all arrays and get the filtered set (same for all)
      set_filtered = nothing
      for i in eachindex(arrs)
        arrs[i], current_set = subset_array(arrs[i], set, fltr)
        if set_filtered === nothing
          set_filtered = current_set  # Store the first filtered set
        end
      end
      set = set_filtered  # Update set to the filtered version
    end
    
    df = to_tidy_dataframe(arrs, set, loc_names)
  else
    @warn "Sets Don't Match Before fltr" 
    if !isempty(fltr)
      for i in eachindex(arrs)
        arrs[i], sets[i] = subset_array(arrs[i], sets[i], fltr)
      end
    end
    if allequal(sets) 
      @warn "Sets match After fltr"
      set = sets[1]
      df = to_tidy_dataframe(arrs, set, loc_names)
    else 
      @warn "Sets still don't match after fltr"
      diagnose_set_mismatch(sets, loc_names; fltr)
      dfs = [to_tidy_dataframe(arrs[i],sets[i]; Value = loc_names[i]) for i in eachindex(arrs)]
      df = join_vars(dfs...)
    end
  end
  
  if diff !== false
    df = add_differences!(df, locs, diff)
  end
  
  if pdiff !== false
    df = add_percent_differences!(df, locs, pdiff)
  end
  
  return df  # Don't forget to return the result!
end

function compare_vars(vnames::Vector{String}, locs::Vector{<:Location}; 
             fltr::Dict{Symbol, <:Any}=Dict{Symbol,Any}(),
             sec::Char=' ',
             diff::Union{Bool, Symbol, Vector{Int}}=false,
             pdiff::Union{Bool, Symbol, Vector{Int}}=false)
  values = Vector{DataFrame}()
for (i, vname) in enumerate(vnames)
    println("\n=== Processing iteration $i ===")
    println("Current vname: '$vname'")
    println("vnames vector: $vnames")
    
    result = var(vname, locs; fltr, sec, diff, pdiff)
    push!(values, result)
    
    println("vnames vector after var(): $vnames")
  end
  variables = Dict(zip(vnames, values))  # Simple Dict instead of NamedTuple
  summary = DataFrame(Name = vnames)
  for loc in locs
    summary[!, Symbol(loc.name)] = [sum(x[:,Symbol(loc.name)]) for x in values]  # or some default value
  end
  if diff != false
    diff_idx = contains.(names(values[1]), "minus")
    for dname ∈ names(values[1])[diff_idx]
      summary[!, Symbol(dname)] = [sum(abs.(x[:,Symbol(dname)])) for x in values]
    end
  end
  if pdiff != false
    pdiff_idx = contains.(names(values[1]), "pdiff")
    for dname ∈ names(values[1])[pdiff_idx]
      summary[!, Symbol(dname)] = [sum(abs.(x[:,Symbol(dname)])) for x in values]
    end
  end
  return(variables,summary)
end

function diff_fast(vname::String, loc1, loc2; 
                   fltr::Dict{Symbol, <:Any}=Dict{Symbol,Any}(),
                   sec::Char=' ')
  vname1, dbname1 = lookup_database(vname, loc1; sec)
  vname2, dbname2 = lookup_database(vname, loc2; sec)
  arr1, set1 = arr_set(vname1, dbname1, loc1)
  arr2, set2 = arr_set(vname2, dbname2, loc2)
  set1 == set2 ? sets = set1 : error("Sets Don't Match")
  arr1, set1 = subset_array(arr1, sets, fltr)
  arr2, set2 = subset_array(arr2, sets, fltr)
  df1 = to_tidy_dataframe(arr1, set1)
  df2 = to_tidy_dataframe(arr2, set2)
  df = diff(df1, df2; name1 = loc1.name, name2 = loc2.name)
end

function join_vars(df1::DataFrame, df2::DataFrame)
  data = names(df1)[end]
  dims = intersect(names(df1), names(df2))
  dims = dims[dims.!="Value"]
  df = outerjoin(df1, df2, on=dims, makeunique=true)
  values = [setdiff(names(df2), dims); setdiff(names(df1), dims)]
  # for c ∈ eachcol(df[!, Symbol.(values)])
  #   replace!(c, missing => 0)
  # end
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

function plot_diff(data; dim="ECC", num=10, 
  title::String = "",
  units::String = "")

  dim = dim isa String ? dim : string(dim)
  df = deepcopy(data)
  ss = tops(df; dim, num)
  others = setdiff(df[:, dim], ss)
  df[in.(df[:, dim], Ref(others)), dim] .= "Other"
  df = @by(df, [Symbol(dim), :Year], :Diff = sum(:Diff))
  cats = categorical(df[:, dim])
  colors = distinguishable_colors(length(unique(cats)))
  fig = Figure()
  ax = Axis(fig[1, 1]; title=title, ylabel = units, xlabel = "Year")
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
  title::String = "",
  units::String = "")
  
  df = deepcopy(data)
  
  # Convert Year from categorical to numeric if needed
  if df.Year isa CategoricalArray
    df.Year = parse.(Int, Vector(df.Year))
  end
  
  if col == ""
    l = last(names(df))
  else 
    l = col
  end
  l_symbol = Symbol(l)
  
  # Filter out zero values - handle categorical columns properly
  if df[!, l_symbol] isa CategoricalArray
    @rsubset! df Vector($l_symbol) != 0
  else
    @rsubset! df $l_symbol != 0
  end
  
  df2 = @by(df, [Symbol(dim)], :Value = sum(abs.(Vector($l_symbol))))
  df2 = @orderby(df2, -:Value)
  ss = first(df2[!, dim], num)
  others = setdiff(df[:, dim], ss)
  df[in.(df[:, dim], Ref(others)), dim] .= "Other"
  
  # Group and sum
  df = @by(df, [Symbol(dim), :Year], $l_symbol = sum(Vector($l_symbol)))
  
  df[!, dim] = droplevels!(df[!, dim])
  cats = categorical(df[:, dim])
  unique_categories = unique(df[:, dim])  # Only categories in final data
  n_categories = length(unique_categories)
  colors = distinguishable_colors(n_categories)
  labels = String.(unique_categories)  
  
  fig = Figure()
  ax = Axis(fig[1, 1]; title=title, ylabel = units, xlabel = "Year")
  
  # Use a mapping from level codes to sequential indices
  level_to_index = Dict(level => i for (i, level) in enumerate(levels(cats)))
  color_indices = [level_to_index[level] for level in cats]
  
  barplot!(ax, df.Year, df[:, l], stack=color_indices, color=colors[color_indices])
  
  elements = [PolyElement(polycolor=colors[i]) for i in 1:length(labels)]
  Legend(fig[1, 2], elements, labels, String(dim))
  display(fig)
end

function plot_lines(data::DataFrame, cols::AbstractVector{<:Union{Symbol,String,Location}};
  title::String = "",
  units::String = "")
  
  if typeof(cols[1]) <: Location
    cols = [Symbol(loc.name) for loc in cols]
  elseif typeof(cols[1]) == String
    cols = [Symbol(col) for col in cols]
  end
  
  df = dropmissing(data)
  dfg = combine(groupby(df, :Year), cols .=> sum, renamecols=false)
  if nrow(dfg) == 1
    println(dfg)
    @error "plot_lines requires multiple years. Only one year of data is present."
    return
  end
  dfg = stack(dfg, cols)
  sort!(dfg, :Year)
  if dfg.Year isa CategoricalArray
    dfg.Year = parse.(Int, Vector(dfg.Year))
  end
  println(first(dfg, 5))
  println(typeof(dfg.Year))

  
  fig = Figure()
  ax = Axis(fig[1, 1], 
           xlabel = "Year",
           ylabel = units, 
           title = title)
  
  # Plot each group separately
  for var in unique(dfg.variable)
    subset_df = filter(row -> row.variable == var, dfg)
    lines!(ax, subset_df.Year, subset_df.value, 
           label = var,
           linewidth = 5,
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
    subset_array(array::AbstractArray, sets::NamedTuple, fltr::Dict{Symbol, <:Any})

Filters a multi-dimensional array based on dimension names and values.

# Arguments
- `array::AbstractArray`: The array to filter
- `sets::NamedTuple`: Named tuple containing the dimension values for each axis
- `fltr::Dict{Symbol, <:Any}`: Dictionary mapping dimension names to filter values

# Returns
- `Tuple{Array, NamedTuple}`: A tuple containing the filtered array and the filtered dimension sets
"""
function subset_array(array::AbstractArray, sets::NamedTuple, fltr::Dict{Symbol, <:Any})
    dim_names = collect(propertynames(sets))
    
    if isempty(fltr)
        return array, sets
    end
    
    indices_by_dim = Vector{Vector{Int}}()
    filtered_sets_values = Vector{Any}()
    
    for dim_name in dim_names
        dim_values = getproperty(sets, dim_name)  # This is categorical
        
        if haskey(fltr, dim_name)
            condition = fltr[dim_name]
            
            # Optimized filtering for categorical arrays
            if condition isa AbstractVector
                # Use Set for O(1) lookup instead of O(n) for each element
                condition_set = Set(condition)
                indices = findall(x -> x in condition_set, dim_values)
            elseif condition isa Function
                indices = findall(condition, dim_values)
            else
                # Single value equality
                indices = findall(==(condition), dim_values)
            end
            
            push!(indices_by_dim, indices)
            # Keep as categorical with reduced levels
            filtered_vals = dim_values[indices]
            push!(filtered_sets_values, categorical(filtered_vals, levels=levels(filtered_vals)))
        else
            # No filter - keep all indices
            indices = collect(1:length(dim_values))
            push!(indices_by_dim, indices)
            push!(filtered_sets_values, dim_values)  # Keep original categorical
        end
    end
    
    # Extract subarray
    indexing = Tuple(indices_by_dim)
    result = array[indexing...]
    
    # Reconstruct sets as NamedTuple
    filtered_sets = NamedTuple{Tuple(dim_names)}(Tuple(filtered_sets_values))
    
    return result, filtered_sets
end

"""
    subset_dataframe(df::DataFrame, fltr::Dict; drop_filtered_cols::Bool=false)

Filter a DataFrame based on conditions specified in a dictionary.

# Arguments
- `df`: The DataFrame to filter
- `fltr`: Dictionary where keys are column names and values are the filter conditions
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
function subset_dataframe(df_in::DataFrame, fltr::Dict{<:Any,<:Any}; drop_filtered_cols::Bool=false)
  df = copy(df_in)
  # Check if the DataFrame is empty
  isempty(df) && return df
  
  # Keep track of columns we've filtered on
  filtered_cols = Symbol[]
  
  # Apply each filter
  for (col_key, condition) in fltr
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
    subset_dataframe!(df::DataFrame, fltr::Dict; drop_filtered_cols::Bool=false)

In-place version of subset_dataframe that modifies the input DataFrame directly.

# Arguments and behavior are the same as subset_dataframe
"""
function subset_dataframe!(df::DataFrame, fltr::Dict; drop_filtered_cols::Bool=false)
  # Check if the DataFrame is empty
  isempty(df) && return df
  
  # Keep track of columns we've filtered on
  filtered_cols = Symbol[]
  
  # Apply each filter
  for (col_key, condition) in fltr
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
    subset_dataframe_optimized(df::DataFrame, fltr::Dict; drop_filtered_cols::Bool=false)

Optimized version of subset_dataframe that uses DataFrames.jl's built-in filtering
capabilities more efficiently. This version is faster for large DataFrames.

# Arguments and behavior are the same as subset_dataframe
"""
function subset_dataframe_optimized(df::DataFrame, fltr::Dict; drop_filtered_cols::Bool=false)
    # Check if the DataFrame is empty
    isempty(df) && return df
    
    result_df = df
    filtered_cols = String[]
    
    for (col, condition) in fltr
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
filtered_array, filtered_sets = subset_array(HDGCCI_j, sets, fltr)
df = to_tidy_dataframe(filtered_array, filtered_sets)
```
"""

function to_tidy_dataframe(array::AbstractArray, sets::NamedTuple; Value::Union{Symbol,String}=:Value)
    Value = Symbol(Value)
    dim_names = collect(propertynames(sets))
    
    if ndims(array) != length(dim_names)
        error("Number of dimensions in array ($(ndims(array))) doesn't match number of sets ($(length(dim_names)))")
    end
    
    # Extract categorical dimension values
    values_by_dim = [getproperty(sets, dim_name) for dim_name in dim_names]
    
    # Pre-calculate total elements
    n_elements = length(array)
    
    # Pre-allocate all columns at once
    dim_columns = [CategoricalVector{eltype(vals)}(undef, n_elements) for vals in values_by_dim]
    value_column = Vector{eltype(array)}(undef, n_elements)
    
    # Single pass through all indices
    idx_counter = 1
    for idx in CartesianIndices(array)
        # Fill dimension columns (categorical)
        for (d, dim_idx) in enumerate(Tuple(idx))
            dim_columns[d][idx_counter] = values_by_dim[d][dim_idx]
        end
        # Fill value column (non-categorical)
        value_column[idx_counter] = array[idx]
        idx_counter += 1
    end
    
    # Create DataFrame with pre-allocated columns
    df = DataFrame()
    for (i, name) in enumerate(dim_names)
        df[!, Symbol(name)] = dim_columns[i]
    end
    df[!, Value] = value_column
    
    # Handle Year parsing if needed
    if :Year ∈ dim_names && eltype(values_by_dim[findfirst(==(Symbol(:Year)), dim_names)]) <: AbstractString
        df.Year = parse.(Int64, df.Year)
    end
    
    return df
end

function to_tidy_dataframe(arrays::Vector{<:AbstractArray}, sets::NamedTuple, value_names::Vector{String})
    if isempty(arrays)
        error("At least one array must be provided")
    end
    
    # Quick validation
    ref_size = size(arrays[1])
    for (i, arr) in enumerate(arrays)
        if size(arr) != ref_size
            error("Array $i has size $(size(arr)), but expected $(ref_size)")
        end
    end
    
    if length(value_names) != length(arrays)
        error("Number of value names ($(length(value_names))) doesn't match number of arrays ($(length(arrays)))")
    end
    
    dim_names = collect(propertynames(sets))
    values_by_dim = [getproperty(sets, dim_name) for dim_name in dim_names]
    
    n_elements = length(arrays[1])
    
    # Pre-allocate: categorical for dimensions, regular vectors for values
    dim_columns = [CategoricalVector{eltype(vals)}(undef, n_elements) for vals in values_by_dim]
    value_columns = [Vector{eltype(arr)}(undef, n_elements) for arr in arrays]
    
    # Single pass fill
    idx_counter = 1
    for idx in CartesianIndices(arrays[1])
        # Fill dimension columns (same for all arrays)
        for (d, dim_idx) in enumerate(Tuple(idx))
            dim_columns[d][idx_counter] = values_by_dim[d][dim_idx]
        end
        # Fill value columns (one per array)
        for (a, arr) in enumerate(arrays)
            value_columns[a][idx_counter] = arr[idx]
        end
        idx_counter += 1
    end
    
    # Build DataFrame
    df = DataFrame()
    for (i, name) in enumerate(dim_names)
        df[!, Symbol(name)] = dim_columns[i]
    end
    for (i, value_name) in enumerate(value_names)
        df[!, Symbol(value_name)] = value_columns[i]
    end
    
    # Handle Year parsing
    if :Year ∈ dim_names && eltype(values_by_dim[findfirst(==(Symbol(:Year)), dim_names)]) <: AbstractString
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
    ReadDiskRaw(db::String, name::String; fltr::Dict = Dict())

Reads the dataset named `name` from the HDF5 file specified by `db` and returns
the raw array with filters applied at read time.

# Arguments
- `db::String`: Path to the HDF5 file
- `name::String`: Name of the dataset to read
- `fltr::Dict=Dict()`: A dictionary mapping dimension names to filter criteria

# Returns
- `Tuple{Array, Dict}`: A tuple containing the filtered array and a dictionary 
  mapping dimension names to their filtered values

# Example
```julia
# Read only data for specific years
data, dim_info = ReadDiskRaw("data.h5", "my_dataset", 
                             fltr = Dict("Year" => [2020, 2030]))
```
"""
function ReadDiskRaw(db::String, name::String; fltr::Dict = Dict())
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
        
        if haskey(fltr, dim)
            filter_criterion = fltr[dim]
            
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

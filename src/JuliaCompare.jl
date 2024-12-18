module JuliaCompare

import CSV
using DataFrames, Chain, DataFramesMeta, HDF5
import PromulaDBA as P
using Makie, CairoMakie
using Colors, CategoricalArrays
import SmallModel: ReadDisk

greet() = print("Hello Randy")

function filter_ns(df_in)
  df = copy(df_in)
  if "Area" ∈ names(df)
    @rsubset! df :Area ∈ ["NS", "PE", "NB"]
  end
  if "Enduse" ∈ names(df)
    @rsubset! df :Enduse == "Heat"
  end
  if "EC" ∈ names(df)
    @rsubset! df :EC == "IronSteel"
  end
  if "Tech" ∈ names(df)
    @rsubset! df :Tech == "Biomass"
  end
  if "Year" ∈ names(df)
    @rsubset! df :Year == 2022
  end
  return (df)
end

function filter_ne(df_in)
  df = copy(df_in)
  if "Area" ∈ names(df)
    @rsubset! df :Area ∈ ["NEng", "MAtl"]
  end
  if "Enduse" ∈ names(df)
    @rsubset! df :Enduse == "Heat"
  end
  if "EC" ∈ names(df)
    @rsubset! df :EC == "Petrochemicals"
  end
  if "Tech" ∈ names(df)
    @rsubset! df :Tech ∈ ["Electric", "Gas", "Coal", "Oil", "LPG"]
  end
  if "Year" ∈ names(df)
    @rsubset! df :Year == 2022
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
  df = df[(e[1]+1):(e[2]-1), :]
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
  vars[i, 1:5]
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
    print(needle, " not found in vars, perhaps your variable is in the list below\n")
    temp2 = findall(occursin.(lowercase.(vars.Variable), lowercase(needle)))
    return (temp2)
  elseif unique(temp, [:Variable, :Dimensions]) == 1
    print("Combining ", nrow(temp), "variables")
    return (map(x -> var_id(x, vars, DATA_FOLDER), temp.RowID))
  else
    print(needle, " could have several values, use var_id to select one of the below\n")
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
    return (ReadDisk(DataFrame, loc.HDF5_path, needle))
  else
    i = findall(lowercase.(vars.Variable) .== lowercase(needle))
    if sum(i) == 1
      needle = vars.Database[i] * "/" * needle
      var(needle, loc)
    else
      print("Multiple variables exist, please specify database")
      df = vars[i, :]
      print(df)
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

end # module JuliaCompare

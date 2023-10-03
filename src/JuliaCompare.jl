module JuliaCompare

import CSV
using DataFrames, Chain, DataFramesMeta
import PromulaDBA as P
using Makie, CairoMakie
using Colors, CategoricalArrays

greet() = print("Hello Randy")

struct Loc
  vars::DataFrame
  DATA_FOLDER::String
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

  if contains(cfilename, "input") | contains(cfilename, "output")
    letter = SubString(cfilename, 1, 1)
    cfile = joinpath(DATA_FOLDER, join([letter, "input.dba"]))
  else
    cfile = e2020db
  end

  makepairs = function (d, cfile, e2020db)
    if d == "Enduse"
      Pair(cfile, "EUDS")
    elseif d in sets
      Pair(cfile, join([d, "DS"]))
    elseif d == "Year"
      Pair(e2020db, "YrvDS")
    elseif d == "FuelEP"
      Pair(e2020db, "FlEPDS")
    elseif d == "Unit"
      Pair(eg, "UnName")
    else
      Pair(e2020db, join([d, "DS"]))
    end
  end

  dim_pairs = [[makepairs(d, cfile, e2020db) for d in ds] for ds in split_dims]

  DataFrame(Variable=vars, Dimensions=dims, Description=desc, Database=cfilename, DPairs=dim_pairs)
end

function get_var(i, vars, DATA_FOLDER)
  fname = joinpath(DATA_FOLDER, string(vars.Database[i], ".dba"))

  vname = string(vars.Variable[i])
  opts = vars.DPairs[i]
  print("fname is: ", fname)
  print("\n")
  # out = P.dataframe(fname, vname, opts..., skip_zeros=false)
  out = P.dataframe(fname, vname, opts...)
  rename!(out,
    Dict(zip(names(out)[1:end-1], split(vars.Dimensions[i], ',')))
  )
  if contains(vars.Dimensions[i], "Year")
    out.Year = parse.(Int64, out.Year)
  end
  return (out)
end

function get_data(i, vars, DATA_FOLDER)
  fname = joinpath(DATA_FOLDER, string(vars.Database[i], ".dba"))
  vname = string(vars.Variable[i])
  P.data(fname, vname)
end

function list_vars(CODE_FOLDER, DATA_FOLDER, db_files)
  dfs = map(x -> list_var(x, CODE_FOLDER, DATA_FOLDER, false), db_files)
  vars = vcat(dfs...)

  vars.RowID = 1:size(vars, 1)
  vars = @chain vars begin
    select!(_, :RowID, Not(:RowID)) # Move ID column to the front
  end
  return (vars)
end

function find_var(needle, vars)
  i = findall(occursin.(lowercase.(vars.Variable), lowercase(needle)))
  vars[i, 1:4]
end

const Canada = ["Ontario", "Quebec", "British Columbia", "Alberta", "Manitoba",
  "Saskatchewan", "New Brunswick", "Nova Scotia", "Newfoundland",
  "Prince Edward Island", "Yukon Territory", "Northwest Territory", "Nunavut"
]

function var(needle, vars, DATA_FOLDER)
  temp = find_var(needle, vars)
  if nrow(temp) == 1
    return (get_var(temp.RowID[1], vars, DATA_FOLDER))
  elseif nrow(temp) == 0
    print(Needle, " not found in vars")
  elseif unique(temp, [:Variable, :Dimensions]) == 1
    print("Combining ", nrow(temp), "variables")
    return (map(x -> get_var(x, vars, DATA_FOLDER), temp.RowID))
  else
    print(needle, " could have several values, use get_var to select one of the below")
    return (temp)
  end
end

function var(needle, loc::Loc)
  (; vars, DATA_FOLDER) = loc
  return (var(needle, vars, DATA_FOLDER))
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

function plot_diff(data; dim="ECC", num=10)
  df = deepcopy(data)
  ss = tops(df; dim, num)
  others = setdiff(df[:, dim], ss)
  df[in.(df[:, dim], Ref(others)), dim] .= "Other"
  df = @by(df, [Symbol(dim), :Year], :Diff = sum(:Diff))
  cats = categorical(df[:, dim])
  colors = distinguishable_colors(length(unique(cats)))
  fig = Figure()
  ax = Axis(fig[1, 1])
  barplot!(ax, df.Year, df.Diff, stack=levelcode.(cats), color=colors[levelcode.(cats)])
  labels = levels(cats)
  # return (labels)
  elements = [PolyElement(polycolor=colors[i]) for i in 1:length(labels)]
  Legend(fig[1, 2], elements, labels, dim)
  display(fig)
end

end # module JuliaCompare

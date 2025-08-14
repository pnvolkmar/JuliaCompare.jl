function list_var(cfilename, CODE_FOLDER, DATA_FOLDER, verbose=false)
  if verbose
    print("cfilename is: ", cfilename, "\n")
  end

  code_path = joinpath(CODE_FOLDER, join([cfilename, ".src"]))
  df = CSV.read(code_path, DataFrame, delim="\t", comment="*", header=[:x])
  
  # Find all "Open Input" lines to track database changes
  open_input_lines = findall(occursin.(r"Open Input.*\.dba", df.x))
  
  # Extract database names from "Open Input" lines
  db_names = String[]
  for line_idx in open_input_lines
    line = df.x[line_idx]
    # Extract database name between quotes
    match_result = match(r"Open Input\s+\"([^\"]+)\"", line)
    if match_result !== nothing
      db_name = replace(match_result.captures[1], ".dba" => "")
      push!(db_names, db_name)
    end
  end
  
  # Find variable definition blocks
  a = occursin.("Define Variable", df.x)
  define_var_lines = findall(a)
  end_define_lines = findall(occursin.("End Define Variable", df.x))
  
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
    elseif contains(d, "TOM")
      Pair(kinput, join([d, "Key"]))
    else
      Pair(e2020db, join([d, "Key"]))
    end
  end
  
  dim_pairs = [[makepairs(d, cfile, e2020db) for d in ds] for ds in split_dims]
  
  DataFrame(Variable=vars, Dimensions=dims, Description=desc, Database=db_assignments, DPairs=dim_pairs)
end

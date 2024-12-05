using HDF5
import SmallModel as M
db = loc3
function print_structure(file::String)
  h5open(file, "r") do fid
    show_tree(fid)
  end
end

# Alternative for more detailed information:
function show_tree(g::Union{HDF5.File,HDF5.Group}, prefix::String="")
  for name in keys(g)
    obj = g[name]
    println(prefix, name)
    if isa(obj, HDF5.Group)
      show_tree(obj, prefix * "  ")
    elseif isa(obj, HDF5.Dataset)
      println(prefix * "  └ Type: ", typeof(read(obj)))
      println(prefix * "  └ Size: ", size(read(obj)))
    end
  end
end

# Usage:
print_structure(loc3)

h5open(loc3, "r") do fid
  for name in keys(fid)
    println(name)
  end
end

using HDF5, DataFrames

function create_hdf5_map(file::String)
  data = Vector{NamedTuple{(:path, :type, :size),Tuple{String,Type,Tuple}}}()

  function traverse(g::Union{HDF5.File,HDF5.Group}, path::String="/")
    for name in keys(g)
      full_path = path == "/" ? "/$name" : "$path/$name"
      obj = g[name]

      if isa(obj, HDF5.Dataset)
        push!(data, (path=full_path, type=typeof(read(obj)), size=size(read(obj))))
      elseif isa(obj, HDF5.Group)
        traverse(obj, full_path)
      end
    end
  end

  h5open(file, "r") do fid
    traverse(fid)
  end

  return DataFrame(data)
end

# Usage:
df = J.create_hdf5_map(loc3)

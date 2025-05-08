"""
    arrays_to_tidy_df(arrays::AbstractArray...; 
                      value_names = ["Value$i" for i in 1:length(arrays)],
                      dim_names = ["Dim$i" for i in 1:ndims(first(arrays))],
                      index_values = nothing)

Convert multiple arrays of the same dimensions into a tidy DataFrame with:
- One row per cell position across all arrays
- Columns for each dimension's index
- One value column per input array

# Arguments
- `arrays`: Two or more arrays with identical dimensions
- `value_names`: Names for the value columns (default: "Value1", "Value2", ...)
- `dim_names`: Names for the dimension columns (default: "Dim1", "Dim2", ...)
- `index_values`: Optional vector of vectors with custom index values for each dimension

# Returns
- A DataFrame in tidy format with dimension indices and array values

# Example
```julia
# Create two arrays
arr1 = [1 2 3; 4 5 6]
arr2 = [0.1 0.2 0.3; 0.4 0.5 0.6]

# Convert to tidy DataFrame with default settings
df = arrays_to_tidy_df(arr1, arr2)

# With custom names
df = arrays_to_tidy_df(arr1, arr2, 
                      value_names=["Count", "Proportion"], 
                      dim_names=["Row", "Col"])

# With custom index values
df = arrays_to_tidy_df(arr1, arr2, 
                      value_names=["Count", "Proportion"], 
                      dim_names=["Year", "Region"],
                      index_values=[["2022", "2023"], ["East", "Central", "West"]])
```
"""
function arrays_to_tidy_df(arrays::AbstractArray...; 
                          value_names = ["Value$i" for i in 1:length(arrays)],
                          dim_names = ["Dim$i" for i in 1:ndims(first(arrays))],
                          index_values = nothing)
    
    # Validate inputs
    isempty(arrays) && throw(ArgumentError("At least one array must be provided"))
    
    first_array = first(arrays)
    dims = size(first_array)
    ndim = ndims(first_array)
    
    for (i, arr) in enumerate(arrays)
        if size(arr) != dims
            throw(DimensionMismatch("All arrays must have the same dimensions. " *
                "Array $i has size $(size(arr)), expected $dims"))
        end
    end
    
    if length(value_names) != length(arrays)
        throw(ArgumentError("Number of value_names ($(length(value_names))) must match " *
            "number of arrays ($(length(arrays)))"))
    end
    
    if length(dim_names) != ndim
        throw(ArgumentError("Number of dim_names ($(length(dim_names))) must match " *
            "number of dimensions ($ndim)"))
    end
    
    # Validate custom index values if provided
    if index_values !== nothing
        if length(index_values) != ndim
            throw(ArgumentError("index_values must have an entry for each dimension"))
        end
        
        for (d, (idxs, dim_size)) in enumerate(zip(index_values, dims))
            if length(idxs) != dim_size
                throw(ArgumentError("index_values for dimension $d has length $(length(idxs)), " *
                    "but the dimension size is $dim_size"))
            end
        end
    end
    
    # Calculate total number of rows in output DataFrame
    n_rows = prod(dims)
    
    # Create the DataFrame with dimension columns
    df = DataFrame()
    
    # Create cartesian indices for all positions
    cart_indices = collect(CartesianIndices(dims))
    
    # Add dimension columns
    for d in 1:ndim
        if index_values !== nothing
            # Use custom index values
            col_values = [index_values[d][cart_idx[d]] for cart_idx in cart_indices]
        else
            # Use numeric indices
            col_values = [cart_idx[d] for cart_idx in cart_indices]
        end
        df[!, dim_names[d]] = col_values
    end
    
    # Add value columns for each array
    for (arr_idx, arr) in enumerate(arrays)
        # Convert to linear indices for faster access
        values = [arr[cart_idx] for cart_idx in cart_indices]
        df[!, value_names[arr_idx]] = values
    end
    
    return df
end

# Alternative implementation optimized for large arrays
"""
    arrays_to_tidy_df_optimized(arrays::AbstractArray...; 
                               value_names = ["Value$i" for i in 1:length(arrays)],
                               dim_names = ["Dim$i" for i in 1:ndims(first(arrays))],
                               index_values = nothing)

Optimized version of arrays_to_tidy_df, more efficient for large arrays.
Uses pre-allocation and avoids collecting CartesianIndices.
"""
function arrays_to_tidy_df_optimized(arrays::AbstractArray...; 
                                    value_names = ["Value$i" for i in 1:length(arrays)],
                                    dim_names = ["Dim$i" for i in 1:ndims(first(arrays))],
                                    index_values = nothing)
    
    # Validation (same as in the previous function)
    isempty(arrays) && throw(ArgumentError("At least one array must be provided"))
    
    first_array = first(arrays)
    dims = size(first_array)
    ndim = ndims(first_array)
    
    for (i, arr) in enumerate(arrays)
        if size(arr) != dims
            throw(DimensionMismatch("All arrays must have the same dimensions. " *
                "Array $i has size $(size(arr)), expected $dims"))
        end
    end
    
    if length(value_names) != length(arrays)
        throw(ArgumentError("Number of value_names ($(length(value_names))) must match " *
            "number of arrays ($(length(arrays)))"))
    end
    
    if length(dim_names) != ndim
        throw(ArgumentError("Number of dim_names ($(length(dim_names))) must match " *
            "number of dimensions ($ndim)"))
    end
    
    if index_values !== nothing
        if length(index_values) != ndim
            throw(ArgumentError("index_values must have an entry for each dimension"))
        end
        
        for (d, (idxs, dim_size)) in enumerate(zip(index_values, dims))
            if length(idxs) != dim_size
                throw(ArgumentError("index_values for dimension $d has length $(length(idxs)), " *
                    "but the dimension size is $dim_size"))
            end
        end
    end
    
    # Calculate total number of rows
    n_rows = prod(dims)
    
    # Pre-allocate columns
    columns = Dict{String, Vector}()
    
    # Pre-allocate dimension columns
    for d in 1:ndim
        if index_values !== nothing
            columns[dim_names[d]] = Vector{eltype(index_values[d])}(undef, n_rows)
        else
            columns[dim_names[d]] = Vector{Int}(undef, n_rows)
        end
    end
    
    # Pre-allocate value columns
    for (arr_idx, arr) in enumerate(arrays)
        columns[value_names[arr_idx]] = Vector{eltype(arr)}(undef, n_rows)
    end
    
    # Fill columns
    row = 1
    # Use nested loops instead of CartesianIndices for better performance
    # This creates the Cartesian product of all indices
    indices = [1:d for d in dims]
    for idx in Base.product(indices...)
        # Fill dimension columns
        for d in 1:ndim
            if index_values !== nothing
                columns[dim_names[d]][row] = index_values[d][idx[d]]
            else
                columns[dim_names[d]][row] = idx[d]
            end
        end
        
        # Fill value columns
        for (arr_idx, arr) in enumerate(arrays)
            columns[value_names[arr_idx]][row] = arr[idx...]
        end
        
        row += 1
    end
    
    # Create DataFrame from columns
    return DataFrame(columns)
end

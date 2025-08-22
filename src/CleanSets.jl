"""
Clean duplicate/empty values from dimension sets and return mapping for array compression
"""
function clean_duplicate_dimensions(sets::NamedTuple)
    dim_names = keys(sets)
    cleaned_values = Vector{Any}()
    index_mapping = Vector{Vector{Int}}()
    
    for dim_name in dim_names
        dim_vals = getproperty(sets, dim_name)
        
        if should_clean_dimension(dim_name, dim_vals)
            unique_vals, mapping = clean_dimension_values(dim_vals)
            push!(cleaned_values, unique_vals)
            push!(index_mapping, mapping)
        else
            # No cleaning needed - keep original
            push!(cleaned_values, dim_vals)
            push!(index_mapping, collect(1:length(dim_vals)))
        end
    end
    
    cleaned_sets = NamedTuple{dim_names}(Tuple(cleaned_values))
    return cleaned_sets, index_mapping
end

function clean_duplicate_dimensions_from_vectors(values::Vector)
    cleaned_values = Vector{Any}()
    index_mapping = Vector{Vector{Int}}()
    
    for vals in values
        # Check if this dimension needs cleaning based on content
        if has_duplicates_or_empties(vals)
            unique_vals, mapping = clean_dimension_values(vals)
            push!(cleaned_values, unique_vals)
            push!(index_mapping, mapping)
        else
            push!(cleaned_values, vals)
            push!(index_mapping, collect(1:length(vals)))
        end
    end
    
    return cleaned_values, index_mapping
end

"""
Determine if a dimension should be cleaned based on its name or content
"""
function should_clean_dimension(dim_name::Symbol, dim_vals)
    # Specific dimensions known to have duplicates/unused slots
    if dim_name in [:OGUnits, :Units, :Unit]
        return true
    end
    
    # Check for empty strings, "Null", or duplicate values
    return has_duplicates_or_empties(dim_vals)
end

function has_duplicates_or_empties(vals)
    # Check for empty strings or "Null" values
    if any(v -> v == "" || v == "Null" || ismissing(v), vals)
        return true
    end
    
    # Check for duplicates
    return length(unique(vals)) < length(vals)
end

"""
Clean a single dimension's values and return unique values plus index mapping
"""
function clean_dimension_values(dim_vals)
    # Find non-empty, non-null, unique values
    valid_indices = findall(v -> v != "" && v != "Null" && !ismissing(v), dim_vals)
    valid_values = dim_vals[valid_indices]
    
    # Get unique values while preserving order
    unique_values = Vector{eltype(valid_values)}()
    old_to_new_mapping = Vector{Int}()
    value_to_new_index = Dict{eltype(valid_values), Int}()
    
    for (i, val) in enumerate(valid_values)
        if !haskey(value_to_new_index, val)
            # First occurrence of this value
            push!(unique_values, val)
            new_idx = length(unique_values)
            value_to_new_index[val] = new_idx
            push!(old_to_new_mapping, new_idx)
        else
            # Duplicate value - map to existing index
            push!(old_to_new_mapping, value_to_new_index[val])
        end
    end
    
    # Create full mapping for all original indices (including invalid ones)
    full_mapping = Vector{Int}(undef, length(dim_vals))
    mapping_idx = 1
    for (orig_idx, val) in enumerate(dim_vals)
        if orig_idx in valid_indices
            full_mapping[orig_idx] = old_to_new_mapping[mapping_idx]
            mapping_idx += 1
        else
            # Invalid values map to 0 (will be filtered out)
            full_mapping[orig_idx] = 0
        end
    end
    
    return unique_values, full_mapping
end

"""
Compress array by removing unused dimension slots
"""
function compress_array(arr::AbstractArray, index_mappings::Vector{Vector{Int}})
    # Find which indices to keep (non-zero mappings)
    valid_indices = Vector{Vector{Int}}()
    
    for mapping in index_mappings
        valid_idx = findall(x -> x > 0, mapping)
        push!(valid_indices, valid_idx)
    end
    
    # Extract the compressed array
    compressed_arr = arr[valid_indices...]
    
    return compressed_arr
end

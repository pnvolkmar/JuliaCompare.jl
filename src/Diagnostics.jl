function diagnose_set_mismatch(sets::Vector{NamedTuple}, loc_names::Vector{String}; 
                              fltr::Dict{Symbol, <:Any}=Dict{Symbol,Any}())
    println("\n=== SET MISMATCH DIAGNOSTIC ===")
    println("Locations: ", join(loc_names, ", "))
    println("Applied filter: ", fltr)
    
    # Check if all sets have the same dimension names
    all_dim_names = [collect(propertynames(s)) for s in sets]
    if !allequal(all_dim_names)
        println("\n❌ DIMENSION NAMES MISMATCH:")
        for (i, (loc, dims)) in enumerate(zip(loc_names, all_dim_names))
            println("  $loc: $(dims)")
        end
        return
    end
    
    dim_names = all_dim_names[1]
    println("\n✅ All sets have same dimensions: $(dim_names)")
    
    # Check each dimension for mismatches
    for dim_name in dim_names
        dim_values = [getproperty(s, dim_name) for s in sets]
        dim_lengths = [length(vals) for vals in dim_values]
        
        println("\n--- Dimension: $dim_name ---")
        println("Lengths: ", Dict(zip(loc_names, dim_lengths)))
        
        if !allequal(dim_lengths)
            println("❌ LENGTH MISMATCH for $dim_name")
            for (loc, vals) in zip(loc_names, dim_values)
                println("  $loc ($(length(vals))): $(vals[1:min(5, length(vals))])$(length(vals) > 5 ? "..." : "")")
            end
        elseif !allequal(dim_values)
            println("❌ VALUE MISMATCH for $dim_name")
            # Find differences
            ref_values = dim_values[1]
            for (i, (loc, vals)) in enumerate(zip(loc_names, dim_values))
                if i == 1 continue end
                differences = setdiff(vals, ref_values)
                missing_vals = setdiff(ref_values, vals)
                if !isempty(differences) || !isempty(missing_vals)
                    println("  $loc vs $(loc_names[1]):")
                    if !isempty(differences)
                        println("    Extra values: $(differences[1:min(3, length(differences))])$(length(differences) > 3 ? "... ($(length(differences)-3) more)" : "")")
                    end
                    if !isempty(missing_vals)
                        println("    Missing values: $(missing_vals[1:min(3, length(missing_vals))])$(length(missing_vals) > 3 ? "... ($(length(missing_vals)-3) more)" : "")")
                    end
                end
            end
        else
            println("✅ $dim_name matches across all locations")
        end
    end
    
    # Check data types
    println("\n--- Data Types ---")
    for dim_name in dim_names
        dim_types = [eltype(getproperty(s, dim_name)) for s in sets]
        if !allequal(dim_types)
            println("❌ TYPE MISMATCH for $dim_name:")
            for (loc, typ) in zip(loc_names, dim_types)
                println("  $loc: $typ")
            end
        end
    end
end

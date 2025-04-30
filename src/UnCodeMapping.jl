# Create a mapping from numeric codes to descriptive names
function create_ecc_mapping()
    descriptions = [
        "SingleFamilyDetached",
        "SingleFamilyAttached",
        "MultiFamily",
        "OtherResidential",
        "Wholesale",
        "Retail",
        "Warehouse",
        "Information",
        "Offices",
        "Education",
        "Health",
        "OtherCommercial",
        "NGDistribution",
        "OilPipeline",
        "NGPipeline",
        "StreetLighting",
        "Food",
        "Textiles",
        "Lumber",
        "Furniture",
        "PulpPaperMills",
        "Petrochemicals",
        "IndustrialGas",
        "OtherChemicals",
        "Fertilizer",
        "Petroleum",
        "Rubber",
        "Cement",
        "Glass",
        "LimeGypsum",
        "OtherNonMetallic",
        "IronSteel",
        "Aluminum",
        "OtherNonferrous",
        "TransportEquipment",
        "OtherManufacturing",
        "IronOreMining",
        "OtherMetalMining",
        "NonMetalMining",
        "LightOilMining",
        "HeavyOilMining",
        "FrontierOilMining",
        "PrimaryOilSands",
        "SAGDOilSands",
        "CSSOilSands",
        "OilSandsMining",
        "OilSandsUpgraders",
        "ConventionalGasProduction",
        "SweetGasProcessing",
        "UnconventionalGasProduction",
        "SourGasProcessing",
        "LNGProduction",
        "CoalMining",
        "Construction",
        "Forestry",
        "OnFarmFuelUse",
        "CropProduction",
        "AnimalProduction",
        "Passenger",
        "Freight",
        "AirPassenger",
        "AirFreight",
        "ForeignPassenger",
        "ForeignFreight",
        "ResidentialOffRoad",
        "CommercialOffRoad",
        "Miscellaneous",
        "H2Production",
        "UtilityGen",
        "BiofuelProduction",
        "Steam",
        "DirectAirCapture",
        "SolidWaste",
        "Wastewater",
        "Incineration",
        "LandUse",
        "RoadDust",
        "OpenSources",
        "ForestFires",
        "Biogenics"
    ]
    
    # Create mapping dictionaries
    numeric_to_desc = Dict("ECC$i" => "ECC$(descriptions[i])" for i in 1:length(descriptions))
    desc_to_numeric = Dict("ECC$(descriptions[i])" => "ECC$i" for i in 1:length(descriptions))
    
    return numeric_to_desc, desc_to_numeric
end

# Function to transform numeric ECC code to descriptive name
function transform_spruce_to_tanoak(code)
    numeric_to_desc, _ = create_ecc_mapping()
    
    # Split the code parts
    parts = split(code, "_")
    
    if length(parts) >= 3 && startswith(parts[3], "ECC")
        # Extract the ECC part
        ecc_part = parts[3]
        
        # Check if this ECC code exists in our mapping
        if haskey(numeric_to_desc, ecc_part)
            # Replace the ECC part with descriptive version
            descriptive_ecc = numeric_to_desc[ecc_part]
            parts[3] = descriptive_ecc
            
            # Rejoin parts and return
            return join(parts, "_")
        end
    end
    
    # Return original code if we can't transform it
    return code
end

# Function to transform descriptive name to numeric ECC code
function transform_tanoak_to_spruce(code)
    _, desc_to_numeric = create_ecc_mapping()
    
    # Split the code parts
    parts = split(code, "_")
    
    if length(parts) >= 3 && startswith(parts[3], "ECC")
        # Find if any descriptive ECC matches
        for (desc_key, num_value) in desc_to_numeric
            if startswith(parts[3], desc_key)
                # Replace with numeric version
                parts[3] = num_value
                return join(parts, "_")
            end
        end
    end
    
    # Return original code if we can't transform it
    return code
end

# Function to check if codes match after transformation
function check_code_match(spruce_code, tanoak_code)
    # Transform Spruce to Tanoak format
    transformed_spruce = transform_spruce_to_tanoak(spruce_code)
    
    # Check if they match after transformation
    return transformed_spruce == tanoak_code
end

# Apply to your dataframe
function reconcile_codes(df)
    # Create a new column with the transformed Spruce codes
    df.TransformedSpruce = [transform_spruce_to_tanoak(code) for code in df.Spruce]
    
    # Check if transformed codes match Tanoak
    df.MatchAfterTransform = df.TransformedSpruce .== df.Tanoak
    
    return df
end


using HDF5

"""
    WriteParameters(filename,basis,mode)
    WriteParameters(file,basis)

write all the parameters to a file
"""
function WriteParameters(filename::String,
                         basis::AbstractAstroBasis,
                         mode::String="r+")

    h5open(filename, mode) do file
        WriteParameters(file,basis)
    end
end

function WriteParameters(file::HDF5.File,
                         basis::AB) where {AB <: AbstractAstroBasis}

    group = create_group(file,"BasisParameters")
    for i = 1:fieldcount(AB)
        # If this field is a string or a number, it is a parameter
        # (Prevent from dumping arrays)
        if fieldtype(AB,i) <: Union{String,Number}
            varname = string(fieldname(AB,i))
            write(group,varname,getfield(basis,i))
        end
    end
end

"""
    GetParameters(basis)

Returns all the basis parameters in a dictionary. The dictionary includes only the fields
that are of type `String` or `Number`, ignoring arrays and other types.

# Arguments
- `basis::AB`: An instance of a type that is a subtype of `AbstractAstroBasis`.

# Returns
- A dictionary where the keys are the field names as strings and the values are the corresponding field values
  if they are of type `String` or `Number`.

# Example
```julia
# Example basis type
struct ExampleBasis <: SphericalBasis
    name::String
    lmax::Int64
    nradial::Int64
    G::Float64
    rb::Float64
    tabPrefU::Array{Float64,2}
    tabPrefD::Array{Float64,2}
    tabUl::Array{Float64,1}
    tabDl::Array{Float64,1}
end

# Creating an instance of ExampleBasis
basis_instance = ExampleBasis("Example", 10, 5, 1.0, 1.0, rand(10, 5), rand(10, 5), rand(10), rand(10))

# Getting parameters
params = GetParameters(basis_instance)
println(params)
"""
function getparameters(basis::AB) where {AB <: AbstractAstroBasis}

    # Initialize an empty dictionary to store the parameters
    paramsdict = Dict{String,Union{String,Number}}()

    # Iterate through each field in the basis type
    for i = 1:fieldcount(AB)
        # Check if the field type is a String or Number
        if fieldtype(AB,i) <: Union{String,Number}
            # Convert the field name to a string
            varname = string(fieldname(AB,i))
            # Get the field value from the basis instance and store it in the dictionary
            paramsdict[varname] = getfield(basis,i)
        end
    end

    # Add the dimensionality parameter
    paramsdict["dimension"] = dimension(basis)

    # Return the dictionary of parameters
    return paramsdict
end
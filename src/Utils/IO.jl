
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

gives all the basis parameters in a dictionnary
"""
function GetParameters(basis::AB) where {AB <: AbstractAstroBasis}

    paramsdict = Dict{String,Union{String,Number}}()
    for i = 1:fieldcount(AB)
        # If this field is a string or a number, it is a parameter
        # (Prevent from dumping arrays)
        if fieldtype(AB,i) <: Union{String,Number}
            varname = string(fieldname(AB,i))
            paramsdict[varname] = getfield(basis,i)
        end
    end
    return paramsdict
end
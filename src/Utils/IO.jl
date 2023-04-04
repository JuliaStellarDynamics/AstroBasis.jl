
using HDF5

"""
    WriteParameters(filename,params,mode)
    WriteParameters(file,params)

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
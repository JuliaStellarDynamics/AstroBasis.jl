

import AstroBasis

UF,DF = AstroBasis.read_and_fill_prefactors(6,20)

println(AstroBasis.UlnpCB73(1,2,0.1,UF))

tabUlnp = zeros(20)

AstroBasis.tabUlnpCB73!(0,0.1,tabUlnp,20,UF)

# run some timing tests
function eval_time(l::Int64,nr::Int64,
                   tabUlnp::Array{Float64,1},
                   nradial::Int64,
                   tabPrefCB73_Ulnp::Matrix{Float64},
                   rb::Float64=1.)
    for x=1:nr
        r = x*0.001
        AstroBasis.tabUlnpCB73!(0,r,tabUlnp,20,UF)
    end
end

@time eval_time(0,1000,tabUlnp,20,UF)

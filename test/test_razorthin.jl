
G, rb = 1., 1.
ltest, nradial = 2, 5


@testset "razorthinbases" begin
    @testset "CB72" begin
        basis = CB72Basis(lmax=ltest,nradial=nradial,G=G, rb=rb)
        @test dimension(basis) == 2
        @test getparameters(basis)["name"] == "CB72"
        @test getDln(basis,ltest,1,rb) == 0.0
        @test getDln(basis,ltest,nradial-1,rb) ≈ 0.33126698841066865 atol=1e-6
        @test getUln(basis,ltest,nradial-1,rb) ≈ -0.3202172114362374 atol=1e-6
        tabUl!(basis,0,rb)
        @test basis.tabUl[1] == -1.0
        tabDl!(basis,0,rb)
        @test basis.tabDl[2] == 0
    end
    @testset "Kalnajs76" begin
        basis = K76Basis(lmax=ltest,nradial=nradial,G=G, rb=rb)
        @test dimension(basis) == 2
        @test getparameters(basis)["name"] == "K76"
        @test getDln(basis,ltest,1,rb) == 0.0
        @test getDln(basis,ltest,nradial-1,rb) == 0.0
        @test getUln(basis,ltest,nradial-1,rb) ≈ -0.3167679231608087 atol=1e-6
        tabUl!(basis,0,rb)
        @test basis.tabUl[1] == -0.8580855308097834
        tabDl!(basis,0,rb)
        @test basis.tabDl[1] == 0
    end
end
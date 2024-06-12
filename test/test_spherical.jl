

G, rb = 1., 1.
ltest, nradial = 2, 5


@testset "sphericalbases" begin
    @testset "CB73" begin
        # build bases with unique cases
        basis = AstroBasis.CB73Basis(lmax=ltest,nradial=1,G=G, rb=rb)
        basis = AstroBasis.CB73Basis(lmax=ltest,nradial=2,G=G, rb=rb)
        # build a standard basis
        basis = AstroBasis.CB73Basis(lmax=ltest,nradial=nradial,G=G, rb=rb)
        @test dimension(basis) == 3
        @test getparameters(basis)["name"] == "CB73"
        # backward compatibility check
        @test AstroBasis.GetParameters(basis)["name"] == "CB73"
        @test getDln(basis,ltest,nradial-1,rb) ≈ 1.6230683210206467 atol=1e-6
        @test getUln(basis,ltest,nradial-1,rb) ≈ -0.418381088294792 atol=1e-6
        tabUl!(basis,0,rb)
        @test basis.tabUl[2] == 0
        tabDl!(basis,0,rb)
        @test basis.tabDl[2] == 0
        # test the basis writing (for coverage!)
        AstroBasis.WriteParameters("tmp.h5",basis,"w")
    end
    @testset "Hernquist" begin
        basis = AstroBasis.HernquistBasis(lmax=ltest,nradial=nradial,G=G, rb=rb)
        @test dimension(basis) == 3
        @test getparameters(basis)["name"] == "Hernquist"
        @test getDln(basis,ltest,nradial-1,rb) ≈ 1.552003244163087 atol=1e-6
        @test getUln(basis,ltest,nradial-1,rb) ≈ -0.8668021315929388 atol=1e-6
        tabUl!(basis,0,rb)
        @test basis.tabUl[2] == 0
        tabDl!(basis,0,rb)
        @test basis.tabDl[2] == 0
    end
end
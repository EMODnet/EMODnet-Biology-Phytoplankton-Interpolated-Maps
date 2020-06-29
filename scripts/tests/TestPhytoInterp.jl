using Test
using Dates
include("../PhytoInterp.jl")

@testset "Get names from file name" begin
    datafile1 = "Pseudo-nitzschia delicatissima-1995-2020.csv"
    speciesname = get_species_name(datafile1)
    speciesslug = get_species_slug(datafile1)
    @test speciesname == "Pseudo-nitzschia delicatissima"
    @test speciesslug == "Pseudo-nitzschia_delicatissima"

    datafile2 = "Plagiogramma brockmanni var. brockmanni-1995-2020.csv"
    speciesname = get_species_slug(datafile2)
    speciesslug = get_species_slug(datafile2)
    @test speciesname == "Plagiogramma_brockmanni_var_brockmanni"
    @test speciesslug == "Plagiogramma_brockmanni_var_brockmanni"
end

@testset "Reading new dataset" begin
    # New dataset (Deltares, greater North Sea)
    testfile1 = "../../data/test/Polykrikos-1995-2020.csv"
    dates, lons, lats, occurs = read_data_phyto(testfile1);
    @test lons[1] == 488813.927651439
    @test lats[2] == 5.8095984857214e6
    @test occurs[10] == 0
    @test dates[1] == Date(1988, 4, 30)
    @test sum(occurs) == 2300
end

@testset "Read old old dataset" begin
    testfile2 = "../../data/test/Cerataulina bergonii-1995-2020.csv";
    dates, lons, lats, occurs = read_data_phyto(testfile2);
    @test lons[11] == 848766.458997052
    @test lats[22] == 6.11699330386247e6
    @test occurs[end] == 0
    @test dates[end] == Date(1999, 12, 14)
    @test sum(occurs) == 2998
end

@testset "Coordinate transformation" begin
    lonp, latp = transform_coords([477850.085718031, 560048.806065514],
                                  [5756212.06783869, 5746827.70084791])

    @test lonp == [2.677667000000006, 3.872167000000008]
    @test latp == [51.956167000000036, 51.86900000000005]
end

@testset "Interpolation" begin
    longrid = 0:1.:100.;
    latgrid = 0:1.:50.;
    xx = [i for i in longrid, j in latgrid]
    yy = [j for i in longrid, j in latgrid];
    function f(x, y)
        return cos.(.3 .* x) .+ sin.(0.2 .* y);
    end
    field = f(xx, yy);

    obslon = [5.3, 4.2]
    obslat = [13.2, 43.1];
    finterp = reinterp_field(longrid, latgrid, field, obslon, obslat);
    @test finterp[1] == 0.4607541402478184
    @test finterp[2] == 1.0234827512907576
    finterp2 = reinterp_field(longrid, latgrid, field, [0.], [0.])
    @test finterp2 == [1.0]
end

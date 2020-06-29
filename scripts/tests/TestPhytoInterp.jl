using Test
include("./PhytoInterp.jl")


speciesname = get_species_name("Archaeperidinium-1995-2020.csv")
@test speciesname == "Archaeperidinium"


testfile1 = "/data/EMODnet/Biology/phytoplankton/GreaterNorthSea/Polykrikos-1995-2020.csv"
dates, lons, lats, occurs = read_data_phyto(testfile1);
@test lons[1] == 488813.927651439
@test lats[2] == 5.8095984857214e6
@test occurs[10] == 0
@test dates[1] == Date(1988, 4, 30)
@test sum(occurs) == 2

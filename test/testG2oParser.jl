using Test

@testset "Test g2o file parsing." begin
    instructions = importG2o("octagon.g2o")
    # Choose random fields to test.
    @test instructions[1][1] == "EDGE_SE2"
    @test instructions[7][1] == "EDGE_SE2"
    @test instructions[6][12] == "6541.252776"
    @test instructions[3][7] == "1211.201664"
    # Input file has 8 lines, with 12 columns.
    @test length(instructions) == 8
    @test length(instructions[2]) == 12
end

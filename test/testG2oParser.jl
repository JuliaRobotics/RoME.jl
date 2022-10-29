using Test
using RoME

@testset "Test g2o file parsing." begin
  
  
  file = @__DIR__()*"/octagon.g2o"
  instructions = importG2o(file)
  # instructions = importG2o("octagon.g2o")
  # Choose random fields to test.
  @test instructions[1][1] == "EDGE_SE2"
  @test instructions[7][1] == "EDGE_SE2"
  @test instructions[6][12] == "6541.252776"
  @test instructions[3][7] == "1211.201664"
  # Input file has 8 lines, with 12 columns.
  @test length(instructions) == 8
  @test length(instructions[2]) == 12

  # build the factor graph from these instructions


  # generate g2o from factor graph and should match orginal file
end


@testset "Test g2o file parsing." begin

  reflines = [
              "EDGE_SE2 0 1 10.0 0.0 1.0471975511965976 10.0 0.0 0.0 10.0 0.0 10.0";
              "LANDMARK 0 2 0.0 20.0 3.162277660168379 0.0 1.0";
              "EDGE_SE2 1 3 10.0 0.0 1.0471975511965976 10.0 0.0 0.0 10.0 0.0 10.0";
              "EDGE_SE2 3 4 10.0 0.0 1.0471975511965976 10.0 0.0 0.0 10.0 0.0 10.0";
              "EDGE_SE2 4 5 10.0 0.0 1.0471975511965976 10.0 0.0 0.0 10.0 0.0 10.0";
              "EDGE_SE2 5 6 10.0 0.0 1.0471975511965976 10.0 0.0 0.0 10.0 0.0 10.0";
              "EDGE_SE2 6 7 10.0 0.0 1.0471975511965976 10.0 0.0 0.0 10.0 0.0 10.0";
              "LANDMARK 7 2 0.0 20.0 3.162277660168379 0.0 1.0"]


  # test file generation
  fg = generateGraph_Hexagonal(graphinit=false)

  filepath = exportG2o(fg)

  fid = open(filepath, "r")
  count = 0
  for line in eachline(fid)
    count += 1
    @test reflines[count] == line
  end
  close(fid)

end

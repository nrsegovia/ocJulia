using CSV
using DataFrames
using Plots
include("lightCurve.jl")
include("fourierFit.jl")
include("ocJulia.jl")

# gaston()
gr()

loc = Base.@__DIR__
testFile = joinpath(loc, "testFiles","OGLE-BLG-RRLYR-13527.dat")#"DataFirstRun", "OGLE-BLG-RRLYR-09998.dat")
plotsDir = joinpath(loc, "Plots")
df = DataFrame(CSV.File(testFile))

testTime = df.HJD
testSignal = df.MAG
testPeriod = 0.43389532

first = LightCurve(testTime, testSignal, testPeriod)

bins = createBins(first, 80)
saveDiagnoseHist(first, bins, joinpath(plotsDir, "histogramTest.png"))
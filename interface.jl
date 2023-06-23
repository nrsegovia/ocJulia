using CSV
using DataFrames
using Plots
include("lightCurve.jl")
include("fourierFit.jl")
include("ocJulia.jl")

# gaston()
gr()

loc = Base.@__DIR__
testFile = joinpath(loc, "testFiles","OGLE-BLG-RRLYR-09998.dat")#"13527", "09998"
plotsDir = joinpath(loc, "Plots")
df = DataFrame(CSV.File(testFile))

testTime = df.HJD
testSignal = df.MAG
testPeriod = 0.3603214 #0.43389532

maskForTemplate = 5290 .< testTime .< 5455
first = LightCurve(testTime, testSignal, testPeriod, median(testTime[maskForTemplate]))
curveForTemplate = LightCurve(testTime[maskForTemplate], testSignal[maskForTemplate], testPeriod, median(testTime[maskForTemplate]))
computePhase(curveForTemplate)

bins = createBins(first, 30.)
# saveDiagnoseHist(first, bins, joinpath(plotsDir, "histogramTest09998.png"))
myFourier = fourierFit(curveForTemplate)
myOC = ocBase(first, myFourier, bins)
allResults = runOC(myOC, 100)

scatter(allResults[1], allResults[2])
gui()
readline()

# O-C computation is working... TODO: make template creation automatic somehow. Implement bootstrap... make it depend on the number of points per bin? 
# And then it should be done! for future: add lightCurve function to "join" LCs from different surveys.
# Not relevant right now as I will use OGLE only...

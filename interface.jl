using CSV
using DataFrames
using Plots
include("lightCurve.jl")
include("fourierFit.jl")

gaston()
loc = Base.Filesystem.pwd()
testFile = joinpath(loc, "OGLE-BLG-RRLYR-13527.dat")#"DataFirstRun", "OGLE-BLG-RRLYR-09998.dat")
df = DataFrame(CSV.File(testFile))

testTime = df.HJD
testSignal = df.MAG
testPeriod = 0.43389532

first = LightCurve(testTime, testSignal, testPeriod)

computePhase(first)


# fourierTest = fourierFit(first, 15, true)
fourierTest = fourierFit(first)
# print(fourierTest.order)
print(fourierTest.coefficients)

x1 = 0:0.01:1
y1 = evalFour(x1, fourierTest)

plot(x1, y1)
plotCurve(first, true)

yflip!(true)

# display(toPlot)
# display(toPlot2)
gui()
readline()
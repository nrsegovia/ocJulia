using CSV
using DataFrames
using Plots
using Printf
include("lightCurve.jl")
include("fourierFit.jl")
include("ocJulia.jl")

# gaston()
gr()

loc = Base.@__DIR__
source = DataFrame(CSV.File(joinpath(loc, "testFiles", "RRcOddList.csv")))
plotsDir = joinpath(loc, "Plots")
outputDataDir = joinpath(loc, "output")

for x in 1:size(source, 1)
    name = source[x, "ID"]
    df = DataFrame(CSV.File(joinpath(loc, "testFiles",name*".dat")))
    time = df.HJD
    signal = df.MAG
    period = source[x, "PeriodOGLE"]
    periodString = @sprintf("%.4f",period)
    maskForTemplate = source[x, "TemplateStart"] .< time .< source[x, "TemplateEnd"]
    lc = LightCurve(time, signal, period, median(time[maskForTemplate]))
    curveForTemplate = LightCurve(time[maskForTemplate], signal[maskForTemplate], period, median(time[maskForTemplate]))
    computePhase(curveForTemplate)

    bins = createBins(lc, 30.) # using 30 days as default...

    # saveDiagnoseHist(lc, bins, joinpath(plotsDir, name+"_histogram.png")) # In case you want to see some details related to the binning process
    
    myFourier = fourierFit(curveForTemplate)
    myOC = ocBase(lc, myFourier, bins)
    allResults = runOC(myOC, 100) # At the moment the bootstrap does nothing

    xInterp = allResults.representativeTimes[1]:1.:allResults.representativeTimes[end]
    yInterp = @. allResults.interpolate(xInterp)

    lcCorrection = @. allResults.interpolate(time)
    correctedTime = @. time - lcCorrection
    newCurve = LightCurve(correctedTime, signal, period)
    computePhase(newCurve)

    scatter(allResults.representativeTimes, allResults.computedOC, label = "O-C")
    plot!(xInterp, yInterp, label = "Interpolation")
    xlabel!("HJD-2,450,000")
    ylabel!("O-C [d]")
    title!(name*", P: "*periodString)
    png(joinpath(plotsDir, name*"_oc.png"))

    scatter(newCurve.phase.phase, newCurve.properties.signal)
    xlabel!("Arbitrary Phase")
    ylabel!("I-Mag")
    yflip!()
    title!(name*", P: "*periodString)
    png(joinpath(plotsDir, name*"_corrected.png"))

    exportOCasCSV(allResults, joinpath(outputDataDir, name*"_oc.dat"))

    df[!, "HJD_New"] = correctedTime
    CSV.write(joinpath(outputDataDir, name*"_new.dat"), df)
end
# add function to store values of corrected LC and O-C values

# gui()
# readline()
# O-C computation is working... TODO: make template creation automatic somehow. Implement bootstrap... make it depend on the number of points per bin? 
# And then it should be done! for future: add lightCurve function to "join" LCs from different surveys.
# Not relevant right now as I will use OGLE only...

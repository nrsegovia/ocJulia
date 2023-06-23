using CSV
using DataFrames
using Plots
using Optim
include("lightCurve.jl")
include("fourierFit.jl")

function phaseCoverage(time, period)
    dataPhased = @. time/period - time÷period
    total = 0
    extent = 0
    leftEdge = 0
    for i in 0.05:0.05:1
        extent = extent + 1
        rightEdge = i
        total = any(leftEdge .< dataPhased .<= rightEdge) ? total + 1 : total
        leftEdge = i
    end
    return total/extent
end

function createBins(data::LightCurve, numberCycles::Int)
    binSize = numberCycles * data.period.value
    binEdges = (data.properties.time[1] - 1):binSize:(data.properties.time[end]+binSize/2)
    return binEdges
end

function saveDiagnoseHist(data::LightCurve, binEdges::AbstractVector{Float64}, savePath::String)
    xP = []
    yP = []
    for p in 1:(length(binEdges) - 1)
        rightEdge = binEdges[p+1]
        leftEdge = binEdges[p]
        push!(xP, (rightEdge + leftEdge)/2.0)
        currentTimes = data.properties.time[leftEdge .< data.properties.time .< rightEdge]
        push!(yP, phaseCoverage(currentTimes, data.period.value))
    end
    stephist(first.properties.time, bins = binEdges, ylabel = "Counts", xlabel = "HJD - 2,450,000")
    scatter!(twinx(), xP, yP, color=:red, xticks=:none, ylabel = "Phase Coverage", label = "Fraction", markersize = 2)
    png(savePath)
end

struct ocBase
    allData::LightCurve
    template::fourierFit
    binEdges::AbstractVector{Float64}
    ocBase(curve::LightCurve, fourier::fourierFit, edges::AbstractVector{Float64}) = new(curve, fourier, edges)
end

function runOC(baseInfo::ocBase)
    # Add bootstrap...
    allEdges = baseInfo.binEdges
    totalEdges = length(allEdges)
    allTimes = baseInfo.allData.properties.time
    period = baseInfo.allData.period.value
    representativeTimeValues = []
    computedOC = []
    uncertaintyOC = []
    for k in 1:(totalEdges[end] -1)
        timeSubSet = allTimes[allEdges[k] .< allTimes .<= allEdges[k + 1]]
        coverage = phaseCoverage(timeSubSet, period)
        if coverage >= 0.5
            push!(representativeTimeValues, median(timeSubSet))
            # Add minimization here...
        end
    end
    return(representativeTimeValues, computedOC, uncertaintyOC)
end
#Use period as scale, to take into account short periods... then go back to days/minutes/whatever


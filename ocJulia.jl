using CSV
using DataFrames
using Plots
using LeastSquaresOptim
include("lightCurve.jl")
include("fourierFit.jl")

function phaseCoverage(time, period)
    dataPhased = @. time/period - time÷period
    dataPhased[dataPhased .< 0] .+= 1
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

function createBins(data::LightCurve, sizeInDays::Float64)
    binEdges = collect((data.properties.time[1] - 1):sizeInDays:(data.properties.time[end]+sizeInDays/2))
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

function binEdgesToPairs(binEdges::AbstractVector{Float64}, periodAtTime::Float64)
    nPairs = length(binEdges) - 1
    pairs = []
    periodIdx = 0
    if periodAtTime > binEdges[end] || periodAtTime < binEdges[1]
        error("Period does not correspond to time baseline")
    end
    check = true
    for k in 1:nPairs
        left = binEdges[k]
        right = binEdges[k+1]
        if check && (left < periodAtTime < right)
            periodIdx = k
            check = false
        end
        push!(pairs, [left, right])
    end
    return pairs, periodIdx
end


function runOC(baseInfo::ocBase, numberBootstrap::Int)
    # Add bootstrap!
    allEdges = baseInfo.binEdges
    pairs, periodIdx = binEdgesToPairs(allEdges, baseInfo.allData.period.atTime)
    totPairs = length(pairs)
    orderedIndxs = periodIdx - totPairs == 0 ? collect(totPairs:-1:1) : vcat(collect(periodIdx:1:totPairs), collect((periodIdx - 1):-1:1))
    allTimes = baseInfo.allData.properties.time
    allMags = baseInfo.allData.properties.signal
    period = baseInfo.allData.period.value
    representativeTimeValues = []
    computedOC = []
    uncertaintyOC = []
    idxPrev = orderedIndxs[1]
    ocFirstResult = 0.0
    ocPrev = 0.0 # Works as initial guess
    firstPass = true
    for k in orderedIndxs
        if abs(k - idxPrev) > 1
            ocPrev = ocFirstResult
        end
        idxPrev = k
        edgeMask = pairs[k][1] .< allTimes .<= pairs[k][2]
        timeSubSet = allTimes[edgeMask]
        nPoints = length(timeSubSet)
        # maybe set 10 as variable instead of fixed value
        if nPoints >= 10
            magSubSet = allMags[edgeMask]
            push!(representativeTimeValues, median(timeSubSet))
            toBeMinimized(oc) = residualsFour((timeSubSet .- oc[1])/period, magSubSet, baseInfo.template)
            minimizing = optimize(toBeMinimized, [ocPrev], LevenbergMarquardt())
            result = minimizing.minimizer[1]
            if firstPass
                ocFirstResult = result
                firstPass = false
            end
            push!(computedOC, result)
            ocPrev = result
        end
    end
    return(representativeTimeValues, computedOC, uncertaintyOC)
end

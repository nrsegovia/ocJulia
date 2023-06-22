using CSV
using DataFrames
include("lightCurve.jl")
include("fourierFit.jl")

function checkPhaseCoverage()
end

function createBins()
end

struct ocBase
    template::fourierFit
    allData::LightCurve
    binEdges::AbstractVector{Float64}

end

function runOC(baseInfo::ocBase)
end
#Use period as scale, to take into account short periods... then go back to days/minutes/whatever
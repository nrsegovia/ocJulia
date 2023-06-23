using CSV
using DataFrames
using Plots

function toPhase(time::AbstractVector{Float64}, period::Float64)
    preProcess = time / period
    return (preProcess - floor.(preProcess))
end

struct Period
    value::Float64
    atTime::Float64
    Period(value,time) = value > 0.0 ? new(value, time) : error("Period must be positive")
    Period(value) = value > 0.0 ? new(value, 0.0) : error("Period must be positive")
    Period() = new(0.0, 0.0)
end

struct LightCurveProperties
    dataPoints::Int64
    time::AbstractVector{Float64}
    signal::AbstractVector{Float64}
    signalError::AbstractVector{Float64}

    LightCurveProperties(dataPoints::Int64, time::AbstractVector{Float64}, signal::AbstractVector{Float64}) = (length(time) == dataPoints & length(signal) == dataPoints) ? 
    new(dataPoints, time, signal, zeros(dataPoints)) : error("Length of vectors does not match")

    LightCurveProperties(dataPoints::Int64, time::AbstractVector{Float64}, signal::AbstractVector{Float64}, signalError::AbstractVector{Float64}) = 
    (length(time) == dataPoints & length(signal) == dataPoints & length(signalError) == dataPoints) ? new(dataPoints, time, signal, signalError) : 
    error("Length of vectors does not match")
end

struct LightCurvePhase
    dataPoints::Int64
    phase::AbstractVector{Float64}
    done::Bool

    LightCurvePhase(dataPoints::Int64, phase::AbstractVector{Float64}) = length(phase) == dataPoints ? new(dataPoints, phase, true) : error("Length of phase vector does not match")
    LightCurvePhase(dataPoints::Int64) = new(dataPoints, zeros(dataPoints), false)

end

mutable struct LightCurve
    properties::LightCurveProperties
    period::Period
    phase::LightCurvePhase
    LightCurve(time::AbstractVector{Float64}, signal::AbstractVector{Float64}) = new(LightCurveProperties(length(time), time, signal), Period(), LightCurvePhase(length(time)))
    LightCurve(time::AbstractVector{Float64}, signal::AbstractVector{Float64}, signalError::AbstractVector{Float64}) = new(LightCurveProperties(length(time), time, signal, signalError), Period(), LightCurvePhase(length(time)))
    LightCurve(time::AbstractVector{Float64}, signal::AbstractVector{Float64}, period::Float64) = new(LightCurveProperties(length(time), time, signal), Period(period), LightCurvePhase(length(time)))
    LightCurve(time::AbstractVector{Float64}, signal::AbstractVector{Float64}, signalError::AbstractVector{Float64}, period::Float64) = new(LightCurveProperties(length(time), time, signal, signalError), Period(period), LightCurvePhase(length(time)))
    LightCurve(time::AbstractVector{Float64}, signal::AbstractVector{Float64}, period::Float64, periodAtTime::Float64) = new(LightCurveProperties(length(time), time, signal), Period(period, periodAtTime), LightCurvePhase(length(time)))
    LightCurve(time::AbstractVector{Float64}, signal::AbstractVector{Float64}, signalError::AbstractVector{Float64}, period::Float64, periodAtTime::Float64) = new(LightCurveProperties(length(time), time, signal, signalError), Period(period, periodAtTime), LightCurvePhase(length(time)))
end


function setPeriod(curve::LightCurve, period::Float64)
     curve.period = Period(period)
end

function computePhase(curve::LightCurve)
    curve.phase = LightCurvePhase(curve.properties.dataPoints, toPhase(curve.properties.time, curve.period.value))
end

function computePhase(curve::LightCurve, period::Float64)
    setPeriod(curve, period)
    curve.phase = LightCurvePhase(curve.properties.dataPoints, toPhase(curve.properties.time, curve.period.value))
end

function plotCurve(curve::LightCurve, phased::Bool)
    if phased
        if curve.phase.done
            plotted = scatter!(curve.phase.phase, curve.properties.signal, label="Folded LC", ms=2, ma=0.5)
        else
            error("Phase not set.")
        end
    else
        plotted = scatter!(curve.properties.time, curve.properties.signal, label="LC", ms=2, ma=0.5)
    end
    return plotted
end

function plotCurve(curve::LightCurve)
    plotted = scatter!(curve.properties.time, curve.properties.signal, label="LC", ms=2, ma=0.5)
    return plotted
end

using CSV
using DataFrames
using GLM
using Lasso
# using Optim
include("lightCurve.jl")


function createExprAndData(lightCurve::LightCurve, N::Int)
    if lightCurve.phase.done
        df = DataFrame(Y = lightCurve.properties.signal)
        strExp = ""
        for i in 1:N
            # print(i)
            A = string("A", i)
            B = string("B", i) 
            df[!,A] = sin.(2 * pi * i * lightCurve.phase.phase)
            df[!,B] = cos.(2 * pi * i * lightCurve.phase.phase)
            if i!=N
                strExp *= string(A, " + ", B, " + ")
            else
                strExp *= string(A, " + ", B)
            end
        end
        return(strExp, df)
    else
        error("Light curve phase has not been computed")
    end
end

function getFourierCoefficientsLasso((strExp, df)::Tuple{String, DataFrame})
    m = fit(LassoModel, term("Y") ~ sum(term.(split(strExp, " + "))), df)
    C = coef(m)
    return(C)
end

function getFourierCoefficients((strExp, df)::Tuple{String, DataFrame})
    ols = lm(term("Y") ~ sum(term.(split(strExp, " + "))), df)
    C = coef(ols.model)
    return(C)
end

# X is phase
function evalFour(X,C,N)
    tot = length(X)
    out = Vector{Float64}(undef,tot)
    out = fill(C[1], tot)
    for i in 1:N
        out += (C[2*i] * sin.(2 * pi * i * X) + C[2*i+1] * cos.(2 * pi * i * X))
    end
    return out
end

function getFourierCoefficients(lightCurve::LightCurve)
    residuals = -99
    outCoeffs = [0]
    N = 0

    for i in 3:15 # defaults to this range... maybe add option to edit it?
        strExp, df = createExprAndData(lightCurve, i)
        ols = lm(term("Y") ~ sum(term.(split(strExp, " + "))), df)
        C = coef(ols.model)
        currentResiduals = sum(abs.(lightCurve.properties.signal - evalFour(lightCurve.phase.phase, C, i)))
        if (residuals < 0) || (residuals > currentResiduals)
            residuals = currentResiduals
            outCoeffs = C
            N = i
        end
    end
    return(N, outCoeffs)
end

mutable struct fourierFit
    order::Int
    coefficients::AbstractVector{Float64}

    fourierFit(lightCurve::LightCurve, N::Int, lass::Bool) = N > 0 ? (lass ? new(N, getFourierCoefficientsLasso(createExprAndData(lightCurve, N))) : new(N, getFourierCoefficients(createExprAndData(lightCurve, N)))) : error("N has to be positive")
    fourierFit(lightCurve::LightCurve) = new(getFourierCoefficients(lightCurve)...)
end

# X is phase
function evalFour(X,F::fourierFit)
    tot = length(X)
    out = Vector{Float64}(undef,tot)
    out = fill(F.coefficients[1], tot)
    for i in 1:F.order
        out += (F.coefficients[2*i] * sin.(2 * pi * i * X) + F.coefficients[2*i+1] * cos.(2 * pi * i * X))
    end
    return out
end

function residualsFour(times::AbstractVector{Float64}, signal::AbstractVector{Float64}, F::fourierFit)
    return sum(abs.(signal - evalFour(times, F)))
end
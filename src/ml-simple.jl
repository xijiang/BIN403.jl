#!/usr/bin/env julia
using Statistics, Distributions, DataFrames, StatsPlots, LaTeXStrings
y = [11, 7, 12]                 # the observations
m = mean(y)
s = std(y)                      # above two functions are from `Statistics`
df = DataFrame(μ = Float64[], σ = Float64[], pd = Float64[])
for μ in m-1:0.01:m+1
    for σ in s-1:.01:s+1
        pd = prod(pdf.(Normal(μ, σ), y)) # `Normal is from Distributions
        push!(df, [μ σ pd])              # note: no `,`, to make it a row vector
    end
end
@df df plot(:μ, :σ, :pd, st=:surface,
            xlabel = L"\mu",    # L"..." is for LaTeX strings
            ylabel = L"\sigma",
            zlable = "Probability density",
            legend = false)
last(sort(df, :pd))             # sort on pd, the last one is where peak is in our grid

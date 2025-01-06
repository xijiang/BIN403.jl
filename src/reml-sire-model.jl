using LinearAlgebra, Distributions, Statistics, Random, Plots
"""
    bhs_simu(ns, nd, nl, σₛ², σₑ²)
---
To simulate a balanced sire model design.
- `ns`: number of sires
- `nd`: number of daughters per sire
- `nl`: number of levels of a fixed effect
- `σₛ²`: sire variance, or `σₐ²/4`.
- `σₑ²`: residual vairance

The function samples the levels of the fixed effects from `N(0, 1)`.
`u` and `e` from `N(0, σₛ²)`, and `N(0, σₑ²)`, respectively.

The returns are:
- `y`: the phenotype observations
- `X`: design matrix of the fixed effects
- `Z`: design matrix of the sires
"""
function bhs_simu(ns, nd, nl, σₛ², σₑ²)
    # Simulation
    nt = ns * nd                       # total number of records
    b = rand(Normal(), nl)             # herd effects ~N(0, 1)
    u = rand(Normal(0, √(σₛ²)), ns)    # sire effects ~N(0, σₛ²)
    e = rand(Normal(0, √(σₑ²)), nt)    # residues ~N(0, σₑ²)
    X = Int.(kron(ones(nt ÷ 2), I(2))) # design matrix for herds
    Z = Int.(kron(I(ns), ones(nd)))    # design matrix for sires
    y = X * b + Z * u + e              # observed data
    @info var(u)
    y, X, Z
end

"""
    function em_reml(y, X, Z, σₛ², σₑ²)
---
One can repeat this function until `|σₛ²hat - σₛ²| < ϵ`.
No A inverse is included here as we assume sires are unrelated.
"""
function em_reml(y, X, Z, σₛ², σₑ²)
    nt, ns = length(y), size(Z)[2]
    iXX = inv(X'X)
    S = I(nt) - X * iXX * X'
    λ = σₑ²/σₛ²
    MMM = Z'S * Z + I(ns) * λ
    iMM = inv(MMM)
    uhat = iMM * Z'S*y
    σₛ²hat = (uhat'uhat + tr(iMM) * σₑ²) / ns
    bhat = iXX * (X'y - X'Z*uhat)
    ehat = y - X * bhat - Z * uhat
    σₑ²hat = (ehat'ehat + tr(iMM * Z'Z) * σₑ²) / nt
    σₛ²hat, σₑ²hat
end


"""
    function ai_reml(y, X, Z, σₛ², σₑ²)
---
Again use true values of variance components as starting values (note: iterative algorithms need starting valuesand good starting values help convergence (especially for AI-REML)).
- slide 8 of `estimation_of_variance_components.ppt` shows calculation of Score functions
- slide 9 shows calculation of second derivatives
- again: no A-inverse is needed here because sires are assumed unrelated
- slide 7 in `Maximum_likelihood_estimation.ppt` gives formula for update in 2nd order algorithms (averageinfo)
- note: `P` is calculated directly below, but it is more efficient to calculate terms like `d'Pd` for any `d` (but treating `d` as records), as `d'Pd = (d-Xh)'inv(V)(d-Xh) = d'inv(R)̂e`, where `h` = estimate of fixed effects IF `d` were actual records, and `̂e` = estimates of residuals if `d` were records.  The latter is only requires estimation of BLUP solutions given `d` as records, not the inversion of `V`.
"""
function ai_reml(y, X, Z, σₛ², σₑ²)
    nt, ns = length(y), size(Z)[2]
    iXX = inv(X'X)
    S = I(nt) - X * iXX * X'
    V = Z * Z' .* σₛ² + I(nt) .* σₑ²
    iV = inv(V)
    P = iV - iV * X * inv(X'iV * X) * X'iV
    λ = σₑ²/σₛ²
    MMM = Z'S * Z + I(ns) * λ
    iMM = inv(MMM)
    uhat = iMM * Z'S * y
    bhat = iXX * (X'y - X'Z * uhat)
    ehat = y - X * bhat - Z * uhat
    Score = [-.5 * (ns/σₛ² - tr(iMM) * σₑ² / σₛ²^2 - uhat'uhat / σₛ²^2),
             -.5 * ((nt -ns) / σₑ² - tr(iMM) * λ / σₑ² - ehat'ehat/σₑ²^2)]
    AI = [uhat'Z'P * Z * uhat/σₛ²^2  ehat'P * Z * uhat / (σₑ² * σₛ²)
          0                          ehat'P * ehat / σₑ²^2 ]
    AI[2, 1] = AI[1, 2]
    [σₛ², σₑ²] - (AI \ Score)
end

function test_remls(n)
    y, X, Z = bhs_simu(10, 10, 2, .1, .9)
    a, b, c, d = .1, .9, .1, .9
    eσ² = zeros(n, 2)
    for i in 1:n
        a, b = em_reml(y, X, Z, a, b)
        c, d = ai_reml(y, X, Z, c, d)
        eσ²[i, :] = [a c]
    end
    plot(1:n, eσ²[:, 1], label = "EM-REML");
    plot!(1:n, eσ²[:, 2], label = "AI_REML");
    savefig("~/Music/tmp/reml.pdf")
end

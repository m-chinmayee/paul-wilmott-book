# Pricing European and American Options via binomial tree
# Author: Chinmayee Mishra 
# Last updated: Jan 10, 2023

using LinearAlgebra

function payoff(A, K, option = "call")
    if option == "put"
        return max(K - A, 0)
    else
        return max(A - K, 0)
    end
end

function parameters(intrate, volatility, expiry, steps)
    # Definition from PW chap3 appendix
    # dt = expiry / steps
    # discount = exp(-intrate * dt)
    # fact1 = exp((intrate + volatility^2) * dt)
    # fact2 = (fact1 + exp(discount + fact1)) / 2.0
    # u = fact2 + sqrt(fact1 * fact2 - 1)
    # d = 1 / u
    # p = (exp(intrate * dt) - d) / (u - d)

    # Standard definitions for discrete time
    dt = expiry / steps
    u = 1 + volatility * sqrt(dt)
    d = 1 / u
    p = 0.5 + intrate * sqrt(dt) / (2.0 * volatility)
    discount = exp(-intrate * dt)

    return u, d, p, discount
end

function assetevolution(s_init, u, d, steps)
    # Initiate arrays and parameters
    s = [zeros(i) for i = 1:steps]

    # Main code block
    s[1][1] = s_init
    for n = 2:steps
        s[n][1] = d * s[n-1][1]
        for j = 2:n
            s[n][j] = (u / d) * s[n][j-1]
        end
    end

    return s
end

function europeanoptionevolution(s, strike, discount, steps, option = "call")
    # Initiate arrays and parameters
    v = [zeros(i) for i = 1:steps]

    # Main code block
    v[steps] = payoff.(s[steps], strike)
    for n = steps-1:-1:1
        for j = 1:n
            v[n][j] = (p * v[n+1][j+1] + (1-p) * v[n+1][j]) * discount
        end
    end

    return v
end

function americanoptionevolution(s, strike, discount, steps, option = "call")
    # Initiate arrays and parameters
    v = [zeros(i) for i = 1:steps]

    # Main code block
    v[steps] = payoff.(s[steps], strike)
    for n = steps-1:-1:1
        for j = 1:n
            v[n][j] = max((p * v[n+1][j+1] + (1-p) * v[n+1][j]) * discount, payoff(s[n][j], strike))
        end
    end

    return v
end

# Global inputs
s0 = 100
r = 0.1
σ = 0.2
K = 100
T = 1
Nt = 5


u, d, p, discount = parameters(r, σ, T, Nt)
s = assetevolution(s0, u, d, Nt)
veu = europeanoptionevolution(s, K, discount, Nt)
vam = americanoptionevolution(s, K, discount, Nt)

println(veu[1][1], " ", vam[1][1])
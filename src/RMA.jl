module RMA

export RMAParameters, vec_f, integrate_ode, critical_r0, e3

using DifferentialEquations
using StaticArrays
using LinearAlgebra

@kwdef struct RMAParameters{T} <: AbstractVector{T}
    c::T = 0.19
    mu::T = 0.03
    nu::T = 0.003
    alpha::T = 800.0
    beta::T = 1.5
    gamma::T = 0.004
    delta::T = 2.2
    r::T
end
Base.size(p::RMAParameters) = (8,)
function Base.getindex(p::RMAParameters, i::Int)
    fieldsym = fieldnames(RMAParameters)[i]
    getfield(p, fieldsym)
end
Base.IndexStyle(::RMAParameters) = IndexLinear()

function vec_f(X, Y; p::RMAParameters)
    (; c, mu, nu, alpha, beta, gamma, delta, r) = p
    [r * X * (1 - c / r * X) * (X - mu) / (X + nu) - alpha * X * Y / (X + beta),
        gamma * (alpha * X * Y) / (X + beta) - delta * Y]
end
function rma_step(u, p, t)
    (; c, mu, nu, alpha, beta, gamma, delta, r) = p
    du1 = r * u[1] * (1 - c / r * u[1]) * ((u[1] - mu) / (nu + u[1])) -
          (alpha * u[1] * u[2]) / (beta + u[1])
    du2 = gamma * (alpha * u[1] * u[2]) / (beta + u[1]) - delta * u[2]
    SA[du1, du2]
end
function integrate_ode(u0, tspan, p)
    ode = ODEProblem(rma_step, u0, tspan, p)
    sol = solve(ode)
end
function critical_r0(p::RMAParameters)
    (; c, mu, nu, alpha, beta, gamma, delta, r) = p
    c * delta * beta / (gamma * alpha - delta)
end
function e3(p::RMAParameters)
    (; c, mu, nu, alpha, beta, gamma, delta, r) = p
    X_star = delta * beta / (gamma * alpha - delta)
    Y_star = r * gamma / delta * X_star * (1 - c / r * X_star) * (X_star - mu) / (nu + X_star)
    SA[X_star, Y_star]
end

end

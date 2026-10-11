using Plots
using LinearAlgebra
using FastGaussQuadrature
using QuadGK
using LaTeXStrings
using ProgressMeter
using Roots

const address = if Sys.iswindows()
    "C:\\Users\\johnf\\Repositories\\Maris-Tandy\\thesis\\images\\"
elseif Sys.isapple()
    "/Users/johnreeg/Documents/Repositories/Maris-Tandy/thesis/images/"
else
    "/home/john-reeg/Documents/Maris-Tandy/thesis/images/"
end

const m = 0.0037
const Lambda2 = 1e6
const epsilon2 = 1e-4
const Lambda_PV = 2e2
const Lambda_QCD = 0.234
const gamma_m = 0.48 # 12/(33 - 2 * Nf{4}) = 0.48
const mu = 19.0
const C_f = 4/3

function QuarkDSE(w::Float64, D::Float64, PV::Bool; radial_steps::Int = 256, angular_steps::Int = 64)
    x, w_x = gausslegendre(radial_steps)
    z, w_z = gausschebyshevu(angular_steps)

    t = 0.5 * (log(Lambda2) - log(epsilon2)) * x .+ 0.5 * (log(Lambda2) + log(epsilon2))
    w_t = 0.5 * (log(Lambda2) - log(epsilon2)) * w_x

    A = 1.6 * ones(length(t))
    B = (m + 1e-3) * ones(length(t))
    Z_2 = 1.0
    Z_4m = 0.7 * m

    k2(p::ComplexF64, q::Float64, z::Float64) = p^2 + q^2 - 2*p*q*z

    alpha_UV(k2::ComplexF64) = 2pi * gamma_m * (1 - exp(-k2)) / (k2 * log(exp(2)-1 + (1 + k2/Lambda_QCD^2)^2))
    alpha_IR(k2::ComplexF64) = D/w^6 * pi * k2 * exp(-k2 / w^2)
    alpha(k2::ComplexF64) = alpha_IR(k2) + alpha_UV(k2)

    function z_intAA(p::ComplexF64, q::Float64)
        if PV == true
            integrand = @. w_z * (p*q*z + 2*(p^2 - p*q*z)*(p*q*z - q^2)/k2(p, q, z)) * alpha(k2(p, q, z)) * Lambda_PV^2 / (k2(p, q, z) + Lambda_PV^2)
        else
            integrand = @. w_z * (p*q*z + 2*(p^2 - p*q*z)*(p*q*z - q^2)/k2(p, q, z)) * alpha(k2(p, q, z))
        end
        return sum(integrand)
    end

    zA_cache = Dict{Tuple{ComplexF64, Float64}, ComplexF64}()

    function z_intA(p::ComplexF64, q::Float64)
        # key = abs.(p) <= q ? (p, q) : (q, p)
        key = (p, q)

        cached = get(zA_cache, key, nothing)
        cached !== nothing && return cached

        result = z_intAA(p, q)
        zA_cache[key] = result

        return result
    end

    function z_intBB(p::ComplexF64, q::Float64)
        if PV == true
            integrand = @. w_z * alpha(k2(p, q, z)) * Lambda_PV^2 / (k2(p, q, z) + Lambda_PV^2)
        else
            integrand = @. w_z * alpha(k2(p, q, z)) 
        end
        return sum(integrand)
    end

    zB_cache = Dict{Tuple{ComplexF64, Float64}, ComplexF64}()

    function z_intB(p::ComplexF64, q::Float64)
        # key = abs.(p) <= q ? (p, q) : (q, p)
        key = (p, q)

        cached = get(zB_cache, key, nothing)
        cached !== nothing && return cached

        result = z_intBB(p, q)
        zB_cache[key] = result

        return result
    end

    function Sigma_A(p2, A, B, Z_2)
        integrand = @. w_t * exp(2*t) * A / (exp(t) * A^2 + B^2) * z_intA(sqrt(Complex(p2)), exp(t/2))
        return Z_2^2 * 16pi/(3*(2pi)^3 * p2) * sum(integrand)
    end

    function Sigma_B(p2, A, B, Z_2)
        integrand = @. w_t * exp(2*t) * B / (exp(t) * A^2 + B^2) * z_intB(sqrt(Complex(p2)), exp(t/2))
        return Z_2^2 * 16pi/(2pi)^3 * sum(integrand)
    end

    function update(A, B, Z_2, Z_4m)
        A_int = Z_2 .+ Sigma_A.(exp.(t), Ref(A), Ref(B), Z_2)
        B_int = Z_4m .+ Sigma_B.(exp.(t), Ref(A), Ref(B), Z_2)
        max_A = maximum(abs.(A .- A_int))
        max_B = maximum(abs.(B .- B_int))
        Z_2  = 1 - Sigma_A(mu^2, A, B, Z_2)
        Z_4m = m - Sigma_B(mu^2, A, B, Z_2)
        return A_int, B_int, Z_2, Z_4m, max(max_A, max_B)
    end

    max_iter = 50
    progress = Progress(41, desc = "Quark DSE:")

    for i in 1:max_iter
        A, B, Z_2, Z_4m, max_error = update(A, B, Z_2, Z_4m)
        next!(progress)
        if max_error < 1e-8
            println()
            break
        end
        if i == max_iter
            throw(error("Failed: max_error = ", max_error))
        end
    end
    A_func(p2) = Z_2 + Sigma_A(p2, A, B, Z_2)
    B_func(p2) = Z_4m + Sigma_B(p2, A, B, Z_2)

    return t, A_func, B_func, Z_2, Z_4m
end

@time t, A, B, Z_2, Z_4m = QuarkDSE(0.4, 0.93312, false)

plot(exp.(t), real.(A.(exp.(t))), xaxis=:log10, xlims = (epsilon2, Lambda2), ylims = (0, 2.0), 
    yticks = 0.4:0.4:2.0)
plot!(exp.(t), real.(B.(exp.(t))))
savefig(address * "realp.pdf")

imagp = -10:0.005:-0.001
plot(imagp, real.(B.(imagp) ./ A.(imagp)), xlims = (-10, 0), ylims = (0, 10.0))
savefig(address * "imagp.pdf")

function EE(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    # k2 = p2 + q2 - 2*pq
    # pk = p2 - pq
    # qk = pq - q2
    # Pk = pP - Pq
    # Delta = pP^2 - p2*P2
    return - 3*Bplus*Bminus - 3*Aplus*Aminus*q2 + 3/4*Aplus*Aminus*P2
end

function EF(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    # k2 = p2 + q2 - 2*pq
    # pk = p2 - pq
    # qk = pq - q2
    # Pk = pP - Pq
    # Delta = pP^2 - p2*P2
    return 3*Aminus*Bplus*Pq - 3/2*Aminus*Bplus*P2 - 3*Aplus*Bminus*Pq - 3/2*Aplus*Bminus*P2
end

function EG(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    # k2 = p2 + q2 - 2*pq
    # pk = p2 - pq
    # qk = pq - q2
    # Pk = pP - Pq
    # Delta = pP^2 - p2*P2
    return - 3/2*Aminus*Bplus*Pq^2 + 3*Aminus*Bplus*q2*Pq - 3/2*Aplus*Bminus*Pq^2 - 3*Aplus*Bminus*q2*Pq
end

function EH(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    # k2 = p2 + q2 - 2*pq
    # pk = p2 - pq
    # qk = pq - q2
    # Pk = pP - Pq
    # Delta = pP^2 - p2*P2
    return - 6*Aplus*Aminus*Pq^2 + 6*Aplus*Aminus*P2*q2
end

function FE(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    k2 = p2 + q2 - 2*pq
    pk = p2 - pq
    qk = pq - q2
    Pk = pP - Pq
    Delta = pP^2 - p2*P2
    return - Aminus*Bplus*Delta^-1*k2^-1*pP*pk*Pk + 2*Aminus*Bplus*Delta^-1*k2^-1*pP*pk*qk + Aminus*Bplus*Delta^-1*pP*pq - 1/2*Aminus*Bplus*Delta^-1*pP^2 + Aminus*Bplus*Delta^-1*p2*k2^-1*Pk^2 - 2*Aminus*Bplus*Delta^-1*p2*k2^-1*qk*Pk - Aminus*Bplus*Delta^-1*p2*Pq + 1/2*Aminus*Bplus*Delta^-1*p2*P2 - Aplus*Bminus*Delta^-1*k2^-1*pP*pk*Pk - 2*Aplus*Bminus*Delta^-1*k2^-1*pP*pk*qk - Aplus*Bminus*Delta^-1*pP*pq - 1/2*Aplus*Bminus*Delta^-1*pP^2 + Aplus*Bminus*Delta^-1*p2*k2^-1*Pk^2 + 2*Aplus*Bminus*Delta^-1*p2*k2^-1*qk*Pk + Aplus*Bminus*Delta^-1*p2*Pq + 1/2*Aplus*Bminus*Delta^-1*p2*P2
end

function FF(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    k2 = p2 + q2 - 2*pq
    pk = p2 - pq
    qk = pq - q2
    Pk = pP - Pq
    Delta = pP^2 - p2*P2
    return 2*Bplus*Bminus*Delta^-1*k2^-1*pP*pk*Pk - 2*Bplus*Bminus*Delta^-1*pP^2 - 2*Bplus*Bminus*Delta^-1*p2*k2^-1*Pk^2 + 2*Bplus*Bminus*Delta^-1*p2*P2 + 3*Bplus*Bminus + 4*Aplus*Aminus*Delta^-1*k2^-1*pP*Pq*pk*qk + 2*Aplus*Aminus*Delta^-1*pP*pq*Pq - 2*Aplus*Aminus*Delta^-1*q2*k2^-1*pP*pk*Pk + 2*Aplus*Aminus*Delta^-1*q2*pP^2 - 1/2*Aplus*Aminus*Delta^-1*P2*k2^-1*pP*pk*Pk - Aplus*Aminus*Delta^-1*P2*pP^2 - 4*Aplus*Aminus*Delta^-1*p2*k2^-1*Pq*qk*Pk - 2*Aplus*Aminus*Delta^-1*p2*Pq^2 + 2*Aplus*Aminus*Delta^-1*p2*q2*k2^-1*Pk^2 + 1/2*Aplus*Aminus*Delta^-1*p2*P2*k2^-1*Pk^2 - 2*Aplus*Aminus*Delta^-1*p2*P2*q2 + Aplus*Aminus*Delta^-1*p2*P2^2 - 3*Aplus*Aminus*q2 + 3/4*Aplus*Aminus*P2
end

function FG(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    k2 = p2 + q2 - 2*pq
    pk = p2 - pq
    qk = pq - q2
    Pk = pP - Pq
    Delta = pP^2 - p2*P2
    return 2*Bplus*Bminus*Delta^-1*k2^-1*pP*Pq*pk*qk + Bplus*Bminus*Delta^-1*pP*pq*Pq - 2*Bplus*Bminus*Delta^-1*p2*k2^-1*Pq*qk*Pk - Bplus*Bminus*Delta^-1*p2*Pq^2 - Aplus*Aminus*Delta^-1*k2^-1*pP*Pq^2*pk*Pk - 1/2*Aplus*Aminus*Delta^-1*pP^2*Pq^2 + 2*Aplus*Aminus*Delta^-1*q2*k2^-1*pP*Pq*pk*qk + Aplus*Aminus*Delta^-1*q2*pP*pq*Pq + 1/2*Aplus*Aminus*Delta^-1*P2*k2^-1*pP*Pq*pk*qk + 1/4*Aplus*Aminus*Delta^-1*P2*pP*pq*Pq + Aplus*Aminus*Delta^-1*p2*k2^-1*Pq^2*Pk^2 - 2*Aplus*Aminus*Delta^-1*p2*q2*k2^-1*Pq*qk*Pk - Aplus*Aminus*Delta^-1*p2*q2*Pq^2 - 1/2*Aplus*Aminus*Delta^-1*p2*P2*k2^-1*Pq*qk*Pk + 1/4*Aplus*Aminus*Delta^-1*p2*P2*Pq^2
end

function FH(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    k2 = p2 + q2 - 2*pq
    pk = p2 - pq
    qk = pq - q2
    Pk = pP - Pq
    Delta = pP^2 - p2*P2
    return 2*Aminus*Bplus*Delta^-1*k2^-1*pP*Pq*pk*Pk + 4*Aminus*Bplus*Delta^-1*k2^-1*pP*Pq*pk*qk + 2*Aminus*Bplus*Delta^-1*pP*pq*Pq - 2*Aminus*Bplus*Delta^-1*pP^2*Pq - 4*Aminus*Bplus*Delta^-1*q2*k2^-1*pP*pk*Pk + 4*Aminus*Bplus*Delta^-1*q2*pP^2 - 2*Aminus*Bplus*Delta^-1*P2*k2^-1*pP*pk*qk - Aminus*Bplus*Delta^-1*P2*pP*pq - 2*Aminus*Bplus*Delta^-1*p2*k2^-1*Pq*Pk^2 - 4*Aminus*Bplus*Delta^-1*p2*k2^-1*Pq*qk*Pk - 2*Aminus*Bplus*Delta^-1*p2*Pq^2 + 4*Aminus*Bplus*Delta^-1*p2*q2*k2^-1*Pk^2 + 2*Aminus*Bplus*Delta^-1*p2*P2*k2^-1*qk*Pk + 3*Aminus*Bplus*Delta^-1*p2*P2*Pq - 4*Aminus*Bplus*Delta^-1*p2*P2*q2 + 3*Aminus*Bplus*Pq - 6*Aminus*Bplus*q2 - 2*Aplus*Bminus*Delta^-1*k2^-1*pP*Pq*pk*Pk + 4*Aplus*Bminus*Delta^-1*k2^-1*pP*Pq*pk*qk + 2*Aplus*Bminus*Delta^-1*pP*pq*Pq + 2*Aplus*Bminus*Delta^-1*pP^2*Pq - 4*Aplus*Bminus*Delta^-1*q2*k2^-1*pP*pk*Pk + 4*Aplus*Bminus*Delta^-1*q2*pP^2 + 2*Aplus*Bminus*Delta^-1*P2*k2^-1*pP*pk*qk + Aplus*Bminus*Delta^-1*P2*pP*pq + 2*Aplus*Bminus*Delta^-1*p2*k2^-1*Pq*Pk^2 - 4*Aplus*Bminus*Delta^-1*p2*k2^-1*Pq*qk*Pk - 2*Aplus*Bminus*Delta^-1*p2*Pq^2 + 4*Aplus*Bminus*Delta^-1*p2*q2*k2^-1*Pk^2 - 2*Aplus*Bminus*Delta^-1*p2*P2*k2^-1*qk*Pk - 3*Aplus*Bminus*Delta^-1*p2*P2*Pq - 4*Aplus*Bminus*Delta^-1*p2*P2*q2 - 3*Aplus*Bminus*Pq - 6*Aplus*Bminus*q2
end

function GE(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    k2 = p2 + q2 - 2*pq
    pk = p2 - pq
    qk = pq - q2
    Pk = pP - Pq
    Delta = pP^2 - p2*P2
    return - Aminus*Bplus*Delta^-1*k2^-1*Pk^2 + 2*Aminus*Bplus*Delta^-1*k2^-1*qk*Pk + Aminus*Bplus*Delta^-1*Pq + Aminus*Bplus*Delta^-1*P2*k2^-1*pP^-1*pk*Pk - 2*Aminus*Bplus*Delta^-1*P2*k2^-1*pP^-1*pk*qk - Aminus*Bplus*Delta^-1*P2*pP^-1*pq - Aplus*Bminus*Delta^-1*k2^-1*Pk^2 - 2*Aplus*Bminus*Delta^-1*k2^-1*qk*Pk - Aplus*Bminus*Delta^-1*Pq + Aplus*Bminus*Delta^-1*P2*k2^-1*pP^-1*pk*Pk + 2*Aplus*Bminus*Delta^-1*P2*k2^-1*pP^-1*pk*qk + Aplus*Bminus*Delta^-1*P2*pP^-1*pq
end

function GF(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    k2 = p2 + q2 - 2*pq
    pk = p2 - pq
    qk = pq - q2
    Pk = pP - Pq
    Delta = pP^2 - p2*P2
    return 2*Bplus*Bminus*Delta^-1*k2^-1*Pk^2 - 2*Bplus*Bminus*Delta^-1*P2*k2^-1*pP^-1*pk*Pk + 4*Aplus*Aminus*Delta^-1*k2^-1*Pq*qk*Pk + 2*Aplus*Aminus*Delta^-1*Pq^2 - 2*Aplus*Aminus*Delta^-1*q2*k2^-1*Pk^2 - 4*Aplus*Aminus*Delta^-1*P2*k2^-1*pP^-1*Pq*pk*qk - 1/2*Aplus*Aminus*Delta^-1*P2*k2^-1*Pk^2 - 2*Aplus*Aminus*Delta^-1*P2*pP^-1*pq*Pq + 2*Aplus*Aminus*Delta^-1*P2*q2*k2^-1*pP^-1*pk*Pk + 1/2*Aplus*Aminus*Delta^-1*P2^2*k2^-1*pP^-1*pk*Pk
end

function GG(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    k2 = p2 + q2 - 2*pq
    pk = p2 - pq
    qk = pq - q2
    Pk = pP - Pq
    Delta = pP^2 - p2*P2
    return 2*Bplus*Bminus*Delta^-1*k2^-1*Pq*qk*Pk + Bplus*Bminus*Delta^-1*Pq^2 - 2*Bplus*Bminus*Delta^-1*P2*k2^-1*pP^-1*Pq*pk*qk - Bplus*Bminus*Delta^-1*P2*pP^-1*pq*Pq - Aplus*Aminus*Delta^-1*k2^-1*Pq^2*Pk^2 + 2*Aplus*Aminus*Delta^-1*q2*k2^-1*Pq*qk*Pk + Aplus*Aminus*Delta^-1*q2*Pq^2 + Aplus*Aminus*Delta^-1*P2*k2^-1*pP^-1*Pq^2*pk*Pk + 1/2*Aplus*Aminus*Delta^-1*P2*k2^-1*Pq*qk*Pk + 1/4*Aplus*Aminus*Delta^-1*P2*Pq^2 - 2*Aplus*Aminus*Delta^-1*P2*q2*k2^-1*pP^-1*Pq*pk*qk - Aplus*Aminus*Delta^-1*P2*q2*pP^-1*pq*Pq - 1/2*Aplus*Aminus*Delta^-1*P2^2*k2^-1*pP^-1*Pq*pk*qk - 1/4*Aplus*Aminus*Delta^-1*P2^2*pP^-1*pq*Pq
end

function GH(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    k2 = p2 + q2 - 2*pq
    pk = p2 - pq
    qk = pq - q2
    Pk = pP - Pq
    Delta = pP^2 - p2*P2
    return 2*Aminus*Bplus*Delta^-1*k2^-1*Pq*Pk^2 + 4*Aminus*Bplus*Delta^-1*k2^-1*Pq*qk*Pk + 2*Aminus*Bplus*Delta^-1*Pq^2 - 4*Aminus*Bplus*Delta^-1*q2*k2^-1*Pk^2 - 2*Aminus*Bplus*Delta^-1*P2*k2^-1*pP^-1*Pq*pk*Pk - 4*Aminus*Bplus*Delta^-1*P2*k2^-1*pP^-1*Pq*pk*qk - 2*Aminus*Bplus*Delta^-1*P2*k2^-1*qk*Pk - 2*Aminus*Bplus*Delta^-1*P2*pP^-1*pq*Pq - Aminus*Bplus*Delta^-1*P2*Pq + 4*Aminus*Bplus*Delta^-1*P2*q2*k2^-1*pP^-1*pk*Pk + 2*Aminus*Bplus*Delta^-1*P2^2*k2^-1*pP^-1*pk*qk + Aminus*Bplus*Delta^-1*P2^2*pP^-1*pq - 2*Aplus*Bminus*Delta^-1*k2^-1*Pq*Pk^2 + 4*Aplus*Bminus*Delta^-1*k2^-1*Pq*qk*Pk + 2*Aplus*Bminus*Delta^-1*Pq^2 - 4*Aplus*Bminus*Delta^-1*q2*k2^-1*Pk^2 + 2*Aplus*Bminus*Delta^-1*P2*k2^-1*pP^-1*Pq*pk*Pk - 4*Aplus*Bminus*Delta^-1*P2*k2^-1*pP^-1*Pq*pk*qk + 2*Aplus*Bminus*Delta^-1*P2*k2^-1*qk*Pk - 2*Aplus*Bminus*Delta^-1*P2*pP^-1*pq*Pq + Aplus*Bminus*Delta^-1*P2*Pq + 4*Aplus*Bminus*Delta^-1*P2*q2*k2^-1*pP^-1*pk*Pk - 2*Aplus*Bminus*Delta^-1*P2^2*k2^-1*pP^-1*pk*qk - Aplus*Bminus*Delta^-1*P2^2*pP^-1*pq
end

function HE(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    k2 = p2 + q2 - 2*pq
    pk = p2 - pq
    qk = pq - q2
    Pk = pP - Pq
    Delta = pP^2 - p2*P2
    return Aplus*Aminus*Delta^-1*k2^-1*Pq*pk*Pk - Aplus*Aminus*Delta^-1*k2^-1*pq*Pk^2 + Aplus*Aminus*Delta^-1*k2^-1*pP*qk*Pk - 1/2*Aplus*Aminus*Delta^-1*pP*Pq - Aplus*Aminus*Delta^-1*P2*k2^-1*pk*qk + 1/2*Aplus*Aminus*Delta^-1*P2*pq
end

function HF(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    k2 = p2 + q2 - 2*pq
    pk = p2 - pq
    qk = pq - q2
    Pk = pP - Pq
    Delta = pP^2 - p2*P2
    return - Aminus*Bplus*Delta^-1*k2^-1*Pq*pk*Pk + Aminus*Bplus*Delta^-1*k2^-1*pq*Pk^2 - Aminus*Bplus*Delta^-1*k2^-1*pP*qk*Pk + 1/2*Aminus*Bplus*Delta^-1*pP*Pq + Aminus*Bplus*Delta^-1*P2*k2^-1*pk*qk - 1/2*Aminus*Bplus*Delta^-1*P2*pq - Aplus*Bminus*Delta^-1*k2^-1*Pq*pk*Pk + Aplus*Bminus*Delta^-1*k2^-1*pq*Pk^2 - Aplus*Bminus*Delta^-1*k2^-1*pP*qk*Pk + 1/2*Aplus*Bminus*Delta^-1*pP*Pq + Aplus*Bminus*Delta^-1*P2*k2^-1*pk*qk - 1/2*Aplus*Bminus*Delta^-1*P2*pq
end

function HG(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    k2 = p2 + q2 - 2*pq
    pk = p2 - pq
    qk = pq - q2
    Pk = pP - Pq
    Delta = pP^2 - p2*P2
    return - 1/2*Aminus*Bplus*Delta^-1*k2^-1*Pq^2*pk*Pk + 1/2*Aminus*Bplus*Delta^-1*k2^-1*pq*Pq*Pk^2 - 1/2*Aminus*Bplus*Delta^-1*k2^-1*pP*Pq*qk*Pk + 1/4*Aminus*Bplus*Delta^-1*pP*Pq^2 + 1/2*Aminus*Bplus*Delta^-1*P2*k2^-1*Pq*pk*qk - 1/4*Aminus*Bplus*Delta^-1*P2*pq*Pq + 1/2*Aplus*Bminus*Delta^-1*k2^-1*Pq^2*pk*Pk - 1/2*Aplus*Bminus*Delta^-1*k2^-1*pq*Pq*Pk^2 + 1/2*Aplus*Bminus*Delta^-1*k2^-1*pP*Pq*qk*Pk - 1/4*Aplus*Bminus*Delta^-1*pP*Pq^2 - 1/2*Aplus*Bminus*Delta^-1*P2*k2^-1*Pq*pk*qk + 1/4*Aplus*Bminus*Delta^-1*P2*pq*Pq
end

function HH(p2, q2, P2, pP, pq, Pq, Aplus, Aminus, Bplus, Bminus)
    k2 = p2 + q2 - 2*pq
    pk = p2 - pq
    qk = pq - q2
    Pk = pP - Pq
    Delta = pP^2 - p2*P2
    return - 2*Bplus*Bminus*Delta^-1*k2^-1*Pq*pk*Pk + 2*Bplus*Bminus*Delta^-1*k2^-1*pq*Pk^2 - 2*Bplus*Bminus*Delta^-1*k2^-1*pP*qk*Pk + Bplus*Bminus*Delta^-1*pP*Pq + 2*Bplus*Bminus*Delta^-1*P2*k2^-1*pk*qk - Bplus*Bminus*Delta^-1*P2*pq + 2*Aplus*Aminus*Delta^-1*q2*k2^-1*Pq*pk*Pk - 2*Aplus*Aminus*Delta^-1*q2*k2^-1*pq*Pk^2 + 2*Aplus*Aminus*Delta^-1*q2*k2^-1*pP*qk*Pk - Aplus*Aminus*Delta^-1*q2*pP*Pq - 1/2*Aplus*Aminus*Delta^-1*P2*k2^-1*Pq*pk*Pk + 1/2*Aplus*Aminus*Delta^-1*P2*k2^-1*pq*Pk^2 - 1/2*Aplus*Aminus*Delta^-1*P2*k2^-1*pP*qk*Pk + 1/4*Aplus*Aminus*Delta^-1*P2*pP*Pq - 2*Aplus*Aminus*Delta^-1*P2*q2*k2^-1*pk*qk + Aplus*Aminus*Delta^-1*P2*q2*pq + 1/2*Aplus*Aminus*Delta^-1*P2^2*k2^-1*pk*qk - 1/4*Aplus*Aminus*Delta^-1*P2^2*pq
end

function build_AB(M, radial_steps, angular_steps)
    x, _ = gausslegendre(radial_steps)
    z, _ = gausslegendre(angular_steps)
    t = 0.5 * (log(Lambda2) - log(epsilon2)) * x .+ 0.5 * (log(Lambda2) + log(epsilon2))
    q2 = exp.(t)

    Ap = Vector{ComplexF64}(undef, radial_steps*angular_steps)
    Am = Vector{ComplexF64}(undef, radial_steps*angular_steps)
    Bp = Vector{ComplexF64}(undef, radial_steps*angular_steps)
    Bm = Vector{ComplexF64}(undef, radial_steps*angular_steps)

    for k in 1:radial_steps
        for l in 1:angular_steps
            kl = (k-1)*angular_steps + l
            qplus2  = q2[k] + im*M*sqrt(q2[k])*z[l] - M^2/4
            qminus2 = q2[k] - im*M*sqrt(q2[k])*z[l] - M^2/4
            Ap[kl] = A(qplus2)
            Am[kl] = A(qminus2)
            Bp[kl] = B(qplus2)
            Bm[kl] = B(qminus2)
        end
    end
    return Ap, Am, Bp, Bm
end

function entries(Func, M, p2, q2, z_q, Aplus, Aminus, Bplus, Bminus; w = 0.4, D = 0.93312, y_steps = 32, beta = pi/4)
    pP = im*M*sqrt(p2)*cos(beta)
    Pq = im*M*sqrt(q2)*z_q
    pq(y) = sqrt(p2*q2)*(sin(beta)*y*sqrt(1-z_q^2) + cos(beta)*z_q)
    qplus2 = q2 + Pq - M^2/4
    qminus2 = q2 - Pq - M^2/4
    alpha_UV(k2) = 2pi * gamma_m * (1 - exp(-k2)) / (k2 * log(exp(2)-1 + (1 + k2/Lambda_QCD^2)^2))
    alpha_IR(k2) = D/w^6 * pi * k2 * exp(-k2 / w^2)
    alpha(k2) = alpha_IR(k2) + alpha_UV(k2)
    k2_val(y) = p2 + q2 - 2*pq(y)
    y, w_y = gausslegendre(y_steps)
    y_integral = sum(@. w_y * alpha(k2_val(y)) * Func(p2, q2, -M^2, pP, pq(y), Pq, Aplus, Aminus, Bplus, Bminus))
    return - C_f * 1/(2pi)^2 * Z_2^2 * q2^2 * sqrt(1-z_q^2)/((qplus2*Aplus^2 + Bplus^2)*(qminus2*Aminus^2 + Bminus^2)) * y_integral
end

function build_kappa(Func, M, t, z, w_t, w_z)
    radial_steps, angular_steps = length(t), length(z)
    N = radial_steps * angular_steps
    q2 = exp.(t)
    Ap, Am, Bp, Bm = build_AB(M, radial_steps, angular_steps)
    Kappa = Matrix{ComplexF64}(undef, N, N)

    # Threads.@threads for k in 1:radial_steps
    for k in 1:radial_steps
        for l in 1:angular_steps
            col = (k - 1) * angular_steps + l
            w = w_t[k] * w_z[l]
            for i in 1:radial_steps
                val = w * entries(Func, M, q2[i], q2[k], z[l], Ap[col], Am[col], Bp[col], Bm[col])
                rows = (i - 1) * angular_steps + 1 : i * angular_steps
                @views Kappa[rows, col] .= val
            end
        end
    end
    return Kappa
end

function lambda(M; radial_steps::Int = 64, angular_steps::Int = 16)
    x, w_x = gausslegendre(radial_steps)
    z, w_z = gausslegendre(angular_steps)

    t = 0.5 * (log(Lambda2) - log(epsilon2)) * x .+ 0.5 * (log(Lambda2) + log(epsilon2))
    w_t = 0.5 * (log(Lambda2) - log(epsilon2)) * w_x

    progress = Progress(16, desc = "BSE, M = $(round(M; digits = 3)):")
    kEE = build_kappa(EE, M, t, z, w_t, w_z)
    next!(progress)
    kEF = build_kappa(EF, M, t, z, w_t, w_z)
    next!(progress)
    kEG = build_kappa(EG, M, t, z, w_t, w_z)
    next!(progress)
    kEH = build_kappa(EH, M, t, z, w_t, w_z)
    next!(progress)
    kFE = build_kappa(FE, M, t, z, w_t, w_z)
    next!(progress)
    kFF = build_kappa(FF, M, t, z, w_t, w_z)
    next!(progress)
    kFG = build_kappa(FG, M, t, z, w_t, w_z)
    next!(progress)
    kFH = build_kappa(FH, M, t, z, w_t, w_z)
    next!(progress)
    kGE = build_kappa(GE, M, t, z, w_t, w_z)
    next!(progress)
    kGF = build_kappa(GF, M, t, z, w_t, w_z)
    next!(progress)
    kGG = build_kappa(GG, M, t, z, w_t, w_z)
    next!(progress)
    kGH = build_kappa(GH, M, t, z, w_t, w_z)
    next!(progress)
    kHE = build_kappa(HE, M, t, z, w_t, w_z)
    next!(progress)
    kHF = build_kappa(HF, M, t, z, w_t, w_z)
    next!(progress)
    kHG = build_kappa(HG, M, t, z, w_t, w_z)
    next!(progress)
    kHH = build_kappa(HH, M, t, z, w_t, w_z)
    next!(progress)
    
    Kappa = [kEE kEF kEG kEH;
             kFE kFF kFG kFH;
             kGE kGF kGG kGH;
             kHE kHF kHG kHH]

    @time states = eigen(Kappa)
    return real(states.values[end]), real.(states.vectors[:, end])
end

@time find_zero(M -> lambda(M)[1] - 1, 0.141) # 0.1408955741336375
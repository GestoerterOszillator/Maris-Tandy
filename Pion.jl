using Plots
using FastGaussQuadrature
using LaTeXStrings
using ProgressMeter
pgfplotsx() # pgfplotsx() or gr()
default(
    fontfamily = "Computer Modern",
    tic(p-q)font   = font(15, "Computer Modern"),
    guidefont  = font(15, "Computer Modern"),
    legendfont = font(15, "Computer Modern"),
    lw = 2,
    grid = :on,
    gridalpha = 0.3,
    size = (500, 400),
    legend = :topright,
    framestyle = :box,
    color_palette = palette(:tab10)
)
# const address = "/Users/johnreeg/Documents/Repositories/Maris-Tandy/thesis/images/"
const address = "/home/john-reeg/Documents/Maris-Tandy/thesis/images/"

# Constants
const m = 0.0037
const Lambda2 = 1e6
const epsilon2 = 1e-4
const Lambda_PV = 2e2
const Lambda_QCD = 0.234
const gamma_m = 0.48 # 12/(33 - 2 * Nf{4}) = 0.48
const mu = 19.0

function Teil_Eins(w::Float64, D::Float64, PV::Bool; radial_steps::Int = 256, angular_steps::Int = 64)
    # Integration Grids
    x, w_x = gausslegendre(radial_steps)
    z, w_z = gausslegendre(angular_steps) # Chebyshev second kind

    t = 0.5 * (log(Lambda2) - log(epsilon2)) * x .+ 0.5 * (log(Lambda2) + log(epsilon2))
    w_t = 0.5 * (log(Lambda2) - log(epsilon2)) * w_x

    # Initial Values
    A = 1.6 * ones(length(t))
    B = (m + 1e-3) * ones(length(t))
    Z_2 = 1
    Z_4m = 0.7 * m

    # Dot Products
    k2(p::ComplexF64, q::Float64, z::Float64) = p^2 + q^2 - 2*p*q*z

    # Gluon stuff
    alpha_UV(k2::ComplexF64) = 2pi * gamma_m * (1 - exp(-k2)) / (k2 * log(exp(2)-1 + (1 + k2/Lambda_QCD^2)^2))
    alpha_IR(k2::ComplexF64) = D/w^6 * pi * k2 * exp(-k2 / w^2)
    alpha(k2::ComplexF64) = alpha_IR(k2) + alpha_UV(k2)

    # Angular Integrals over z = cos(psi)
    function z_intAA(p::ComplexF64, q::Float64)
        if PV == true
            integrand = @. w_z * sqrt(1-z^2) * (p*q*z + 2*(p^2 - p*q*z)*(p*q*z - q^2)/k2(p, q, z)) * alpha(k2(p, q, z)) * Lambda_PV^2 / (k2(p, q, z) + Lambda_PV^2)
        else
            integrand = @. w_z * sqrt(1-z^2) * (p*q*z + 2*(p^2 - p*q*z)*(p*q*z - q^2)/k2(p, q, z)) * alpha(k2(p, q, z))
        end
        return sum(integrand)
    end

    zA_cache = Dict{Tuple{ComplexF64,Float64}, Float64}()

    function z_intA(p::ComplexF64, q::Float64)
        key = p <= q ? (p, q) : (q, p)

        cached = get(zA_cache, key, nothing)
        cached !== nothing && return cached

        result = z_intAA(p, q)
        zA_cache[key] = result

        return result
    end

    function z_intBB(p::ComplexF64, q::Float64)
        if PV == true
            integrand = @. w_z * sqrt(1-z^2) * alpha(k2(p, q, z)) * Lambda_PV^2 / (k2(p, q, z) + Lambda_PV^2)
        else
            integrand = @. w_z * sqrt(1-z^2) * alpha(k2(p, q, z)) 
        end
        return sum(integrand)
    end

    zB_cache = Dict{Tuple{ComplexF64,Float64}, Float64}()

    function z_intB(p::ComplexF64, q::Float64)
        key = p <= q ? (p, q) : (q, p)

        cached = get(zB_cache, key, nothing)
        cached !== nothing && return cached

        result = z_intBB(p, q)
        zB_cache[key] = result

        return result
    end

    # Final Radial Integrals
    function Sigma_A(p2, A, B, Z_2)
        integrand = @. w_t * exp(2*t) * A / (exp(t) * A^2 + B^2) * z_intA(sqrt(p2), exp(t/2))
        return Z_2^2 * 16pi/(3*(2pi)^3 * p2) * sum(integrand)
    end

    function Sigma_B(p2, A, B, Z_2)
        integrand = @. w_t * exp(2*t) * B / (exp(t) * A^2 + B^2) * z_intB(sqrt(p2), exp(t/2))
        return Z_2^2 * 16pi/(2pi)^3 * sum(integrand)
    end

    # Iteration Step Function
    function update(A, B, Z_2, Z_4m)
        A_int = Z_2 .+ Sigma_A.(exp.(t), Ref(A), Ref(B), Z_2)
        B_int = Z_4m .+ Sigma_B.(exp.(t), Ref(A), Ref(B), Z_2)
        max_A = maximum(abs.(A .- A_int))
        max_B = maximum(abs.(B .- B_int))
        Z_2  = 1 - Sigma_A(mu^2, A, B, Z_2)
        Z_4m = m - Sigma_B(mu^2, A, B, Z_2)
        return A_int, B_int, Z_2, Z_4m, max(max_A, max_B)
    end

    # Iteration Loop
    max_iter = 100
    p = Progress(max_iter, desc = "Berechne...")

    for i in 1:max_iter
        A, B, Z_2, Z_4m, max_error = update(A, B, Z_2, Z_4m)
        next!(p)
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

@time t, A, B, Z_2, Z_4m = Teil_Eins(0.4, 1.0, true)

function EE(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return  - 3*Bplus*Bminus - 3*Aplus*Aminus*q2 + 3/4*Aplus*Aminus*P2
end

function EF(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return 3*Aminus*Bplus*Pq - 3/2*Aminus*Bplus*P2 - 3*Aplus*Bminus*Pq - 3/2*
      Aplus*Bminus*P2
end

function EG(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return  - 3/2*Aminus*Bplus*Pq^2 + 3*Aminus*Bplus*q2*Pq - 3/2*Aplus*Bminus*
      Pq^2 - 3*Aplus*Bminus*q2*Pq
end

function EH(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return  - 6*Aplus*Aminus*Pq^2 + 6*Aplus*Aminus*P2*q2
end

function FE(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return  - Aminus*Bplus*Delta^-1*k2^-1*pP*pk*Pk + 2*Aminus*Bplus*Delta^-1*
      k2^-1*pP*pk*qk + Aminus*Bplus*Delta^-1*pP*pq - 1/2*Aminus*Bplus*Delta^-1*
      pP^2 + Aminus*Bplus*Delta^-1*p2*k2^-1*Pk^2 - 2*Aminus*Bplus*Delta^-1*p2*
      k2^-1*qk*Pk - Aminus*Bplus*Delta^-1*p2*Pq + 1/2*Aminus*Bplus*Delta^-1*
      p2*P2 - Aplus*Bminus*Delta^-1*k2^-1*pP*pk*Pk - 2*Aplus*Bminus*Delta^-1*
      k2^-1*pP*pk*qk - Aplus*Bminus*Delta^-1*pP*pq - 1/2*Aplus*Bminus*Delta^-1*
      pP^2 + Aplus*Bminus*Delta^-1*p2*k2^-1*Pk^2 + 2*Aplus*Bminus*Delta^-1*p2*
      k2^-1*qk*Pk + Aplus*Bminus*Delta^-1*p2*Pq + 1/2*Aplus*Bminus*Delta^-1*
      p2*P2
end

function FF(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return 2*Bplus*Bminus*Delta^-1*k2^-1*pP*pk*Pk - 2*Bplus*Bminus*Delta^-1*
      pP^2 - 2*Bplus*Bminus*Delta^-1*p2*k2^-1*Pk^2 + 2*Bplus*Bminus*Delta^-1*
      p2*P2 + 3*Bplus*Bminus + 4*Aplus*Aminus*Delta^-1*k2^-1*pP*Pq*pk*qk + 2*
      Aplus*Aminus*Delta^-1*pP*pq*Pq - 2*Aplus*Aminus*Delta^-1*q2*k2^-1*pP*pk*
      Pk + 2*Aplus*Aminus*Delta^-1*q2*pP^2 - 1/2*Aplus*Aminus*Delta^-1*P2*
      k2^-1*pP*pk*Pk - Aplus*Aminus*Delta^-1*P2*pP^2 - 4*Aplus*Aminus*Delta^-1*
      p2*k2^-1*Pq*qk*Pk - 2*Aplus*Aminus*Delta^-1*p2*Pq^2 + 2*Aplus*Aminus*
      Delta^-1*p2*q2*k2^-1*Pk^2 + 1/2*Aplus*Aminus*Delta^-1*p2*P2*k2^-1*Pk^2
       - 2*Aplus*Aminus*Delta^-1*p2*P2*q2 + Aplus*Aminus*Delta^-1*p2*P2^2 - 3*
      Aplus*Aminus*q2 + 3/4*Aplus*Aminus*P2
end

function FG(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return 2*Bplus*Bminus*Delta^-1*k2^-1*pP*Pq*pk*qk + Bplus*Bminus*Delta^-1*
      pP*pq*Pq - 2*Bplus*Bminus*Delta^-1*p2*k2^-1*Pq*qk*Pk - Bplus*Bminus*
      Delta^-1*p2*Pq^2 - Aplus*Aminus*Delta^-1*k2^-1*pP*Pq^2*pk*Pk - 1/2*Aplus*
      Aminus*Delta^-1*pP^2*Pq^2 + 2*Aplus*Aminus*Delta^-1*q2*k2^-1*pP*Pq*pk*
      qk + Aplus*Aminus*Delta^-1*q2*pP*pq*Pq + 1/2*Aplus*Aminus*Delta^-1*P2*
      k2^-1*pP*Pq*pk*qk + 1/4*Aplus*Aminus*Delta^-1*P2*pP*pq*Pq + Aplus*Aminus*
      Delta^-1*p2*k2^-1*Pq^2*Pk^2 - 2*Aplus*Aminus*Delta^-1*p2*q2*k2^-1*Pq*qk*
      Pk - Aplus*Aminus*Delta^-1*p2*q2*Pq^2 - 1/2*Aplus*Aminus*Delta^-1*p2*P2*
      k2^-1*Pq*qk*Pk + 1/4*Aplus*Aminus*Delta^-1*p2*P2*Pq^2
end

function FH(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return 2*Aminus*Bplus*Delta^-1*k2^-1*pP*Pq*pk*Pk + 4*Aminus*Bplus*Delta^-1*
    k2^-1*pP*Pq*pk*qk + 2*Aminus*Bplus*Delta^-1*pP*pq*Pq - 2*Aminus*Bplus*
      Delta^-1*pP^2*Pq - 4*Aminus*Bplus*Delta^-1*q2*k2^-1*pP*pk*Pk + 4*Aminus*
      Bplus*Delta^-1*q2*pP^2 - 2*Aminus*Bplus*Delta^-1*P2*k2^-1*pP*pk*qk - 
      Aminus*Bplus*Delta^-1*P2*pP*pq - 2*Aminus*Bplus*Delta^-1*p2*k2^-1*Pq*
      Pk^2 - 4*Aminus*Bplus*Delta^-1*p2*k2^-1*Pq*qk*Pk - 2*Aminus*Bplus*
      Delta^-1*p2*Pq^2 + 4*Aminus*Bplus*Delta^-1*p2*q2*k2^-1*Pk^2 + 2*Aminus*
      Bplus*Delta^-1*p2*P2*k2^-1*qk*Pk + 3*Aminus*Bplus*Delta^-1*p2*P2*Pq - 4*
      Aminus*Bplus*Delta^-1*p2*P2*q2 + 3*Aminus*Bplus*Pq - 6*Aminus*Bplus*q2
       - 2*Aplus*Bminus*Delta^-1*k2^-1*pP*Pq*pk*Pk + 4*Aplus*Bminus*Delta^-1*
      k2^-1*pP*Pq*pk*qk + 2*Aplus*Bminus*Delta^-1*pP*pq*Pq + 2*Aplus*Bminus*
      Delta^-1*pP^2*Pq - 4*Aplus*Bminus*Delta^-1*q2*k2^-1*pP*pk*Pk + 4*Aplus*
      Bminus*Delta^-1*q2*pP^2 + 2*Aplus*Bminus*Delta^-1*P2*k2^-1*pP*pk*qk + 
      Aplus*Bminus*Delta^-1*P2*pP*pq + 2*Aplus*Bminus*Delta^-1*p2*k2^-1*Pq*
      Pk^2 - 4*Aplus*Bminus*Delta^-1*p2*k2^-1*Pq*qk*Pk - 2*Aplus*Bminus*
      Delta^-1*p2*Pq^2 + 4*Aplus*Bminus*Delta^-1*p2*q2*k2^-1*Pk^2 - 2*Aplus*
      Bminus*Delta^-1*p2*P2*k2^-1*qk*Pk - 3*Aplus*Bminus*Delta^-1*p2*P2*Pq - 4*
      Aplus*Bminus*Delta^-1*p2*P2*q2 - 3*Aplus*Bminus*Pq - 6*Aplus*Bminus*q2
end

function GE(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return  - Aminus*Bplus*Delta^-1*k2^-1*Pk^2 + 2*Aminus*Bplus*Delta^-1*k2^-1*
      qk*Pk + Aminus*Bplus*Delta^-1*Pq + Aminus*Bplus*Delta^-1*P2*k2^-1*pP^-1*
      pk*Pk - 2*Aminus*Bplus*Delta^-1*P2*k2^-1*pP^-1*pk*qk - Aminus*Bplus*
      Delta^-1*P2*pP^-1*pq - Aplus*Bminus*Delta^-1*k2^-1*Pk^2 - 2*Aplus*Bminus*
      Delta^-1*k2^-1*qk*Pk - Aplus*Bminus*Delta^-1*Pq + Aplus*Bminus*Delta^-1*
      P2*k2^-1*pP^-1*pk*Pk + 2*Aplus*Bminus*Delta^-1*P2*k2^-1*pP^-1*pk*qk + 
      Aplus*Bminus*Delta^-1*P2*pP^-1*pq
end

function GF(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return 2*Bplus*Bminus*Delta^-1*k2^-1*Pk^2 - 2*Bplus*Bminus*Delta^-1*P2*
      k2^-1*pP^-1*pk*Pk + 4*Aplus*Aminus*Delta^-1*k2^-1*Pq*qk*Pk + 2*Aplus*
      Aminus*Delta^-1*Pq^2 - 2*Aplus*Aminus*Delta^-1*q2*k2^-1*Pk^2 - 4*Aplus*
      Aminus*Delta^-1*P2*k2^-1*pP^-1*Pq*pk*qk - 1/2*Aplus*Aminus*Delta^-1*P2*
      k2^-1*Pk^2 - 2*Aplus*Aminus*Delta^-1*P2*pP^-1*pq*Pq + 2*Aplus*Aminus*
      Delta^-1*P2*q2*k2^-1*pP^-1*pk*Pk + 1/2*Aplus*Aminus*Delta^-1*P2^2*k2^-1*
      pP^-1*pk*Pk
end

function GG(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return 2*Bplus*Bminus*Delta^-1*k2^-1*Pq*qk*Pk + Bplus*Bminus*Delta^-1*Pq^2
       - 2*Bplus*Bminus*Delta^-1*P2*k2^-1*pP^-1*Pq*pk*qk - Bplus*Bminus*
      Delta^-1*P2*pP^-1*pq*Pq - Aplus*Aminus*Delta^-1*k2^-1*Pq^2*Pk^2 + 2*
      Aplus*Aminus*Delta^-1*q2*k2^-1*Pq*qk*Pk + Aplus*Aminus*Delta^-1*q2*Pq^2
       + Aplus*Aminus*Delta^-1*P2*k2^-1*pP^-1*Pq^2*pk*Pk + 1/2*Aplus*Aminus*
      Delta^-1*P2*k2^-1*Pq*qk*Pk + 1/4*Aplus*Aminus*Delta^-1*P2*Pq^2 - 2*Aplus*
      Aminus*Delta^-1*P2*q2*k2^-1*pP^-1*Pq*pk*qk - Aplus*Aminus*Delta^-1*P2*
      q2*pP^-1*pq*Pq - 1/2*Aplus*Aminus*Delta^-1*P2^2*k2^-1*pP^-1*Pq*pk*qk - 1/
      4*Aplus*Aminus*Delta^-1*P2^2*pP^-1*pq*Pq
end

function GH(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return 2*Aminus*Bplus*Delta^-1*k2^-1*Pq*Pk^2 + 4*Aminus*Bplus*Delta^-1*
      k2^-1*Pq*qk*Pk + 2*Aminus*Bplus*Delta^-1*Pq^2 - 4*Aminus*Bplus*Delta^-1*
      q2*k2^-1*Pk^2 - 2*Aminus*Bplus*Delta^-1*P2*k2^-1*pP^-1*Pq*pk*Pk - 4*
      Aminus*Bplus*Delta^-1*P2*k2^-1*pP^-1*Pq*pk*qk - 2*Aminus*Bplus*Delta^-1*
      P2*k2^-1*qk*Pk - 2*Aminus*Bplus*Delta^-1*P2*pP^-1*pq*Pq - Aminus*Bplus*
      Delta^-1*P2*Pq + 4*Aminus*Bplus*Delta^-1*P2*q2*k2^-1*pP^-1*pk*Pk + 2*
      Aminus*Bplus*Delta^-1*P2^2*k2^-1*pP^-1*pk*qk + Aminus*Bplus*Delta^-1*
      P2^2*pP^-1*pq - 2*Aplus*Bminus*Delta^-1*k2^-1*Pq*Pk^2 + 4*Aplus*Bminus*
      Delta^-1*k2^-1*Pq*qk*Pk + 2*Aplus*Bminus*Delta^-1*Pq^2 - 4*Aplus*Bminus*
      Delta^-1*q2*k2^-1*Pk^2 + 2*Aplus*Bminus*Delta^-1*P2*k2^-1*pP^-1*Pq*pk*Pk
       - 4*Aplus*Bminus*Delta^-1*P2*k2^-1*pP^-1*Pq*pk*qk + 2*Aplus*Bminus*
      Delta^-1*P2*k2^-1*qk*Pk - 2*Aplus*Bminus*Delta^-1*P2*pP^-1*pq*Pq + Aplus*
      Bminus*Delta^-1*P2*Pq + 4*Aplus*Bminus*Delta^-1*P2*q2*k2^-1*pP^-1*pk*Pk
       - 2*Aplus*Bminus*Delta^-1*P2^2*k2^-1*pP^-1*pk*qk - Aplus*Bminus*
      Delta^-1*P2^2*pP^-1*pq
end

function HE(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return Aplus*Aminus*Delta^-1*k2^-1*Pq*pk*Pk - Aplus*Aminus*Delta^-1*k2^-1*
      pq*Pk^2 + Aplus*Aminus*Delta^-1*k2^-1*pP*qk*Pk - 1/2*Aplus*Aminus*
      Delta^-1*pP*Pq - Aplus*Aminus*Delta^-1*P2*k2^-1*pk*qk + 1/2*Aplus*Aminus*
      Delta^-1*P2*pq
end

function HF(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return  - Aminus*Bplus*Delta^-1*k2^-1*Pq*pk*Pk + Aminus*Bplus*Delta^-1*
      k2^-1*pq*Pk^2 - Aminus*Bplus*Delta^-1*k2^-1*pP*qk*Pk + 1/2*Aminus*Bplus*
      Delta^-1*pP*Pq + Aminus*Bplus*Delta^-1*P2*k2^-1*pk*qk - 1/2*Aminus*Bplus*
      Delta^-1*P2*pq - Aplus*Bminus*Delta^-1*k2^-1*Pq*pk*Pk + Aplus*Bminus*
      Delta^-1*k2^-1*pq*Pk^2 - Aplus*Bminus*Delta^-1*k2^-1*pP*qk*Pk + 1/2*
      Aplus*Bminus*Delta^-1*pP*Pq + Aplus*Bminus*Delta^-1*P2*k2^-1*pk*qk - 1/2*
      Aplus*Bminus*Delta^-1*P2*pq
end

function HG(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return  - 1/2*Aminus*Bplus*Delta^-1*k2^-1*Pq^2*pk*Pk + 1/2*Aminus*Bplus*
      Delta^-1*k2^-1*pq*Pq*Pk^2 - 1/2*Aminus*Bplus*Delta^-1*k2^-1*pP*Pq*qk*Pk
       + 1/4*Aminus*Bplus*Delta^-1*pP*Pq^2 + 1/2*Aminus*Bplus*Delta^-1*P2*
      k2^-1*Pq*pk*qk - 1/4*Aminus*Bplus*Delta^-1*P2*pq*Pq + 1/2*Aplus*Bminus*
      Delta^-1*k2^-1*Pq^2*pk*Pk - 1/2*Aplus*Bminus*Delta^-1*k2^-1*pq*Pq*Pk^2
       + 1/2*Aplus*Bminus*Delta^-1*k2^-1*pP*Pq*qk*Pk - 1/4*Aplus*Bminus*
      Delta^-1*pP*Pq^2 - 1/2*Aplus*Bminus*Delta^-1*P2*k2^-1*Pq*pk*qk + 1/4*
      Aplus*Bminus*Delta^-1*P2*pq*Pq
end

function HH(p2, q2, P2, pP, pq, Pq, A, B)
    Aplus = A(q2 + Pq + P2/4); Aminus = A(q2 - Pq + P2/4); Bplus = B(q2 + Pq + 
P2/4); Bminus = B(q2 - Pq + P2/4); k2 = p2 + q2 - 2*pq; pk = p2 - pq; qk = pq -
 q2; Pk = pP - Pq; Delta = pP^2 - p2*P2
    return  - 2*Bplus*Bminus*Delta^-1*k2^-1*Pq*pk*Pk + 2*Bplus*Bminus*Delta^-1*
      k2^-1*pq*Pk^2 - 2*Bplus*Bminus*Delta^-1*k2^-1*pP*qk*Pk + Bplus*Bminus*
      Delta^-1*pP*Pq + 2*Bplus*Bminus*Delta^-1*P2*k2^-1*pk*qk - Bplus*Bminus*
      Delta^-1*P2*pq + 2*Aplus*Aminus*Delta^-1*q2*k2^-1*Pq*pk*Pk - 2*Aplus*
      Aminus*Delta^-1*q2*k2^-1*pq*Pk^2 + 2*Aplus*Aminus*Delta^-1*q2*k2^-1*pP*
      qk*Pk - Aplus*Aminus*Delta^-1*q2*pP*Pq - 1/2*Aplus*Aminus*Delta^-1*P2*
      k2^-1*Pq*pk*Pk + 1/2*Aplus*Aminus*Delta^-1*P2*k2^-1*pq*Pk^2 - 1/2*Aplus*
      Aminus*Delta^-1*P2*k2^-1*pP*qk*Pk + 1/4*Aplus*Aminus*Delta^-1*P2*pP*Pq
       - 2*Aplus*Aminus*Delta^-1*P2*q2*k2^-1*pk*qk + Aplus*Aminus*Delta^-1*P2*
      q2*pq + 1/2*Aplus*Aminus*Delta^-1*P2^2*k2^-1*pk*qk - 1/4*Aplus*Aminus*
      Delta^-1*P2^2*pq
end


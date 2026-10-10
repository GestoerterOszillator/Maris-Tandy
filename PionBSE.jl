using Plots
using LinearAlgebra
using FastGaussQuadrature
using QuadGK
using LaTeXStrings
using ProgressMeter
using Roots

const address = "C:\\Users\\johnf\\Repositories\\Maris-Tandy\\thesis\\images\\"
# const address = "/Users/johnreeg/Documents/Repositories/Maris-Tandy/thesis/images/"
# const address = "/home/john-reeg/Documents/Maris-Tandy/thesis/images/"

const m = 0.0037
const Lambda2 = 1e6
const epsilon2 = 1e-4
const Lambda_PV = 2e2
const Lambda_QCD = 0.234
const gamma_m = 0.48 # 12/(33 - 2 * Nf{4}) = 0.48
const mu = 19.0

function Teil_Eins(w::Float64, D::Float64, PV::Bool; radial_steps::Int = 512, angular_steps::Int = 64)
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

@time t, A, B, Z_2, Z_4m = Teil_Eins(0.4, 0.93312, false)

plot(exp.(t), real.(A.(exp.(t))), xaxis=:log10, xlims = (epsilon2, Lambda2), ylims = (0, 2.0), 
    yticks = 0.4:0.4:2.0)
plot!(exp.(t), real.(B.(exp.(t))))
savefig(address * "realp.pdf")

imagp = -10:0.005:-0.001
plot(imagp, real.(B.(imagp) ./ A.(imagp)), xlims = (-10, 0), ylims = (0, 10.0))
savefig(address * "imagp.pdf")

function entries(M, p2, q2, z_p, z_q; w = 0.4, D = 0.93312, y_steps = 128) # Es fehlen w's und E(q2, Pq)
    Pq = im*M*sqrt(q2)*z_q
    qplus2 = q2 + Pq - M^2/4
    qminus2 = q2 - Pq - M^2/4
    Aplus = A(qplus2)
    Aminus = A(qminus2)
    Bplus = B(qplus2)
    Bminus = B(qminus2)
    alpha_UV(k2) = 2pi * gamma_m * (1 - exp(-k2)) / (k2 * log(exp(2)-1 + (1 + k2/Lambda_QCD^2)^2))
    alpha_IR(k2) = D/w^6 * pi * k2 * exp(-k2 / w^2)
    alpha(k2) = alpha_IR(k2) + alpha_UV(k2)
    k2_val(y) = p2 + q2 - sqrt(2*p2*q2)*(y*sqrt(1-z_q^2) + z_q)
    y, w_y = gausslegendre(y_steps)
    y_integral = sum(@. w_y * alpha(k2_val(y)))
    return 1/pi^2 * Z_2^2 * q2^2 * sqrt(1-z_q^2) * (Aplus*Aminus*(q2 + M^2/4) + Bplus*Bminus)/((qplus2*Aplus^2 + Bplus^2)*(qminus2*Aminus^2 + Bminus^2)) * y_integral
end

function build_kappa(M, t, z, w_t, w_z)
    radial_steps, angular_steps = length(t), length(z)
    N = radial_steps * angular_steps
    q2 = exp.(t)                       # exp nur einmal pro Stützstelle
    Kappa = Matrix{Float64}(undef, N, N)   # alles wird überschrieben, zeros unnötig
    progress = Progress(radial_steps)

    Threads.@threads for k in 1:radial_steps
        for l in 1:angular_steps
            col = (k - 1) * angular_steps + l
            w = w_t[k] * w_z[l]
            for i in 1:radial_steps
                # entries hängt nicht von z[j] ab -> z[1] ist nur Platzhalter
                val = w * entries(M, q2[i], q2[k], z[1], z[l])
                rows = (i - 1) * angular_steps + 1 : i * angular_steps
                @views Kappa[rows, col] .= val   # alle j auf einmal, zusammenhängend im Speicher
            end
        end
        next!(progress)
    end
    return Kappa
end

function lambda(M; radial_steps::Int = 64, angular_steps::Int = 16)
    x, w_x = gausslegendre(radial_steps)
    z, w_z = gausslegendre(angular_steps)

    t = 0.5 * (log(Lambda2) - log(epsilon2)) * x .+ 0.5 * (log(Lambda2) + log(epsilon2))
    w_t = 0.5 * (log(Lambda2) - log(epsilon2)) * w_x
    
    Eigens = eigen(build_kappa(M, t, z, w_t, w_z))
    return real(Eigens.values[end]), real.(Eigens.vectors[:, end])
end

find_zero(M -> lambda(M)[1] - 1, 0.122)
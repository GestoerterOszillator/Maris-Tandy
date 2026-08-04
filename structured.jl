using Plots
using FastGaussQuadrature
using LaTeXStrings
using ProgressMeter
pgfplotsx() # pgfplotsx() or gr()
default(
    fontfamily = "Computer Modern",
    tickfont   = font(15, "Computer Modern"),
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
const address = "/Users/johnreeg/Documents/Repositories/Maris-Tandy/thesis/images/"
# const address = "/home/john-reeg/Documents/Maris-Tandy/thesis/images/"

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
    k2(p::Float64, q::Float64, z::Float64) = p^2 + q^2 - 2*p*q*z

    # Gluon stuff
    alpha_UV(k2::Float64) = 2pi * gamma_m * (1 - exp(-k2)) / (k2 * log(exp(2)-1 + (1 + k2/Lambda_QCD^2)^2))
    alpha_IR(k2::Float64) = D/w^6 * pi * k2 * exp(-k2 / w^2)
    alpha(k2::Float64) = alpha_IR(k2) + alpha_UV(k2)

    # Angular Integrals over z = cos(psi)
    function z_intAA(p::Float64, q::Float64)
        if PV == true
            integrand = @. w_z * sqrt(1-z^2) * (p*q*z + 2*(p^2 - p*q*z)*(p*q*z - q^2)/k2(p, q, z)) * alpha(k2(p, q, z)) * Lambda_PV^2 / (k2(p, q, z) + Lambda_PV^2)
        else
            integrand = @. w_z * sqrt(1-z^2) * (p*q*z + 2*(p^2 - p*q*z)*(p*q*z - q^2)/k2(p, q, z)) * alpha(k2(p, q, z))
        end
        return sum(integrand)
    end

    zA_cache = Dict{Tuple{Float64,Float64}, Float64}()

    function z_intA(p::Float64, q::Float64)
        key = (p, q)

        if haskey(zA_cache, key)
            return zA_cache[key]
        end

        result = z_intAA(p, q)
        zA_cache[key] = result

        return result
    end

    function z_intBB(p::Float64, q::Float64)
        if PV == true
            integrand = @. w_z * sqrt(1-z^2) * alpha(k2(p, q, z)) * Lambda_PV^2 / (k2(p, q, z) + Lambda_PV^2)
        else
            integrand = @. w_z * sqrt(1-z^2) * alpha(k2(p, q, z)) 
        end
        return sum(integrand)
    end

    zB_cache = Dict{Tuple{Float64,Float64}, Float64}()

    function z_intB(p::Float64, q::Float64)
        key = (p, q)

        if haskey(zB_cache, key)
            return zB_cache[key]
        end

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
    return t, A, B, Z_2, Z_4m
end

@time t, A, B, _, _ = Teil_Eins(0.4, 1.0, true; angular_steps = 64)

# Plots
realp = plot(exp.(t), A, xaxis=:log10, xlims = (epsilon2, Lambda2), ylims = (0, 2.0), 
    yticks = 0.4:0.4:2.0, label = L"A(p^2)", xlabel = L"$p^2 \; \left[\mathrm{GeV^2}\right]$")
realp = plot!(exp.(t), B, label = L"$B(p^2) \; \left[\mathrm{GeV}\right]$")
realp = plot!(exp.(t), B./A, label = L"$M(p^2) \; \left[\mathrm{GeV}\right]$")
savefig(realp, address * "realp.pdf")

function Teil_Zwei(t, A, B, Z_2, Z_4m; w = 0.16, D = 0.93)
    k2(p::ComplexF64, q::Float64, z::Float64) = p^2 + q^2 - 2*p*q*z

    alpha_UV(k2) = 2pi * gamma_m * (1 - exp(-k2)) / (k2 * log(exp(2)-1 + (1 + k2/Lambda_QCD^2)^2))
    alpha_IR(k2) = D/w^6 * pi * k2 * exp(-k2 / w^2)
    alpha(k2) = alpha_IR(k2) + alpha_UV(k2)

    z, w_z = gausslegendre(16)
    
    function z_intA(p::ComplexF64, q::Float64)
        K = k2.(p, q, z)
        integrand = @. w_z * sqrt(1-z^2) * (p*q*z + 2*(p^2 - p*q*z)*(p*q*z - q^2)/K) * alpha(K)
        return sum(integrand)
    end

    function z_intB(p::ComplexF64, q::Float64)
        K = k2.(p, q, z)
        integrand = @. w_z * sqrt(1-z^2) * alpha(K)
        return sum(integrand)
    end

    _, w_x = gausslegendre(length(t))
    w_t = 0.5 * (log(Lambda2) - log(epsilon2)) * w_x

    function Sigma_A(p2)
        integrand = @. w_t * exp(2*t) * A / (exp(t) * A^2 + B^2) * z_intA(sqrt(p2), exp(t/2))
        return Z_2^2 * 16pi/(3*(2pi)^3 * p2) * sum(integrand)
    end

    function Sigma_B(p2)
        integrand = @. w_t * exp(2*t) * B / (exp(t) * A^2 + B^2) * z_intB(sqrt(p2), exp(t/2))
        return Z_2^2 * 16pi/(2pi)^3 * sum(integrand)
    end

    p2 = Complex.(-10 : 0.04 : -0.01)
    AA = Z_2 .+ Sigma_A.(p2)
    BB = Z_4m .+ Sigma_B.(p2)
    for (AAA, BBB) in zip(AA, BB)
        if abs(imag(AAA)) > abs(real(AAA)*1e-3)
            throw(error("Schlimm bei A: ", AAA))
        end
        if abs(imag(BBB)) > abs(real(BBB)*1e-3)
            throw(error("Schlimm bei B: ", BBB))
        end
    end
    MM = abs.(BB./AA)
    mask = isfinite.(MM) .& (abs.(MM) .< 1e2)
    return plot(real.(p2)[mask], MM[mask], xlims = (-10, 0), ylims = (0, 10), label = L"M(p^2)", xlabel = L"p^2")
end

imagp = Teil_Zwei(Teil_Eins(0.16, 0.93, false; angular_steps = 16)...)
savefig(imagp, address * "imagp.pdf")
imagp
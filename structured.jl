using LinearAlgebra
using Plots
using QuadGK
using FastGaussQuadrature
using LaTeXStrings
using ProgressMeter
using Dierckx
using Interpolations
using LsqFit
using SpecialFunctions
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
# const address = "/Users/johnreeg/Documents/Repositories/Maris-Tandy/images/"
const address = "/home/john-reeg/Documents/Maris-Tandy/images/";

# Constants
const m = 0.0037
const Lambda2 = 1e6
const epsilon2 = 1e-4
const Lambda_PV = 2e2
const Lambda_QCD = 0.234
const gamma_m = 0.48 # 12/(33 - 2 * Nf{4}) = 0.48
const mu = 19.0

# Dot Products
pk(p, q, z) = p^2 - p*q*z
qk(p, q, z) = p*q*z - q^2
k2(p, q, z) = p^2 + q^2 - 2*p*q*z

function Teil_Eins(w::Float64, D::Float64, PV::Bool; radial_steps::Int = 100, angular_steps::Int = 20)
    # Integration Grids
    x, w_x = gausslegendre(radial_steps)
    z, w_z = gausschebyshevu(angular_steps) # Chebyshev second kind

    t = 0.5 * (log(Lambda2) - log(epsilon2)) * x .+ 0.5 * (log(Lambda2) + log(epsilon2))
    w_t = 0.5 * (log(Lambda2) - log(epsilon2)) * w_x

    # Initial Values
    A = 1.6 * ones(length(t))
    B = (m + 1e-3) * ones(length(t))
    Z_2 = 1
    Z_4m = 0.7 * m

    # Gluon stuff
    alpha_UV(k2) = 2pi * gamma_m * (1 - exp(-k2)) / log(exp(2) - 1 + (1 + k2 / Lambda_QCD^2)^2)
    alpha_IR(k2) = D/w^6 * pi * k2^2 * exp(-k2 / w^2)
    alpha(k2) = alpha_IR(k2) + alpha_UV(k2)

    # Angular Integrals over z = cos(psi)
    function z_intA(p, q)
        K2 = k2.(p, q, z)
        if PV == true
            integrand = @. w_z * ( p*q*z + 2*pk(p, q, z)*qk(p, q, z)/K2 ) * alpha(K2) / K2 * Lambda_PV^2 / (K2 + Lambda_PV^2)
        else
            integrand = @. w_z * ( p*q*z + 2*pk(p, q, z)*qk(p, q, z)/K2 ) * alpha(K2) / K2
        end
        return sum(integrand)
    end

    function z_intB(p, q)
        K2 = k2.(p, q, z)
        if PV == true
            integrand = @. w_z * alpha(K2) / K2 * Lambda_PV^2 / (K2 + Lambda_PV^2)
        else
            integrand = @. w_z * alpha(K2) / K2
        end
        return sum(integrand)
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
    max_iter = 200
    p = Progress(max_iter, desc = "Berechne...")

    for i in 1:max_iter
        A, B, Z_2, Z_4m, max_error = update(A, B, Z_2, Z_4m)
        next!(p)
        if max_error < 1e-8
            println("\nFinished: i = ", i)
            break
        end
        if i == max_iter
            throw(error("Failed: max_error = ", max_error))
        end
    end

    return t, A, B, Z_2, Z_4m
end

t, A, B, Z_2, Z_4m = Teil_Eins(0.16, 0.93, false)

# Plots
# realp = plot(exp.(t), A, xaxis=:log10, xlims = (epsilon2, Lambda2), ylims = (0, 2), 
#     yticks = 0.4:0.4:2.0, label = L"A(p^2)", xlabel = L"p^2")
# realp = plot!(exp.(t), B, label = L"B(p^2)")
# realp = plot!(exp.(t), B./A, label = L"M(p^2)")
# savefig(realp, address * "D1.pdf")

function Teil_Zwei(t, A, B)
    # fitter(x, p) = @. p[1] .+ (p[2] - p[1]) ./ (1 .+ exp.(-p[3] .* (x .- p[4]))).^(1 / p[5]) - p[6] * x
    # par_A = curve_fit(fitter, t, A, [2.24, 1.,     1.27, -0.01, 1.1,  0.001  ]).param
    # par_B = curve_fit(fitter, t, B, [1.28, 0.0035, 2.5,  -0.58, 2.54, 0.00026]).param

    # println("max. abs. err. A: ", maximum(abs.(A - fitter(t, par_A))))
    # println("max. rel. err. A: ", maximum(abs.(A - fitter(t, par_A)) ./ A))
    # println("max. abs. err. B: ", maximum(abs.(B - fitter(t, par_B))))
    # println("max. rel. err. B: ", maximum(abs.(B - fitter(t, par_B)) ./ B))

    # A_itp(x) = @. par_A[1] .+ (par_A[2] - par_A[1]) ./ (1 .+ exp.(-par_A[3] .* (x .- par_A[4]))).^(1 / par_A[5]) - par_A[6] * x
    # B_itp(x) = @. par_B[1] .+ (par_B[2] - par_B[1]) ./ (1 .+ exp.(-par_B[3] .* (x .- par_B[4]))).^(1 / par_B[5]) - par_B[6] * x
    return Spline1D(t, A, k = 3), Spline1D(t, B, k = 3)
end

function Teil_Drei(A_itp, B_itp, Z_2, Z_4m; D = 0.93, w = 0.16)
    # Gluon stuff
    alpha_UV(k2) = 2pi * gamma_m * (1 - exp(-k2)) / log(exp(2) - 1 + (1 + k2 / Lambda_QCD^2)^2)
    alpha_IR(k2) = D/w^6 * pi * k2^2 * exp(-k2 / w^2)
    alpha(k2) = alpha_IR(k2) + alpha_UV(k2)

    function z_intA(p, q; delta = 0.01)
        integrand(z) = sqrt(1 - z^2) * (p*q*z + 2*pk(p,q,z)*qk(p,q,z)/k2(p,q,z)) *
                    alpha(k2(p,q,z)) / k2(p,q,z)

        z0 = (p^2 + q^2) / (2p*q)   # exact pole location, k2(p,q,z0) = 0

        if abs(imag(z0)) < 0.5 && abs(real(z0)) < 1.5
            # pole is close to the real contour -> route around it
            zc = real(z0) + im*delta*sign(imag(z0) == 0 ? 1.0 : -imag(z0))
            # clamp the real part into a sane range so the waypoint stays near the path
            zc = clamp(real(zc), -1.0, 1.0) + im*imag(zc)
            return quadgk(integrand, -1.0, zc, 1.0)[1]
        else
            return quadgk(integrand, -1.0, 1.0)[1]
        end
    end

    z, w_z = gausslegendre(20)
    z1  = 0.5 * z .- 0.5
    z2  = 0.5 * z .+ 0.5
    w_z = 0.5 * w_z;

    function z_intAA(p, q)
        K1 = k2.(p, q, z1)
        K2 = k2.(p, q, z2)
        integrand1 = @. w_z * sqrt(1-z1^2) * (p*q*z1 + 2*(p^2 - p*q*z1)*(p*q*z1 - q^2)/(K1 + 1e-3im)) * alpha(K1) / K1
        integrand2 = @. w_z * sqrt(1-z2^2) * (p*q*z2 + 2*(p^2 - p*q*z2)*(p*q*z2 - q^2)/(K2 + 1e-3im)) * alpha(K2) / K2
        return sum(integrand1) + sum(integrand2)
    end

    function z_intB(p, q)
        K1 = k2.(p, q, z1)
        K2 = k2.(p, q, z2)
        integrand1 = @. w_z * sqrt(1-z1^2) * alpha_UV(K1) / K1
        integrand2 = @. w_z * sqrt(1-z2^2) * alpha_UV(K2) / K2
        return sum(integrand1) + sum(integrand2) + D * pi^2 / w^4 * exp(-(p^2 + q^2)/w^2) * ( (p^2 + q^2 + 2*w^2)/(2p*q) * besseli(1, 2p*q/w^2) - besseli(0, 2p*q/w^2) )
    end

    function Sigma_A(p2)
        integrand(t) = exp(2*t) * A_itp(t) / (exp(t) * A_itp(t)^2 + B_itp(t)^2) * z_intA(sqrt(p2), exp(t/2))
        return Z_2^2 * 16pi/(3*(2pi)^3 * p2) * quadgk(integrand, log(epsilon2), log(abs(p2)), log(Lambda2), maxevals=10^3)[1]
    end

    function Sigma_B(p2)
        integrand(t) = exp(2*t) * B_itp(t) / (exp(t) * A_itp(t)^2 + B_itp(t)^2) * z_intB(sqrt(p2), exp(t/2))
        return Z_2^2 * 16pi/(2pi)^3 * quadgk(integrand, log(epsilon2), log(abs(p2)), log(Lambda2))[1]
    end

    p2 = Complex.(-10 : 0.1 : -0.001)
    AA = Z_2 .+ Sigma_A.(p2)
    BB = Z_4m .+ Sigma_B.(p2)
    MM = real.(BB./AA)
    mask = isfinite.(MM) .& (abs.(MM) .< 1e4)
    return plot(real.(p2)[mask], MM[mask], xlims = (-10, 0), ylims = (0, 600), label = L"M(p^2)", xlabel = L"p^2")
    # return (Z_4m .+ Sigma_B.(Complex(-8)))/(Z_2 + Sigma_A(Complex(-8)))
end

@time Teil_Drei(Teil_Zwei(t, A, B)..., Z_2, Z_4m)
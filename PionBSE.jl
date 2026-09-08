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

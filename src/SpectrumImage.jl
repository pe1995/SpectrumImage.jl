module SpectrumImage

using Reexport
using PythonPlot
using DelimitedFiles
using PythonCall
using ArgParse
using Images
using Plots
using Interpolations
using StatsBase
using KernelDensity
using LsqFit
using NPZ
@reexport using LaTeXStrings

export spectrum
export read_spectrum, read_line_indicators
export spectrum_gif
export normalize, continuum

# constants
const KBoltzmann = 1.380658E-16                # Boltzman's cst. [erg/K]
const CLight     = 2.99792458E+10              # Speed of light [cm/s]
const HPlanck    = 6.6260755E-27               # Planck's constant [erg s]
const ev_to_erg  = 1.60218e-12                 # conversion
const twohc2     = 2.0e0 *HPlanck*CLight^2
const hc_k       = HPlanck*CLight/KBoltzmann
const aa_to_cm   = 1.0e-8
const σ_S        = 5.6704e-5

#= Python modules =#
const numpy = PythonCall.pynew()

__init__() = begin 
    PythonCall.pycopy!(numpy, pyimport("numpy"))
end

# Write your package code here.
include("_spectrum_image.jl")

end

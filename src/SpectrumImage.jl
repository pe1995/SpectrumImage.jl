module SpectrumImage

using PythonPlot
using DelimitedFiles
using PythonCall
using LaTeXStrings

plt = matplotlib.pyplot
matplotlib.style.use("dark_background")

export spectrum

# Write your package code here.
include("_spectrum_image.jl")

end

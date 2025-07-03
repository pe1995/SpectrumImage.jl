using Pkg; Pkg.activate(".")
using SpectrumImage

SpectrumImage.ArgParse.parse_item(::Type{Vector{A}}, x::AbstractString) where {A<:AbstractString} = split(x, ",", keepempty=false)
SpectrumImage.ArgParse.parse_item(::Type{Vector{A}}, x::AbstractString) where {A<:AbstractFloat} = parse.(Float64, split(x, ",", keepempty=false))

# available command line arguments
s = SpectrumImage.ArgParseSettings()
SpectrumImage.@add_arg_table s begin
    "--file", "-f"
        help = "CSV file containing the spectra."
        arg_type = String
        required=true
    "--outfile", "-o"
        help = "Path to store the spectrum image."
        arg_type = String
        required=true
    "--figsize"
        help = "Figsize in inches. Pass x,y."
        arg_type = Vector{Float64}
        default=[9, 6]    
    "--dpi"
        help = "DPI of the final image. Higher resolution for e.g. png images."
        arg_type = Int
        default=-1
    "--header"
        help = "Number of lines in the CSV file that are header. Those will be skipped."
        arg_type = Int
        default = 0
    "--separator"
        help = "Separator between columns in the CSV file."
        arg_type = Char
        default = ','
    "--rows"
        help = "Number of rows in the spectrum will be split in."
        arg_type = Int
        default=30
    "--separator_width"
        help = "Width of the row separator in the final image (the black lines between rows)."
        arg_type = Float64
        default=1.5
    "--show_lambda_range"
        help = "Show the spectral range at the top of the figure"
        action = :store_true
    "--lambda_UV"
        help = "Wavelength below which every point is shown with the UV color of the colormap."
        arg_type = Float64
        default=-1.0
    "--lambda_IR"
        help = "Wavelength above which every point is shown with the IR color of the colormap."
        arg_type = Float64
        default = -1.0
    "--F_low"
        help = "Set the minimum brightness in the final image to this value. Points with this flux will be black."
        arg_type = Float64
        default=-1.0
    "--F_high"
        help = "Set the maximum brightness in the final image to this value. Points with this flux will be the brightest."
        arg_type = Float64
        default=-1.0
    "--line_indicators"
        help = "Plot white line indicators at these positions."
        arg_type = Vector{Float64}
        default=Float64[]
    "--indicator_fontsize"
        help = "Size of the line indicator labels."
        arg_type = String
        default="small"
    "--units"
        help = "Units to add behind the wavelengths. Note that the string will be passed through LaTeX."
        arg_type = String
        default= ""
    "--colormap"
        help = """
        Colormap to use. All matplotlib colormaps that are available. Please note the the spectrum will be shown from red (large lambda) to blue (small lambda). 
        If you use a spectral colormap make sure that the beginning is red and the end blue. If reverse, consider the `_r` option, e.g. `rainbow_r`.
        """
        arg_type = String
        default="gist_rainbow"
end

# read and parse command line arguments
arguments = SpectrumImage.parse_args(ARGS, s)
data = read_spectrum(arguments["file"], arguments["separator"], skipstart=arguments["header"])

if !(eltype(data) <: AbstractFloat)
    error("The given data set could not be read. Is there a header? Are columns separated by ','? Consider --separator and --header.")
end

u = arguments["units"]
f, ax = spectrum(
    data[:,1], data[:,2], 
    rows=arguments["rows"], 
    figsize=(arguments["figsize"][1], arguments["figsize"][2]), 
    λ_IR=arguments["lambda_IR"],
    λ_UV=arguments["lambda_UV"],
    F_low=arguments["F_low"],
    F_high=arguments["F_high"],
    separator_width=arguments["separator_width"],
    indicator_fontsize=arguments["indicator_fontsize"],
    show_lambda_range=arguments["show_lambda_range"], 
    line_indicators=arguments["line_indicators"], 
    units=length(u) == 0 ? u : latexstring(u),
    colormap=arguments["colormap"]
)

if arguments["dpi"] > 0
    f.savefig(arguments["outfile"], dpi=arguments["dpi"])
else
    f.savefig(arguments["outfile"])
end
@info "Image created and stored at $(arguments["outfile"])."
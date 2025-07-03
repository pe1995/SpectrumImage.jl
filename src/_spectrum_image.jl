"""
    spectrum(λ, F; colormap="gist_rainbow", rows=30, separator_width=1.5, show_lambda_range=false, λ_UV=nothing, λ_IR=nothing, F_low=minimum(F), F_high=maximum(F), line_indicators=[], indicator_fontsize="small", units="", kwargs...)

Create a 2D spectrum image from wavelength and Flux arrays. 
Wavelength will be used as the indicator for color, and should be chosen from `red` to `blue`.
`λ` and `F` will be sorted automatically. More `rows` will increase the rows of the image. 

Kwargs will be passed to the plt.figure constructor. `plt.figure` and `ax` will be returned and can be used to add anything you want.

By default, the entire colormap will be used from the red end of `λ` to the blue end. You can optically specify where the red part should stop (`λ_IR`) and
where the violet part should be begin (`λ_UV`). Space before and after will be filled by the respective color.

You can specify `line_indicators` in wavelength. At the corresponding positions there will be white line indicators shown in the final image.

The `units` string will be pasted directly behind the line indicators, if wanted. `show_lambda_range` will add the lambda range to the top of the image.
You can also specify `F_low` and `F_high` to adjust the maximum and minimum for the normalization.

# Example:
```julia
# read a spectrum from file with `,` separator, uses `readdlm`. Skip e.g. the first line in this example.
data = read_spectrum("my_spectrum.csv", ',', skipstart=1)

# create the spectrum image
f, ax = spectrum(data[:, 1], data[:, 2]; rows=30, figsize=(9, 6), show_lambda_range=true, λ_IR=5500, line_indicators=[5500, 5400]);
```
"""
function spectrum(λ, F; colormap="gist_rainbow", rows=30, separator_width=1.5, show_lambda_range=false, λ_UV=-1, λ_IR=-1, F_low=-1, F_high=-1, line_indicators=[], indicator_fontsize="small", units="", kwargs...)
    plt = matplotlib.pyplot
    matplotlib.style.use("dark_background")

    λ_sort = sortperm(λ, rev=true)

    λ = λ[λ_sort]
    F = F[λ_sort]

    notgiven(a) = a < 0

    min_l, max_l = notgiven(λ_IR) ? first(λ) : λ_IR, notgiven(λ_UV) ? last(λ) : λ_UV
    λ_norm = abs.((λ .- min_l) ./ (max_l - min_l))
    λ_norm[λ .> min_l] .= 0.0
    λ_norm[λ .< max_l] .= 1.0

    min_F, max_F = notgiven(F_low) ? minimum(F) : F_low, notgiven(F_high) ? maximum(F) : F_high
    F_norm = (F .- min_F) ./ (max_F - min_F)
    
    columns = floor(Int, length(λ) / rows)
    image_matrix = zeros(rows, columns, 3)

    cmap = plt.get_cmap(colormap)
    colors = pyconvert.(Array, cmap(λ_norm))

    line_indicator_index = Dict()
    for line in line_indicators
        if (line > first(λ)) | (line < last(λ))
            @warn "Line $(line) out or range for the given spectrum ."
        else
            line_indicator_index[argmin(abs.(λ .- line))] = line
        end
    end


    c = 1
    for i in axes(image_matrix, 1)
        for j in axes(image_matrix, 2)
            image_matrix[i, j, :] = colors[c, 1:3] * F_norm[c]
            c += 1
        end
    end

    f, ax = plt.subplots(1, 1; layout="tight", kwargs...)
    ax.imshow(image_matrix, aspect="auto", origin="upper", interpolation="none")
    
    c = 1
    jsep = size(image_matrix, 2) / 200
    for i in axes(image_matrix, 1)
        for j in axes(image_matrix, 2)
            if c in keys(line_indicator_index)
                ax.vlines(j, i -1 -0.5, i - 1 +0.5, color="w", lw=separator_width)
                ax.text(j+jsep, i -1, "$(line_indicator_index[c])"*units, color="w", ha="left", va="center", fontsize=indicator_fontsize)
            end
            c += 1
        end
    end

    for i in 0:rows
        x = range(0, size(image_matrix,2)-1, length=200)|> collect
        y = fill!(similar(x), i-0.5)
        ax.plot(x, y, ls="-", color="k", lw=separator_width)
    end
    

    if show_lambda_range
        ax.text(0.0, 1.01, L"\rm \lambda: "*"$(round(maximum(λ), sigdigits=5))"*units*L"\rm\, - \,"*"$(round(minimum(λ), sigdigits=5))"*units, transform=ax.transAxes)
    end

    ax.axis("off")
    f, ax
end


read_spectrum = readdlm
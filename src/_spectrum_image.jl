"""
    spectrum(λ, F; colormap="gist_rainbow", rows=50, separator_width=2, show_lambda_range=false, λ_UV=nothing, λ_IR=nothing, kwargs...)

Create a 2D spectrum image from wavelength and Flux arrays. 
Wavelength will be used as the indicator for color, and should be chosen from `red` to `blue`.
`λ` and `F` will be sorted automatically. More `rows` will increase the rows of the image. 

Kwargs will be passed to the plt.figure constructor. `plt.figure` and `ax` will be returned and can be used to add anything you want.

By default, the entire colormap will be used from the red end of `λ` to the blue end. You can optically specify where the red part should stop (`λ_IR`) and
where the violet part should be begin (`λ_UV`). Space before and after will be filled by the respective color.
"""
function spectrum(λ, F; colormap="gist_rainbow", rows=50, separator_width=2, show_lambda_range=false, λ_UV=nothing, λ_IR=nothing, kwargs...)
    plt = matplotlib.pyplot
    matplotlib.style.use("dark_background")

    λ_sort = sortperm(λ, rev=true)

    λ = λ[λ_sort]
    F = F[λ_sort]

    min_l, max_l = isnothing(λ_IR) ? first(λ) : λ_IR, isnothing(λ_UV) ? last(λ) : λ_UV
    λ_norm = abs.((λ .- min_l) ./ (max_l - min_l))
    λ_norm[λ .> min_l] .= 0.0
    λ_norm[λ .< max_l] .= 1.0

    min_F, max_F = minimum(F), maximum(F)
    F_norm = (F .- min_F) ./ (max_F - min_F)
    
    columns = floor(Int, length(λ) / rows)
    image_matrix = zeros(rows, columns, 3)

    cmap = plt.get_cmap(colormap)
    colors = pyconvert.(Array, cmap(λ_norm))

    c = 1
    for i in axes(image_matrix, 1)
        for j in axes(image_matrix, 2)
            image_matrix[i, j, :] = colors[c, 1:3] * F_norm[c]
            c += 1
        end
    end

    f, ax = plt.subplots(1, 1; layout="tight", kwargs...)
    ax.imshow(image_matrix, aspect="auto", origin="upper", interpolation="none")
    
    for i in 0:rows
        x = range(0, size(image_matrix,2)-1, length=200)|> collect
        y = fill!(similar(x), i-0.5)
        ax.plot(x, y, ls="-", color="k", lw=separator_width)
    end

    if show_lambda_range
        ax.text(0.0, 1.01, L"\rm \lambda: "*"$(round(maximum(λ), sigdigits=5))"*L"\rm -"*"$(round(minimum(λ), sigdigits=5))", transform=ax.transAxes)
    end

    ax.axis("off")
    f, ax
end
"""
    spectrum(λ, F; colormap="gist_rainbow", rows=30, figsize=(9,6), separator_width=1.5, show_lambda_range=false, λ_UV=nothing, λ_IR=nothing, F_low=minimum(F), F_high=maximum(F), line_indicators=[], indicator_fontsize="small", units="", kwargs...)

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
function spectrum(λ, F; colormap="gist_rainbow", figsize=(9,6), rows=30, separator_width=1.5, color_spacing="index", unify_spacing=false, show_lambda_range=true, λ_UV=-1, λ_IR=-1, F_low=-1, F_high=-1, line_indicators=[], windows=nothing, window_gap_fraction=0.05, window_gap_color=[0,0,0], indicator_fontsize="small", units=L"\rm \,\AA", header=nothing, kwargs...)
    plt = matplotlib.pyplot
    matplotlib.style.use("dark_background")

    λ_sort = sortperm(λ, rev=true)
    λ = λ[λ_sort]
    F = F[λ_sort]

    ll = length(λ)

    notgiven(a) = a < 0

    windows = if isnothing(windows)
        [eachindex(λ)|>collect]
    elseif windows == "4MOST"
        @info "Cutting to 4MOST windows 3926Å - 4355Å, 5160Å - 5730Å, and 6100Å - 6790Å."
        [findall(3926 .<= λ .<= 4355), findall(5160 .<= λ .<= 5730), findall(6100 .<= λ .<= 6790)]
    else
        [findall(wmin .<= λ .<= wmax) for (wmin, wmax) in windows]
    end

    red_limits = λ[first.(windows)]
    blue_limits = λ[last.(windows)]
    window_sort = sortperm(red_limits, rev=true)
    windows = windows[window_sort]
    red_limits = red_limits[window_sort]
    blue_limits = blue_limits[window_sort]

    if !unify_spacing
        for w in windows
            d = abs.(diff(λ[w]))
            if !(all(d .≈ d[1]))
                @warn "Your input spectrum is not uniform in λ. This means that each pixel in the final image may correspond to a different wavelength spacing! Consider `unify_spacing=true`."
                @warn "min:$(minimum(d)), max:$(maximum(d)),  mean:$(sum(d)/length(d))"
            end
        end
    end

    # now we construct a new wavelength araay based on these windows
    λ, F = if !unify_spacing
        vcat([λ[w] for w in windows]...), vcat([F[w] for w in windows]...)
    else
        @info "Linear interpolating fluxes to uniform λ grid for each window..."
        lam_all, f_all = [], []
        for w in windows
            # sort increasing
            λ_sort = sortperm(λ[w])
            l = λ[w][λ_sort]
            fl = F[w][λ_sort]

            # interpolate
            l_new = range(first(l), last(l), length=length(l)) |> collect
            f_new = linear_interpolation(Interpolations.deduplicate_knots!(l), fl).(l_new)
            λ_sort = sortperm(l_new, rev=true)

            append!(lam_all, [l_new[λ_sort]])
            append!(f_all, [f_new[λ_sort]])
        end

        vcat(lam_all...), vcat(f_all...)
    end

    min_l, max_l = notgiven(λ_IR) ? first(λ) : λ_IR, notgiven(λ_UV) ? last(λ) : λ_UV
    λ_norm = if color_spacing == "index" 
        @info "Colors computed based on index in λ array."
        @info "Note: If λ spacing is not uniform, this will result in a stretching of spectral features in color, as every pixel gets a new color."
        abs.((collect(eachindex(λ)) .- 1.0) ./ (length(λ) - 1.0)) 
    else
        @info "Colors computed based on wavelength values."
        @info "Note that this means that points closer together will have a more similar color. This also means that gaps in the spectrum will cause gaps in the color."
        # this scales wavelength based on the numbers --> non-continous when there gaps
        # however, this will respect non-uniform spacing
        abs.((λ .- min_l) ./ (max_l - min_l))
    end
    
    λ_norm[λ .> min_l] .= 0.0
    λ_norm[λ .< max_l] .= 1.0

    min_F, max_F = notgiven(F_low) ? minimum(F) : F_low, notgiven(F_high) ? maximum(F) : F_high
    F_norm = (F .- min_F) ./ (max_F - min_F)

    cmap = plt.get_cmap(colormap)
    colors = pyconvert.(Array, cmap(λ_norm))
    colors = [colors[i, 1:3] for i in axes(colors, 1)]
    
    # now we add black limiters in between the windows
    if length(windows) > 1
        length_sep = floor(Int, length(λ_norm) / rows * window_gap_fraction)
        for i in 1:length(windows)-1
            w = windows[i]
            insert_at_index = findfirst(x -> blue_limits[i] > x, λ) -1
            for j in 1:length_sep
                insert!(λ, insert_at_index, blue_limits[i])
                insert!(F, insert_at_index, 0)
                insert!(F_norm, insert_at_index, 1)
                insert!(λ_norm, insert_at_index, λ_norm[insert_at_index])
                insert!(colors, insert_at_index, floor.(Int, window_gap_color))
            end
        end
    end

    columns = floor(Int, length(λ_norm) / rows)
    image_matrix = zeros(rows, columns, 3)

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
            image_matrix[i, j, :] = colors[c] * F_norm[c]
            c += 1
        end
    end

    f, ax = plt.subplots(1, 1; layout="tight", figsize=figsize)
    ax.imshow(image_matrix, aspect="auto", origin="upper", interpolation="none", rasterized=true)
    
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
        ax.text(0.0, 1.01, L"\rm \lambda: "*"$(round(minimum(λ), sigdigits=5))"*units*L"\rm\, - \,"*"$(round(maximum(λ), sigdigits=5))"*units, transform=ax.transAxes, ha="left", va="bottom")
    end

    if !isnothing(header)
        ax.text(1.0, 1.01, header, transform=ax.transAxes, ha="right", va="bottom", color="w")
    end

    ax.axis("off")
    f, ax
end


read_spectrum = readdlm


"""
    spectrum_gif(out_path, λ, Fs; store_at=mktempdir(), fps=5, dpi=300, figsize=(9,6), kwargs...)

Create an animation from multiple fluxes and store it at `out_path`. Parameters passed as kwargs are passed on to `spectrum`.
`Fs` is assumed to be an array of size (nLambda,nFluxes).
"""
function spectrum_gif(out_path, λ, Fs; store_at=mktempdir(), fps=5, dpi=300, figsize=(9,6), header=[nothing for _ in size(Fs, 2)], kwargs...)
    @info "Building Animation..."
    
    fls = []
    for i in axes(Fs, 2)
        F = Fs[:, i]
        matplotlib.pyplot.close()
        f, ax = spectrum(λ, F; figsize=figsize, header=header[i], kwargs...)
        
        f.savefig(joinpath(store_at,"cube_$(i).png"), bbox_inches="tight", dpi=dpi)
        append!(fls, [joinpath(store_at,"cube_$(i).png")])
    end
    
    v_images = Images.load.(fls)
    anim = @animate for i ∈ eachindex(v_images)
        Plots.plot(
            v_images[i], 
            axis=([], false), 
            background_color=:black,
            size=(figsize[1]*dpi, figsize[2]*dpi)
        )
    end every 1
    g = gif(anim, out_path, fps=fps)
    rm.(fls)

    @info "...Animation built."

    g
end
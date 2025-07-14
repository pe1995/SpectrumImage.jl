"""
    spectrum(λ, F; colormap="gist_rainbow", figsize=(9,6), rows=30, separator_width=1.5, color_spacing="index", unify_spacing=true, show_lambda_range=true, λ_UV=-1, λ_IR=-1, F_low=-1, F_high=-1, line_indicators=[], windows=nothing, window_gap_fraction=0.05, window_gap_color=[0,0,0], indicator_fontsize="small", units="AA", header=nothing, kwargs...)

Create a 2D spectrum image from wavelength and Flux arrays. 
Wavelength will be used as the indicator for color, and should be chosen from `red` to `blue`.
`λ` and `F` will be sorted automatically. More `rows` will increase the rows of the image. 

Kwargs will be passed to the plt.figure constructor. `plt.figure` and `ax` will be returned and can be used to add anything you want.

By default, the entire colormap will be used from the red end of `λ` to the blue end. You can optically specify where the red part should stop (`λ_IR`) and
where the violet part should be begin (`λ_UV`). Space before and after will be filled by the respective color.

You can specify `line_indicators` in wavelength. At the corresponding positions there will be white line indicators shown in the final image. 
You can also pass a vector of pairs, where the second entry corresponds to the text you want to put in the image, e.g. [6562.8=>"Hα", 3968.5=>"Ca II - H", 3933.7=>"Ca II - K"].

The `units` string will be pasted directly behind the line indicators, if wanted. `show_lambda_range` will add the lambda range to the top of the image.
You can also specify `F_low` and `F_high` to adjust the maximum and minimum for the normalization. 
You can also specify `F_low` and `F_high` per window, i.e. for 2 windows you can do `F_low=[-1, 0.5], F_high=[1.1, -1]` to set the limits. -1 will put the limit automatically to the max and min in that window.

You can specify `windows` to cut your input spectrum accordingly and insert window separators in between (set `window_gap_color` for RGB color of the line, and `window_gap_fraction` for the line width in units of row fraction.)
Note that the colormap will be uniform across all windows if `color_spacing="index"`, and respect the actual wavelength difference if you use `color_spacing="wavelength"`. There will be jumps in color if there are jumps in windows.

If your input spectrum is not equally spaced in wavelength, this means that in the final image pixels directly next to each other are not corresponding to the same wavelength step.
If you have more points in spectral lines for example this causes the lines to be spread. If you set `unify_spacing=true` the code will interpolate the flux in each window to avoid this.

# Example:
```julia
# read a spectrum from file with `,` separator, uses `readdlm`. Skip e.g. the first line in this example.
data = read_spectrum("my_spectrum.csv", ',', skipstart=1)

# create the spectrum image
f, ax = spectrum(data[:, 1], data[:, 2]; rows=30, figsize=(9, 6), show_lambda_range=true, λ_IR=5500, line_indicators=[5500, 5400]);
```
"""
function spectrum(λ, F; colormap="gist_rainbow", 
    figsize=(9,6), rows=30, separator_width=1.5, color_spacing="index", 
    unify_spacing=true, show_lambda_range=true, 
    λ_UV=-1, λ_IR=-1, F_low=-1, F_high=-1,
    line_indicators=[], windows=nothing, window_gap_fraction=0.05, window_gap_color=[0,0,0],
    indicator_fontsize="small", units=L"\rm \,\AA", header=nothing, header_location="right", kwargs...)
    # matplotlib for plotting
    plt = matplotlib.pyplot
    matplotlib.style.use("dark_background")

    λ_sort = sortperm(λ, rev=true)
    λ = λ[λ_sort]
    F = F[λ_sort]

    notgiven(a) = a < 0 

    windows = if isnothing(windows)
        [eachindex(λ)|>collect]
    elseif windows == "4MOST"
        @info "Cutting to 4MOST windows 3926Å - 4355Å, 5160Å - 5730Å, and 6100Å - 6790Å."
        [findall(3926 .<= λ .<= 4355), findall(5160 .<= λ .<= 5730), findall(6100 .<= λ .<= 6790)]
    else
        [findall(wmin .<= λ .<= wmax) for (wmin, wmax) in windows]
    end

    F_low = (length(F_low) == 1) ? [F_low] : F_low
    F_high = (length(F_high) == 1) ? [F_high] : F_high
    F_by_window = (length(F_high) > 1) | (length(F_low) > 1)
    F_low, F_high = if F_by_window
        fl = (length(F_low) == 1) ? [F_low[1] for _ in windows] : F_low
        fh = (length(F_high) == 1) ? [F_high[1] for _ in windows] : F_high
        fl, fh
    else
        F_low[1], F_high[1]
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

    # now we construct a new wavelength array based on these windows
    λ, F = if !unify_spacing
        vcat([λ[w] for w in windows]...), vcat([F[w] for w in windows]...)
    else
        @info "Linear interpolating fluxes to uniform λ grid for each window..."
        lam_all, f_all = [], []
        for (wi, w) in enumerate(windows)
            # sort increasing
            λ_sort = sortperm(λ[w])
            l = λ[w][λ_sort]
            fl = F[w][λ_sort]

            # interpolate
            l_new = range(first(l), last(l), length=length(l)) |> collect
            f_new = linear_interpolation(Interpolations.deduplicate_knots!(l), fl).(l_new)
            λ_sort = sortperm(l_new, rev=true)

            # normalize by window if wanted
            if F_by_window
                min_F, max_F = notgiven(F_low[wi]) ? minimum(f_new) : F_low[wi], notgiven(F_high[wi]) ? maximum(f_new) : F_high[wi]
                f_new .= (f_new .- min_F) ./ (max_F - min_F)
            end

            append!(lam_all, [l_new[λ_sort]])
            append!(f_all, [f_new[λ_sort]])
        end

        vcat(lam_all...), vcat(f_all...)
    end

    min_l, max_l = notgiven(λ_IR) ? first(λ) : λ_IR, notgiven(λ_UV) ? last(λ) : λ_UV
    λ_norm = if color_spacing == "index" 
        @info "Colors computed based on index in λ array."
        if !unify_spacing
            @info "Note: If λ spacing is not uniform, this will result in a stretching of spectral features in color, as every pixel gets a new color."
        end
        iminl = argmin(abs.(λ .- min_l))
        imaxl = argmin(abs.(λ .- max_l))
        abs.((collect(eachindex(λ)) .- iminl) ./ (imaxl - iminl)) 
    elseif color_spacing == "wavelength" 
        @info "Colors computed based on wavelength values."
        @info "Note that this means that points closer together will have a more similar color. This also means that gaps in the spectrum will cause gaps in the color."
        # this scales wavelength based on the numbers --> non-continous when there gaps
        # however, this will respect non-uniform spacing
        abs.((λ .- min_l) ./ (max_l - min_l))
    else
        error("Given colorspacing not available. Please put `index` or `wavelength`.")
    end
    
    λ_norm[λ .> min_l] .= 0.0
    λ_norm[λ .< max_l] .= 1.0

    F_norm = if F_by_window
        F
    else
        min_F, max_F = notgiven(F_low) ? minimum(F) : F_low, notgiven(F_high) ? maximum(F) : F_high
        (F .- min_F) ./ (max_F - min_F)
    end

    F_norm[F_norm.<0.0] .= 0.0
    F_norm[F_norm.>1.0] .= 1.0

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
    line_indicator_text = Dict()
    for lineind in line_indicators
        line, text = if typeof(lineind) <: Pair
            first(lineind), last(lineind)
        else
            lineind, "$(lineind)"*units
        end
        if (line > first(λ)) | (line < last(λ))
            @warn "Line $(line) out or range for the given spectrum ."
        else
            line_indicator_index[argmin(abs.(λ .- line))] = line
            line_indicator_text[argmin(abs.(λ .- line))] = text
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
                if j<columns-floor(Int, columns*0.05)
                    ax.text(j+jsep, i -1, line_indicator_text[c], color="w", ha="left", va="center", fontsize=indicator_fontsize)
                else
                    ax.text(j-jsep, i -1, line_indicator_text[c], color="w", ha="right", va="center", fontsize=indicator_fontsize)
                end
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
        hx, hy, ha, va = header_location=="right" ? (1.0, 1.01, "right", "bottom") : (0.0, 1.01, "left", "bottom")
        ax.text(hx, hy, header, transform=ax.transAxes, ha=ha, va=va, color="w")
    end

    ax.axis("off")
    f, ax
end

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





#= Continuum normalization =#

function binned_maxima(x, y; nbins=max(floor(Int, 3*length(x)/1000), 10))
    counts, bins = numpy.histogram(x, weights=y, bins=nbins)
    counts = floor.(Int, SpectrumImage.pyconvert(Array, counts))
    bins = SpectrumImage.pyconvert(Array, bins)

    # compute maxima for continuum estimate
    bin_centers = (bins[1:end-1] .+ bins[2:end]) ./ 2.0
    bin_maxima = [maximum(y[bins[i-1] .< x .< bins[i]]) for i in 2:length(bins)]

    bin_centers, bin_maxima
end

"""
    continuum(x, y; bins=floor(Int, length(x)/1000), bandwidth=length(bins)/10)

Normalize the spectrum approximatelly to remore brightness gradient from image.
This procedure splits the spectrum in a number of `nbins` bins and makes a Kernel Density estimate 
of the resulting weighted wavelength distribution. A scaling factor is fitted to the maxima. The continuum is returned.
"""
continuum(x, y; kwargs...) = continuum_function(x, y; kwargs...)(x)

function continuum_function(x, y; nbins=max(floor(Int, 3*length(x)/1000), 10), bandwidth=std(x, Weights(y))/10)
    # bin and get maxima
    bin_centers, bin_maxima = binned_maxima(x, y; nbins=nbins)

    # perform kernel density estimate and fit scaling factor
    k = kde(bin_centers, weights=bin_maxima, bandwidth=bandwidth)
    kde_scaled(x, p) = pdf(k, x) .* p[1]
    r = curve_fit(kde_scaled, bin_centers, bin_maxima, [1.0])

    xi -> kde_scaled(xi, r.param)
end

function iterate_continuum(x, y; nbins=max(floor(Int, 3*length(x)/1000), 10))
    # bin and get maxima
    bin_centers, bin_maxima = binned_maxima(x, y; nbins=nbins)
    ci(xi, p) = begin
        f = continuum_function(x, y; bandwidth=abs.(p[1]), nbins=nbins)
        f(xi)
    end
    r = curve_fit(ci, bin_centers, bin_maxima, [(maximum(x) - minimum(x))/60])
    ci(x, r.param)
end

normalize(x, y; kwargs...) = y ./ continuum(x, y; kwargs...)






#= Read spectra =#

function read_line_indicators(f, args...; as_latex=true, kwargs...)
    wavelname = SpectrumImage.readdlm(f, args...; kwargs...)
    wavel, name = wavelname[:, 1], wavelname[:, 2]
    name = as_latex ? latexstring.(name) : name
    [w=>n for (w, n) in zip(wavel, name)]
end

"""
    read_spectrum(args...; format=:readdlm, kwargs...)

Read the spectrum from a file. Calls the function specified as `:func` with all arguments in args and kwargs.
Use e.g. `:readdlm` or `:npzread` for `.npy` files.
"""
read_spectrum(args...; format=:readdlm, kwargs...) = getfield(@__MODULE__, format)(args...; kwargs...)






#= Aux functions =#

"""
    planck_λ(λ, T)

Compute the Planck function in cgs units.
"""
planck_λ(λ, T) = twohc2 /λ^5 /(exp(hc_k / (λ*T)) - 1.0)

"""
    planck_ν(ν, T)

Compute the Planck function in cgs units.
"""
planck_ν(ν, T) = 2 * HPlanck * ν^3 / CLight^2 / (exp(HPlanck * ν / (KBoltzmann * T)) - 1.0)
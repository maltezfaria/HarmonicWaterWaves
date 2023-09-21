# adapted from MakiePublication.jl
# TEXTWIDTH = 469.75499 # in pt
function makie_theme(;
    width=469.75499/2,
    height=width*0.68,
)
    axis_theme = (
        xlabelsize=10,
        ylabelsize=10,
        spinewidth=1.1,
        xticklabelsize=8,
        yticklabelsize=8,
        # xgridstyle=:dash, ygridstyle=:dash,
        xgridvisible=false,
        ygridvisible=false,
        xtickalign=1,
        ytickalign=1,
        xticksize=5,
        yticksize=5,
        xtickwidth=0.8,
        ytickwidth=0.8,
        xminorticksvisible=true,
        yminorticksvisible=true,
        xminortickalign=1,
        yminortickalign=1,
        xminorticks=IntervalsBetween(5),
        yminorticks=IntervalsBetween(5),
        xminorticksize=3,
        yminorticksize=3,
        xminortickwidth=0.75,
        yminortickwidth=0.75,
        xlabelpadding=-2,
        ylabelpadding=2,
    )

    line_theme = (
        linewidth=2.0,
    )

    scatter_theme = (
        markersize=7,
    )

    legend_theme = (
        nbanks=1,
        framecolor=(:grey, 0.5),
        framevisible=true,
        bgcolor=(:white,0.5),
        labelsize=7.5,
        padding=(1,1,-4,-4), # usually too much white space on top of legends
        margin=(0, 0, 0, 0),
        # position=:rt, # l=left, r=right, c=center; b=bottom, t=top, c=center
        rowgap=-10,
        colgap=4,
    )

    heatmap_theme = (
        colormap=:inferno,
        interpolate=true
    )

    colorbar_theme = (
        ticklabelsize=8,
        labelsize=8,
        labelrotation=0,
        colormap=:inferno,
    )

    return Theme(figure_padding=8,
        resolution=(width, height),
        # font="Helvetica",
        Axis=axis_theme,
        Lines=line_theme,
        Scatter=scatter_theme,
        Legend=legend_theme,
        Heatmap=heatmap_theme,
        Colorbar=colorbar_theme,
    )
end

function set_makie_theme!(;kwargs...)
    theme = makie_theme(;kwargs...)
    set_theme!(theme)
    return theme
end

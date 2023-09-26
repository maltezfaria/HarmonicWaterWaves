using Pkg
using Logging

dir_paper = @__DIR__
Pkg.activate(joinpath(dir_paper, "."))
Pkg.instantiate()

dir_scripts = joinpath(dir_paper, "scripts")

isdir(joinpath(dir_paper, "logs")) || mkdir(joinpath(dir_paper, "logs"))
isdir(joinpath(dir_paper, "figures")) || mkdir(joinpath(dir_paper, "figures"))

dir_logs = joinpath(dir_paper, "logs")
dir_figures = joinpath(dir_paper, "figures")
dir_animations = joinpath(dir_paper, "animations")

# cleanup figures, logs, and animations
for file in readdir(dir_figures)
    if file != "WW_problem.pdf"
        rm(joinpath(dir_figures, file))
    end
end
foreach(file -> rm(joinpath(dir_logs, file)), readdir(dir_logs))
foreach(file -> rm(joinpath(dir_animations, file)), readdir(dir_logs))

# map script file to figure file
script2fig = Dict(
    "wavemaker_solution.jl" => "fig2a",
    "wavemaker_pml_convergence.jl" => "fig2b",
    "modes_vary_impedance.jl" => "fig3",
    "pml_real_and_imag.jl" => "fig4a",
    "wavemaker_stretching_convergence.jl" => "fig4b",
    "wavemaker_mesh_convergence.jl" => "fig5",
    "jellyfish_scattering.jl" => "fig6",
    "jellyfish_scattering_mesh_convergence.jl" => "fig7",
    "double_piercing_scattering.jl" => "fig8",
    "double_piercing_scattering_mesh_convergence.jl" => "fig9",
    "step_scattering.jl" => "fig10",
    "double_piercing_eigenvalue.jl" => "fig11",
)

# generate all figures
for (file, fig) in script2fig
    @info "Executing $file to generate $fig"
    full_path = joinpath(dir_scripts, file)
    log_file  = joinpath(dir_logs, "$(fig).log")
    io        = open(log_file, "w+")
    logger    = ConsoleLogger(io)
    # wrap each include around a module to avoid namespace pollution and
    # conflicts
    with_logger(logger) do
        @eval module $(gensym(file))
        include($full_path)
        end
    end
    flush(io)
    close(io)
end

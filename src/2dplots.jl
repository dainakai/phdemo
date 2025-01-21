using JSON
using Plots

function twodimplots()
    date = "20250121"
    scenenum = 10
    data_dir = joinpath(date, "06_analysis", "C000H001S"*lpad(scenenum, 4, '0'))
    output_dir = joinpath(date, "06_analysis", "C000H001S"*lpad(scenenum, 4, '0'))

    trajectory_stats = JSON.parsefile(joinpath(data_dir, "trajectory_stats.json"))

    # Convert to DataFrame for visualization
    df = DataFrame(
        particle_id=collect(keys(trajectory_stats)),
        avg_diameter=[x[1] for x in values(trajectory_stats)],
        avg_velocity=[x[2] for x in values(trajectory_stats)]
    )

    # Visualization - Scatter plot
    scatter(df.avg_diameter, df.avg_velocity,
        xlabel="Average Diameter",
        ylabel="Average Velocity",
        title="Particle Size vs Velocity",
        legend=false,
        markersize=3,
        markercolor=:blue,
        ylims=(-300, 300),
        dpi=300,
    )
    savefig(joinpath(output_dir, "size_velocity_scatter_x.png"))

    # Visualization - Heatmap
    histogram2d(df.avg_diameter, df.avg_velocity,
        xlabel="Average Diameter",
        ylabel="Average Velocity",
        title="Particle Size-Velocity Distribution",
        bins = (200, 200),
        color=:heat,
        # xlims=(0, 100),
        ylims=(-300, 300),
        dpi=300,
    )
    savefig(joinpath(output_dir, "size_velocity_heatmap_x.png"))

    # Convert to DataFrame for visualization
    df = DataFrame(
        particle_id=collect(keys(trajectory_stats)),
        avg_diameter=[x[1] for x in values(trajectory_stats)],
        avg_velocity=[x[3] for x in values(trajectory_stats)]
    )

    # Visualization - Scatter plot
    scatter(df.avg_diameter, df.avg_velocity,
        xlabel="Average Diameter",
        ylabel="Average Velocity",
        title="Particle Size vs Velocity",
        legend=false,
        markersize=3,
        markercolor=:blue,
        ylims=(-300, 300),
        dpi=300,
    )
    savefig(joinpath(output_dir, "size_velocity_scatter_y.png"))

    # Visualization - Heatmap
    histogram2d(df.avg_diameter, df.avg_velocity,
        xlabel="Average Diameter",
        ylabel="Average Velocity",
        title="Particle Size-Velocity Distribution",
        bins = (200, 200),
        color=:heat,
        # xlims=(0, 100),
        ylims=(-300, 300),
        dpi=300,
    )
    savefig(joinpath(output_dir, "size_velocity_heatmap_y.png"))

    # Convert to DataFrame for visualization
    df = DataFrame(
        particle_id=collect(keys(trajectory_stats)),
        avg_diameter=[x[1] for x in values(trajectory_stats)],
        avg_velocity=[x[4] for x in values(trajectory_stats)]
    )

    # Visualization - Scatter plot
    scatter(df.avg_diameter, df.avg_velocity,
        xlabel="Average Diameter",
        ylabel="Average Velocity",
        title="Particle Size vs Velocity",
        legend=false,
        markersize=3,
        markercolor=:blue,
        ylims=(-300, 300),
        dpi=300,
    )
    savefig(joinpath(output_dir, "size_velocity_scatter_z.png"))

    # Visualization - Heatmap
    histogram2d(df.avg_diameter, df.avg_velocity,
        xlabel="Average Diameter",
        ylabel="Average Velocity",
        title="Particle Size-Velocity Distribution",
        bins = (200, 200),
        color=:heat,
        # xlims=(0, 100),
        ylims=(-300, 300),
        dpi=300,
    )
    savefig(joinpath(output_dir, "size_velocity_heatmap_z.png"))
end

twodimplots()
using JSON
using DataFrames
using Plots
using Statistics
using ParticleHolography
using MetaGraphsNext
using Graphs
using UUIDs
using Glob
using YAML
using CSV
using LinearAlgebra
using DelimitedFiles

using SavitzkyGolay
using ProgressMeter

using DataStructures

function enum_paths(edges::Vector{Vector{Int}})
    adj_list = DefaultDict{Int, Set{Int}}(() -> Set{Int}())
    for edge in edges
        u, v = edge
        push!(adj_list[u], v)
        push!(adj_list[v], u)
    end

    visited = Set{Int}()
    paths = []

    function dfs(node, path)
        push!(visited, node)
        push!(path, node)
        for neighbor in adj_list[node]
            if neighbor ∉ visited
                dfs(neighbor, path)
            end
        end
    end

    for node in keys(adj_list)
        if node ∉ visited
            path = []
            dfs(node, path)
            push!(paths, sort(path))
        end
    end

    return unique(paths)
end

function remove_duplicate!(paths)
    for path in paths
        path .= unique(path)
    end
end

function smoothed_traj(paths, fulldict, windowsize)
    smootheddict = Dict()
    @showprogress for path in paths
        data = reduce(hcat, [fulldict[pt] for pt in path])
        xdata = data[2,:]
        ydata = data[3,:]
        zdata = data[4,:]

        xdata_smooth = savitzky_golay(xdata, windowsize, 3)
        ydata_smooth = savitzky_golay(ydata, windowsize, 3)
        zdata_smooth = savitzky_golay(zdata, windowsize, 3)
        vxdata_smooth = savitzky_golay(xdata, windowsize, 3, deriv=1)
        vydata_smooth = savitzky_golay(ydata, windowsize, 3, deriv=1)
        vzdata_smooth = savitzky_golay(zdata, windowsize, 3, deriv=1)

        for (i,pt) in enumerate(path)
            if haskey(smootheddict, pt)
                @warn "Overwriting existing data for $pt"
            end
            smootheddict[pt] = [fulldict[pt][1],  xdata_smooth.y[i], ydata_smooth.y[i], zdata_smooth.y[i], vxdata_smooth.y[i], vydata_smooth.y[i], vzdata_smooth.y[i]]
        end
    end
    return smootheddict
end

function reconnect_cand_bool(path1end, path2top, missingpoints, smootheddict, inspectradius, dim3weight)
    @assert length(first(smootheddict)[2]) >= 7 "dict should have at least 7 values: [framenum, x, y, z, vx, vy, vz]"
    rend = deepcopy(smootheddict[path1end][2:4])
    vend = deepcopy(smootheddict[path1end][5:7])
    rtop = deepcopy(smootheddict[path2top][2:4])
    vtop = deepcopy(smootheddict[path2top][5:7])

    for _ in 1:missingpoints+1
        rend += vend/2
        rtop -= vtop/2 
    end

    if node_distance(rend, rtop, dim3weight) < inspectradius
        return true, node_distance(rend, rtop, dim3weight)
    else
        return false, 0
    end
end

function reconnect_traj(paths, smootheddict, max_time_interval=3, inspectradius=10, dim3weight=1)
    pathheads = [path[1] for path in paths]
    pathtails = [path[end] for path in paths]

    pathconnectgraph = MetaGraph(
        DiGraph(),
        label_type=Int,
        vertex_data_type=NTuple{2,UUID},
    )

    for idx in eachindex(paths)
        add_vertex!(pathconnectgraph, idx, (pathheads[idx], pathtails[idx]))
    end

    for (tailidx, pathtail) in enumerate(pathtails)
        candinateheads = []
        candinatevolediffs = []
        for (headidx, pathhead) in enumerate(pathheads)
            if pathhead == pathtail
                continue
            end

            time_interval = smootheddict[pathhead][1] - smootheddict[pathtail][1]
            if time_interval > 0 && time_interval <= max_time_interval
                flag, velodiff = reconnect_cand_bool(pathtail, pathhead, time_interval-1, smootheddict, inspectradius, dim3weight)
                if flag
                    push!(candinateheads, headidx)
                    push!(candinatevolediffs, velodiff)
                end
            end
        end
        if length(candinateheads) > 0
            bestheadidx = candinateheads[argmin(candinatevolediffs)]
            add_edge!(pathconnectgraph, tailidx, bestheadidx)
        end
    end

    connectedpaths = []
    binconnections = enum_edge(pathconnectgraph, Int)

    connections = enum_paths(binconnections)

    for connection in connections
        append!(connectedpaths, [vcat(paths[connection]...)])
    end

    isolatedpaths = filter(v -> degree(pathconnectgraph, v) == 0, vertices(pathconnectgraph))
    display(isolatedpaths)
    for isolatedpath in isolatedpaths
        append!(connectedpaths, [paths[isolatedpath]])
    end

    return connectedpaths
end

function analysis()
    # Configuration
    date = "20250121"
    scenenum = 10
    data_dir = joinpath(date, "03_particles", "C000H001S"*lpad(scenenum, 4, '0'))
    output_dir = joinpath(date, "06_analysis", "C000H001S"*lpad(scenenum, 4, '0'))
    mkpath(output_dir)

    variables = YAML.load_file(joinpath(date, "variables.yaml"))

    # Load JSON files and create trajectory graphs
    json_files = glob(joinpath(data_dir, "*.json"))
    dicts = ParticleHolography.dictload.(json_files)

    # Create connection graphs between consecutive frames
    graphs = [Labonte(dict1, dict2, dim3weight=1) for (dict1, dict2) in zip(dicts[1:end-1], dicts[2:end])]

    # Initialize paths from first graph
    paths = ParticleHolography.enum_edge(graphs[1])

    # Connect paths across all graphs
    for graph in graphs[2:end]
        append_path!(paths, graph)
    end

    # Filter short trajectories
    filter!(path -> length(path) > 5, paths)

    # Create full dictionary of all particle positions
    fulldict = gen_fulldict(json_files)

    smoothedtraj = smoothed_traj(paths, fulldict, 5)
    connectedpaths = reconnect_traj(paths, smoothedtraj, 3, 20, 1)

    # Initialize data structures
    trajectory_stats = Dict{UUID,Tuple{Float64,Float64,Float64, Float64}}() # avg_diameter, avg_vx, avg_vy, avg_vz

    # Calculate statistics for connected paths
    for path in paths
        vxdata = [smoothedtraj[id][5] for id in path] .* variables["dx"]
        vydata = [smoothedtraj[id][6] for id in path] .* variables["dx"]
        vzdata = [smoothedtraj[id][7] for id in path] .* variables["dz"]

        diamdata = [fulldict[id][end] for id in path] .* variables["dx"]

        # Use first UUID in path as identifier
        trajectory_stats[path[1]] = (mean(diamdata), mean(vxdata), mean(vydata), mean(vzdata))
    end

    # save trajectory stats to JSON
    JSON.open(joinpath(output_dir, "trajectory_stats.json"), "w") do io
        JSON.print(io, trajectory_stats)
    end

    trajsavedir = joinpath(date, "07_trajectories", "C000H001S"*lpad(scenenum, 4, '0'))
    !isdir(trajsavedir) && mkpath(trajsavedir)
    for path in connectedpaths
        # 軌跡ごとに情報を保存。smoothedtrajを参照して、位置、速度、粒径をCSVで保存。columns: frame, x, y, z, vx, vy, vz, diameter
        frames = [smoothedtraj[id][1] for id in path]
        xdata = [smoothedtraj[id][2] for id in path]
        ydata = [smoothedtraj[id][3] for id in path]
        zdata = [smoothedtraj[id][4] for id in path]
        vxdata = [smoothedtraj[id][5] for id in path]
        vydata = [smoothedtraj[id][6] for id in path]
        vzdata = [smoothedtraj[id][7] for id in path]
        diameterdata = [fulldict[id][end] for id in path]
        df = DataFrame([frames, xdata, ydata, zdata, vxdata, vydata, vzdata, diameterdata], [:frame, :x, :y, :z, :vx, :vy, :vz, :diameter])
        CSV.write(joinpath(trajsavedir, "paths_$(path[1]).csv"), df)
    end

end

analysis()
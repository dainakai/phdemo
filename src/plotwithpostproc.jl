using CSV
using DataFrames
using UUIDs
using ParticleHolography
using Plots

date = "20250121"
scenenum = 10
traj_dir = joinpath(date, "07_trajectories", "C000H001S"*lpad(scenenum, 4, '0'))
output_dir = joinpath(date, "06_analysis", "C000H001S"*lpad(scenenum, 4, '0'))
csv_files = readdir(traj_dir, join=true)

trajectories = Dict{Int, DataFrame}()
global traj_counter = 1

for file in csv_files
    df = CSV.read(file, DataFrame)
    
    df[!, :traj_id] .= traj_counter
    
    trajectories[traj_counter] = df
    
    global traj_counter += 1
end

combined_df = vcat(values(trajectories)...)

result = combine(groupby(combined_df, :frame)) do group
    (traj_ids = group.traj_id,
     x_coords = group.x,
     y_coords = group.y,
     z_coords = group.z)
end

CSV.write(joinpath(output_dir, "trajectory_results.csv"), result)

colors = palette(:tab20)

plot()
anim = @animate for idx in 1:100
    particleplot(Dict(uuid4() => Float32.([-100, -100, -1000000000])), legend = false, scaling=(10.0, 10.0, -100.0), shift=(0.0, 0.0, 1e5), xlabel="x [µm]", ylabel="z [µm]", zlabel="y [µm]", xlim=(0,10240), ylim=(0,1e5), zlim=(0,10240), dpi=300)

    for idx in findall(result[:,1] .== idx)
        particleplot!(Dict(uuid4() => Float32.([result[idx,3], result[idx,4], result[idx,5]])), color=colors[result[idx,2] % length(colors) + 1], scaling=(10.0, 10.0, -100.0), shift=(0.0, 0.0, 1e5), camera=(30, 30))
    end
end

gif(anim, joinpath(output_dir, "smoothedplot.gif"), fps = 10)


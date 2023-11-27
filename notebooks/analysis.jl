### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 3bd56622-8ae6-11ee-2261-e39faf579a7e
using DrWatson

# ╔═╡ 3e4df9ac-014f-4c4f-bbf2-82c94050b31d
@quickactivate

# ╔═╡ cc3d4f10-b13b-4381-8bbb-c5a5292ef8ea
begin
  using OrdinaryDiffEq
  using StaticArrays
  using LinearAlgebra
  BLAS.set_num_threads(6)
  using Optimization

  using GLMakie
  using HDF5
  using DataFrames
end

# ╔═╡ cb6dea58-b3e3-4ab5-82b5-dcee93452d6c
include(srcdir("RMA.jl"))

# ╔═╡ 14cce842-94a4-4e2d-a444-fde8b54ae0d6
df = collect_results!(datadir("basins_dense"))

# ╔═╡ 6ef38bfd-2bdf-4745-bc12-364c5a87b4e8
begin
  row = df[10, :]
  nt = struct2dict(row.p) |> values |> first
  p = RMAParameters(nt...)
end

# ╔═╡ f7f83bfc-df1a-437c-927f-3c3285316814
p.r

# ╔═╡ 6c9b84e3-71bc-4981-9961-92f9435c552e
md"""Basin of attraction"""

# ╔═╡ bdce918d-dbae-4a1e-89b4-a9b950b176e9
let xs = row.xs, ys = row.ys, zs = row.zs, eq = e3(p)
  fig, ax, hm = heatmap(
    xs, ys, zs,
    colormap=Makie.Categorical(:viridis),
    axis=(; xlabel="predator", ylabel="prey")
  )
  scatter!(eq[1], eq[2],
    color=:red,
    #marker=:cross,
    markersize=10)
  scatter!(p.r / p.c, 0,
    color=:red,
    #marker=:cross,
    markersize=10)
  scatter!(p.mu, 0,
    color=:red,
    #marker=:cross,
    markersize=10)
  #Makie.save(plotsdir("basin.png"), fig, px_per_unit = 6)
  current_figure()
end

# ╔═╡ 5218dfe1-84d9-4669-a882-672bf2ee07f5
begin
  zs(i) = df[i, :].zs
  function p_fn(i)
    row = df[i, :]
    nt = struct2dict(row.p) |> values |> first
    p = RMAParameters(nt...)
  end
  unstab_zs = zs(10) - zs(1)
end

# ╔═╡ bc95b6b1-4d17-4368-a2a2-8d0d58e6c94b
md"""Basin unstable regions"""

# ╔═╡ 02e83fc3-33a7-41b2-8b13-6be5f9da846e
let xs = row.xs, ys = row.ys, zs = unstab_zs
  fig, ax, hm = heatmap(
    xs, ys, zs,
    colormap=Makie.Categorical(:viridis),
    axis=(; xlabel="predator", ylabel="prey")
  )
  eq = e3(p_fn(10))
  scatter!(eq[1], eq[2],
    color=:red,
    #marker=:cross,
    markersize=10)
  scatter!(p_fn(10).r / p.c, 0,
    color=:red,
    #marker=:cross,
    markersize=10)
  scatter!(p_fn(10).mu, 0,
    color=:red,
    #marker=:cross,
    markersize=10)

  eq = e3(p_fn(1))
  scatter!(eq[1], eq[2],
    color=:black,
    #marker=:cross,
    markersize=10)
  scatter!(p_fn(1).r / p.c, 0,
    color=:black,
    #marker=:cross,
    markersize=10)
  scatter!(p_fn(1).mu, 0,
    color=:black,
    #marker=:cross,
    markersize=10)
  #Makie.save(plotsdir("basin_instab.png"), fig, px_per_unit=6)
  current_figure()
end

# ╔═╡ fbc3ef99-cde3-4764-bbcd-3151d278eada
function plot_basin!(f, xs, ys, zs, p; zs_orig)
	diff_zs = zs + zs_orig
  heatmap!(
    f,
    xs, ys, diff_zs,
    colormap=Makie.Categorical(:viridis)
  )
  eq = e3(p)
  scatter!(eq[1], eq[2],
    color=:red,
    #marker=:cross,
    markersize=10)
  scatter!(p.r / p.c, 0,
    color=:red,
    #marker=:cross,
    markersize=10)
  scatter!(p.mu, 0,
    color=:red,
    #marker=:cross,
    markersize=10)
  #Makie.save(plotsdir("basin.png"), fig, px_per_unit = 6)
end

# ╔═╡ 65a0fd01-dd0d-45ea-ba61-1058e4b02045
zs_iter = [zs(i) for i in 1:10]

# ╔═╡ 2784c9ad-3b0e-419f-9289-fd05b84a9a1a
# testing
let
  vid_scene = Figure()
  vid_ax = Axis(vid_scene[1, 1]; xlabel="predator", ylabel="prey")
  plot_basin!(vid_ax, row.xs, row.ys, zs_iter[1], p_fn(1); zs_orig=zs_iter[10])
  current_figure()
end

# ╔═╡ a2547a36-5b4f-49dd-a1e0-1ef312f59aa1
begin
  vid_scene = Figure()
  vid_ax = Axis(vid_scene[1, 1]; xlabel="predator", ylabel="prey")
  record(vid_scene, "ptip.mp4", reverse(1:10);
    framerate=5) do i
    plot_basin!(vid_ax, row.xs, row.ys, zs(i), p_fn(i); zs_orig = zs(10))
  end
end

# ╔═╡ Cell order:
# ╠═3bd56622-8ae6-11ee-2261-e39faf579a7e
# ╠═3e4df9ac-014f-4c4f-bbf2-82c94050b31d
# ╠═cc3d4f10-b13b-4381-8bbb-c5a5292ef8ea
# ╠═cb6dea58-b3e3-4ab5-82b5-dcee93452d6c
# ╠═14cce842-94a4-4e2d-a444-fde8b54ae0d6
# ╠═6ef38bfd-2bdf-4745-bc12-364c5a87b4e8
# ╠═f7f83bfc-df1a-437c-927f-3c3285316814
# ╠═6c9b84e3-71bc-4981-9961-92f9435c552e
# ╠═bdce918d-dbae-4a1e-89b4-a9b950b176e9
# ╠═5218dfe1-84d9-4669-a882-672bf2ee07f5
# ╠═bc95b6b1-4d17-4368-a2a2-8d0d58e6c94b
# ╠═02e83fc3-33a7-41b2-8b13-6be5f9da846e
# ╠═fbc3ef99-cde3-4764-bbcd-3151d278eada
# ╠═65a0fd01-dd0d-45ea-ba61-1058e4b02045
# ╠═2784c9ad-3b0e-419f-9289-fd05b84a9a1a
# ╠═a2547a36-5b4f-49dd-a1e0-1ef312f59aa1

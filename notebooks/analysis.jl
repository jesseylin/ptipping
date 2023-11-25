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

  using CairoMakie
  using HDF5
  using DataFrames
end

# ╔═╡ cb6dea58-b3e3-4ab5-82b5-dcee93452d6c
include(srcdir("RMA.jl"))

# ╔═╡ 14cce842-94a4-4e2d-a444-fde8b54ae0d6
df = collect_results!(datadir("basins"))

# ╔═╡ 6ef38bfd-2bdf-4745-bc12-364c5a87b4e8
begin
  row = df[10, :]
  nt = struct2dict(row.p) |> values |> first
  p = RMAParameters(nt...)
end

# ╔═╡ ad0371f9-c0ad-4dfb-af4d-5b1fbc8dedc8
p.r

# ╔═╡ bdce918d-dbae-4a1e-89b4-a9b950b176e9
let xs = row.xs, ys = row.ys, zs = row.zs, eq = e3(p)
  fig, ax, hm = heatmap(
    xs, ys, zs,
    colormap=Makie.Categorical(:viridis),
    axis=(; xlabel="X", ylabel="Y")
  )
  scatter!(eq[1], eq[2],
    color=:red,
    marker=:cross,
    markersize=20)
  scatter!(p.r / p.c, 0,
    color=:red,
    marker=:cross,
    markersize=20)
  scatter!(p.mu, 0,
    color=:red,
    marker=:cross,
    markersize=20)
  current_figure()
end

# ╔═╡ 5218dfe1-84d9-4669-a882-672bf2ee07f5
begin
  zs(i) = df[i, :].zs
  unstab_zs = zs(10) - zs(1)
end

# ╔═╡ 02e83fc3-33a7-41b2-8b13-6be5f9da846e
let xs = row.xs, ys = row.ys, zs = unstab_zs, eq = e3(p)
  fig, ax, hm = heatmap(
    xs, ys, zs,
    colormap=Makie.Categorical(:viridis),
    axis=(; xlabel="X", ylabel="Y")
  )
  scatter!(eq[1], eq[2],
    color=:red,
    marker=:cross,
    markersize=20)
  scatter!(p.r / p.c, 0,
    color=:red,
    marker=:cross,
    markersize=20)
  scatter!(p.mu, 0,
    color=:red,
    marker=:cross,
    markersize=20)
  current_figure()
end

# ╔═╡ Cell order:
# ╠═3bd56622-8ae6-11ee-2261-e39faf579a7e
# ╠═3e4df9ac-014f-4c4f-bbf2-82c94050b31d
# ╠═cc3d4f10-b13b-4381-8bbb-c5a5292ef8ea
# ╠═cb6dea58-b3e3-4ab5-82b5-dcee93452d6c
# ╠═14cce842-94a4-4e2d-a444-fde8b54ae0d6
# ╠═6ef38bfd-2bdf-4745-bc12-364c5a87b4e8
# ╠═ad0371f9-c0ad-4dfb-af4d-5b1fbc8dedc8
# ╠═bdce918d-dbae-4a1e-89b4-a9b950b176e9
# ╠═5218dfe1-84d9-4669-a882-672bf2ee07f5
# ╠═02e83fc3-33a7-41b2-8b13-6be5f9da846e

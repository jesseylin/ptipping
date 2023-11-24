using DrWatson
@quickactivate

include(srcdir("RMA.jl"))

using DifferentialEquations
using StaticArrays
using LinearAlgebra
BLAS.set_num_threads(6)
using ProgressMeter


function poincare_iterate(
    u0, p;
    section_transverse=0, section_int=(0, -1e-4), tmax=1000
)
    @assert section_int[1] < section_int[2]
    is_in_section(u) = section_int[1] < u[2] < section_int[2]
    if !is_in_section(u0)
        error("""
              Must choose initial guess u0 inside the section.
              Chose u0=$(u0), but
              section_transverse=$(section_transverse)
              and
              section_int=$(section_int)
              """)
    end

    condition(u, t, integrator) = u[1] - section_transverse


    function affect!(int)
        if is_in_section(int.u)
            terminate!(int)
        end
    end

    cb = ContinuousCallback(condition, affect!)
    ode = ODEProblem(RMA.rma_step, u0, (0.0, tmax), p)
    sol = solve(ode, AutoVern7(Rodas5()); callback=cb, save_everystep=false)

    SA[section_transverse, sol.u[end][2]], sol.t[end]
end
function integrate_to_limit(u0, p; tequil=1000)
    ode = ODEProblem(RMA.rma_step, u0, (0, tequil), p)
    # solve using a more sophisticated stiff solver
    sol = solve(ode, AutoVern7(Rodas5()); save_everystep=false)
    p_lc = sol[end]
    p_lc
end
function find_fixed_point(u0, p; maxtries=1e5, reltol=1e-4)
    # ensure we're limiting
    u0 = integrate_to_limit(u0, p)

    # ensure we're on the section
    eq = RMA.e3(p)
    u0, _ = poincare_iterate(u0, p; section_transverse=eq[1], section_int=(0, 1))

    # determine size of section
    d = abs(u0[2] - eq[2])
    section_int = (u0[2] - d, u0[2] + d)

    # iterate
    err = Inf
    num_try = 0
    T = 0.0
    while err > reltol
        if num_try > maxtries
            println("exceeded maxtries")
            @show err
            break
        end
        u0_old = copy(u0)
        u0, T = poincare_iterate(u0_old, p;
            section_transverse=eq[1],
            section_int=section_int)
        err = abs((u0[2] - u0_old[2]) / u0_old[2])
        num_try = num_try + 1
    end
    @show err
    return u0, T
end
function find_limit_cycle(r; reltol=1e-8, maxtries=1e5)
    p = RMA.RMAParameters(r=r)
    eq = RMA.e3(p)

    # start within the limit cycle basin
    u0 = eq + SA[0, -1e-6]

    # do iteration
    find_fixed_point(u0, p; maxtries=maxtries, reltol=reltol)
end

# ╔═╡ 45ea640c-0ef9-461e-82f9-6452724be1b8
function brute_force_discriminator(X, Y; p, tequil=1000, abstol=1e-3)
    u0 = SA[X, Y]
    u0_lc = integrate_to_limit(u0, p; tequil=tequil)
    if (u0_lc[1] < -1) | (u0_lc[2] < -1)
        return 0
    else
        return norm(u0_lc) < abstol ? 0 : 1
    end
end

# ╔═╡ 62bd9bc4-03c5-44ef-95d0-7fdbf200cec4
function basin(r; num_samples=100)
    num_samples = num_samples
    p = RMA.RMAParameters(r=r)
    xs = LinRange(0, p.r / p.c + 8, num_samples)
    ys = LinRange(0, 0.03, num_samples)
    zs = [brute_force_discriminator(x, y; p=p) for x in xs, y in ys]
    xs, ys, zs
end

function main()
    r__range = 1.6:0.1:2.5
    num_samples = 300
    @showprogress for r in r_range
        p = RMA.RMAParameters(r=r)
        # guard to not redo experiment
        if isfile(savename(p, "jld2"))
            println("Already did r=$(r), skipping...")
            continue
        end
        xs, ys, zs = basin(r; num_samples=num_samples)
        output = @strdict xs ys zs p
        @tagsave(datadir("basins", savename(p, "jld2")), output)
    end
end

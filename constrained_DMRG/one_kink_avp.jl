#!/usr/bin/env julia

# N=100 one-kink constrained-DMRG/AVP reference pipeline.
# The solver source is bundled in this directory; no Rivanna path is needed.
# H_eff = H_XXZ - sum_i lambda[i] * Sz[i].

include(joinpath(@__DIR__, "DMRGOpt", "DMRGOpt.jl"))
using .DMRGOpt
using LinearAlgebra
using Printf
using Random
using Serialization
using Statistics

envint(key, default) = parse(Int, get(ENV, key, string(default)))
envfloat(key, default) = parse(Float64, get(ENV, key, string(default)))

function settings()
    n = envint("AVP_N", 100)
    n % 4 == 0 || error("AVP_N must be divisible by four so the one-kink target has total Sz=0")
    (; N=n, delta=envfloat("AVP_DELTA", 3.0), chi=envint("AVP_CHI", 64),
       amplitude=envfloat("AVP_AMPLITUDE", 0.439364),
       stages=envint("AVP_STAGES", 250), max_sweeps=envint("AVP_MAX_SWEEPS", 100),
       local_max_iter=envint("AVP_LOCAL_MAX_ITER", 500),
       local_tol=envfloat("AVP_LOCAL_TOL", 1e-2),
       global_tol=envfloat("AVP_GLOBAL_TOL", 1e-2),
       lambda_tol=envfloat("AVP_LAMBDA_TOL", 1e-2),
       alpha=envfloat("AVP_ALPHA", 0.005),
       dt=envfloat("AVP_DT", 0.05), gamma=envfloat("AVP_GAMMA", 0.05),
       seed=envint("AVP_SEED", 1), bias=envfloat("AVP_MPS_BIAS", 0.03),
       output=get(ENV, "AVP_OUTPUT", joinpath(@__DIR__, "output", "N100_Delta3_one_kink")))
end

function kink_target(p)
    wall = p.N ÷ 2
    [p.amplitude * (isodd(i) ? 1.0 : -1.0) * (i <= wall ? 1.0 : -1.0)
     for i in 1:p.N]
end

function biased_mps(p, target)
    Random.seed!(p.seed)
    psi = random_mps(p.N, 2, p.chi)
    for i in 1:p.N
        disfavored = target[i] > 0 ? 2 : 1
        @views psi.tensors[i][:, disfavored, :] .*= p.bias
    end
    psi.tensors[1] ./= sqrt(real(measure_norm(psi)))
    psi
end

function save_checkpoint(path, state, stage, target, lambda, measured, accepted)
    serialize(path, (; stage, target=copy(target), lambda=copy(lambda),
                     measured=copy(measured), accepted,
                     tensors=deepcopy(state.mps.tensors)))
end

function save_t0_exports(output, state, target, lambda, measured, accepted)
    checkpoint = joinpath(output, "t0_mps.jls")
    save_checkpoint(checkpoint, state, 0, target, lambda, measured, accepted)
    open(joinpath(output, "t0_profile.csv"), "w") do io
        println(io, "site,target_Sz,measured_Sz,lambda")
        for i in eachindex(target)
            println(io, join((i, target[i], measured[i], lambda[i]), ','))
        end
    end
    open(joinpath(output, "t0_lambda.csv"), "w") do io
        println(io, "site,lambda")
        for i in eachindex(lambda)
            println(io, join((i, lambda[i]), ','))
        end
    end
end

function main()
    p = settings()
    p.N >= 4 && p.chi >= 1 && 0 < p.amplitude < 0.5 || error("Invalid N, chi, or amplitude")
    p.stages >= 0 && p.max_sweeps > 0 && p.local_max_iter > 0 || error("Invalid iteration count")
    0 < p.bias <= 1 && p.alpha > 0 && p.dt > 0 && p.gamma > 0 || error("Invalid bias or step size")
    isdir(p.output) && !isempty(readdir(p.output)) && error("Output directory is not empty: $(p.output)")
    mkpath(p.output)
    BLAS.set_num_threads(envint("AVP_BLAS_THREADS", 8))

    target = kink_target(p)
    abs(sum(target)) < 1e-12 || error("Target must sum to zero")
    lambda = 2.0 .* target
    lambda .-= mean(lambda)
    psi = biased_mps(p, target)
    state = MPSState(psi, build_xxz_mpo(p.N, 1.0, p.delta; d=2, h=-lambda))
    solver = LanczosSolver(3, 1) # Original local Krylov setting.
    options = DMRGOptions(p.chi, 1e-12, 2)

    @printf("N=%d Delta=%.2f chi=%d kink=%d|%d stages=%d\n",
            p.N, p.delta, p.chi, p.N÷2, p.N÷2+1, p.stages)
    @printf("H_eff=H_XXZ-sum(lambda_i Sz_i); local alpha=%.4g; dt=%.4g gamma=%.4g\n",
            p.alpha, p.dt, p.gamma)
    println("output: ", p.output)

    open(joinpath(p.output, "sweeps.csv"), "w") do sweeps
        open(joinpath(p.output, "profiles.csv"), "w") do profiles
            open(joinpath(p.output, "stages.csv"), "w") do stages
                println(sweeps, "stage,time,sweep,max_abs_g,max_abs_delta_lambda,right_failures,left_failures,Eeff,E_XXZ,total_Sz,accepted")
                println(profiles, "stage,time,sweep,site,target_Sz,measured_Sz,lambda")
                println(stages, "stage,time,sweeps_used,accepted,max_abs_g,max_abs_delta_lambda,E_XXZ,core_left_Sz,core_right_Sz")

                for stage in 0:p.stages
                    stable = 0
                    accepted = false
                    maxg = Inf
                    dlambda = Inf
                    measured = sz_profile(state.mps)
                    bare_energy = NaN
                    sweeps_used = 0

                    for sweep in 1:p.max_sweeps
                        sweeps_used = sweep
                        lambda_before = copy(lambda)
                        # At every site: solve local H_eff, update local lambda,
                        # re-solve, then QR/LQ moves to the next site.
                        right = dmrg_sweep_one_site_constrained!(
                            state, solver, options, :right, lambda, target;
                            local_max_iter=p.local_max_iter, local_tol=p.local_tol,
                            alpha_init=p.alpha, alpha_min=p.alpha, alpha_max=p.alpha,
                            adaptive=false, sweep_number=sweep)
                        left = dmrg_sweep_one_site_constrained!(
                            state, solver, options, :left, lambda, target;
                            local_max_iter=p.local_max_iter, local_tol=p.local_tol,
                            alpha_init=p.alpha, alpha_min=p.alpha, alpha_max=p.alpha,
                            adaptive=false, sweep_number=sweep)

                        measured = sz_profile(state.mps)
                        maxg = maximum(abs, measured .- target)
                        dlambda = maximum(abs, lambda .- lambda_before)
                        e_eff = real(measure_energy(state)) / real(measure_norm(state.mps))
                        bare_energy = e_eff + dot(lambda, measured)
                        stable = maxg <= p.global_tol && dlambda <= p.lambda_tol ? stable + 1 : 0
                        accepted = stable >= 2

                        println(sweeps, join((stage, stage*p.dt, sweep, maxg, dlambda,
                            right.failed_sites, left.failed_sites, e_eff, bare_energy,
                            sum(measured), accepted), ','))
                        for i in 1:p.N
                            println(profiles, join((stage, stage*p.dt, sweep, i,
                                target[i], measured[i], lambda[i]), ','))
                        end
                        flush(sweeps)
                        flush(profiles)
                        @printf("stage=%3d sweep=%3d max|g|=%.4e max|dλ|=%.4e E_XXZ=%.9f core=(%+.5f,%+.5f) accepted=%s\n",
                            stage, sweep, maxg, dlambda, bare_energy,
                            measured[p.N÷2], measured[p.N÷2+1], accepted)
                        accepted && break
                    end

                    println(stages, join((stage, stage*p.dt, sweeps_used, accepted,
                        maxg, dlambda, bare_energy, measured[p.N÷2], measured[p.N÷2+1]), ','))
                    flush(stages)
                    save_checkpoint(joinpath(p.output, @sprintf("stage_%04d.jls", stage)),
                        state, stage, target, lambda, measured, accepted)
                    if stage == 0
                        save_t0_exports(p.output, state, target, lambda, measured, accepted)
                    end
                    accepted || @warn "Stage did not meet the convergence criterion; next target uses this approximate state" stage maxg dlambda
                    stage == p.stages && break

                    # AVP target evolution. Projection removes the constant
                    # multiplier mode and keeps sum_i target_i = 0.
                    projected_lambda = lambda .- mean(lambda)
                    target .-= (p.dt * p.gamma) .* projected_lambda
                    maximum(abs, target) < 0.5 || error("Evolved target exceeded spin-1/2 range")
                end
            end
        end
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end

# Algorithms/dmrg.jl

"""
    dmrg_sweep_one_site(state, solver, options, direction)

Single-site DMRG sweep. Bond-preserving (no SVD/truncation), so `chi_max` and
`cutoff` in `options` are ignored. Caller is responsible for ensuring the
initial MPS already has the working bond dimension — single-site DMRG cannot
grow χ from a product state.

Right sweep solves at sites 1..N-1 (orthog center at `state.center` must be 1
at entry; ends at N). Left sweep solves at sites N..2 (entry center N; ends 1).
"""
function dmrg_sweep_one_site(state::MPSState, solver::LanczosSolver, options::DMRGOptions, direction::Symbol)
    N = length(state.mps.tensors)
    last_energy = 0.0
    max_chi = 0

    if direction == :right
        for site in 1:N-1
            chi_l, d, chi_r = size(state.mps.tensors[site])

            left_env  = site == 1 ? state.environment.tensors[N+1] : state.environment.tensors[site-1]
            right_env = state.environment.tensors[site+1]
            Heff = OneSiteEffectiveHamiltonian(left_env, state.mpo.tensors[site], right_env)

            evec, eval = _solve(solver, Heff, vec(state.mps.tensors[site]))
            last_energy = real(eval)

            F = qr(reshape(evec, (chi_l*d, chi_r)))
            new_chi = size(F.R, 1)
            Q_full = Matrix(F.Q)
            Q = size(Q_full, 2) > new_chi ? Q_full[:, 1:new_chi] : Q_full
            R = Matrix(F.R)

            state.mps.tensors[site] = reshape(Q, (chi_l, d, new_chi))
            _update_left_environment(state, site)

            @tensoropt new_next[-1,-2,-3] := R[-1,4] * state.mps.tensors[site+1][4,-2,-3]
            state.mps.tensors[site+1] = new_next
            state.environment.tensors[site+1] = nothing
            state.center = site + 1
            max_chi = max(max_chi, new_chi)
        end
    else
        for site in N:-1:2
            chi_l, d, chi_r = size(state.mps.tensors[site])

            left_env  = state.environment.tensors[site-1]
            right_env = state.environment.tensors[site+1]
            Heff = OneSiteEffectiveHamiltonian(left_env, state.mpo.tensors[site], right_env)

            evec, eval = _solve(solver, Heff, vec(state.mps.tensors[site]))
            last_energy = real(eval)

            F = lq(reshape(evec, (chi_l, d*chi_r)))
            new_chi = size(F.L, 2)
            L = Matrix(F.L)
            Q_full = Matrix(F.Q)
            Q = size(Q_full, 1) > new_chi ? Q_full[1:new_chi, :] : Q_full

            state.mps.tensors[site] = reshape(Q, (new_chi, d, chi_r))
            _update_right_environment(state, site)

            @tensoropt new_prev[-1,-2,-3] := state.mps.tensors[site-1][-1,-2,4] * L[4,-3]
            state.mps.tensors[site-1] = new_prev
            state.environment.tensors[site-1] = nothing
            state.center = site - 1
            max_chi = max(max_chi, new_chi)
        end
    end

    return (E=last_energy, max_chi=max_chi)
end

function dmrg_sweep(state::MPSState,solver::LanczosSolver,options::DMRGOptions,direction::Symbol)
    N = length(state.mps.tensors)
    last_energy = 0.0
    if direction == :right sites = 1:N-1 else sites=N:-1:2 end
    if direction == :right
        for site in sites
            @tensoropt theta[-1,-2,-3,-4] := state.mps.tensors[site][-1,-2,5]*state.mps.tensors[site+1][5,-3,-4]
            chi_l,d1,d2,chi_r = size(theta)

            left_env = site == 1 ? state.environment.tensors[N+1] : state.environment.tensors[site-1]
            right_env = state.environment.tensors[site+2]
            Heff = TwoSiteEffectiveHamiltonian(left_env,state.mpo.tensors[site],state.mpo.tensors[site+1],right_env)

            evec,eval = _solve(solver,Heff,vec(theta))
            last_energy = real(eval)
            U,S,V = _svd_truncate(reshape(evec,(chi_l*d1,d2*chi_r)),options.chi_max,options.cutoff)
            state.mps.tensors[site] = reshape(U,(chi_l,d1,:))
            @tensoropt state.mps.tensors[site+1][-1,-2,-3] := Diagonal(S)[-1,4]*reshape(V,(:,d2,chi_r))[4,-2,-3]
            state.center = site+1
            _update_left_environment(state,site)
        end
    else
        for site in sites
            @tensoropt theta[-1,-2,-3,-4] := state.mps.tensors[site-1][-1,-2,5]*state.mps.tensors[site][5,-3,-4]
            chi_l,d1,d2,chi_r = size(theta)

            left_env = site == 2 ? state.environment.tensors[N+1] : state.environment.tensors[site-2]
            right_env = state.environment.tensors[site+1]
            Heff = TwoSiteEffectiveHamiltonian(left_env,state.mpo.tensors[site-1],state.mpo.tensors[site],right_env)

            evec,eval = _solve(solver,Heff,vec(theta))
            last_energy = real(eval)
            U,S,V = _svd_truncate(reshape(evec,(chi_l*d1,d2*chi_r)),options.chi_max,options.cutoff)
            state.mps.tensors[site] = reshape(V,(:,d2,chi_r))
            @tensoropt state.mps.tensors[site-1][-1,-2,-3] := reshape(U,(chi_l,d1,:))[-1,-2,4]*Diagonal(S)[4,-3]
            state.center = site-1
            _update_right_environment(state,site)
        end
    end
    return last_energy
end

# ============================================================
# Constrained local DMRG
#
# At each one-site / two-site visit:
#   1. freeze left/right environments
#   2. solve local effective GS
#   3. measure local <Sz>
#   4. update local lambda(s)
#   5. re-solve the SAME local problem
#   6. only after convergence, move orthogonality center
#
# Convention:
# H_eff = H - sum_i lambda_i Sz_i
#
# Therefore:
# lambda_i <- lambda_i + alpha * (target_i - measured_i)
# ============================================================

function _local_sz_one(A, Sz)
    d = size(A, 2)
    X = reshape(permutedims(A, (2,1,3)), d, :)
    den = sum(abs2, X)
    return real(sum(conj.(X) .* (Sz * X)) / den)
end

function _local_sz_two(theta, Sz)
    d1 = size(theta, 2)
    d2 = size(theta, 3)

    X1 = reshape(permutedims(theta, (2,1,3,4)), d1, :)
    X2 = reshape(permutedims(theta, (3,1,2,4)), d2, :)

    m1 = real(sum(conj.(X1) .* (Sz * X1)) / sum(abs2, X1))
    m2 = real(sum(conj.(X2) .* (Sz * X2)) / sum(abs2, X2))

    return m1, m2
end

function _adapt_local_alpha(alpha, residual, previous_residual;
                            alpha_min=1e-6,
                            alpha_max=0.05,
                            shrink=0.5,
                            grow=1.05,
                            good_ratio=0.8)

    previous_residual === nothing && return alpha, "initial"

    a = abs(residual)
    ap = abs(previous_residual)

    # Strong sign reversal = overshooting.
    if residual * previous_residual < 0 && a > good_ratio * ap
        return max(alpha_min, alpha * shrink), "sign_flip"
    end

    # Residual became larger.
    if a > ap
        return max(alpha_min, alpha * shrink), "residual_grew"
    end

    # Clearly improving.
    if a < good_ratio * ap
        return min(alpha_max, alpha * grow), "good_decrease"
    end

    return alpha, "keep"
end

function _adapt_local_alpha_vec(alpha, residual, previous_residual;
                                alpha_min=1e-6,
                                alpha_max=0.05,
                                shrink=0.5,
                                grow=1.05,
                                good_ratio=0.8)

    previous_residual === nothing && return alpha, "initial"

    nr = norm(residual)
    np = norm(previous_residual)

    if dot(residual, previous_residual) < 0 && nr > good_ratio * np
        return max(alpha_min, alpha * shrink), "sign_flip"
    end

    if nr > np
        return max(alpha_min, alpha * shrink), "residual_grew"
    end

    if nr < good_ratio * np
        return min(alpha_max, alpha * grow), "good_decrease"
    end

    return alpha, "keep"
end

# ============================================================
# ONE-SITE constrained sweep
# ============================================================

function dmrg_sweep_one_site_constrained!(
    state::MPSState,
    solver::LanczosSolver,
    options::DMRGOptions,
    direction::Symbol,
    lambda::Vector{Float64},
    target::Vector{Float64};
    local_max_iter::Int=500,
    local_tol::Float64=1e-4,
    alpha_init::Float64=0.05,
    alpha_min::Float64=1e-6,
    alpha_max::Float64=0.05,
    adaptive::Bool=true,
    sweep_number::Int=0,
    trace_io=nothing,
    verbose::Bool=false,
    site_visit_callback=nothing,
    local_iteration_callback=nothing
)
    N = length(state.mps.tensors)
    Sz = Matrix{Float64}(real.(spin_ops(size(state.mps.tensors[1],2))[:Z]))

    last_energy = NaN
    max_chi = 0
    max_local_used = 0
    failed_sites = 0

    sites = direction == :right ? (1:N-1) : (N:-1:2)

    for site in sites
        chi_l, d, chi_r = size(state.mps.tensors[site])

        # IMPORTANT:
        # These environments remain fixed throughout the local lambda solve.
        if direction == :right
            left_env = site == 1 ?
                state.environment.tensors[N+1] :
                state.environment.tensors[site-1]
            right_env = state.environment.tensors[site+1]
        else
            left_env = state.environment.tensors[site-1]
            right_env = state.environment.tensors[site+1]
        end

        alpha = alpha_init
        previous_residual = nothing
        local_converged = false
        local_used = 0
        local_tensor = state.mps.tensors[site]
        final_local_measured = NaN
        final_local_residual = NaN

        for local_iter in 1:local_max_iter
            local_used = local_iter

            # Rebuild Heff every iteration because lambda_i changed,
            # hence the MPO tensor at this site changed.
            Heff = OneSiteEffectiveHamiltonian(
                left_env,
                state.mpo.tensors[site],
                right_env
            )

            evec, eval = _solve(
                solver,
                Heff,
                vec(local_tensor)
            )

            last_energy = real(eval)
            local_tensor = reshape(evec, (chi_l,d,chi_r))

            measured = _local_sz_one(local_tensor, Sz)
            residual = target[site] - measured
            final_local_measured = measured
            final_local_residual = residual

            alpha_action = "keep"

            if adaptive
                alpha, alpha_action = _adapt_local_alpha(
                    alpha,
                    residual,
                    previous_residual;
                    alpha_min=alpha_min,
                    alpha_max=alpha_max
                )
            end

            if trace_io !== nothing
                println(trace_io, join((
                    sweep_number,
                    direction,
                    site,
                    local_iter,
                    target[site],
                    measured,
                    residual,
                    lambda[site],
                    alpha,
                    alpha_action,
                    last_energy
                ), ','))
                flush(trace_io)
            end

            # Optional hook for diagnosing the nonlocal response caused by
            # every individual lambda_i update. `local_tensor` is the state
            # optimized at the current lambda_i; the callback may temporarily
            # insert it into the MPS to measure global observables, but must
            # restore the stored tensor before returning.
            if local_iteration_callback !== nothing
                is_converged = abs(residual) <= local_tol
                will_update_lambda = !is_converged && local_iter < local_max_iter
                lambda_next = will_update_lambda ?
                    lambda[site] + alpha * residual : lambda[site]
                local_iteration_callback(state, (
                    direction=direction,
                    site=site,
                    local_iteration=local_iter,
                    target=target[site],
                    local_measured=measured,
                    local_residual=residual,
                    lambda_used=lambda[site],
                    lambda_next=lambda_next,
                    alpha=alpha,
                    alpha_action=alpha_action,
                    local_energy=last_energy,
                    local_converged=is_converged,
                    will_update_lambda=will_update_lambda,
                ), local_tensor, lambda, target)
            end

            if verbose && (
                local_iter == 1 ||
                local_iter % 50 == 0 ||
                abs(residual) <= local_tol
            )
                println(
                    "sweep=$sweep_number $direction site=$site ",
                    "local=$local_iter ",
                    "target=$(round(target[site],digits=6)) ",
                    "Sz=$(round(measured,digits=6)) ",
                    "r=$(round(residual,digits=6)) ",
                    "lambda=$(round(lambda[site],digits=6)) ",
                    "alpha=$(round(alpha,sigdigits=4))"
                )
            end

            # Local constraint solved.
            if abs(residual) <= local_tol
                local_converged = true
                break
            end

            # Do NOT change lambda after the final allowed solve,
            # otherwise lambda and local_tensor would be inconsistent.
            if local_iter == local_max_iter
                break
            end

            previous_residual = residual

            # Update only lambda_i.
            lambda[site] += alpha * residual

            # Only the current on-site field changed. Keep the captured
            # left/right environments frozen during this local solve.
            set_lambda_site!(state, site, lambda[site])
        end

        max_local_used = max(max_local_used, local_used)

        if !local_converged
            failed_sites += 1
        end

        # Use the locally optimized tensor.
        state.mps.tensors[site] = local_tensor

        # ----------------------------------------------------
        # Only NOW move the orthogonality center.
        # ----------------------------------------------------

        if direction == :right
            F = qr(reshape(local_tensor, (chi_l*d, chi_r)))

            new_chi = size(F.R,1)
            Q_full = Matrix(F.Q)
            Q = size(Q_full,2) > new_chi ?
                Q_full[:,1:new_chi] :
                Q_full
            R = Matrix(F.R)

            state.mps.tensors[site] =
                reshape(Q, (chi_l,d,new_chi))

            _update_left_environment(state, site)

            @tensoropt new_next[-1,-2,-3] :=
                R[-1,4] *
                state.mps.tensors[site+1][4,-2,-3]

            state.mps.tensors[site+1] = new_next
            state.environment.tensors[site+1] = nothing
            state.center = site + 1

            max_chi = max(max_chi,new_chi)

        else
            F = lq(reshape(local_tensor, (chi_l,d*chi_r)))

            new_chi = size(F.L,2)
            L = Matrix(F.L)
            Q_full = Matrix(F.Q)
            Q = size(Q_full,1) > new_chi ?
                Q_full[1:new_chi,:] :
                Q_full

            state.mps.tensors[site] =
                reshape(Q, (new_chi,d,chi_r))

            _update_right_environment(state, site)

            @tensoropt new_prev[-1,-2,-3] :=
                state.mps.tensors[site-1][-1,-2,4] *
                L[4,-3]

            state.mps.tensors[site-1] = new_prev
            state.environment.tensors[site-1] = nothing
            state.center = site - 1

            max_chi = max(max_chi,new_chi)
        end

        # Optional diagnostic hook. At this point the optimized tensor has
        # been installed, the orthogonality center has moved, and the MPS
        # represents the post-visit physical state. Ordinary runs leave this
        # as `nothing` and therefore pay no measurement/storage overhead.
        if site_visit_callback !== nothing
            site_visit_callback(state, (
                direction=direction,
                site=site,
                local_iterations=local_used,
                local_converged=local_converged,
                target=target[site],
                local_measured=final_local_measured,
                local_residual=final_local_residual,
                lambda=lambda[site],
                local_energy=last_energy,
            ), lambda, target)
        end
    end

    return (
        E=last_energy,
        max_chi=max_chi,
        max_local_iterations=max_local_used,
        failed_sites=failed_sites
    )
end

# ============================================================
# TWO-SITE constrained sweep
#
# At each two-site block, optimize BOTH local lambdas until
#
# <Sz_left>  = target_left
# <Sz_right> = target_right
#
# or local_max_iter is reached.
# ============================================================

function dmrg_sweep_two_site_constrained!(
    state::MPSState,
    solver::LanczosSolver,
    options::DMRGOptions,
    direction::Symbol,
    lambda::Vector{Float64},
    target::Vector{Float64};
    local_max_iter::Int=500,
    local_tol::Float64=1e-4,
    alpha_init::Float64=0.05,
    alpha_min::Float64=1e-6,
    alpha_max::Float64=0.05,
    adaptive::Bool=true,
    sweep_number::Int=0,
    trace_io=nothing,
    verbose::Bool=false
)
    N = length(state.mps.tensors)
    Sz = Matrix{Float64}(real.(spin_ops(size(state.mps.tensors[1],2))[:Z]))

    last_energy = NaN
    max_chi = 0
    max_local_used = 0
    failed_blocks = 0

    sites = direction == :right ? (1:N-1) : (N:-1:2)

    for site in sites

        # Define physical two-site block.
        if direction == :right
            i = site
            j = site + 1
        else
            i = site - 1
            j = site
        end

        @tensoropt theta[-1,-2,-3,-4] :=
            state.mps.tensors[i][-1,-2,5] *
            state.mps.tensors[j][5,-3,-4]

        chi_l,d1,d2,chi_r = size(theta)

        # Freeze environments outside this two-site block.
        if direction == :right
            left_env = i == 1 ?
                state.environment.tensors[N+1] :
                state.environment.tensors[i-1]

            right_env = state.environment.tensors[j+1]
        else
            left_env = i == 1 ?
                state.environment.tensors[N+1] :
                state.environment.tensors[i-1]

            right_env = state.environment.tensors[j+1]
        end

        alpha = alpha_init
        previous_residual = nothing
        local_converged = false
        local_used = 0
        local_theta = theta

        for local_iter in 1:local_max_iter
            local_used = local_iter

            # Rebuild because lambda_i and/or lambda_j may have changed.
            Heff = TwoSiteEffectiveHamiltonian(
                left_env,
                state.mpo.tensors[i],
                state.mpo.tensors[j],
                right_env
            )

            evec,eval = _solve(
                solver,
                Heff,
                vec(local_theta)
            )

            last_energy = real(eval)
            local_theta =
                reshape(evec,(chi_l,d1,d2,chi_r))

            m1,m2 = _local_sz_two(local_theta,Sz)

            residual = [
                target[i] - m1,
                target[j] - m2
            ]

            alpha_action = "keep"

            if adaptive
                alpha,alpha_action = _adapt_local_alpha_vec(
                    alpha,
                    residual,
                    previous_residual;
                    alpha_min=alpha_min,
                    alpha_max=alpha_max
                )
            end

            if trace_io !== nothing
                println(trace_io, join((
                    sweep_number,
                    direction,
                    i,
                    j,
                    local_iter,
                    target[i],
                    m1,
                    residual[1],
                    lambda[i],
                    target[j],
                    m2,
                    residual[2],
                    lambda[j],
                    alpha,
                    alpha_action,
                    last_energy
                ), ','))
                flush(trace_io)
            end

            maxres = maximum(abs,residual)

            if verbose && (
                local_iter == 1 ||
                local_iter % 50 == 0 ||
                maxres <= local_tol
            )
                println(
                    "sweep=$sweep_number $direction block=($i,$j) ",
                    "local=$local_iter ",
                    "Sz=($(round(m1,digits=6)),",
                    "$(round(m2,digits=6))) ",
                    "r=($(round(residual[1],digits=6)),",
                    "$(round(residual[2],digits=6))) ",
                    "lambda=($(round(lambda[i],digits=6)),",
                    "$(round(lambda[j],digits=6))) ",
                    "alpha=$(round(alpha,sigdigits=4))"
                )
            end

            if maxres <= local_tol
                local_converged = true
                break
            end

            if local_iter == local_max_iter
                break
            end

            previous_residual = copy(residual)

            # Update the TWO local fields simultaneously.
            lambda[i] += alpha * residual[1]
            lambda[j] += alpha * residual[2]

            set_lambda!(state,lambda)
        end

        max_local_used = max(max_local_used,local_used)

        if !local_converged
            failed_blocks += 1
        end

        # ----------------------------------------------------
        # Only after local lambda optimization: SVD / move OC
        # ----------------------------------------------------

        U,S,V = _svd_truncate(
            reshape(local_theta,(chi_l*d1,d2*chi_r)),
            options.chi_max,
            options.cutoff
        )

        new_chi = length(S)
        max_chi = max(max_chi,new_chi)

        if direction == :right
            state.mps.tensors[i] =
                reshape(U,(chi_l,d1,new_chi))

            @tensoropt state.mps.tensors[j][-1,-2,-3] :=
                Diagonal(S)[-1,4] *
                reshape(V,(new_chi,d2,chi_r))[4,-2,-3]

            state.center = j
            _update_left_environment(state,i)

        else
            state.mps.tensors[j] =
                reshape(V,(new_chi,d2,chi_r))

            @tensoropt state.mps.tensors[i][-1,-2,-3] :=
                reshape(U,(chi_l,d1,new_chi))[-1,-2,4] *
                Diagonal(S)[4,-3]

            state.center = i
            _update_right_environment(state,j)
        end
    end

    return (
        E=last_energy,
        max_chi=max_chi,
        max_local_iterations=max_local_used,
        failed_blocks=failed_blocks
    )
end

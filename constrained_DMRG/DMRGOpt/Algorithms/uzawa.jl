# Algorithms/uzawa.jl
#
# Uzawa (per-sweep) updates of Lagrange multipliers for
#   min <H>  s.t.  <Sz_i> = m_i.
#
# H_eff = H - Σ λ_i Sz_i is injected as Zeeman fields h = -λ.
#
# Global Uzawa:
#   λ <- λ + η (m - <Sz>)
#
# For local constrained DMRG:
#   set_lambda_site! changes ONLY the current site's MPO field
#   and does NOT rebuild environments.

"""
    rebuild_environment!(state)

Rebuild left/right environments from the current MPS, MPO, and center.
Call after a global MPO-field update when cached environments are stale.
"""
function rebuild_environment!(state::MPSState)
    state.environment = _build_environment(state.mps, state.mpo, state.center)
    return state
end

"""
    set_lambda!(state, λ)

Write h = -λ into the entire XXZ MPO and rebuild environments.

Use this for GLOBAL Uzawa updates between complete DMRG/TDVP sweeps.
Do NOT use this inside the local constrained-DMRG inner loop.
"""
function set_lambda!(state::MPSState, λ::AbstractVector{<:Real})
    update_fields!(state.mpo, Vector{Float64}(-λ))
    rebuild_environment!(state)
    return state
end

"""
    set_lambda_site!(state, site, λi)

Update ONLY the local field at one site.

Convention:
    H_eff = H - Σ_i λ_i Sz_i
so the XXZ MPO field is
    h_i = -λ_i.

No environments are rebuilt.
No other MPO tensor is changed.

This is the function to use inside local constrained DMRG.
"""
function set_lambda_site!(
    state::MPSState,
    site::Int,
    λi::Real
)
    N = length(state.mpo.tensors)
    @assert 1 <= site <= N "site must satisfy 1 <= site <= $N"

    W = state.mpo.tensors[site]
    d = size(W, 3)
    Sz = Matrix{Float64}(spin_ops(d)[:Z])

    h = -Float64(λi)

    if site == 1
        W[1,1,:,:] = h * Sz
    elseif site == N
        W[end,1,:,:] = h * Sz
    else
        W[end,1,:,:] = h * Sz
    end

    return state
end

"""
    sz_profile(mps; d=2)

Return <Sz_i> at all sites.
"""
function sz_profile(mps::MPS; d::Int=2)
    Sz = Matrix{Float64}(real.(spin_ops(d)[:Z]))
    nrm = real(measure_norm(mps))
    N = length(mps.tensors)

    return [
        real(measure_local_observable(mps, Sz, i)) / nrm
        for i in 1:N
    ]
end

"""
    physical_energy(state, λ)

Recover the physical XXZ energy from

    H_eff = H - Σ_i λ_i Sz_i.

Therefore

    E_phys = E_eff + Σ_i λ_i <Sz_i>.
"""
function physical_energy(
    state::MPSState,
    λ::AbstractVector{<:Real}
)
    nrm = real(measure_norm(state.mps))
    E_eff = real(measure_energy(state)) / nrm
    sz = sz_profile(state.mps)

    return E_eff + dot(λ, sz)
end

"""
    uzawa_update!(λ, sz, m, η)

In-place dual step

    λ <- λ + η (m - sz)

Returns ||sz - m||_2.
"""
function uzawa_update!(
    λ::Vector{Float64},
    sz::Vector{Float64},
    m::Vector{Float64},
    η::Float64
)
    g = sz .- m
    λ .+= η .* (m .- sz)
    return norm(g)
end

"""
    set_delta_and_lambda!(state, Jxy, Δ, λ; d=2)

Replace the XXZ MPO with a new Δ while keeping h = -λ,
then rebuild environments.
"""
function set_delta_and_lambda!(
    state::MPSState,
    Jxy::Real,
    Δ::Real,
    λ::AbstractVector{<:Real};
    d::Int=2
)
    N = length(state.mpo.tensors)

    mpo = build_xxz_mpo(
        N,
        Jxy,
        Δ;
        d=d,
        h=Vector{Float64}(-λ)
    )

    for i in 1:N
        state.mpo.tensors[i] = mpo.tensors[i]
    end

    rebuild_environment!(state)
    return state
end

"""
    dmrg_pair!(state, solver, opts)

One complete right + left two-site DMRG pair.
"""
function dmrg_pair!(
    state::MPSState,
    solver::LanczosSolver,
    opts::DMRGOptions
)
    dmrg_sweep(state, solver, opts, :right)
    dmrg_sweep(state, solver, opts, :left)
    return state
end

"""
    tdvp_pair!(state, solver, opts)

One complete right + left two-site TDVP pair.
"""
function tdvp_pair!(
    state::MPSState,
    solver::KrylovExponential,
    opts::TDVPOptions
)
    tdvp_sweep(state, solver, opts, :right)
    tdvp_sweep(state, solver, opts, :left)
    return state
end

"""
    tdvp_pair_one_site!(state, solver, opts)

One complete right + left one-site TDVP pair.
"""
function tdvp_pair_one_site!(
    state::MPSState,
    solver::KrylovExponential,
    opts::TDVPOptions
)
    tdvp_sweep_one_site(state, solver, opts, :right)
    tdvp_sweep_one_site(state, solver, opts, :left)
    return state
end

"""
    uzawa_until!(state, λ, φ, sweep!;
                 η=0.25,
                 maxiter=40,
                 tol=5e-3,
                 n_inner=1)

Repeat n_inner complete sweeps, then perform one global Uzawa update:

    λ <- λ + η/sqrt(k) * (φ - <Sz>)

After the global λ update, the MPO fields and environments are rebuilt.

Stops when

    ||<Sz> - φ||_2 < tol.
"""
function uzawa_until!(
    state::MPSState,
    λ::Vector{Float64},
    φ::Vector{Float64},
    sweep!;
    η::Float64=0.25,
    maxiter::Int=40,
    tol::Float64=5e-3,
    n_inner::Int=1
)
    g = Inf

    for k in 1:maxiter
        for _ in 1:n_inner
            sweep!(state)
        end

        sz = sz_profile(state.mps)

        g = uzawa_update!(
            λ,
            sz,
            φ,
            η / sqrt(k)
        )

        set_lambda!(state, λ)

        g < tol && return g
    end

    return g
end

"""
    m_stag(sz)

Staggered magnetization

    m_stag = (1/N) Σ_i (-1)^(i+1) Sz_i.
"""
m_stag(sz::AbstractVector) =
    sum((-1.0)^(i+1) * sz[i] for i in eachindex(sz)) / length(sz)

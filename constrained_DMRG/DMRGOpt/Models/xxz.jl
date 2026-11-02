# src/xxz.jl
#
# Hand-built MPO for the spin-S XXZ chain:
#
#   H = Σᵢ [ (Jxy/2)(S⁺ᵢS⁻ᵢ₊₁ + S⁻ᵢS⁺ᵢ₊₁) + Jz Sᶻᵢ Sᶻᵢ₊₁ ] + Σᵢ hᵢ Sᶻᵢ
#
# MPO bond dimension D = 5.  Tensor indices: W[α_L, α_R, σ, σ']
#
# Bulk W matrix (rows = left bond α, cols = right bond β):
#
#   α\β │  1       2          3          4       5
#   ────┼──────────────────────────────────────────────
#    1  │  I       0          0          0       0
#    2  │  S⁺      0          0          0       0
#    3  │  S⁻      0          0          0       0
#    4  │  Sᶻ      0          0          0       0
#    5  │  hᵢSᶻ   (Jxy/2)S⁻  (Jxy/2)S⁺  Jz·Sᶻ   I
#
# Left boundary (site 1): row α=5 → shape (1, D, d, d)
# Right boundary (site N): col β=1 → shape (D, 1, d, d)

"""
    build_xxz_mpo(N, Jxy, Jz; d=2, h=zeros(N))

Build the MPO for an N-site XXZ chain with local dimension `d = 2S+1`.
The vector `h` provides site-local Zeeman fields — inject Lagrange
multiplier fields here for constrained DMRG: set h[i] = −λᵢ.
"""
function build_xxz_mpo(N::Int, Jxy::Real, Jz::Real;
                       d::Int=2, h::Vector{Float64}=zeros(Float64, N))
    @assert length(h) == N "h must have length N=$N"
    D = 5

    ops = spin_ops(d)
    Sz = Matrix{Float64}(ops[:Z])
    Sp = Matrix{Float64}(ops[:Sp])
    Sm = Matrix{Float64}(ops[:Sm])
    Id = Matrix{Float64}(ops[:I])

    tensors = Vector{Array{Float64,4}}(undef, N)

    for i in 1:N
        if i == 1
            # Left boundary: (1, D, d, d)
            W = zeros(Float64, 1, D, d, d)
            W[1, 1, :, :] = h[i] * Sz
            W[1, 2, :, :] = (Jxy / 2) * Sm
            W[1, 3, :, :] = (Jxy / 2) * Sp
            W[1, 4, :, :] = Jz * Sz
            W[1, 5, :, :] = Id
        elseif i == N
            # Right boundary: (D, 1, d, d)
            W = zeros(Float64, D, 1, d, d)
            W[1, 1, :, :] = Id
            W[2, 1, :, :] = Sp
            W[3, 1, :, :] = Sm
            W[4, 1, :, :] = Sz
            W[5, 1, :, :] = h[i] * Sz
        else
            # Bulk: (D, D, d, d)
            W = zeros(Float64, D, D, d, d)
            W[1, 1, :, :] = Id
            W[2, 1, :, :] = Sp
            W[3, 1, :, :] = Sm
            W[4, 1, :, :] = Sz
            W[5, 1, :, :] = h[i] * Sz
            W[5, 2, :, :] = (Jxy / 2) * Sm
            W[5, 3, :, :] = (Jxy / 2) * Sp
            W[5, 4, :, :] = Jz * Sz
            W[5, 5, :, :] = Id
        end
        tensors[i] = W
    end

    return MPO(tensors)
end

"""
    update_fields!(mpo::MPO, h::Vector{Float64})

Replace all on-site Zeeman fields in an XXZ MPO in-place.
For constrained DMRG, call this with h[i] = −λᵢ after each
multiplier update — no need to rebuild the full MPO.
"""
function update_fields!(mpo::MPO, h::Vector{Float64})
    N = length(mpo.tensors)
    d = size(mpo.tensors[1], 3)
    Sz = Matrix{Float64}(spin_ops(d)[:Z])

    for i in 1:N
        if i == 1
            mpo.tensors[1][1, 1, :, :] = h[i] * Sz
        elseif i == N
            mpo.tensors[N][end, 1, :, :] = h[i] * Sz
        else
            mpo.tensors[i][end, 1, :, :] = h[i] * Sz
        end
    end
end

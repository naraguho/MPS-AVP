# src/random_mps.jl
#
# Random MPS initialization for DMRG.
# Bond dimensions respect exact Hilbert space growth from both edges,
# capped at the target χ.

"""
    random_mps(N, d, chi; T=Float64)

Create a normalized random MPS on `N` sites with local dimension `d`
and maximum bond dimension `chi`.

Bond dimensions grow as dⁱ from the left and d^(N−i) from the right,
capped at `chi`, so the MPS never wastes parameters on linearly
dependent components.

Returns an `MPS{T}` ready for `MPSState(mps, mpo)`.
"""
function random_mps(N::Int, d::Int, chi::Int; T::Type=Float64)
    @assert N ≥ 2 "Need at least 2 sites"
    @assert d ≥ 2 "Local dimension must be ≥ 2"
    @assert chi ≥ 1 "Bond dimension must be ≥ 1"

    # Compute bond dimensions: bonds[i] sits between site i-1 and site i
    # bonds[1] = left boundary = 1, bonds[N+1] = right boundary = 1
    bonds = ones(Int, N + 1)
    for i in 2:N
        bonds[i] = min(chi, bonds[i-1] * d)
    end
    for i in N:-1:2
        bonds[i] = min(bonds[i], bonds[i+1] * d)
    end

    tensors = Vector{Array{T,3}}(undef, N)
    for i in 1:N
        tensors[i] = randn(T, bonds[i], d, bonds[i+1])
    end

    mps = MPS(tensors)

    # Normalize
    nrm = sqrt(abs(measure_norm(mps)))
    if nrm > 0
        tensors[1] ./= nrm
    end

    return mps
end

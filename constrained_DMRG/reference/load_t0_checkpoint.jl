"""Load the bundled accepted t=0 checkpoint and its DMRGOpt MPS type.

Usage from the repository root:
    include("constrained_DMRG/reference/load_t0_checkpoint.jl")
    data = load_t0_checkpoint()
    data.mps, data.lambda, data.target, data.measured
"""

include(joinpath(@__DIR__, "..", "DMRGOpt", "DMRGOpt.jl"))
using .DMRGOpt
using Serialization

function load_t0_checkpoint(path=joinpath(@__DIR__, "one_kink_t0_mps.jls"))
    checkpoint = deserialize(path)
    # Earlier run checkpoints used `converged`; the bundled pipeline writes `accepted`.
    accepted = hasproperty(checkpoint, :accepted) ? checkpoint.accepted : checkpoint.converged
    return merge(checkpoint, (; accepted, mps=MPS(checkpoint.tensors)))
end

# Core/states.jl
# MPSState struct is defined in Core/types.jl

function MPSState(mps::MPS{Tmps}, mpo::MPO{Tmpo}; center=1) where {Tmps,Tmpo}
    # Environment type is the promotion of MPS and MPO types
    Tenv = promote_type(Tmps, Tmpo)

    # NO CONVERSIONS OR COPIES! Just use the inputs directly
    make_canonical(mps, center)  # Modifies in-place

    # Build environment with natural type promotion
    env = _build_environment(mps, mpo, center)

    return MPSState{Tmps,Tmpo,Tenv}(mps, mpo, env, center)
end

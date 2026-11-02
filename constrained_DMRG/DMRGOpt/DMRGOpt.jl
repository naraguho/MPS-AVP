module DMRGOpt

using LinearAlgebra
using TensorOperations

export
    # Core types
    TensorNetwork,
    MPS,
    MPO,
    MPDO,
    Environment,
    MPSState,
    DMRGOptions,
    TDVPOptions,
    SweepSchedule,

    # Sites
    SpinSite,
    SuperSite,
    spin_ops,

    # Measurements
    measure_energy,
    measure_norm,
    measure_local_observable,
    measure_correlation,
    mps_to_vector,
    mps_overlap,
    pad_mps,

    # Canonicalization
    make_canonical,
    is_left_orthogonal,
    is_right_orthogonal,
    is_orthogonal,

    # Decomposition utilities
    entropy,
    truncation_error,

    # Solvers
    EffectiveHamiltonian,
    OneSiteEffectiveHamiltonian,
    TwoSiteEffectiveHamiltonian,
    ZeroSiteEffectiveHamiltonian,
    LanczosSolver,
    KrylovExponential,

    # DMRG
    dmrg_sweep,
    dmrg_sweep_one_site,
    dmrg_sweep_one_site_constrained!,
    dmrg_sweep_two_site_constrained!,

    # TDVP
    tdvp_sweep,
    tdvp_sweep_one_site,

    # Uzawa / constraints
    rebuild_environment!,
    set_lambda!,
    set_lambda_site!,
    sz_profile,
    physical_energy,
    uzawa_update!,
    set_delta_and_lambda!,
    dmrg_pair!,
    tdvp_pair!,
    tdvp_pair_one_site!,
    uzawa_until!,
    m_stag,

    # Models
    build_xxz_mpo,
    update_fields!,
    random_mps,

    # ED
    embed_operator,
    embed_two_site,
    embed_multi_site,
    basis_state,
    neel_state,
    build_xxz_hamiltonian,
    EDEigensystem,
    diagonalize,
    ground_state,
    low_energy_states,
    ed_expectation,
    ed_local_expectation,
    ed_local_profile,
    ed_time_evolve

# Include files in dependency order
include("Core/types.jl")
include("Core/site.jl")
include("TensorOps/canonicalization.jl")
include("TensorOps/decomposition.jl")
include("TensorOps/environment.jl")
include("TensorOps/measurements.jl")
include("TensorOps/padding.jl")
include("Core/states.jl")
include("Models/random_mps.jl")
include("Models/xxz.jl")
include("Algorithms/solvers.jl")
include("Algorithms/dmrg.jl")
include("Algorithms/tdvp.jl")
include("Algorithms/uzawa.jl")
include("ED/operators.jl")
include("ED/hamiltonian.jl")
include("ED/eigensolver.jl")
include("ED/observables.jl")

end # module DMRGOpt

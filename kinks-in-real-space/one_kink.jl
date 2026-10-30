#!/usr/bin/env julia

# Zero-field one-kink trial at bond 50|51. Run from the repository root:
# julia --startup-file=no --threads=1 --project=. kinks-in-real-space/one_kink.jl
# This file includes only the QuenchDyn primitives it needs from ../src.

using KrylovKit, LinearAlgebra, Printf, Serialization, Statistics, TensorOperations

module QD
using LinearAlgebra, TensorOperations
const SOURCE = normpath(joinpath(@__DIR__, "..", "src"))
for file in ("Core/types.jl", "Core/site.jl", "Core/states.jl",
             "TensorOps/decomposition.jl", "TensorOps/canonicalization.jl",
             "TensorOps/environment.jl", "TensorOps/measurements.jl",
             "Algorithms/solvers.jl", "Algorithms/dmrg.jl",
             "Builders/xxz_builder.jl")
    include(joinpath(SOURCE, file))
end
end

const N = 100
const DELTA = 3.0
const CHI = 64
const KINK_BOND = 50
const CORE_L, CORE_R = 45, 56
const BACKGROUND_SWEEPS = parse(Int, get(ENV, "KINK_BACKGROUND_SWEEPS", "6"))
const MAX_ROUNDS = parse(Int, get(ENV, "KINK_MAX_ROUNDS", "30"))
const OUT = get(ENV, "KINK_OUTPUT", joinpath(@__DIR__, "outputs", "one-kink"))
const LOCAL_ALG = KrylovKit.Lanczos(krylovdim=32,maxiter=30,tol=1e-10,verbosity=0)

function energy(psi, mpo)
    left = ones(Float64, 1, 1, 1)
    for i in 1:N
        left = QD._contract_left_environment(left, psi.tensors[i], mpo.tensors[i])
    end
    return real(left[1,1,1]) / real(QD.measure_norm(psi))
end

function profile(psi)
    # Transfer environments give the entire profile in two passes.
    right = Vector{Matrix{Float64}}(undef, N+1)
    right[N+1] = ones(1,1)
    for i in N:-1:1
        a = psi.tensors[i]
        right[i] = zeros(size(a,1), size(a,1))
        for s in 1:2
            block = @view a[:,s,:]
            right[i] .+= block*right[i+1]*transpose(block)
        end
    end
    left = ones(1,1)
    m = zeros(N)
    for i in 1:N
        a = psi.tensors[i]
        up, down = @view(a[:,1,:]), @view(a[:,2,:])
        lu, ld = transpose(up)*left*up, transpose(down)*left*down
        m[i] = 0.5*sum((lu-ld).*right[i+1])
        left = lu+ld
    end
    return m ./ left[1,1]
end

staggered(m, sites) = mean((isodd(i) ? 1.0 : -1.0)*m[i] for i in sites)

function background(mpo)
    tensors = [zeros(Float64,1,2,1) for _ in 1:N]
    for i in 1:N
        tensors[i][1,isodd(i) ? 1 : 2,1] = 1.0
    end
    # The product seed chooses branch A. Every sweep uses zero field.
    state = QD.MPSState(QD.MPS{Float64}(tensors),mpo;center=1)
    opts = QD.DMRGOptions(CHI,1e-10,2)
    solver = QD.LanczosSolver(3,12)
    for sweep in 1:BACKGROUND_SWEEPS
        QD.dmrg_sweep(state,solver,opts,isodd(sweep) ? :right : :left)
        @printf("background sweep=%2d E_XXZ=%+.12f\n", sweep, energy(state.mps,mpo))
        flush(stdout)
    end
    m = profile(state.mps)
    staggered(m,10:35)>0.3 && staggered(m,65:90)>0.3 ||
        error("The zero-field background did not retain one Néel branch")
    return state.mps
end

function solve_site(state, i)
    left = state.environment.tensors[i-1]
    right = state.environment.tensors[i+1]
    heff = QD.OneSiteEffectiveHamiltonian(left,state.mpo.tensors[i],right)
    dims = size(state.mps.tensors[i])
    values, vectors, _ = KrylovKit.eigsolve(
        v -> QD._apply(heff,v),vec(state.mps.tensors[i]),1,:SR,LOCAL_ALG)
    x = vectors[1]
    residual = norm(QD._apply(heff,x)-real(values[1])*x)/max(1,abs(real(values[1])))
    residual < 1e-8 || error("Local eigensolver residual $residual at site $i")
    return reshape(x,dims)
end

function core_round!(state)
    state.center == CORE_L || error("Core round must begin at site $CORE_L")
    for i in CORE_L:CORE_R-1
        tensor = solve_site(state,i)
        dl,d,dr = size(tensor)
        F = qr(reshape(tensor,dl*d,dr))
        chi = size(F.R,1)
        state.mps.tensors[i] = reshape(Matrix(F.Q)[:,1:chi],dl,d,chi)
        QD._update_left_environment(state,i)
        next = state.mps.tensors[i+1]
        state.mps.tensors[i+1] = reshape(Matrix(F.R)*reshape(next,dr,:),
                                         chi,size(next,2),size(next,3))
        state.environment.tensors[i+1] = nothing
        state.center = i+1
    end
    state.mps.tensors[CORE_R] = solve_site(state,CORE_R)
    for i in CORE_R:-1:CORE_L+1
        tensor = solve_site(state,i)
        dl,d,dr = size(tensor)
        F = lq(reshape(tensor,dl,d*dr))
        chi = size(F.L,2)
        state.mps.tensors[i] = reshape(Matrix(F.Q)[1:chi,:],chi,d,dr)
        QD._update_right_environment(state,i)
        previous = state.mps.tensors[i-1]
        state.mps.tensors[i-1] = reshape(reshape(previous,:,dl)*Matrix(F.L),
                                         size(previous,1),size(previous,2),chi)
        state.environment.tensors[i-1] = nothing
        state.center = i-1
    end
    state.mps.tensors[CORE_L] = solve_site(state,CORE_L)
end

function save_profile(io,round,m)
    for i in 1:N
        @printf(io,"%d,%d,%.17g\n",round,i,m[i])
    end
    flush(io)
end

const PLOT = raw"""
import csv, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows = list(csv.DictReader(open(sys.argv[1])))
round_id = max(int(r['round']) for r in rows)
last = sorted((r for r in rows if int(r['round']) == round_id),
              key=lambda r: int(r['site']))
x = [int(r['site']) for r in last]
m = [float(r['Sz']) for r in last]
fig, ax = plt.subplots(figsize=(11,4), layout='constrained')
ax.plot(x,m,'o-',color='#de6b16',ms=2.5,lw=1.4,label=r'$\langle S_i^z\rangle$')
ax.axvspan(45,56,color='#ffd468',alpha=.23,label='optimized core')
ax.axvline(50.5,color='#444',ls='--',lw=1)
ax.axhline(0,color='#777',lw=.7)
ax.set(xlim=(1,100),ylim=(-.52,.52),xlabel='site $i$',
       ylabel=r'$\langle S_i^z\rangle$',
       title=r'One kink: $N=100$, $\Delta=3$, zero field')
ax.legend(loc='upper right',frameon=False)
fig.savefig(sys.argv[2],dpi=190)
"""

function main()
    BLAS.set_num_threads(8)
    isdir(OUT) && !isempty(readdir(OUT)) && error("Output exists: $OUT")
    mkpath(OUT)
    mpo = QD.build_xxz_mpo(N,1.0,DELTA)
    A = background(mpo)
    psi = QD.MPS{Float64}(deepcopy(A.tensors))
    QD.make_canonical(psi,CORE_L)
    for i in KINK_BOND+1:N
        psi.tensors[i] = psi.tensors[i][:,2:-1:1,:]
    end
    state = QD.MPSState(psi,mpo;center=CORE_L)
    frozen = deepcopy(state.mps.tensors)
    old_e,old_m = energy(state.mps,mpo),profile(state.mps)
    stable = 0
    converged = false
    open(joinpath(OUT,"profiles.csv"),"w") do profiles
        println(profiles,"round,site,Sz")
        save_profile(profiles,0,old_m)
        open(joinpath(OUT,"round_summary.csv"),"w") do summary
            println(summary,"round,energy,dE,max_dSz,total_Sz,left_staggered,right_staggered,Sz50,Sz51")
            @printf(summary,"0,%.17g,NaN,NaN,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    old_e,sum(old_m),staggered(old_m,10:35),
                    staggered(old_m,65:90),old_m[50],old_m[51])
            for round in 1:MAX_ROUNDS
                core_round!(state)
                all(state.mps.tensors[i] == frozen[i] for i in 1:N
                    if i<CORE_L || i>CORE_R) || error("A frozen tensor changed")
                m = profile(state.mps)
                e = energy(state.mps,mpo)
                de,dm = e-old_e,maximum(abs.(m.-old_m))
                left,right = staggered(m,10:35),staggered(m,65:90)
                @printf("round=%2d E=%+.12f dE=%+.2e max|dSz|=%.2e LR=(%+.4f,%+.4f) core=(%+.5f,%+.5f)\n",
                        round,e,de,dm,left,right,m[50],m[51])
                @printf(summary,"%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                        round,e,de,dm,sum(m),left,right,m[50],m[51])
                flush(stdout); flush(summary)
                save_profile(profiles,round,m)
                left*right < 0 && min(abs(left),abs(right))>0.2 ||
                    error("The optimized state lost its one-kink domains")
                stable = abs(de)<1e-10 && dm<1e-8 ? stable+1 : 0
                old_e,old_m = e,m
                if stable>=3
                    converged = true
                    break
                end
            end
        end
    end
    serialize(joinpath(OUT,"one_kink_mps.jls"),
              (;N,Delta=DELTA,chi=CHI,core=(CORE_L,CORE_R),
                energy=old_e,converged,tensors=deepcopy(state.mps.tensors)))
    @printf("one-kink result: converged=%s E=%+.12f M=%+.6e\n",
            string(converged),old_e,sum(old_m))
    csv=joinpath(OUT,"profiles.csv"); png=joinpath(OUT,"one_kink.png")
    if success(`python3 -c "import matplotlib"`)
        run(`python3 -c $PLOT $csv $png`)
        println("Saved $png")
    else
        println("Install Python matplotlib to generate the PNG; CSV and MPS were saved")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

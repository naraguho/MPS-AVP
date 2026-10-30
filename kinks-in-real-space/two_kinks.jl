#!/usr/bin/env julia

# Zero-field A|B|A trial with two well-separated defects at 36|37 and 64|65.
# Run from the repository root:
# julia --startup-file=no --threads=1 --project=. kinks-in-real-space/two_kinks.jl
# Only the two defect blocks are optimized after the Neel background is ready.

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
const DEFECT_CHI = 128
const LEFT_BOND, RIGHT_BOND = 36, 64
const BACKGROUND_SWEEPS = parse(Int, get(ENV,"TWO_KINK_BACKGROUND_SWEEPS","6"))
const MAX_ROUNDS = parse(Int, get(ENV,"TWO_KINK_MAX_ROUNDS","80"))
const NROOTS = parse(Int, get(ENV,"TWO_KINK_NROOTS","8"))
const OUT = get(ENV,"TWO_KINK_OUTPUT",joinpath(@__DIR__,"outputs","two-kinks"))
const SZ = Matrix{Float64}(real.(QD.spin_ops(2)[:Z]))

function energy(psi,mpo)
    left = ones(Float64,1,1,1)
    for i in 1:N
        left = QD._contract_left_environment(left,psi.tensors[i],mpo.tensors[i])
    end
    return real(left[1,1,1]) / real(QD.measure_norm(psi))
end

function profile(psi)
    right = Vector{Matrix{Float64}}(undef,N+1)
    right[N+1] = ones(1,1)
    for i in N:-1:1
        a = psi.tensors[i]
        right[i] = zeros(size(a,1),size(a,1))
        for s in 1:2
            block = @view a[:,s,:]
            right[i] .+= block*right[i+1]*transpose(block)
        end
    end
    left = ones(1,1)
    m = zeros(N)
    for i in 1:N
        a = psi.tensors[i]
        up,down = @view(a[:,1,:]),@view(a[:,2,:])
        lu,ld = transpose(up)*left*up,transpose(down)*left*down
        m[i] = 0.5*sum((lu-ld).*right[i+1])
        left = lu+ld
    end
    return m ./ left[1,1]
end

staggered(m,sites) = mean((isodd(i) ? 1.0 : -1.0)*m[i] for i in sites)
function walls(m)
    q = [(isodd(i) ? 1.0 : -1.0)*m[i] for i in 1:N]
    return [i for i in 1:N-1 if q[i]*q[i+1]<0]
end

function same_spin(psi,j)
    corr = QD.measure_correlation(psi,SZ,j,SZ,j+1) / QD.measure_norm(psi)
    return 0.5 + 2*real(corr)
end

function observables(psi,mpo)
    m = profile(psi)
    return (;e=energy(psi,mpo),m,M=sum(m),wall_bonds=walls(m),
            left=staggered(m,10:25),middle=staggered(m,45:55),
            right=staggered(m,75:90),
            pleft=same_spin(psi,LEFT_BOND),pright=same_spin(psi,RIGHT_BOND))
end

function allowed(o)
    # This guard defines the retained A|B|A variational branch. It tests
    # domain order, the actual wall positions, and same-spin wall weight.
    w = o.wall_bonds
    return o.left>0.25 && o.middle< -0.25 && o.right>0.25 &&
           length(w)==2 && abs(w[1]-LEFT_BOND)<=5 &&
           abs(w[2]-RIGHT_BOND)<=5 && o.pleft>=0.09 && o.pright>=0.09
end

function background(mpo)
    tensors = [zeros(Float64,1,2,1) for _ in 1:N]
    for i in 1:N
        tensors[i][1,isodd(i) ? 1 : 2,1] = 1.0
    end
    # The product seed chooses branch A. Every sweep uses zero field.
    state = QD.MPSState(QD.MPS{Float64}(tensors),mpo;center=1)
    solver = QD.LanczosSolver(3,12)
    opts = QD.DMRGOptions(CHI,1e-10,2)
    for sweep in 1:BACKGROUND_SWEEPS
        QD.dmrg_sweep(state,solver,opts,isodd(sweep) ? :right : :left)
        @printf("background sweep=%2d E_XXZ=%+.12f\n",sweep,energy(state.mps,mpo))
        flush(stdout)
    end
    m = profile(state.mps)
    staggered(m,10:25)>0.3 && staggered(m,75:90)>0.3 ||
        error("The background failed to select a consistent Néel branch")
    return state.mps
end

function seed(A)
    psi = QD.MPS{Float64}(deepcopy(A.tensors))
    for i in LEFT_BOND+1:RIGHT_BOND
        psi.tensors[i] = psi.tensors[i][:,2:-1:1,:]
    end
    psi.tensors[1] ./= sqrt(real(QD.measure_norm(psi)))
    return psi
end

function candidate(canonical,j,x,dims)
    dl,d,dr = dims
    F = svd(reshape(x,dl*d,d*dr))
    keep = min(DEFECT_CHI,length(F.S))
    discarded = sum(abs2,F.S[keep+1:end])/sum(abs2,F.S)
    discarded<=1e-8 || return nothing
    psi = QD.MPS{Float64}(deepcopy(canonical.tensors))
    psi.tensors[j] = reshape(F.U[:,1:keep],dl,d,keep)
    psi.tensors[j+1] = reshape(Diagonal(F.S[1:keep])*F.Vt[1:keep,:],keep,d,dr)
    psi.tensors[j+1] ./= sqrt(real(QD.measure_norm(psi)))
    return psi
end

function update_block(psi,mpo,j)
    # Canonicalization changes MPS gauges but leaves its physical state fixed.
    state = QD.MPSState(QD.MPS{Float64}(deepcopy(psi.tensors)),mpo;center=j)
    a,b = state.mps.tensors[j],state.mps.tensors[j+1]
    dl,d,mid = size(a); dr = size(b,3)
    theta = reshape(reshape(a,dl*d,mid)*reshape(b,mid,d*dr),dl,d,d,dr)
    heff = QD.TwoSiteEffectiveHamiltonian(state.environment.tensors[j-1],
            mpo.tensors[j],mpo.tensors[j+1],state.environment.tensors[j+2])
    solver = KrylovKit.Lanczos(krylovdim=max(80,10NROOTS),maxiter=60,
                                tol=1e-9,verbosity=0)
    vals,vecs,info = KrylovKit.eigsolve(v->QD._apply(heff,v),
                                       vec(theta),NROOTS,:SR,solver)
    old = observables(state.mps,mpo)
    allowed(old) || error("The current state left the A|B|A branch before bond $j")
    best_psi,best_obs,best_root,best_res = state.mps,old,0,NaN
    for k in 1:min(NROOTS,length(vals))
        x,ev = vecs[k],real(vals[k])
        residual = norm(QD._apply(heff,x)-ev*x)/max(1.0,abs(ev))
        residual<=1e-7 || continue
        trial = candidate(state.mps,j,x,(dl,d,dr))
        trial===nothing && continue
        o = observables(trial,mpo)
        if allowed(o) && o.e<best_obs.e-1e-10
            best_psi,best_obs,best_root,best_res = trial,o,k,residual
        end
    end
    @printf("  block %d|%d: requested roots=%d returned=%d selected=%d residual=%.2e E=%+.12f\n",
            j,j+1,NROOTS,length(vals),best_root,best_res,best_obs.e)
    flush(stdout)
    return best_psi,best_obs,best_root
end

function save_profile(io,round,m)
    for i in 1:N
        @printf(io,"%d,%d,%.17g\n",round,i,m[i])
    end
    flush(io)
end

const PLOT = raw"""
import csv,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows=list(csv.DictReader(open(sys.argv[1])))
round_id=max(int(r['round']) for r in rows)
last=sorted((r for r in rows if int(r['round'])==round_id),
            key=lambda r:int(r['site']))
x=[int(r['site']) for r in last]; m=[float(r['Sz']) for r in last]
q=[(-1)**(i+1)*v for i,v in enumerate(m)]
walls=[i+.5 for i in range(1,len(q)) if q[i-1]*q[i]<0]
fig,ax=plt.subplots(figsize=(11,4),layout='constrained')
ax.plot(x,m,'o-',color='#de6b16',ms=2.5,lw=1.4,
        label=r'$\langle S_i^z\rangle$')
for k,w in enumerate(walls):
    ax.axvline(w,color='#555',ls='--',lw=1,
               label='observed wall' if k==0 else None)
ax.axhline(0,color='#777',lw=.7)
ax.set(xlim=(1,100),ylim=(-.52,.52),xlabel='site $i$',
       ylabel=r'$\langle S_i^z\rangle$',
       title=r'Two kinks: $N=100$, $\Delta=3$, zero field')
ax.legend(loc='upper right',frameon=False)
fig.savefig(sys.argv[2],dpi=190)
"""

function main()
    BLAS.set_num_threads(8)
    isdir(OUT) && !isempty(readdir(OUT)) && error("Output exists: $OUT")
    mkpath(OUT)
    mpo = QD.build_xxz_mpo(N,1.0,DELTA)
    psi = seed(background(mpo))
    old = observables(psi,mpo)
    allowed(old) || error("The initial two-kink seed failed its domain checks")
    stable=0; stalled=0; converged=false
    open(joinpath(OUT,"profiles.csv"),"w") do profiles
        println(profiles,"round,site,Sz")
        save_profile(profiles,0,old.m)
        open(joinpath(OUT,"round_summary.csv"),"w") do summary
            println(summary,"round,energy,dE,max_dSz,total_Sz,left_order,middle_order,right_order,P_left,P_right,wall_left,wall_right,root_left,root_right")
            @printf(summary,"0,%.17g,NaN,NaN,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,0,0\n",
                    old.e,old.M,old.left,old.middle,old.right,old.pleft,
                    old.pright,old.wall_bonds[1],old.wall_bonds[2])
            for round in 1:MAX_ROUNDS
                roots=zeros(Int,2)
                for (k,j) in enumerate((LEFT_BOND,RIGHT_BOND))
                    psi,_,roots[k] = update_block(psi,mpo,j)
                end
                o = observables(psi,mpo)
                de,dm = o.e-old.e,maximum(abs.(o.m.-old.m))
                @printf("round=%2d E=%+.12f dE=%+.2e max|dSz|=%.2e M=%+.3e walls=%s P=(%.3f,%.3f)\n",
                        round,o.e,de,dm,o.M,string(o.wall_bonds),o.pleft,o.pright)
                @printf(summary,"%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%d,%d\n",
                        round,o.e,de,dm,o.M,o.left,o.middle,o.right,
                        o.pleft,o.pright,o.wall_bonds[1],o.wall_bonds[2],roots...)
                flush(stdout);flush(summary)
                save_profile(profiles,round,o.m)
                allowed(o) || error("A|B|A structure was lost")
                stable = abs(de)<1e-7 && dm<1e-4 ? stable+1 : 0
                stalled = roots==[0,0] ? stalled+1 : 0
                old=o
                if stable>=2
                    converged=true
                    break
                end
                if stalled>=2
                    println("No admissible energy-lowering local root in two rounds")
                    break
                end
            end
        end
    end
    serialize(joinpath(OUT,"two_kink_mps.jls"),
              (;N,Delta=DELTA,chi=CHI,requested_bonds=(LEFT_BOND,RIGHT_BOND),
                observed_bonds=old.wall_bonds,energy=old.e,converged,
                tensors=deepcopy(psi.tensors)))
    @printf("two-kink result: converged=%s E=%+.12f M=%+.6e walls=%s\n",
            string(converged),old.e,old.M,string(old.wall_bonds))
    csv=joinpath(OUT,"profiles.csv");png=joinpath(OUT,"two_kinks.png")
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

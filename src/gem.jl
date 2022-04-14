#-----------------------------functions related to gaussian basises-----------------------------------
mutable struct GemBasis
    nmax::Int8
    r1::Float64
    rnmax::Float64
end

function gembasisn(gb::GemBasis,n)::Float64
    nmax,r1, rnmax=gb.nmax,gb.r1,gb.rnmax
    return 1 / r1^2 * (r1 / rnmax)^((2n - 2) / (nmax - 1))
end

#---------------------------functions related to Hamiltonian stuff-------------------------------------
mutable struct GemChannel
    id::Int8
    μ::Float64
    ΔE::Float64
    l::Int8
end

mutable struct GemHamiltonian
    N::Int8
    T::Vector{GemChannel}
    V 
    function GemHamiltonian(vμ,vΔE,vl,Vf)
        if !(length(vμ)==length(vΔE)&&length(vμ)==length(vl))
            error("The input vectors should be in the same length!")
        end
        len=length(vμ)
        tmp=Vector{GemChannel}(undef,len)
        for i in 1:len
            tmp[i]=GemChannel(i,Float64(vμ[i]),Float64(vΔE[i]),Int8(vl[i]))
        end
        new(len,tmp,Vf)
    end
end

mutable struct GemModel 
    Tmat::Matrix{Float64}
    Vmat::Matrix{Float64}
    ΔEmat::Matrix{Float64}
    Nmat::Matrix{Float64}

    evals::Vector{ComplexF64}
    evecs::Matrix{ComplexF64}
    
    function GemModel(gh::GemHamiltonian,gb::GemBasis)
        nc,nb=gh.N,gb.rnmax
        new(zeros(Float64,nc*nb,nc*nb),zeros(Float64,nc*nb,nc*nb),zeros(Float64,nc*nb,nc*nb)
        ,zeros(Float64,nc*nb,nc*nb),zeros(ComplexF64,nc*nb),zeros(ComplexF64,nc*nb,nc*nb))
    end
end



function GemSolve!(mod::GemModel,hal::GemHamiltonian)
    # defined by user
end
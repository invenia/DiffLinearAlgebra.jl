module DiffLinAlg

    using FDM

    # Some aliases used repeatedly throughout the package.
    const AA, AM, AVM = AbstractArray, AbstractMatrix, AbstractVecOrMat
    const AS, ASVM = Union{Real, AA}, Union{Real, AVM}
    const Arg1, Arg2 = Type{Val{1}}, Type{Val{2}}


    # Some aliases used repeatedly throughout the package.
    export ∇

    include("generic.jl")
    # include("blas.jl")
    # include("strided.jl")
    # include("diagonal.jl")
    # include("triangular.jl")
    # include("uniformscaling.jl")
    # include("factorization/cholesky.jl")

end # module

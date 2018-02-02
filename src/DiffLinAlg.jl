module DiffLinAlg

    using FDM

    # Some aliases used repeatedly throughout the package.
    const AA, AM, AVM = AbstractArray, AbstractMatrix, AbstractVecOrMat
    const AS, ASVM = Union{Real, AA}, Union{Real, AVM}
    const Arg1, Arg2, Arg3 = Type{Val{1}}, Type{Val{2}}, Type{Val{3}}
    const Arg4, Arg5, Arg6 = Type{Val{4}}, Type{Val{5}}, Type{Val{6}}
    const SV, SM, SVM, SA = StridedVector, StridedMatrix, StridedVecOrMat, StridedArray

    # Some aliases used repeatedly throughout the package.
    export ∇

    include("generic.jl")
    include("blas.jl")
    # include("strided.jl")
    # include("diagonal.jl")
    # include("triangular.jl")
    # include("uniformscaling.jl")
    # include("factorization/cholesky.jl")

end # module

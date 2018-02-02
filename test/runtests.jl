using DiffLinAlg, Base.Test

import FDM: assert_approx_equal, central_fdm
import DiffLinAlg: AA, AM, AVM, AS, ASVM, Arg1, Arg2

import BenchmarkTools: @benchmark

@testset "DiffLinearAlgebra" begin
    include("test_util.jl")
    include("generic.jl")
    include("blas.jl")
    include("diagonal.jl")
    include("triangular.jl")
    include("uniformscaling.jl")
    include("factorization/cholesky.jl")
end

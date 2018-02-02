@testset "Triangular" begin

    for f in [LowerTriangular, UpperTriangular]
        let P = 10, rng = MersenneTwister(123456), N = 10

            sc = ()->randn(rng)
            mPP = ()->abs.(randn(rng, P, P))
            dPP = ()->abs.(randn(rng, P, P))
            tPP = ()->f(abs.(randn(rng, P, P)))
            m0 = fill!(similar(mPP()), zero(eltype(mPP())))

            # Construction.
            @test check_errs(N, unary_ȲD(f)..., tPP, mPP, mPP)
            @test check_errs(N, unary_ȲD_inplace(f, m0)..., tPP, mPP, mPP)

            # Determinant.
            @test check_errs(N, unary_ȲD(det)..., sc, tPP, tPP)
            @test check_errs(N, unary_ȲD_inplace(det, zero(tPP()))..., sc, tPP, tPP)
            @test check_errs(N, unary_ȲD_inplace(det, zero(dPP()))..., sc, tPP, tPP)

            # Log Determinant.
            @test check_errs(N, unary_ȲD(logdet)..., sc, tPP, tPP)
            @test check_errs(N, unary_ȲD_inplace(logdet, zero(tPP()))..., sc, tPP, tPP)
            @test check_errs(N, unary_ȲD_inplace(logdet, zero(dPP()))..., sc, tPP, tPP)
        end
    end

    # let rng = MersenneTwister(123456), N = 10
    #     for _ in 1:10
    #         A = LowerTriangular(randn(rng, N, N))
    #         VA = LowerTriangular(randn(rng, N, N))
    #         @test check_errs(det, 10.0, A, VA)
    #     end
    # end
    # let rng = MersenneTwister(123456), N = 3
    #     for _ in 1:10
    #         A = UpperTriangular(randn(rng, N, N))
    #         VA = UpperTriangular(randn(rng, N, N))
    #         @test check_errs(det, 10.0, A, VA)
    #     end
    # end
    # let rng = MersenneTwister(123456), N = 10
    #     for _ in 1:10
    #         A = LowerTriangular(exp.(randn(rng, N, N)))
    #         VA = LowerTriangular(randn(rng, N, N))
    #         @test check_errs(logdet, 10.0, A, VA)
    #     end
    # end
    # let rng = MersenneTwister(123456), N = 3
    #     for _ in 1:10
    #         A = UpperTriangular(exp.(randn(rng, N, N)))
    #         VA = UpperTriangular(randn(rng, N, N))
    #         @test check_errs(logdet, 10.0, A, VA)
    #     end
    # end

    # # Check that the optimisations occur correctly and produce the required types when
    # # everything is Diagonal.
    # let rng = MersenneTwister(123456)
    #     A = UpperTriangular(exp.(randn(rng, 10, 10)))

    #     @test ∇(det)(A)[1] isa Diagonal
    #     @test ∇(A->det(A) + det(A))(A)[1] isa Diagonal
    #     @test ∇(logdet)(A)[1] isa Diagonal
    #     @test ∇(A->logdet(A) + det(A))(A)[1] isa Diagonal
    # end
end

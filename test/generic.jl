@testset "Generic" begin

    let P = 5, Q = 4, rng = MersenneTwister(123456), N = 100

        # Utility to create the closures required for unit testing.
        unary_ȲD(f) = (Y, Ȳ, X)->∇(f, Val{1}, (), Y, Ȳ, X)
        function binary_ȲD(f, arg::Int, G)
            G_ = G()
            if arg == 1
                return X->f(X, G_), (Y, Ȳ, X)->∇(f, Val{1}, (), Y, Ȳ, X, G_)
            else
                return X->f(G_, X), (Y, Ȳ, X)->∇(f, Val{2}, (), Y, Ȳ, G_, X)
            end
        end

        # Utility for generating square matrices, vectors, non-square matrices and scalars.
        mPP, mQQ = ()->randn(rng, P, P), ()->randn(rng, Q, Q)
        mPQ, mQP = ()->randn(rng, P, Q), ()->randn(rng, Q, P)
        v, sc = ()->randn(rng, P), ()->randn(rng)
        psd = ()->(A = randn(rng, P, P); transpose(A) * A + 1e-3I)

        # Test all of the unary sensitivities.
        @test repeatedly_check_errs(N, -, unary_ȲD(-), mPQ, mPQ, mPQ)
        @test repeatedly_check_errs(N, trace, unary_ȲD(trace), sc, mPP, mPP)
        @test repeatedly_check_errs(N, inv, unary_ȲD(inv), mPP, mPP, mPP)
        @test repeatedly_check_errs(N, det, unary_ȲD(det), sc, mPP, mPP)
        @test repeatedly_check_errs(N, logdet, unary_ȲD(logdet), sc, psd, psd)
        @test repeatedly_check_errs(N, transpose, unary_ȲD(transpose), mQP, mPQ, mPQ)
        @test repeatedly_check_errs(N, ctranspose, unary_ȲD(ctranspose), mQP, mPQ, mPQ)
        @test repeatedly_check_errs(N, vecnorm, unary_ȲD(vecnorm), sc, mPQ, mPQ)
        @test repeatedly_check_errs(N, vecnorm, unary_ȲD(vecnorm), sc, sc, sc)

        # Test all of the binary sensitivities.
        @test repeatedly_check_errs(N, binary_ȲD(*, 1, mQP)..., mPP, mPQ, mPQ)
        @test repeatedly_check_errs(N, binary_ȲD(*, 2, mPQ)..., mPP, mQP, mQP)
        @test repeatedly_check_errs(N, binary_ȲD(At_mul_B, 1, mQP)..., mPP, mQP, mQP)
        @test repeatedly_check_errs(N, binary_ȲD(At_mul_B, 2, mQP)..., mPP, mQP, mQP)
        @test repeatedly_check_errs(N, binary_ȲD(A_mul_Bt, 1, mPQ)..., mPP, mPQ, mPQ)
        @test repeatedly_check_errs(N, binary_ȲD(A_mul_Bt, 2, mPQ)..., mPP, mPQ, mPQ)
        @test repeatedly_check_errs(N, binary_ȲD(At_mul_Bt, 1, mPQ)..., mPP, mQP, mQP)
        @test repeatedly_check_errs(N, binary_ȲD(At_mul_Bt, 2, mQP)..., mPP, mPQ, mPQ)
        @test repeatedly_check_errs(N, binary_ȲD(Ac_mul_B, 1, mQP)..., mPP, mQP, mQP)
        @test repeatedly_check_errs(N, binary_ȲD(Ac_mul_B, 2, mQP)..., mPP, mQP, mQP)
        @test repeatedly_check_errs(N, binary_ȲD(A_mul_Bc, 1, mPQ)..., mPP, mPQ, mPQ)
        @test repeatedly_check_errs(N, binary_ȲD(A_mul_Bc, 2, mPQ)..., mPP, mPQ, mPQ)
        @test repeatedly_check_errs(N, binary_ȲD(Ac_mul_Bc, 1, mPQ)..., mPP, mQP, mQP)
        @test repeatedly_check_errs(N, binary_ȲD(Ac_mul_Bc, 2, mQP)..., mPP, mPQ, mPQ)
        @test repeatedly_check_errs(N, binary_ȲD(/, 1, mQQ)..., mPQ, mPQ, mPQ)
        @test repeatedly_check_errs(N, binary_ȲD(/, 2, mPQ)..., mPQ, mQQ, mQQ)
        @test repeatedly_check_errs(N, binary_ȲD(At_rdiv_B, 1, mQQ)..., mPQ, mQP, mQP)
        @test repeatedly_check_errs(N, binary_ȲD(At_rdiv_B, 2, mQP)..., mPQ, mQQ, mQQ)
        @test repeatedly_check_errs(N, binary_ȲD(A_rdiv_Bt, 1, mQQ)..., mPQ, mPQ, mPQ)
        @test repeatedly_check_errs(N, binary_ȲD(A_rdiv_Bt, 2, mPQ)..., mPQ, mQQ, mQQ)
        @test repeatedly_check_errs(N, binary_ȲD(At_rdiv_Bt, 1, mQQ)..., mPQ, mQP, mQP)
        @test repeatedly_check_errs(N, binary_ȲD(At_rdiv_Bt, 2, mQP)..., mPQ, mQQ, mQQ)
        @test repeatedly_check_errs(N, binary_ȲD(Ac_rdiv_B, 1, mQQ)..., mPQ, mQP, mQP)
        @test repeatedly_check_errs(N, binary_ȲD(Ac_rdiv_B, 2, mQP)..., mPQ, mQQ, mQQ)
        @test repeatedly_check_errs(N, binary_ȲD(A_rdiv_Bc, 1, mQQ)..., mPQ, mPQ, mPQ)
        @test repeatedly_check_errs(N, binary_ȲD(A_rdiv_Bc, 2, mPQ)..., mPQ, mQQ, mQQ)
        @test repeatedly_check_errs(N, binary_ȲD(Ac_rdiv_Bc, 1, mQQ)..., mPQ, mQP, mQP)
        @test repeatedly_check_errs(N, binary_ȲD(Ac_rdiv_Bc, 2, mQP)..., mPQ, mQQ, mQQ)
        @test repeatedly_check_errs(N, binary_ȲD(\, 1, mQP)..., mQP, mQQ, mQQ)
        @test repeatedly_check_errs(N, binary_ȲD(\, 2, mQQ)..., mQP, mQP, mQP)
        @test repeatedly_check_errs(N, binary_ȲD(At_ldiv_B, 1, mQP)..., mQP, mQQ, mQQ)
        @test repeatedly_check_errs(N, binary_ȲD(At_ldiv_B, 2, mQQ)..., mQP, mQP, mQP)
        @test repeatedly_check_errs(N, binary_ȲD(A_ldiv_Bt, 1, mPQ)..., mQP, mQQ, mQQ)
        @test repeatedly_check_errs(N, binary_ȲD(A_ldiv_Bt, 2, mQQ)..., mQP, mPQ, mPQ)
        @test repeatedly_check_errs(N, binary_ȲD(At_ldiv_Bt, 1, mPQ)..., mQP, mQQ, mQQ)
        @test repeatedly_check_errs(N, binary_ȲD(At_ldiv_Bt, 2, mQQ)..., mQP, mPQ, mPQ)
        @test repeatedly_check_errs(N, binary_ȲD(Ac_ldiv_B, 1, mQP)..., mQP, mQQ, mQQ)
        @test repeatedly_check_errs(N, binary_ȲD(Ac_ldiv_B, 2, mQQ)..., mQP, mQP, mQP)
        @test repeatedly_check_errs(N, binary_ȲD(A_ldiv_Bc, 1, mPQ)..., mQP, mQQ, mQQ)
        @test repeatedly_check_errs(N, binary_ȲD(A_ldiv_Bc, 2, mQQ)..., mQP, mPQ, mPQ)
        @test repeatedly_check_errs(N, binary_ȲD(At_ldiv_Bt, 1, mPQ)..., mQP, mQQ, mQQ)
        @test repeatedly_check_errs(N, binary_ȲD(At_ldiv_Bt, 2, mQQ)..., mQP, mPQ, mPQ)
        @test repeatedly_check_errs(N, binary_ȲD(vecnorm, 1, sc)..., sc, mPQ, mPQ)
        @test repeatedly_check_errs(N, binary_ȲD(vecnorm, 2, mPQ)..., sc, sc, sc)
        @test repeatedly_check_errs(N, binary_ȲD(vecnorm, 1, sc)..., sc, sc, sc)
        @test repeatedly_check_errs(N, binary_ȲD(vecnorm, 2, sc)..., sc, sc, sc)
    end
end

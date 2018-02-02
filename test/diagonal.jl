@testset "Diagonal" begin

    # diag:
    let P = 10, rng = MersenneTwister(123456), N = 10, k = 2

        # Some rngs.
        vP, mPP, vPk = ()->randn(rng, P), ()->randn(rng, P, P), ()->randn(rng, P - k)

        # Test on-central-diagonal `diag`.
        X0 = fill!(similar(mPP()), zero(eltype(mPP())))
        @test check_errs(N, unary_ȲD(diag)..., vP, mPP, mPP)
        @test check_errs(N, unary_ȲD_inplace(diag, X0)..., vP, mPP, mPP)

        # Test off-central-diagonal `diag`.
        _diag = X->diag(X, k)
        _∇diag = (y, ȳ, X)->∇(diag, Val{1}, (), y, ȳ, X, k)
        _∇diag_inp = (y, ȳ, X)->∇(fill!(similar(X), 0), diag, Val{1}, (), y, ȳ, X, k)
        @test check_errs(N, _diag, _∇diag, vPk, mPP, mPP)
        @test check_errs(N, _diag, _∇diag_inp, vPk, mPP, mPP)
    end

    # diagm:
    let P = 10, rng = MersenneTwister(123456), N = 10, k = 3

        # Some rngs.
        sc, vP = ()->randn(rng), ()->randn(rng, P)
        mPP, mPPk = ()->randn(rng, P, P), ()->randn(rng, P + k, P + k)

        # Test on-central-diagonal `diagm`.
        x0 = fill!(similar(vP()), zero(eltype(vP())))
        @test check_errs(N, unary_ȲD(diagm)..., mPP, vP, vP)
        @test check_errs(N, unary_ȲD_inplace(diagm, x0)..., mPP, vP, vP)

        # Test off-central-diagonal `diagm`.
        _diagm = x->diagm(x, k)
        _∇diagm = (Y, Ȳ, x)->∇(diagm, Val{1}, (), Y, Ȳ, x, k)
        _∇diagm_inp = (Y, Ȳ, x)->∇(fill!(similar(x), 0), diagm, Val{1}, (), Y, Ȳ, x, k)
        @test check_errs(N, _diagm, _∇diagm, mPPk, vP, vP)
        @test check_errs(N, _diagm, _∇diagm_inp, mPPk, vP, vP)

        # Test scalar input to `diagm`.
        @test check_errs(N, unary_ȲD(diagm)..., ()->randn(rng, 1, 1), sc, sc)
    end

    # Diagonal:
    let P = 10, rng = MersenneTwister(123456), N = 10
        sc, vP = ()->abs(randn(rng)), ()->abs.(randn(rng, P))
        mPP, dPP = ()->abs.(randn(rng, P, P)), ()->Diagonal(abs.(randn(rng, P)))
        D0 = Diagonal(fill!(similar(vP()), 0))

        # Construction.
        @test check_errs(N, unary_ȲD(Diagonal)..., dPP, vP, vP)
        @test check_errs(N, unary_ȲD_inplace(Diagonal, zeros(P))..., dPP, vP, vP)
        @test check_errs(N, unary_ȲD(Diagonal)..., dPP, mPP, mPP)
        @test check_errs(N, unary_ȲD_inplace(Diagonal, zeros(P, P))..., dPP, mPP, mPP)

        # Determinant.
        @test check_errs(N, unary_ȲD(det)..., sc, dPP, dPP)
        @test check_errs(N, unary_ȲD_inplace(det, D0)..., sc, dPP, dPP)

        # Log Determinant.
        @test check_errs(N, unary_ȲD(det)..., sc, dPP, dPP)
        @test check_errs(N, unary_ȲD_inplace(det, D0)..., sc, dPP, dPP)
    end

    # let rng = MersenneTwister(123456), N = 10
    #     λ_2 = x->diagm(x, 2)
    #     λ_m3 = x->diagm(x, -3)
    #     λ_0 = x->diagm(x, 0)
    #     λ_false = x->diagm(x, false)
    #     λ_true = x->diagm(x, true)
    #     for _ in 1:10

    #         # Test vector case.
    #         x, vx = randn.(rng, [N, N])
    #         @test check_errs(diagm, diagm(randn(rng, N)), x, vx)

    #         @test check_errs(λ_2, λ_2(randn(rng, N)), x, vx)
    #         @test check_errs(λ_m3, λ_m3(randn(rng, N)), x, vx)
    #         @test check_errs(λ_0, λ_0(randn(rng, N)), x, vx)
    #         @test check_errs(λ_false, λ_false(randn(rng, N)), x, vx)
    #         @test check_errs(λ_true, λ_true(randn(rng, N)), x, vx)

    #         # Test scalar case.
    #         x, vx = randn(rng), randn(rng)
    #         @test check_errs(diagm, diagm(randn(rng)), x, vx)
    #     end
    # end
    # let rng = MersenneTwister(123456), N = 10
    #     for _ in 1:10
    #         A = randn(rng, N, N)
    #         VA = randn(rng, N, N)
    #         @test check_errs(x -> diag(x), randn(rng, N), A, VA)

    #         # Check various diagonals.
    #         for k = -3:3
    #             @test check_errs(x -> diag(x, k), randn(rng, N - abs(k)), A, VA)
    #         end
    #     end
    # end
    # let rng = MersenneTwister(123456), N = 10
    #     for _ in 1:10
    #         A = randn(rng, N)
    #         VA = randn(rng, N)
    #         @test check_errs(Diagonal, Diagonal(randn(rng, N)), A, VA)
    #     end
    # end
    # let rng = MersenneTwister(123456), N = 10
    #     for _ in 1:10
    #         A = randn(rng, N, N)
    #         VA = randn(rng, N, N)
    #         @test check_errs(Diagonal, Diagonal(randn(rng, N)), A, VA)
    #     end
    # end
    # let rng = MersenneTwister(123456), N = 10
    #     for _ in 1:10
    #         A = Diagonal(randn(rng, N))
    #         VA = Diagonal(randn(rng, N))
    #         @test check_errs(det, 10.0, A, VA)
    #     end
    # end
    # let rng = MersenneTwister(123456), N = 10
    #     for _ in 1:10
    #         A = Diagonal(exp.(randn(rng, N)))
    #         VA = Diagonal(randn(rng, N))
    #         @test check_errs(logdet, 10.0, A, VA)
    #     end
    # end
end

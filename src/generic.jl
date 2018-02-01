import Base: -, trace, inv, det, logdet, transpose, ctranspose, vecnorm

# Unary sensitivities.
∇(::typeof(-), ::Arg1, p, Y::AA, Ȳ::AA, X::AA) = map(-, Ȳ)
∇(::typeof(trace), ::Arg1, p, Y::Real, Ȳ::Real, X::AM) = Diagonal(fill!(similar(X), Ȳ))
∇(::typeof(inv), ::Arg1, p, Y::AM, Ȳ::AM, X::AM) = -transpose(Y) * Ȳ * transpose(Y)
∇(::typeof(det), ::Arg1, p, Y::Real, Ȳ::Real, X::AM) = Y * Ȳ * transpose(inv(X))
∇(::typeof(logdet), ::Arg1, p, Y::Real, Ȳ::Real, X::AM) = Ȳ * transpose(inv(X))
∇(::typeof(transpose), ::Arg1, p, Y::AVM, Ȳ::AVM, X::AVM) = transpose(Ȳ)
∇(::typeof(ctranspose), ::Arg1, p, Y::AVM, Ȳ::AVM, X::AVM) = ctranspose(Ȳ)
∇(::typeof(vecnorm), ::Arg1, p, Y::Real, Ȳ::Real, X::AA) = Ȳ ./ Y .* abs2.(X) ./ X
∇(::typeof(vecnorm), ::Arg1, p, Y::Real, Ȳ::Real, X::Real) = Ȳ * sign(X)

# Binary sensitivities.
∇(::typeof(*), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::AVM, B::AVM) = Ȳ * B'
∇(::typeof(*), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::AVM, B::AVM) = A' * Ȳ
∇(::typeof(At_mul_B), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::AVM, B::AVM) = B * Ȳ'
∇(::typeof(At_mul_B), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::AVM, B::AVM) = A * Ȳ
∇(::typeof(A_mul_Bt), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::AVM, B::AVM) = Ȳ * B
∇(::typeof(A_mul_Bt), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::AVM, B::AVM) = Ȳ' * A
∇(::typeof(At_mul_Bt), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = B' * Ȳ'
∇(::typeof(At_mul_Bt), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = Ȳ' * A'
∇(::typeof(Ac_mul_B), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = B * Ȳ'
∇(::typeof(Ac_mul_B), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = A * Ȳ
∇(::typeof(A_mul_Bc), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = Ȳ * B
∇(::typeof(A_mul_Bc), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = Ȳ' * A
∇(::typeof(Ac_mul_Bc), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = B' * Ȳ'
∇(::typeof(Ac_mul_Bc), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = Ȳ' * A'
∇(::typeof(/), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = Ȳ / B'
∇(::typeof(/), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = -Y' * (Ȳ / B')
∇(::typeof(At_rdiv_B), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = B \ Ȳ'
∇(::typeof(At_rdiv_B), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = -Y' * (B \ Ȳ')'
∇(::typeof(A_rdiv_Bt), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = Ȳ / B
∇(::typeof(A_rdiv_Bt), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = -(B' \ Ȳ') * Y
∇(::typeof(At_rdiv_Bt), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = B' \ Ȳ'
∇(::typeof(At_rdiv_Bt), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = -(B' \ Ȳ') * Y
∇(::typeof(Ac_rdiv_B), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = B \ Ȳ'
∇(::typeof(Ac_rdiv_B), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = -Y' * (B \ Ȳ')'
∇(::typeof(A_rdiv_Bc), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = Ȳ / B
∇(::typeof(A_rdiv_Bc), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = -(B' \ Ȳ') * Y
∇(::typeof(Ac_rdiv_Bc), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = B' \ Ȳ'
∇(::typeof(Ac_rdiv_Bc), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = -(B' \ Ȳ') * Y
∇(::typeof(\), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = -(A' \ Ȳ) * Y'
∇(::typeof(\), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = A' \ Ȳ
∇(::typeof(At_ldiv_B), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = -Y *(A \ Ȳ)'
∇(::typeof(At_ldiv_B), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = A \ Ȳ
∇(::typeof(A_ldiv_Bt), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = -(Ȳ' / A)' * Y'
∇(::typeof(A_ldiv_Bt), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = Ȳ' / A
∇(::typeof(At_ldiv_Bt), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = -Y * (Ȳ' / A')
∇(::typeof(At_ldiv_Bt), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = Ȳ' / A'
∇(::typeof(Ac_ldiv_B), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = -Y * (A \ Ȳ)'
∇(::typeof(Ac_ldiv_B), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = A \ Ȳ
∇(::typeof(A_ldiv_Bc), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = -(Ȳ' / A)' * Y'
∇(::typeof(A_ldiv_Bc), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = Ȳ' / A
∇(::typeof(Ac_ldiv_Bc), ::Arg1, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = -Y * (Ȳ' / A')
∇(::typeof(Ac_ldiv_Bc), ::Arg2, p, Y::ASVM, Ȳ::ASVM, A::ASVM, B::ASVM) = Ȳ' / A'
∇(::typeof(vecnorm), ::Arg1, p, Y::Real, Ȳ::Real, A::AA, B::Real) =
    Ȳ .* Y^(1 - B) .* abs.(A).^B ./ A
∇(::typeof(vecnorm), ::Arg2, p, Y::Real, Ȳ::Real, A::AA, B::Real) =
    Ȳ * (Y^(1 - B) * sum(abs.(A).^B .* log.(abs.(A))) - Y * log(Y)) / B
∇(::typeof(vecnorm), ::Arg1, p, Y::Real, Ȳ::Real, A::Real, B::Real) = Ȳ * sign(A)
∇(::typeof(vecnorm), ::Arg2, p, Y::Real, Ȳ::Real, A::Real, B::Real) = 0

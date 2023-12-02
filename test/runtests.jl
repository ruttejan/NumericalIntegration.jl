using NumericalIntegration
using Test

@testset "NumericalIntegration.jl" begin
    @testset "numerical_int()" begin
        f(x) = -x^2 + 1
        a = -1
        b = 1
        n = 2
        @test numerical_int(Midpoint(), f, a, b, n) == 0.75
        @test numerical_int(Trapezoid(), f, a, b, n) == 1
        @test numerical_int(Simpson(), f, a, b, n) == 4/3
    end

    @testset "error_int()" begin
        f(x) = -x^2 + 1
        a = -1
        b = 1
        n = 2
        @test false
        @test false
        @test false
    end

    @testset "integrate @ integrate_compare" begin
        f(x) = -x^2 + 1
        a = -1
        b = 1
        @test isapprox(integrate(Midpoint(), f, a, b)[1], 4/3; atol=1e-6)
        @test isapprox(integrate(Trapezoid(), f, a, b)[1], 4/3; atol=1e-6)
        @test isapprox(integrate(Simpson(), f, a, b)[1], 4/3; atol=1e-6)
        @test isapprox(integrate_compare(Midpoint(), f, a, b)[1], 4/3; atol=1e-6)
        @test isapprox(integrate_compare(Trapezoid(), f, a, b)[1], 4/3; atol=1e-6)
        @test isapprox(integrate_compare(Simpson(), f, a, b)[1], 4/3; atol=1e-6)
    end
end

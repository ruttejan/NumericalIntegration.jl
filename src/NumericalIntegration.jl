module NumericalIntegration
using Plots, ForwardDiff, QuadGK

export Trapezoid, Midpoint, Simpson, numerical_int, error_int, plot_int, integrate, integrate_compare

abstract type IntMethod end
struct Trapezoid <: IntMethod end
struct Midpoint <: IntMethod end
struct Simpson <: IntMethod end

### ------------------------------------------------------------------------------------------------
# Trapezoid method

"""
    numerical_int(method, f, a, b, n)

Calculates the numerical integral of the function f on interval [a,b] 
using a selected method and certain number of split points.

# Arguments
- `method::IntMethod`: The method to use for numerical integration [Trapezoid(), Midpoint(), Simpson()].
- `f::Function`: The function to integrate.
- `a::Real`: The lower bound of the integration interval.
- `b::Real`: The upper bound of the integration interval.
- `n::Int`: The number of split points.

# Returns
- `integral::Real`: The numerical integral of the function f on interval [a,b].

# Example
```julia-repl
julia> f(x) = x^2
julia> numerical_int(Trapezoid(), f, 0.0, 1.0, 10)
"""
function numerical_int(::Trapezoid,
    f::Function,
    a::Real,
    b::Real,
    n::Int
)
    h = (b - a) / n
    x_values = collect((a+h):h:(b-h))
    integral = (f(a) + f(b))
    for i = eachindex(x_values)
        integral += 2 * f(x_values[i])
    end
    return (h / 2) * integral
end

# Error upper bound of the trapezoid method

"""

    error_int(method, f, a, b, n)

Calculates the error upper bound of the numerical integral of the function f on interval [a,b]
using a selected method and certain number of split points.

# Arguments
- `method::IntMethod`: The method to use for numerical integration [Trapezoid(), Midpoint(), Simpson()].
- `f::Function`: The function to integrate.
- `a::Real`: The lower bound of the integration interval.
- `b::Real`: The upper bound of the integration interval.
- `n::Int`: The number of split points.

# Returns
- `error::Real`: The error upper bound of the numerical integral of the function f on interval [a,b].

# Example
```julia-repl
julia> f(x) = x^2
julia> error_int(Trapezoid(), f, 0.0, 1.0, 10)
"""
function error_int(::Trapezoid,
    f::Function,
    a::Real,
    b::Real,
    n::Int
)
    fd1(x) = ForwardDiff.derivative(f, x)
    fd2(x) = ForwardDiff.derivative(fd1, x)
    fx_max = findmax(x -> abs(fd2(x)), a:b)[1]
    return (b - a)^3 / 12 * fx_max / n^2
end

## Plotting the trapezoid method

"""

    plot_int(method, f, a, b, n)

Plots the function f and the numerical integral of the function f on interval [a,b]
using a selected method and certain number of split points.

# Arguments
- `method::IntMethod`: The method to use for numerical integration [Trapezoid(), Midpoint(), Simpson()].
- `f::Function`: The function to integrate.
- `a::Real`: The lower bound of the integration interval.
- `b::Real`: The upper bound of the integration interval.
- `n::Int`: The number of split points.

# Example
```julia-repl
julia> f(x) = x^2
julia> plot_int(Trapezoid(), f, 0.0, 1.0, 10)
"""
function plot_int(::Trapezoid,
    f::Function,
    a::Real,
    b::Real,
    n::Int
)
    h = (b - a) / n
    x_values = collect(a:h:b)
    y_values = collect(f.(x_values))
    plot(f, a, b, label="Function", legend=true)
    for i in 2:n+1
        x1 = x_values[i-1]
        x2 = x_values[i]
        y1 = y_values[i-1]
        y2 = y_values[i]
        plot!([x1, x2, x2, x1, x1], [0.0, 0.0, y2, y1, 0.0],
            seriestype=:shape,
            lw=1,
            color=:red,
            fillalpha=0.2,
            label=nothing)
    end
    scatter!(x_values, zeros(size(x_values)), color=:black, markersize=4, label="n-points")
    scatter!(x_values, y_values, color=:blue, markersize=4, label="Data Points")

end

### ------------------------------------------------------------------------------------------------
# Midpoint method
function numerical_int(::Midpoint,
    f::Function,
    a::Real,
    b::Real,
    n::Int
)
    h = (b - a) / n
    x_values = collect(a:h:b)
    integral = 0
    for i = 2:n
        integral += f((x_values[i-1] + x_values[i]) / 2)
    end
    return h * integral
end

#Error upper bound of the midpoint method

function error_int(::Midpoint, 
    f :: Function,
    a :: Real, 
    b :: Real, 
    n :: Int)
    fd1(x) = ForwardDiff.derivative(f, x)
    fd2(x) = ForwardDiff.derivative(fd1, x)
    fx_max = findmax(x -> abs(fd2(x)), a:b)[1]
    return (((b - a)^3) * fx_max) / (24 * n^2)
end

## Plotting the midpoint method

function plot_int(::Midpoint,
    f::Function,
    a::Real,
    b::Real,
    n::Int
)
    h = (b - a) / n
    x_values = collect(a:h:b)
    x_mids = []
    y_mids = []
    plot(f, a, b, label="Function", legend=true)#, size=(800, 300))
    for i in 2:n+1
        x1 = x_values[i-1]
        x2 = x_values[i]
        mid = (x1 + x2) / 2
        push!(x_mids, mid)
        y = f(mid)
        push!(y_mids, y)
        plot!([x1, x2, x2, x1, x1], [0.0, 0.0, y, y, 0.0],
            seriestype=:shape,
            lw=1,
            color=:red,
            fillalpha=0.2,
            label=nothing)
    end

    scatter!(x_mids, y_mids, color=:blue, markersize=4, label="Midpoints")
    scatter!(x_values, zeros(size(x_values)), color=:black, markersize=4, label="n-points")

end

### ------------------------------------------------------------------------------------------------
#Simpsons method
function numerical_int(::Simpson,
    f::Function,
    a::Real,
    b::Real,
    n::Int
)
    h = (b - a) / n
    h2 = h / 2
    x_edge = collect((a+h):h:(b-h))
    x_mids = collect((a+h2):h:(b-h/2))
    integral = f(a) + f(b)
    for i in eachindex(x_edge)
        integral += 2 * f(x_edge[i])
    end

    for i in eachindex(x_mids)
        integral += 4 * f(x_mids[i])
    end
    return integral * (h / 6)
end

#Error upper bound of the Simpson's method

function error_int(::Simpson,
    f::Function,
    a::Real,
    b::Real,
    n::Int
)

    fd1(x) = ForwardDiff.derivative(f, x)
    fd2(x) = ForwardDiff.derivative(fd1, x)
    fd3(x) = ForwardDiff.derivative(fd2, x)
    fd4(x) = ForwardDiff.derivative(fd3, x)


    fx_max = findmax(x -> abs(fd4(x)), a:b)[1]
    return (((b - a)^5) * fx_max) / (180 * n^4)
end

## Plotting the Simpson's method
function plot_int(::Simpson,
    f::Function,
    a::Real,
    b::Real,
    n::Int
)
    h = (b - a) / n
    x_values = collect(a:h:b)
    x_mids = []
    y_mids = []
    plot(f, a, b, label="Function", legend=true)
    for i in 2:n+1
        x1 = x_values[i-1]
        x2 = x_values[i]
        mid = (x1 + x2) / 2
        push!(x_mids, mid)
        y = f(mid)
        push!(y_mids, y)
        px = [x1, mid, x2]
        py = [f(x1), f(mid), f(x2)]
        A = hcat([px[i]^2 for i in 1:3], px, [1, 1, 1])
        b = py
        a, b, c = A \ b
        parabola(x, a, b, c) = a * x^2 + b * x + c
        points_x = vcat([x2, x1], collect(range(x1, x2, length=100)))
        points_y = vcat([0.0, 0.0], [parabola(x, a, b, c) for x in points_x[3:length(points_x)]])
        plot!(points_x, points_y;
            seriestype=:shape,
            lw=1,
            color=:red,
            fillalpha=0.2,
            label=nothing)
    end

    scatter!(x_mids, y_mids, color=:blue, markersize=4, label="Midpoints")
    scatter!(x_values, zeros(size(x_values)), color=:black, markersize=4, label="n-points")

end

"""

    integrate(method, f, a, b; threshold, numIter)

Calculates the numerical integral of the function f on interval [a,b]
using a selected method.
Iterates until the difference between two consecutive approximations of the integral
is less than the threshold or the maximum number of iterations is reached.

# Arguments
- `method::IntMethod`: The method to use for numerical integration [Trapezoid(), Midpoint(), Simpson()].
- `f::Function`: The function to integrate.
- `a::Real`: The lower bound of the integration interval.
- `b::Real`: The upper bound of the integration interval.
- `threshold::Real`: The threshold for the difference between two consecutive approximations of the integral.
- `numIter::Int`: The maximum number of iterations.

# Returns
- `integral::Real`: The numerical integral of the function f on interval [a,b].
- `error::Real`: The error upper bound of the numerical integral of the function f on interval [a,b].
- `n::Int`: The number of split points.

# Example
```julia-repl
julia> f(x) = x^2
julia> integrate(Trapezoid(), f, 0.0, 1.0, threshold = 1e-15, numIter = 25)
```
"""
function integrate(
    method::IntMethod,
    f::Function,
    a::Real,
    b::Real;
    threshold::Real=1e-8,
    numIter::Int=20
)

    n = 2
    prev_integral = numerical_int(method, f, a, b, n)
    i = 0
    while i < numIter
        n *= 2
        integral = numerical_int(method, f, a, b, n)
        difference = abs(integral - prev_integral)

        if difference < threshold
            error = error_int(method, f, a, b, n)
            return integral, error, n
        end

        prev_integral = integral
        i += 1
    end
    error = error_int(method, f, a, b, n)
    return prev_integral, error, n
end


"""

    integrate_compare(method, f, a, b; threshold, numIter)

Calculates the numerical integral of the function f on interval [a,b]
using a selected method.
Iterates until the difference between an approximations of the integral
and the reference result is less than the threshold or the maximum number of iterations is reached.

# Arguments
- `method::IntMethod`: The method to use for numerical integration [Trapezoid(), Midpoint(), Simpson()].
- `f::Function`: The function to integrate.
- `a::Real`: The lower bound of the integration interval.
- `b::Real`: The upper bound of the integration interval.
- `threshold::Real`: The threshold for the difference between two consecutive approximations of the integral.
- `numIter::Int`: The maximum number of iterations.

# Returns
- `integral::Real`: The numerical integral of the function f on interval [a,b].
- `error::Real`: The error upper bound of the numerical integral of the function f on interval [a,b].
- `n::Int`: The number of split points.

# Example
```julia-repl
julia> f(x) = x^2
julia> integrate_compare(Trapezoid(), f, 0.0, 1.0, threshold = 1e-15, numIter = 25)
```
"""
function integrate_compare(
    method::IntMethod,
    f::Function,
    a::Real,
    b::Real;
    threshold::Real=1e-8,
    numIter::Int=20
)

    ref_result, _ = quadgk(f, a, b)
    n = 1
    i = 0
    integral = numerical_int(method, f, a, b, n)
    while i < numIter
        n *= 2
        integral = numerical_int(method, f, a, b, n)
        difference = abs(ref_result - integral)

        if difference < threshold
            error = error_int(method, f, a, b, n)
            return integral, error, n
        end
        i += 1
    end
    error = error_int(method, f, a, b, n)
    return integral, error, n

end

end

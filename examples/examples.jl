using Revise
using NumericalIntegration, QuadGK

# initialization of function, bounds and number of subintervals 
f(x) = (x - 1)^5 + x^4 - 4*x + 1
a = 0.0
b = 2.0
n = 5

# Integral calculation
intf_Midpoint = numerical_int(Midpoint(), f, a, b, n)
intf_Trapezoid = numerical_int(Trapezoid(), f, a, b, n)
intf_Simpson = numerical_int(Simpson(), f, a, b, n)

# Vizualization of previous calculation
plot_int(Midpoint(), f, a, b, n)
plot_int(Trapezoid(), f, a, b, n)
plot_int(Simpson(), f, a ,b, n)

# error upper bound of previous calculations
errf_Midpoint = error_int(Midpoint(), f, a, b, n)
errf_Trapezoid = error_int(Trapezoid(), f, a, b, n)
errf_Simpson = error_int(Simpson(), f, a, b, n)


# Integral calculation with certain precision 

# Method #1
result, error = quadgk(f, a, b);
result1, error1, n1 = integrate(Trapezoid(), f, a, b);
result2, error2, n2 = integrate(Midpoint(), f, a, b);
result3, error3, n3 = integrate(Simpson(), f, a, b);

println("Reference Integral using quadgk: ", result, " with error upper bound: ", error)
println("Approximated Integral using Trapezoid method: ", result1, ", with error upper bound: ", error1, ", n: ", n1)
println("Approximated Integral using Midpoint method: ", result2, ", with error upper bound: ", error2, ", n: ", n2)
println("Approximated Integral using Simpson's method: ", result3, ", with error upper bound: ", error3, ", n: ", n3)

# Method #2
result21, error21, n21 = integrate_compare(Trapezoid(), f, a, b);
result22, error22, n22 = integrate_compare(Midpoint(), f, a, b);
result23, error23, n23 = integrate_compare(Simpson(), f, a, b);

println("Reference Integral using quadgk: ", result, " with error upper bound: ", error)
println("Approximated Integral using Trapezoid method: ", result21, ", with error upper bound: ", error21, ", n: ", n21)
println("Approximated Integral using Midpoint method: ", result22, ", with error upper bound: ", error22, ", n: ", n22)
println("Approximated Integral using Simpson's method: ", result23, ", with error upper bound: ", error23, ", n: ", n23)

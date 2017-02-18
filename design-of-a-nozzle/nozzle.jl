using NLsolve

F = 1e3
R0 = 8.314
Tc = 1500
M = 14.4e-3
γ = 1.4
pa = 1e5
pe = 3e5
pc = 100e5

function A2diam(A)
    return 2 * sqrt(A/pi)
end

ε = sqrt((γ - 1)/2) * (2 / (γ+1))^((γ+1) / 2 / (γ-1)) / (pe/pc)^(1/γ) / sqrt(1 - (pe/pc)^((γ - 1)/γ))
Vₑ = sqrt( 2*γ*R0*Tc / (γ-1) / M * (1 - (pe/pc)^((γ-1)/γ)))
Γ = sqrt(γ) * (2 / (γ+1))^((γ+1) / 2 / (γ-1));
Cs = sqrt(R0 / M * Tc) / Γ;

function eqs(x, eq)
    eq[1] = -ε^2 + 1/x[2]^2 * (2 / (γ+1) * (1 + (γ-1) / 2 * x[2]^2 ))^((γ+1) / (γ-1))
    eq[2] = -F + x[1] * Vₑ + x[3] * ε * (pe - pa)
    eq[3] = -x[1] + pc * x[3] / Cs
end

x = nlsolve(eqs, [0.4838, 2.9355, 9.2e-3]).zero
println("ε =", ε)
println("Vₑ =", Vₑ)
println("C⋆ =", Cs)
println("ṁ =", x[1]) 
println("Mₑ =", x[2])
println("⌀t =", A2diam(x[3]))
println("⌀e =", A2diam(ε * x[3]))

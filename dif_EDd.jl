include("ED_module.jl")
using .ED_dist

function projective_parametr(A::Matrix{Int})
    # gives a parametrising map of the projective variety
    (n, m) = size(A)
    Matrix{Int64}(vcat(hcat(A, zeros(n)), transpose(ones(m+1))))
end


### Example 6.6 of "The Euclidean distance degree of an algebraic variety"
A = [1 0 1;
     0 1 1]
EdA = EDdeg(A)

B = projective_parametr(A)
EdB = EDdeg(B)

EdA == EdB


### Jukes Cantor variety - quartet trees
A=[[0 1 0 0 1 1 0 1 1 1 1 1]
   [0 1 1 1 0 0 1 0 1 1 1 1]
   [1 0 0 1 1 0 1 1 0 1 1 1]
   [1 0 1 0 0 1 1 1 1 0 1 1]
   [0 0 1 1 1 1 1 1 1 1 0 1]]
EdA = EDdeg(A)
gEdA = EDdeg(A, true)

B = projective_parametr(A)
EdB = EDdeg(B)
gEdB = EDdeg(B, true)

EdA == EdB
gEdA == gEdB

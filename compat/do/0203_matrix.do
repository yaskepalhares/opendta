* Fase 2: matrix.
matrix drop _all
matrix A = (1,2\3,4)
matrix list A
matrix B = A'
matrix C = A*B
matrix list C
matrix D = A + 2*I(2)
matrix list D
matrix S = (4,2\2,3)
matrix list S
matrix Si = inv(S)
matrix list Si
matrix list Si, format(%9.4f)
matrix rownames A = a b
matrix colnames A = x y
matrix list A
display rowsof(A) " " colsof(A) " " el(A,2,1) " " A[1,2] " " trace(S) " " det(S)
local r : rownames A
display "`r'"
matrix E = A[1..2, "y"]
matrix list E
matrix F = nullmat(F) \ (1,2)
matrix F = F \ (3,4)
matrix list F
matrix v = vecdiag(S)
matrix list v
matrix J1 = J(2,3,0)
matrix list J1
capture noisily matrix X = A + J(3,3,1)
capture noisily matrix list nada
matrix dir

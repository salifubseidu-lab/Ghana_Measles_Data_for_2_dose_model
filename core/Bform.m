function v = Bform(x, a, b, P)
%BFORM Symmetric bilinear form B(a,b) = sum_jk d2F/dx_j dx_k a_j b_k (polarisation).
v = (dir_deriv(x, a + b, P, 2) - dir_deriv(x, a - b, P, 2))/4;
end

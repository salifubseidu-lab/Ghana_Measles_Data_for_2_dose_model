function v = Cform(x, a, b, c, P)
%CFORM Symmetric trilinear form C(a,b,c) by the polarisation identity
%   C(a,b,c) = (1/48) sum_{e in {+-1}^3} e1 e2 e3 D^3F(e1 a + e2 b + e3 c).
v = zeros(6,1);
for e1 = [-1 1]
    for e2 = [-1 1]
        for e3 = [-1 1]
            v = v + e1*e2*e3*dir_deriv(x, e1*a + e2*b + e3*c, P, 3);
        end
    end
end
v = v/48;
end

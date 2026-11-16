function U = oc_unpack(zv, S)
%OC_UNPACK Map scaled optimisation variables to the control matrix U ((n+1)-by-3).
U = zeros(S.n + 1, 3);
U(:, S.idx) = reshape(zv, S.n + 1, numel(S.idx)).*(ones(S.n+1,1)*S.umax(S.idx));
end

function [H, nties] = perm_entropy(s, m, tau)
% Normalised permutation entropy (Bandt & Pompe 2002) of one signal.
% nties: embedding vectors with two or more equal values (ranked in their original order).
s = s(:);
n = numel(s) - (m - 1) * tau;
E = zeros(n, m);
for i = 1:m
    E(:, i) = s((i-1)*tau + (1:n));
end
[Es, idx] = sort(E, 2);
nties = nnz(any(diff(Es, 1, 2) == 0, 2));
code = (idx - 1) * (m .^ (0:m-1))';
p = accumarray(code + 1, 1, [m^m 1]);
p = p(p > 0) / n;
H = -sum(p .* log(p)) / log(factorial(m));
end

function r = dirichletrnd3x(alpha, n)

A = repmat(alpha, [1, 1, n]);
r = randg(A);
sr = sum(r, 2);
r = r./sr;


end

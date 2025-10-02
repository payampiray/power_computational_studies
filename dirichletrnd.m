
function r = dirichletrnd(num_sim, alf)

A = repmat(alf, [num_sim, 1]);
r = randg(A);
sr = sum(r, 2);
r = r./sr;

end
function U = vuCSReconIST (F, K, Idim, rofwgt, roftol, rectol)
%%MRISTA Iterative shrinkage-threshold MR reconstruction.

n = prod(Idim);

Phi = @(x) reshape(roffgp(abs(reshape(x,Idim)), rofwgt, roftol), n, []);

f = F(:);
G = G_nuf(K, Idim);

%cond = max(eig(reshape(G' * G * eye(Idim(1)), Idim)));
cond = norm(reshape(G' * G * ones(Idim), Idim));  % <-- new
fprintf('nufft condition number: %g\n', cond);

up = zeros(n,1);
u = zeros(n,1);

%t = 1e-2 / cond;
t = 0.99 / cond;    % <-- new
maxiter = 1000;

for i = 1:maxiter
  u = u + t * G' * (f - G * u); 
  u = Phi(u);
  chg = norm(u - up) / norm(u + up);
  if chg <= rectol
    break;
  end
  fprintf('chg:%g\n', chg);
  up = u;
end
if i == maxiter
  disp('MRISTA has reached the maximum allowed number of steps.');
end
U = abs(reshape(u, Idim));

end
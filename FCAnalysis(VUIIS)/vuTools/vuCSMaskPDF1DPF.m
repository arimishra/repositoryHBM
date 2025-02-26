function P = vuCSMaskPDF1DPF (n, norm, q, pf)
%%VUCSMASKPDF Symmetric array of sampling probabilities.

if nargin < 4, pf = 1; end

ks = (1:n) - ceil(n/2) - 1;
kmax = floor(n/2);
npf = round(pf*n);
klo = ks(n-npf+1);

for kw = 1:kmax
  P = pdf(ks,kw,klo,q);
  if sum(P) >= norm, break; end
end

P = fftshift(P);
if mod(n,2), P = [1; P]; end

function p = pdf (k, kw, klo, q)
p = (abs(k)/kw).^(-q);
p(k==0) = 0;
p(abs(k) <= kw) = 1;
p(k < klo) = 0;
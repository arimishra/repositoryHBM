function P = vuCSMaskPDF1D (n, norm, q)
%%VUCSMASKPDF Symmetric array of sampling probabilities.


kmax = floor(n / 2);
if mod(n,2), norm = norm - 1; end
norm = norm / 2;

if q == 1
  a = norm / log(kmax);
else
  a = norm*(1-q) / (kmax^(1-q) - 1);
end

if a < 1  % no window required
  P = a * (1:kmax)'.^-q;
else
  if q == 1
    f = @(x) x*(1 + log(kmax/x)) - norm;
    %f = @(x) x + (norm-x)*log(kmax/x)/log(kmax) - norm;
  else
    a = q / kmax^(1-q);
    b = norm * (1-q) / kmax^(1-q);
    f = @(x) x^q - a * x - b;
  end
  
  %options = optimset('Display','iter');
  kw = round(fzero(f, norm));
  P = ones(kmax,1);
  P(kw:kmax) = ((kw:kmax) / kw) .^ (-q);
  
end

P = [P; P(end:-1:1)];
if mod(n,2), P = [1; P]; end

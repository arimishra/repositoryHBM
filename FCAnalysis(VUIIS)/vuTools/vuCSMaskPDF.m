function P = vuCSMaskPDF (dim, norm, q)
%%VUCSMASKPDF Symmetric array of sampling probabilities.

if nargin < 3, q = 1; end

n = dim(1);
if mod(n,2), norm = norm - 1; end
kmax = floor(n / 2);
norm = norm / 2;

if q == 1
	s = @(x) x*(1 + log(kmax/x)) - norm;
else
	a = q / kmax^(1-q);
	b = norm * (1-q) / kmax^(1-q);
	s = @(x) x^q - a * x - b;
end

P = ones(kmax,1);
for kw = 1:kmax
	if s(kw) >= 0, break; end
end
P(kw:kmax) = ((kw:kmax) / kw) .^ (-q);

P = [P; P(end:-1:1)];
if mod(n,2), P = [1; P]; end





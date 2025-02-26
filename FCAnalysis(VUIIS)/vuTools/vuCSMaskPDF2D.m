function P = vuCSMaskPDF2D (dims, norm, q)
%%VUCSMASKPDF Symmetric array of sampling probabilities.

if nargin < 3, q = 1; end

nz = dims(2);
ny = dims(1);
yc = round(ny/2);
zc = round(nz/2);
rmax = sqrt((ny-yc)^2 + (nz-zc)^2);

[Z Y] = meshgrid(1:nz,1:ny);
R = sqrt((Y-yc).^2+(Z-zc).^2);

for rw = 1:round(rmax)
  P = ones(ny,nz);
  P(R >= rw) = (R(R >= rw) / rw) .^ (-q);
  if sum(P(:)) >= norm, break; end
end



function P = vuCSMaskPDF2DPI (dims, norm, q, sense_factor)
%%VUCSMASKPDFPI Symmetric array of sampling probabilities.

if nargin < 4
  sense_factor = 1;
  if nargin < 3
    q = 1; 
  end
end

nz = dims(2);
ny = dims(1);
yc = round(ny/2);
zc = round(nz/2);
rmax = sqrt((ny-yc)^2 + (nz-zc)^2);

[Z Y] = meshgrid(1:nz,1:ny);

R = sqrt((Y-yc).^2+(Z-zc).^2);

Z = abs(Z - zc);
Y = abs(Y - yc);

for rw = 1:round(rmax)
  P = ones(ny,nz)/sense_factor;
  C = Z <= rw & Y <= rw;
  W = Z > rw | Y > rw;
  P(W) = (R(W) / rw) .^ (-q);
  if sum(P(:)) >= norm, break; end
end
if sense_factor  > 1
  D = C;
  for s = 2:sense_factor
    D(s:sense_factor:end,:) = 0;
  end
  P(C) = 0;
  P(D) = 1;
end






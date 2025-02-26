function M = vuCSRandMask (dim, accel, q, seed, pf)
%%VUCSRANDMASK  Random binary phase encode mask.

%  HISTORY
%    20101018 :: dss :: Initial version
%    20110210 :: dss :: Bug and robustness fixes.
%    20110404 :: dss :: Added 3D.
%    20110606 :: dss :: Added partial Fourier for 2D and SENSE for 3D.

if nargin < 3, q = 1; end
if nargin < 4, seed = 11235; end
if nargin < 5, pf = 1; end

rs = RandStream.create('mt19937ar','Seed',seed);
RandStream.setDefaultStream(rs);

if numel(dim) == 2   % 2D mask
  
  nro = dim(1);
  npe = dim(2);
  nacq = round(npe/accel);  % requested # of lines in acquisition
  
  P = vuCSMaskPDF1DPF(npe, nacq, q, pf)';
  %P = vuCSMaskPDF1D(npe, nacq, q);
  
  while true
    M = rand(npe,1) <= P;
    if sum(M) == nacq, break; end
  end
  
  % remove partial Fourier plane and compensate sampling density
  
  
  M = logical(M);
  M = repmat(M,[1 nro]);
  
elseif numel(dim) == 3  % 3D mask
  
  nro = dim(1);
  npe1 = dim(2);
  npe2 = dim(3);
  
  nacq = round(npe1*npe2/accel);
  %P = vuCSMaskPDF2D([npe1 npe2], nacq, q);
  P = vuCSMaskPDF2DPI([npe1 npe2], nacq, q, pf);
  R = rand(npe1,npe2);
  M = R <= P;

  nchosen = sum(M(:));  % correct for inexact number chosen
  if nchosen > nacq
    outerOn = find(M & P ~= 1);
    numToFlip = nchosen - nacq;
    idxs = randperm(numel(outerOn));
    M(outerOn(idxs(1:numToFlip))) = false;
  elseif nchosen < nacq
    outerOff = find(~M);
    idxs = randperm(numel(outerOff));
    numToFlip = nacq - nchosen;
    M(outerOff(idxs(1:numToFlip))) = true;
  end
  
  
  M = shiftdim(repmat(M, [1 1 nro]), 2);
  M = ifftshift(M);
else 
  error('VUCSRANDMASK> Mask dimension must be 2 or 3.');
end

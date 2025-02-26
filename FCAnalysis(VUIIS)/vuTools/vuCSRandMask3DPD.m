function [M C] = vuCSRandMask3DPD (dim, accel, q, seed, pf)
%%VUCSRANDMASK3DPD  Random binary phase encode mask.

%  HISTORY
%    20101018 :: dss :: Initial version
%    20110210 :: dss :: Bug and robustness fixes.
%    20110404 :: dss :: Added 3D.
%    20110606 :: dss :: Added partial Fourier for 2D and SENSE for 3D.
%    20110815 :: dss :: Added Poisson disk-type pattern for 3D.
%
%
%  USAGE
%    For a Poisson disk pattern, pass the DIMs as usual.  Pass the accel
%    as sqrt of the desired accel (e.g. ACCEL=2 if you want 4X).  
%    Finally, Q is the width of the fully sampled center window.  PF is not
%    used, but SEED can be set.


if nargin < 4, seed = 11235; end
if nargin < 5, pf = 1; end

rs = RandStream.create('mt19937ar','Seed',seed);
RandStream.setDefaultStream(rs);

if numel(dim) == 2   % 2D mask
  if nargin < 3, q = 1; end
  nro = dim(1);
  npe = dim(2);
  nacq = round(npe/accel);  % requested # of lines in acquisition
  
  P = vuCSMaskPDF1DPF(npe, nacq, q, pf)';
  C =  [];  % hackity-hack
  %P = vuCSMaskPDF1D(npe, nacq, q);

  while true
    M = rand(npe,1) <= P;
    if sum(M) == nacq, break; end
  end
  
  % remove partial Fourier plane and compensate sampling density
  
  
  M = logical(M);
  M = repmat(M,[1 nro]);
  
elseif numel(dim) == 3  % 3D mask
  
  if nargin < 3, q = 10; end
  nro = dim(2);
  npe1 = dim(1);
  npe2 = dim(3);
  
  nacq = round(npe1*npe2/accel);
  %P = vuCSMaskPDF2D([npe1 npe2], nacq, q);
  twiddle = accel / 2;
  M = false(npe1,npe2);
  M(1:accel:end,1:accel:end) = true;
  
  [jnz knz] = find(M);
  
  jnz = jnz + randi(3,numel(jnz),1) - 2;
  knz = knz + randi(3,numel(knz),1) - 2;
  jnz = max(jnz,1);
  knz = max(knz,1);
  jnz = min(jnz,npe1);
  knz = min(knz,npe2);
  
  idxs = sub2ind([npe1 npe2],jnz,knz);
  M = false(npe1,npe2);
  M(idxs) = true;
  M([1:q end-q+1:end],[1:q end-q+1:end]) = true;
  
  % ky x kz > kz x ky > kz x ky x kx > ky x kx x kz
  M = shiftdim(repmat(shiftdim(M,1), [1 1 nro]), 1);
  M = ifftshift(M);
  if nargout > 1
    C = shiftdim(repmat(shiftdim(C,1), [1 1 nro]), 1);
    C = ifftshift(C);
  end
else 
  error('VUCSRANDMASK> Mask dimension must be 2 or 3.');
end

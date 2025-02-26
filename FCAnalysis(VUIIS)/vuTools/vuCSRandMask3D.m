function [M P] = vuCSRandMask3D (dim, accel, q, seed)
%%VUCSRANDMASK3D  Random binary mask for 3D scans.
%
%  HISTORY
%    20110401 :: dss :: Initial version, based on vuCSRandMask2
%
%  Strategy is to build up 2D masks in slice direction that each have a
%  number of phase encodes that are dependent on the slice number, so that
%  the outer slices are more undersampled than the inner slices.

npe1 = dim(1);
nro = dim(2);
npe2 = dim(3);

if nargin < 3, q = 1; end
if nargin < 4, seed = 11235; end

rs = RandStream.create('mt19937ar','Seed',seed);
RandStream.setDefaultStream(rs);

nacq = round(npe1*npe2/accel);
P = vuCSMaskPDF2D([npe1 npe2], nacq, q);
M = rand(npe1,npe2) <= P;
M = logical(M);
M = repmat(M, [1 1 nro]);
M = shiftdim(M,2);

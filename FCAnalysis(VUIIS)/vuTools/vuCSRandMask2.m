function M = vuCSRandMask2 (dim, accel, q, seed)
%%VUCSRANDMASK  Random binary mask.
%
%  vuCSRandMask(dim,fu,factor) selects a fraction fu of k-space to create
%  the mask.  k-space is divided into lines, of which fu are kept.
%
%  GPU: needs testing

%  HISTORY
%    20101018 :: dss :: Initial version
%    20110210 :: dss :: Bug and robustness fixes.
%
%  TODO:
%  make work with accel = 1 and high accel, asymptote to DC only.

% problem is when no central window but sum(P) is still much larger than
% nacq.  What to do?

npe = dim(1);
nro = dim(2);

nacq = round(npe/accel);  % requested # of lines in acquisition

if nargin < 3, q = 1; end
if nargin < 4, seed = 11235; end

rs = RandStream.create('mt19937ar','Seed',seed);
RandStream.setDefaultStream(rs);

P = vuCSMaskPDF(npe, nacq, q);

while true
	M = rand(npe,1) <= P;
	if sum(M) == nacq, break; end
end

M = logical(M);
M = repmat(M,[1 nro]);
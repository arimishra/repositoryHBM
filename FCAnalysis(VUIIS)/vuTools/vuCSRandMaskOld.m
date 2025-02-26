function M = vuCSRandMaskOld (dim, accel, q, seed)
%%VUCSRANDMASK  2-D random binary mask.
%
%  HISTORY
%    20101018 :: dss :: Initial version
%    20110405 :: dss :: Changed undersample fraction 'fu' to acceleration
%                       'accel'.

if nargin < 3, q = 1; end
if nargin < 4, seed = 11235; end

rs = RandStream.create('mt19937ar','Seed',seed);
RandStream.setDefaultStream(rs)

npe = dim(1);
nro = dim(2);
K = npe / 2;

P = zeros(K,1);
P(end) = 2; % ensure highest k is sampled

for k = 1:K-1
  P(k) = 2;
  P(k+1:(K-1)) = 2 * (((k+1):(K-1)) / (k+1)) .^ (-q);
  if sum(P) > round(npe/accel), break; end
end
P = [P(end:-1:1)' P']';
P = P / 2;

while true
  M = zeros(npe,1);
  r = rand(npe,1) <= P;
  M(r) = 1;
  if sum(M,1) == round(npe/accel), break; end
end

M = repmat(M,[1 nro]);

M = logical(ifftshift(M));
% figure(1);
% imagesc(M);
% figure(2);
% stem(P);


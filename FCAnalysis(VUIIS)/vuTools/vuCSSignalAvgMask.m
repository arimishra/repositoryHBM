function [M krepeat] = vuCSSignalAvgMask (dim, nctr, nsa, seed)
%%VUCSSIGNALAVGMASK  Mask which samples center 2x instead of random outer k lines.
%
%  vuCSSignalAvgMask(dim,nctr,seed) drops nctr k lines from the group of lines 
%  outside the centermost nctr lines.
%
%  GPU: needs testing

%  HISTORY
%    20101031 :: dss :: Initial version
%    20110107 :: dss :: Added high freq multisampling

if nargin < 3, nsa = 2; end
if nargin < 4, seed = 11235; end

rs = RandStream.create('mt19937ar','Seed',seed);
RandStream.setDefaultStream(rs)

rows = dim(1);
cols = dim(2);

N = rows - nctr;
fu = (N - nctr*(nsa-1)) / N;

M = vuCSRandMask([N N], fu, 0.5, seed);
M = M(:,1);
M = padarray(M, nctr/2, true);
M = repmat(M, [1 cols]);

krepeat = [1:round(nctr/2)  rows-nctr+round(nctr/2)+1:rows];
function [A B] = vuCSNormScale (A, B)
%%VUCSNORMSCALE Scale arrays to be roughly normal and have the minimum norm
%%difference.

nz = size(A,3);

% normalize A
for k = 1:nz
  T = A(:,:,k);
  A(:,:,k) = (T - min(T(:))) / range(T(:));
  
  [e a b] = vuCSMinNormDiff(B(:,:,k), A(:,:,k));  % find scaling for min norm
  fprintf('err = %g\n', e);
  B(:,:,k) = a*B(:,:,k) + b; % perform the scaling
end


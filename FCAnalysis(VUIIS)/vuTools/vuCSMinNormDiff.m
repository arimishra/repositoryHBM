function [e x] = vuCSMinNormDiff (I1, I2, p)
%%VUCSMINNORMDIFF  Least norm of image difference, allowing free scaling.
%
% [NORM COEFS] = VUCSMINNORMDIFF(IMG1,IMG2)  Scales IMG1 until the minimum
% L2 norm of the difference between IMG1 and IMG2 is found, then returns
% that norm and the scaling coefficients.
%
% Not GPU compatible yet.

% HISTORY
%   20101018 :: dss :: Initial revision. 
%   20101105 :: dss :: Removed cvx for L2 norms.
%
% TODO

if nargin < 3, p = 2; end

I1 = double(I1(:));
I2 = double(I2(:));

if p == 2
  A = ones(numel(I1),2);
  A(:,1) = I1;
  x = A \ I2;
else
  cvx_quiet true;
  cvx_begin
    variable x(2);
    minimize norm(x(1)*I1 - I2 + x(2), p);
    subject to
    x(1) >= 0;
  cvx_end
end

e = norm(x(1)*I1 - I2 + x(2), p)^p / numel(I1)^(1/p) / range(I2(:));
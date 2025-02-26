function [D e] = vuCSDiffImage (I1, I2, mask)
%%VUCSDIFFIMAGE Returns the RMS minimum difference image.
%
%  D = VUCSDIFFIMAGE(IMAGE1,IMAGE2) returns the difference image that
%  minimizes the L2 norm of the difference between the images.  IMAGE2 is
%  considered the reference and is scaled to the [0,1] interval before the
%  best scaling for IMAGE1 is found.
%
%  D = VUCSDIFFIMAGE(IMAGE1,IMAGE2,MASK) only considers the area in MASK in
%  minimizing the norm.
%
%  [D err] = VUCSDIFFIMAGE(IMAGE1,IMAGE2) returns the difference image and
%  the RMS difference.
%
%  Written by David S. Smith, 12 Oct 2010.
%
%  History
%    v1.0 - Original release

D = zeros(size(I1));
if nargin == 3
  I1 = I1(mask);
  I2 = I2(mask);
end
I2 = normalize(I2);
[e x] = vuCSMinNormDiff(I1, I2);
if exist('mask')
  D(mask) = x(1)*I1 + x(2) - I2;
else
  D = x(1)*I1 + x(2) - I2;
end
end

function I = normalize (I)
I = (I - min(I(:))) / range(I(:));
end
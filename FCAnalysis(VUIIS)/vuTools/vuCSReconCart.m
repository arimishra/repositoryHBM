function [U time] = vuCSReconCart (F, M, user_opts)
%%CSRECON2DCART Recon Cartesian data from partial Fourier measurements.
%
% IMG = VUCSRECONCART(DATA,MASK,OPTS) reconstructs IMG from DATA given
% partial Fourier MASK and recon options OPTS.
%
% GPU compatible
%
% Copyright (c) 2010 - Vanderbilt University Institute of Imaging Science

% HISTORY
%   20101018 :: dss :: Initial revision.
%   20110404 :: dss :: Cleanup, added 3D.

global opts

opts = vuCSReconDefaults;

if nargin == 3 && ~strcmp(class(user_opts), 'double')
  user_fields = fieldnames(user_opts);
  for k = 1:length(user_fields)  % override defaults with any specified
    fld = user_fields(k);
    opts.(fld{1}) = user_opts.(fld{1});
  end
end

if numel(M) == 0  % guess a mask
	M = abs(F) ~= 0;
end

time = -1;
if numel(size(F)) == 2
		[U time] = vuCSReconSB2D(F, M);
elseif numel(size(F)) == 3
	U = vuCSReconSB3D(F, M);
else
	error('vuCSReconCart> Collected Fourier data must have 2 or 3 dimensions.');
end

if ~opts.complex_image
  U = abs(U);
end

function opts = vuCSReconDefaults
%%vuCSDefaultOpts Default option structure.
%
% VUCSDEFAULTOPTS returns a structure filled with default options for the
% CS recon. Users can modify this file to change default recon parameters.

% HISTORY
%   20101018 :: dss :: Initial revision.


opts = [];

opts.verbose = false;

opts.complex_image = false;
opts.nu = 1;
opts.mu = 20;  % controls smoothing; higher reduces smoothing
opts.lambda = 1;

% quality = 0
%opts.ninner = 3;
%opts.nouter = 10;

% quality = 1   
opts.ninner = 5;
opts.nouter = 10; 

% quality = 2
% opts.ninner 7;
% opts.nouter = 15;

% quality = 3
% opts.ninner = 20;
% opts.nouter = 20;

opts.gamma = 1e-6;
opts.norm_p = 1; % L1 norm is the default


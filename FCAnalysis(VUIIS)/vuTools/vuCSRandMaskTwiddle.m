function M = vuCSRandMaskTwiddle (dims, yaccel, zaccel, w, omitCorners)
%%vuCSRandMaskTwiddle  Regularly spaced points randomly twiddled by -1,0,1
% in one direction.
%
% M = vuCSRandMaskTwiddle(dims,yaccel,zaccel,w,omitCorners) returns a
% binary acquisition mask in Matlab's Fourier ordering.  
%
%  dims: [ny nx nz] data size. ny is the 1st PE dimension, nx is the
%  readout dimension, and nz is the 2nd PE dimension.
%
%  yaccel, zaccel: acceleration factors in each PE direction, similar to
%  SENSE factors
%
%  w: OPTIONAL radius of central fully sampled disk.  DEFAULT: 10.
%  
%  omitCorners: OPTIONAL boolean option to trim off corners from mask, leading 
%  to an inscribed circle pattern that could be useful with homodyning in the
%  future.  DEFAULT: false.

nx = dims(1);
ny = dims(2);
nz = dims(3);

if nargin < 4
  w = 10;
end
if nargin < 5;
  omitCorners = false;
end

pad = 2;
ny = ny + 2*pad;
nz = nz + 2*pad;

% Requested acceleration factors
totalAccel = yaccel*zaccel;
fprintf('Requested acceleration: %g\n', yaccel*zaccel);

% correct for center window's effect
effAccel = (ny*nz - round(pi*w^2)) / (ny*nz / totalAccel - round(pi*w^2));
fprintf('Effective outer acceleration: %g\n', effAccel);

%% generate a regular sampling grid
[Z Y] = meshgrid(round(1:zaccel:nz),round(1:yaccel:ny));
idxs = sub2ind([ny nz], Y, Z);
M0 = false(ny,nz);
M0(idxs) = true;

%% add central window
M0(ny/2-w+1:ny/2+w,nz/2-w+1:nz/2+w) = true;

%% generate displacements in a random direction by -1, 0, or 1.
Ty = randi([-1,1],size(Y));
Tz = randi([-1,1],size(Z));
D = logical(randi([0,1], size(Y)));
Ty(D) = 0;
Tz(~D) = 0;
Yt = Y + Ty;
Zt = Z + Tz;
Yt(Yt < 1) = 1;
Yt(Yt > ny) = ny;
Zt(Zt < 1) = 1;
Zt(Zt > nz) = nz;

idxsTwiddled = sub2ind([ny nz], Yt, Zt);
M = false(ny,nz);
M(idxsTwiddled) = true;

%% add central window and zero corners
[Zfull Yfull] = meshgrid(1:nz,1:ny);
R2 = (Yfull - ny/2).^2 + (Zfull - nz/2).^2;
%M(ny/2-w+1:ny/2+w,nz/2-w+1:nz/2+w) = true;
M(R2 <= w^2) = true;
if omitCorners
  fprintf('>> corners omitted\n');
  M(R2 >= max(ny,nz)^2/4) = false;
end

%% crop back down to original size
M = M(pad+1:end-pad,pad+1:end-pad);
M0 = M0(pad+1:end-pad,pad+1:end-pad);
ny = ny - 2*pad;
nz = nz - 2*pad;

%% correct for inexact number of PEs
PEdiff = round(numel(M) / totalAccel) - nnz(M);
[Zfull Yfull] = meshgrid(1:nz,1:ny);
R2 = (Yfull - ny/2).^2 + (Zfull - nz/2).^2;
OuterMask = R2 > w^2;
if PEdiff < 0                     % placed too many
  idxs = find(OuterMask & M);
  r = randperm(numel(idxs));
  M(idxs(r(1:-PEdiff))) = false;  % turn off some randomly
elseif PEdiff > 0                 % placed too few
  idxs = find(OuterMask & ~M);
  r = randperm(numel(idxs));
  M(idxs(r(1:PEdiff))) = true;   % turn on some randomly
end


%% plot it
%imagesc([ M]);
%axis image;
%colormap(gray);

% ky x kz > kz x ky > kz x ky x kx > ky x kx x kz
%M = shiftdim(repmat(shiftdim(M,1), [1 1 nx]), 1);
M = shiftdim(repmat(M, [1 1 nx]),2);
M = ifftshift(M);

  
fprintf('Delivered regular acceleration: %g\n', numel(M0) / nnz(M0));
fprintf('Delivered twiddled acceleration: %g\n', numel(M) / nnz(M));





function U = vuCSReconSB3D (F, M)
%%VUCSRECONSB3D 3-D TV-regularized partial Fourier recon.


global opts;

lambda = opts.lambda;
mu = opts.mu;
gamma = opts.gamma;
nu = opts.nu;

[rows cols slices] = size(F);

% normalize the data so that standard parameter values work
normFactor = getNormalizationFactor(F,M);
M = glogical(M);
F = gsingle(F);
F = normFactor * F;

% Reserve memory for the auxillary variables
F0 = F;
U = gzeros(rows,cols,slices,'single');
X = gzeros(rows,cols,slices,'single');
Bx = gzeros(rows,cols,slices,'single');
Y = gzeros(rows,cols,slices,'single');
By = gzeros(rows,cols,slices,'single');
Z = gzeros(rows,cols,slices,'single');
Bz = gzeros(rows,cols,slices,'single');

% Build Kernels
scale = sqrt(rows*cols);
muMF = ifftn(mu * M .* F) * scale;
% 
% uker = gzeros(rows,cols,'single');
% uker(1,1) = 4;
% uker(1,2) = -1;
% uker(2,1) = -1;
% uker(rows,1) = -1;
% uker(1,cols) = -1;
% uker = 1 ./ (mu * R + lambda * fft2(uker) + gamma);
K = gzeros(rows,cols,slices,'single');
K(1,1,1) = 8;
K(2,1,1) = -1;
K(1,2,1) = -1;
K(1,1,2) = -1;
K(rows,1,1) = -1;
K(1,cols,1) = -1;
K(1,1,slices) = -1;
K = 1 ./ (mu * M + lambda * fftn(K) + gamma);

%  Do the reconstruction
for outer = 1:opts.nouter
  for inner = 1:opts.ninner
    % update u
    if opts.verbose
      fprintf('%d/%d\n', inner+(outer-1)*opts.ninner, ...
        opts.nouter*opts.ninner);
    end
    RHS = muMF + lambda*Dxt(X-Bx) + lambda*Dyt(Y-By) + ...
      lambda*Dzt(Z-Bz)/nu + gamma*U;
    U = ifftn(fftn(RHS) .* K);
    
    % update x and y
    Ax = Dx(U) + Bx;
    Ay = Dy(U) + By;
    Az = nu * Dz(U) + Bz;
    [X Y Z] = shrink3(Ax, Ay, Az, 1/lambda, opts.norm_p);
    
    % update bregman parameters
    Bx = Ax - X;
    By = Ay - Y;
    Bz = Az - Z; 
  end
  
  F = F + F0 - M .* fftn(U) / scale;
  muMF = ifftn(mu * M .* F) * scale;
end

% undo the normalization so that results are scaled properly
U = single(U / normFactor / scale);
end

function normFactor = getNormalizationFactor(F,M)
normFactor = 1/norm(F(:)/size(M==1,1));
end

function d = Dx(u)
cols = size(u,2);
d = u - u(:,[cols 1:cols-1],:);
end

function d = Dxt(u)
cols = size(u,2);
d = u - u(:,[2:cols 1],:);
end

function d = Dy(u)
rows = size(u,1);
d = u - u([rows 1:rows-1],:,:);
end

function d = Dyt(u)
rows = size(u,1);
d = u - u([2:rows 1],:,:);  
end

function d = Dz(u)
nt = size(u,3);
d = u - u(:,:,[nt 1:nt-1]);
%d = zeros(size(u));
end

function d = Dzt(u)
nt = size(u,3);
d = u - u(:,:,[2:nt 1]);  
%d = zeros(size(u));
end


function [xs ys zs] = shrink3(x, y, z, lambda, p)
s = sqrt(x.*conj(x) + y.*conj(y) + z.*conj(z));
t = lambda ./ (s .^ (1 - p));
%t = lambda ./ sqrt(s);
ss = s - t;
ss = ss .* (ss > 0);
s = s + (s < t);
ss = ss ./ s;
xs = ss .* x;
ys = ss .* y;
zs = ss .* z;
end

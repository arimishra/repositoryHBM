function [u time] = vuCSReconSB2D (f, R)
%%VUCSRECONSBWLP GPU-based split Bregman Lp-regularized least squares
%  reconstructor with TV and Wavelet Lp-norm penalties.
%
% Original Help:
%           mrics.m  by Tom Goldstein (TomGoldstein1@gmail.com)
%     This file contains methods for performing compressed sensing
%  recontructions of images from k-space data using the Split Bregman
%  method.
%     To use the method, simply add this "m" file to your current directory,
%  and then call the following method:
%
%              u = mrics(R,F, mu, lambda, gamma, nInner, nOuter);
%
%  The inputs to this function are described below:
%
%  R - This is a matrix which determines which elements of K-Space
%           are known.  Take R(m,n)=0 if Fourier mode (m,n) is unknown,
%           and R(m,n)=1 is the corresponding mode is known.
%  F - This is the K-Space (Fourier) data.  In other words, F is the
%           Fourier transform of the image you wish to recover.  If a
%           Fourier mode is known, then it should have a non-zero value.
%           If a Fourier mode is unknown, then simple set the corresponding
%           entry in this matrix to zero.  If you have set the values
%           in this matrix properly, then you should have (R.*F==F).
%           The method will automatically scale this data so that moderate
%           values of the function paramters (see below) will be
%           appropriate.
%  mu- The parameter on the fidelity term in the Split Bregman method.
%      The method automatically scales the input data so that extreme value
%      of this variable are not necessary.  For this reason, mu=1 should
%      work for many applications.
% lambda - The coefficient of the constraint term in the Split Bregman
%      model.  Becuase of the way that the data is scaled, I suggest using
%      lambda=mu=1.
% gamma - This is a regularization parameter.  I suggest that you take
%      gamma = mu/100.
% nInner - This determines how many "inner" loops the Split Bregman method
%      performs (i.e. loop to enforce the constraint term).  I suggest
%      using nInner = 30 to be safe.  This will usually guarantee good
%      convergence, but will make things a bit slow.  You may find that you
%      can get away with nInner = 5-10
% nOuter - The number of outer (fidelity term) Bregman Iterations.  This
%      parameter depends on how noisy your data is, but I find that
%      nOuter=5 is usually about right.

global opts;

t = tic;

lambda = opts.lambda;
mu = opts.mu;
gamma = opts.gamma;
nInner = opts.ninner;
nBreg = opts.nouter;

[rows cols] = size(f);


% normalize the data so that standard parameter values work
normFactor = getNormalizationFactor(R,f);

R = gsingle(R);
f = gsingle(f);
f = normFactor * f;

% Reserve memory for the auxillary variables
f0 = f;
u = gzeros(rows,cols,'single');
x = gzeros(rows,cols,'single');
bx = gzeros(rows,cols,'single');
y = gzeros(rows,cols,'single');
by = gzeros(rows,cols,'single');

% Build Kernels
scale = sqrt(rows*cols);
murf = ifft2(mu * R .* f) * scale;

uker = gzeros(rows,cols,'single');
uker(1,1) = 4;
uker(1,2) = -1;
uker(2,1) = -1;
uker(rows,1) = -1;
uker(1,cols) = -1;
uker = 1 ./ (mu * R + lambda * fft2(uker) + gamma);


%  Do the reconstruction
for outer = 1:nBreg;
  for inner = 1:nInner;
    % update u
    rhs = murf + lambda*Dxt(x-bx) + lambda*Dyt(y-by) + gamma*u;
    u = ifft2(fft2(rhs) .* uker);
    
    % update x and y
    ax = Dx(u) + bx;
    ay = Dy(u) + by;
    [x y] = shrink2(ax, ay, 1/lambda);
    
    % update bregman parameters
    bx = ax - x;
    by = ay - y;
  end
  
  f = f + f0 - R .* fft2(u) / scale;
  murf = ifft2(mu * R .* f) * scale;
end

% undo the normalization so that results are scaled properly
u = single(u / normFactor / scale);

time = toc(t);

end


function normFactor = getNormalizationFactor(R,f)
normFactor = 1/norm(f(:)/size(R==1,1));
end

function d = Dx(u)
cols = size(u,2);
d = u - u(:,[cols 1:cols-1]);
end

function d = Dxt(u)
cols = size(u,2);
d = u - u(:,[2:cols 1]);
end

function d = Dy(u)
rows = size(u,1);
d = u - u([rows 1:rows-1],:);
end

function d = Dyt(u)
rows = size(u,1);
d = u - u([2:rows 1],:);  
end


function [xs ys] = shrink2(x, y, lambda)
global opts
s = sqrt(x.*conj(x) + y.*conj(y));
t = lambda ./ (s .^ (1 - opts.norm_p));
%t = lambda ./ sqrt(s);
ss = s - t;
ss = ss .* (ss > 0);
s = s + (s < t);
ss = ss ./ s;
xs = ss .* x;
ys = ss .* y;
end

function xs = shrink(x, gamma)
global opts
s = abs(x);
t = gamma ./ (s .^ (1 - opts.norm_p ));
%t = gamma ./ sqrt(s);
ss = s - t;
ss = ss .* (ss > 0);
s = s + (s < t);
ss = ss ./ s;
xs = ss .* x;
end

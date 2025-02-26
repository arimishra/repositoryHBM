function M = vuCSEPIMask (dims, accel, nctr, blipMax, seed)
%%VUCSEPIMASK  Random acquisition with maximum gap between adjacent phase encodes.

if numel(dims) > 2
  error('3D not supported yet.');
end

if nargin < 5, seed = 11235; end

rs = RandStream.create('mt19937ar','Seed',seed);
RandStream.setDefaultStream(rs);

n = dims(1);

fprintf('High freq default gap is %d.\n', accel);

M = false(n);
idxs0 = nctr+1:accel:n-nctr;

% want 2*twiddleMax + accel <= blipMax

twiddleMax = 0.5*(blipMax - accel);
fprintf('Max twiddled gap possible is %d.\n', accel+2*twiddleMax);

twiddles = randi(2*twiddleMax+1,1,numel(idxs0)) - floor(twiddleMax) - 1;
idxs = idxs0 + twiddles;

M(idxs,:) = true;
M([1:nctr n-nctr+1:n],:) = true;

I = vuCSReconCart(M.*fft2(phantom(n)), M);
imagesc([M; I]);

fprintf('Effective acceleration: %g\n', n / nnz(M(:,1)));
Msense = M;
Msense([2:2:nctr n-nctr+2:2:n],:) = false;
fprintf('Effective acceleration with central SENSE 2x: %g\n', n / nnz(Msense(:,1)));
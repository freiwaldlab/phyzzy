%   RBFipol.m
%   =========
%   Radial basis function interpolation with optional normalization.
%   This function takes the pre-calculated interpolation coefficients
%   lambda for the RBFs centered at the DIM-dimensional coordinates x
%   to interpolate the function fip over the M interpolation points
%   xip, using Gaussians of width sigma.
%
%   Usage: fip = RBFipol(lambda, x, xip, sigm, normalize);
%       lambda(N, 1)            : interpolation coefficients from RBFcoef.m
%       x(DIM, N)               : sample coordinates
%       xip(DIM, M)             : test point coordinates
%       sigm                    : width of Gaussians
%       normalize               : (optional, default=0) normalize to yield
%                                 a 'better' interpolation
%
%       fip(N, 1)               : interpolation at test points
%

function fip = RBFipol(lambda, x, xip, sigm, normalize)


% process optional arguments
if ~exist('normalize') normalize = 0; end


% check for correct argument dimensionality
s = size(x);
sip = size(xip);
dim = s(1);
if (dim ~= sip(1))
	error('x and xip must have same dimensionality')
end
sl = size(lambda);
n = s(2);
m = sip(2);


% interpolate the data points
norm = 1/(0.5*sigm.^2);
idx_m = 1:m;
onesn = ones([1 n]);
idx_n = 1:n;
fip = zeros([m 1]);
%fprintf('Interpolating %i data points: ', m);
for i=idx_m
	%fprintf('.');
	d = sum((x(:, idx_n)-xip(:, i*onesn)).^2, 1);
	e = exp(-norm*d);
	if (normalize == 1)
		fip(i) = sum(lambda'.*e)./sum(e);
	else
		fip(i) = sum(lambda'.*e);
	end
end
%fprintf('\n');


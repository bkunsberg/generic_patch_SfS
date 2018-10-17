%This function creates a test surface from which to test my DjN method.


function [image, Nf, f, N_center] = buildTestImg(resolution, gridSize, L)

[X, Y] = meshgrid(linspace(-gridSize, gridSize, resolution), linspace(-gridSize, gridSize, resolution));
sz = size(X);

%Various test functions below:

%f = 3*(1-X).^2.*exp(-(X.^2) - (Y+1).^2) - 10*(X/5 - X.^3 - Y.^5).*exp(-X.^2-Y.^2) - 1/3*exp(-(X+1).^2 - Y.^2); %peaks
%R = sqrt((2*(X)).^2 + (2*(Y).^2)) + eps; f = sin(R)./R;
%f = sin(3*X) + cos(2.2*Y); %the good one to use for testing 4/27 %done a
f = .2*X + .4*Y + sin(3*X) + 2*cos(2*Y); %CURRENT SURFACE FOR POSTER
%gif, r3_norm_change
%f = X + 1.2*Y + sin(3*X) + cos(2.2*Y); 
%f = -(sin(1.7*X) + cos(2.2*Y));
%f = X.^4 + Y.^2;
%f = X.^4 + (Y-1.1).^3;
%f = sqrt(1 - X.^2);
%f = sqrt(1 - (X-0.3).^2); %good cylinder one
%f = sqrt(1 - (X).^2); %good cylinder one
%f = Y.^2 + X.^5 + (X.^3).*(10*(Y - 0.2).^5); %normal_field_to_height_patch tests don't work
%f = exp(-((X-1).^2 + 1.5*Y.^2)/2); %excellent one %done a gif - also, it's additional fig_1
%r3_norm_cyl
%f = exp(-((X).^2 + 1.5*Y.^2)/2); %good one
%f = exp(-X.^2 - Y.^2);
%f = cos(X/2); %good one
%f = cos(X); %good one
%f = cos(2*X); %bad one
%f = sin(2*X) + Y.^2; %bad one
%f = cos(4*X) + cos(1.5*Y);
%f = X.^3 + Y.^2 + 1./((X - 2) + (10^(-5)).^2); %good one, interesting
%saddle
%f = sqrt(1 - X.^2 - Y.^2);
%I think singularities in the image (critical points) cause a lot of
%problems.  (Maybe multiple singularities of opposite sign?)  Why?


if ~exist('L', 'var')
    L = [0, 0, 1]; %a generic lightsource
end

h_spacing = 2*gridSize/sz(1); v_spacing = 2*gridSize/sz(2);
[fx, fy] = gradient(f, h_spacing, v_spacing); 
Nf = zeros(sz(1), sz(2), 3);
Nf(:, :, 1) = -fx./sqrt(1 + fx.^2 + fy.^2);
Nf(:, :, 2) = -fy./sqrt(1 + fx.^2 + fy.^2);
Nf(:, :, 3) = ones(sz(1), sz(2))./sqrt(1 + fx.^2 + fy.^2);

tmp = Nf(:, :, 1)*L(1) + Nf(:, :, 2)*L(2) + Nf(:, :, 3)*L(3);
image = tmp; %CAN'T USE IMAGE VALUES IN CALCULATION OF LT if you've scaled the image.
%image = (tmp - min(tmp(:) + epsilon))./(max(tmp(:)) - (min(tmp(:)) + epsilon));

sz = size(image);
center = [round(sz(1)/2), round(sz(2)/2)];
N_center = squeeze(Nf(center(1), center(2), :));

end
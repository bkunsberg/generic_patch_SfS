
%4/15/2015
%Neccessary input:
%1.
%a.  an image
%b.  a structure of necessary parameters: n - the Taylor order, the nx2
%generic params matrix,
%c.  possible additional params for testing: L - light source, Nf - true
%normal field, N0 - central normal

%Output:
%An option to save to a file: the DjN, the normal field, and the
%reconstructed image (given true parameters)

function [reconstructed_image, adj_DjN, recon_Nf, height_patch, DjI_center] = main_DjN_from_image(fig, image, taylor_order, generic_params, true_params, method_type)
%% Input Control

if ~exist('generic_params', 'var')
    generic_params = zeros(taylor_order, 2);
end

if ~exist('true_params', 'var')
    % true_params contains both the true light source and true central
    % normal (and can contain the true heights also)
    tmp = [0.3; 0.2; 1];
    true_params = struct('Nf', [0; 0; 1], 'L', tmp/norm(tmp));
end

if ~exist('fig', 'var')
    fig = 0;
end

assert(length(size(image)) == 2, 'Image is the wrong size');
assert(taylor_order == round(taylor_order), 'Taylor Order is not an integer');
assert(length(generic_params) == taylor_order, 'Wrong size of generic_params array');

%% Initializations
n = taylor_order;
gridSize = 1/2;
[~, DjI_center] = calc_DjI_from_image_2(image, n, gridSize); % Note that DjI center's first index is intensity.
DjN = cell(n, 1); % The jth cell is a 3 x 2^j matrix.
epsilon = (10^(-15));
sz = size(image);
resolution = sz(1);
center = [round(sz(1)/2), round(sz(2)/2)];
image_center = image(center(1), center(2));

% Do not run for the D1N case, as the normalization terms don't help.
D1N = zeros(3, 2);
D1N(2, :) = DjI_center{1}/(sqrt(1 - image_center^2)+ epsilon);
if 0
    disp('NOTE: Using boundary for D1N!');
    D1N(3, :) = D1N(2, :);
    D1N(3, :) = [0, 0];
else
    D1N(3, :) = generic_params(1, :);
end
DjN{1} = D1N;

%% Actual Algorithm for DjN, j > 2
DjN = DjN_from_DjI(image_center, DjI_center, DjN, generic_params, method_type); %this takes in the cell arrays (only filled in up to i) and outputs the i+1 entry



%% Check reconstructed image
% This code portion is correct and can be tested by latter portion of
% "testing_reconstruction_from_DjI.m"
N = true_params.Nf;
L = true_params.L;

%CONVINCE YOURSELF THAT YOU NEED TO MULTIPLY BY P and NOT PT
lt = (transpose(L) - image_center*transpose(N))/sqrt(1 - image_center^2); %not norm 1!  Fix!
lt = lt/norm(lt);
P = horzcat(N, transpose(lt), cross(N, transpose(lt)));
adj_DjN = pre_multiply_P(P, DjN);

%% Cheat Enabler
%For testing purposes: inherit some true levels from true_DjN
if isfield(true_params, 'true_DjN') && isfield(true_params, 'true_levels')
    true_DjN = true_params.true_DjN;
    switch_levels = true_params.true_levels;
    disp(['Cheat Enabled: Inheriting true_DjN Levels ...', num2str(switch_levels)]);
    for i = 1:length(switch_levels)
        adj_DjN{switch_levels(i)} = true_DjN{switch_levels(i)};
    end
end

%% Rendering the DjN to recover normals
%scaled_true_image = (image - min(image(:) + epsilon))./(max(image(:)) - (min(image(:)) + epsilon));
rendered_image = render_DjI(DjI_center, image, resolution, gridSize); %from image derivatives
[unscaled_image, recon_Nf] = render_DjN(L, N, adj_DjN, resolution, gridSize, 0);

%NOTE:  recon_Nf is not a unit length fields, in fact it blows up in all
%components!
field_norm = sqrt(recon_Nf(:, :, 1).^2 + recon_Nf(:, :, 2).^2 + recon_Nf(:, :, 3).^2);
recon_Nf = recon_Nf./cat(3, field_norm, field_norm, field_norm); %renormalize

reconstructed_image = (unscaled_image - min(unscaled_image(:)))/(max(unscaled_image(:)) - min(unscaled_image(:)));

%% Recover heights from normals
f11 = 0;
fx = - recon_Nf(:, :, 1)./recon_Nf(:, :, 3);
fy = - recon_Nf(:, :, 2)./recon_Nf(:, :, 3);
spacing = 2*gridSize/resolution;
height_patch = intgrad2(fx,fy,spacing,spacing,f11);


%height_patch = normal_field_to_height_patch(recon_Nf, gridSize, 0, 0);
%%old way



%% Figures

if fig
    %figure; imshow(field_norm, [0, 2]); title(['Reconstructed normals norm, Taylor Order: ', num2str(taylor_order)]); xlabel('X axis'); ylabel('Y axis');
    figure('units','normalized','outerposition',[0 0 0.75 0.5]);
    subplot(2, 3, 1); imshow(image); title('Original Unscaled Image'); xlabel('X axis'); ylabel('Y axis');
    subplot(2, 3, 2); imshow(rendered_image); title(['Re-rendered Image using Image Derivatives, Taylor Order: ', num2str(taylor_order)]); xlabel('X axis'); ylabel('Y axis');
    subplot(2, 3, 3); imshow(unscaled_image); title(['Reconstructed Unscaled Image, Taylor Order: ', num2str(taylor_order)]); xlabel('X axis'); ylabel('Y axis');
    subplot(2, 3, 4); imshow(abs(image - unscaled_image), [0, 0.2]); colorbar; title(['Difference between Reconstruction and Original, Taylor Order: ', num2str(taylor_order)]); xlabel('X axis'); ylabel('Y axis');
    subplot(2, 3, 5); contour(image, 20, 'red'); title('Original Isophotes'); xlabel('X axis'); ylabel('Y axis');
    subplot(2, 3, 6); contour(unscaled_image, 20, 'red'); title('Reconstructed Isophotes'); xlabel('X axis'); ylabel('Y axis');
    
    figure;
    l = surf(height_patch);
    set(l, 'edgecolor','none');
    title('Reconstructed Heights'); xlabel('X axis'); ylabel('Y axis');
    
    %If given a true heights, plot them also
    if isfield(true_params, 'f')
        figure;
        h = surf(true_params.f);
        set(h, 'edgecolor','none');
        title('True Heights'); xlabel('X axis'); ylabel('Y axis');
    end
    
    if 0 %for now, plot the curve where e3.N = 0
        e3 = cross(N, transpose(lt))
        imshow(abs(e3(1)*recon_Nf(:, :, 1) + e3(2)*recon_Nf(:, :, 2) + e3(3)*recon_Nf(:, :, 3)) < 10^(-2), [])
        
        %plot reconstructed normals to check
        figure; imshow(field_norm, [0, 2]); title(['Reconstructed normals norm, Taylor Order: ', num2str(taylor_order)]); xlabel('X axis'); ylabel('Y axis');
    end
end
end



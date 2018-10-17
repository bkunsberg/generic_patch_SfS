%This script runs an example.  
%The goal is to take a polynomial shaded patch (up to some 'taylor order')
%and output a height fields that can create that image assuming Lambertian shading.

%


clc; close all;
addpath(genpath(pwd));

%% Construct a random polynomial image as the input
resolution = 500;
gridSize = .4;
L = [-0.2, 0.1, 1];
[image, Nf, f, N_center] = buildTestImg(resolution, gridSize, L);
figure; imshow(image); title('Random Testing Image');

%% Find a possible generic surface with that image and a chosen 2xn vector of 'generic params'

%This is the approximation order desired for the image patch.
taylor_order = 7;

%Choose some generic parameters for visualization.  For any choice of the
%below parameters, you can find a generic surface with those parameters.
generic_params = [-1*ones(taylor_order, 1), -0.2*rand(taylor_order, 1)];

%The true params are recorded (but not used) for comparison purposes.
true_params = struct('Nf', N_center, 'L', transpose(L), 'f', f);
[reconstructed_image, adj_DjN, recon_Nf, height_patch, DjI_center] = ...
    main_DjN_from_image(1, image, taylor_order, generic_params, true_params, 'generic');



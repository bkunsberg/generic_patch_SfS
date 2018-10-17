%This function calculates, from a 2D image, the DjI for j <= n.
%Basically problems when image derivatives get really large > 10^(5);


function [Iderivs, DjI_center] = calc_DjI_from_image_2(image, n, gridSize)
%%
%Create blur filters
sz = size(image);
%sigma = round(sz(1)/100);
sigma = 3;
window_size = round(3*sigma);
gauss_filter = fspecial('gaussian', window_size, sigma);
center = [round(sz(1)/2), round(sz(2)/2)];
h_spacing = 2*gridSize/sz(1); v_spacing = 2*gridSize/sz(2);  


%%
%Construct the different differential levels of intensity

%First build the structure that houses the I_{xxyy}, etc.
%Each cell is a row vector, where the index corresponds to 1 + (num of y
%derivs)

Iderivs = cell(n, 1); %Each cell contains a n x n x m array, where m is the number of derivatives.
blur_img = filter2(gauss_filter, image);
[I_x, I_y] = gradient(blur_img, h_spacing, v_spacing);
tmp = zeros(sz(1), sz(2), 2); tmp(:, :, 1) = I_x; tmp(:, :, 2) = I_y;
Iderivs{1} = tmp;


for i = 1:n-1
    len = size(Iderivs{i}, 3);
    new_Iderivs = zeros(sz(1), sz(2), len + 1);
    for j = 1:2:len
        
        curr = filter2(gauss_filter, Iderivs{i}(:, :, j)); %this line and the next take up 3/4 of running time
        [curr_x, curr_y] = gradient(curr, h_spacing, v_spacing);
        new_Iderivs(:, :, j) = curr_x;
        new_Iderivs(:, :, j+1) = curr_y;
        
    end
    if ~mod(len, 2) %even
        %need to add in last one I_{yy.. y}
        curr = filter2(gauss_filter, Iderivs{i}(:, :, len));
        [~, curr_y] = gradient(curr, h_spacing, v_spacing);
        new_Iderivs(:, :, len + 1) = curr_y;
    end
    Iderivs{i+1} = new_Iderivs;
    
end



%And then map it into DjI via (# of 1's) of (binary(n-1))
DjI_center = cell(n, 1);

%Actually, can just map it directly into a DjI_center, no need to store the
%DjIs.

%This part has been tested and seems to work. (4/20/15)

for i = 1:n
    %Build each DjI{i} cell as  vector
    tmp = zeros(1, 2^i);
    for k = 1:length(tmp)
        bin_rep = dec2bin(k - 1);
        digits_bin = sscanf(strrep(num2str(bin_rep,8),'.',''),'%1d');
        num_y_derivs = sum(digits_bin);
        Iderivs_index = num_y_derivs + 1;
        tmp(k) = Iderivs{i}(center(1), center(2), Iderivs_index);
    end
    DjI_center{i} = tmp;
end



end

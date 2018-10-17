%This function takes as input a collection of DjN from i = 1 to m.
%Note that each DjN is a 3x2^j matrix.  If the DjN is the output of your
%algorithm, remember to premultiply each by P.

%The collection should be a cell array of 3x2^j matrices.
%Note that the first cell will be D0N = N. (a 3x1 vector).
%Also, need a light source (3x1).

%The output is a rendered image.

%NOTE: April 20, 2015.  This function is taking about half of the
%total running time.




function [rendered_image, recon_Nf] = render_DjN(L, trueN, DjN, resolution, gridSize, chop) %resolution is an odd number
if ~exist('chop', 'var')
    chop = 0;
end

if ~exist('resolution', 'var')
    resolution = 500;
end

[X, Y] = meshgrid(linspace(-gridSize, gridSize, resolution), linspace(-gridSize, gridSize, resolution)); 
center_index = round(resolution/2);
p0 = [X(center_index, center_index), Y(center_index, center_index)];
%Use Taylor Expansion around (0, 0)
sz = size(X);
component_level = zeros(sz(1), sz(2), 3);

%need to get the order of the for loops right
%need to add in the starting normal

for k = 1:3 %for each component of the normal
    polyn = zeros(sz);
    for i = 1:length(DjN) %for each differential power
        add_term = zeros(size(X));
        diff_tensor = DjN{i};
        unfolded_DjN = tenmat(diff_tensor, 1);
        diff_level = unfolded_DjN(k, :); %grab the proper row
        
        %possibly chop due to noise
        mn = mean(abs(diff_level(:)));
        stdev = std(abs(diff_level(:)));
        threshold = mn + 1*stdev;
        
        for j = 1:i+1; %for each unique element in the row
            x_exp = (i+1-(j));
            y_exp = (j-1);
 
            coeff = diff_level(2^(j-1)); %I'm actually grabbing the unique (I_{xxy}, etc.) terms and just summing those.  Good!
            if chop
                fixed_coeff = sign(coeff)*min(threshold, abs(coeff));
            else 
                fixed_coeff = coeff;
            end
            partial_term = (1/(factorial(x_exp)*factorial(y_exp)))*fixed_coeff*((X - p0(1)).^(x_exp).*((Y - p0(2)).^(y_exp)));
            add_term = add_term + partial_term;
            
        end
        polyn = polyn + add_term;
    end
    component_level(:, :, k) = polyn;
end
recon_Nf = zeros(size(component_level));

recon_Nf(:, :, 1) = trueN(1) + component_level(:, :, 1); 
recon_Nf(:, :, 2) = trueN(2) + component_level(:, :, 2);
recon_Nf(:, :, 3) = trueN(3) + component_level(:, :, 3);

rendered_image = L(1)*recon_Nf(:, :, 1) + L(2)*recon_Nf(:, :, 2) + L(3)*recon_Nf(:, :, 3);
end


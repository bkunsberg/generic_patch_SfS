%This function takes as input a collection of DjI from i = 1 to m.
%The collection should be a cell array.
%The output is a rendered image.


%Example cell array
% m = 3;
% DjI = cell(m+1, 1); %j starts from 0 (i.e. the intensity)
% DjI{1} = 0.3;
% DjI{2} = [0.2, 0.5];
% DjI{3} = [2, 0, 1];
% DjI{4} = [3, 2, 2, 1];


function rendered_image = render_DjI(DjI, image, resolution, gridSize) %resolution should be an odd number

if ~exist('resolution', 'var')
    resolution = 500;
end

[X, Y] = meshgrid(linspace(-gridSize, gridSize, resolution), linspace(-gridSize, gridSize, resolution)); 
center_index = round(resolution/2);
p0 = [X(center_index, center_index), Y(center_index, center_index)];

%Use Taylor Expansion around (0, 0)
polyn = image(center_index, center_index)*ones(size(X)); %0th order term

for i = 1:length(DjI)
    add_term = zeros(size(X));
    diff_level = DjI{i}; 
    for j = 1:i+1; 
        
        %The main problem is that we have to go through all the terms in
        %DjI{k}, but the sum of the exponents is going to me much lower and
        %we'll have to give the same exponent pair to many elements of
        %DjI{k}.  E.g. {xxyy = xyxy = xyyx = yyxx, etc.}.  Fixed this by
        %using 2^(j-1) indexing.
        
        
        x_exp = (i+1-(j));
        y_exp = (j-1);
        coeff = diff_level(2^(j-1));
        partial_term = (1/(factorial(x_exp)*factorial(y_exp)))*coeff*((X - p0(1)).^(x_exp).*((Y - p0(2)).^(y_exp)));
        add_term = add_term + partial_term;
    end
    polyn = polyn + add_term;
end

rendered_image = polyn;

if 0
    figure; imshow(rendered_image, []);
end
end
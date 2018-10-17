%This function builds r2 from a mix of the image, DjI and r1.

%Formula is: r2 = (DjI + I r1 )/Sqrt( 1 - I^2)

function r2 = tensor_unfold_DjI(image, DjI, r1, k)

curr_DjI = DjI{k};
row_DjI = reshape(curr_DjI, [1, 2^k]);
epsilon = (10^(-15));

r2 = (row_DjI - image*r1)/(epsilon + sqrt(1 - image.^2)); 

end
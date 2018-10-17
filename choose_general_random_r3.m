%This function is for generating a random r3, to build non-generic
%surfaces.  It should be called from DjN_from_DjI.
%May 14, 2015
%
%Input: r1 (just to get appropriate length), norm to give correct normed
%values
function r3 = choose_general_random_r3(r1, norm)

len = length(r1);
num_choices = log2(len) + 1;
rand_choices = norm*rand(num_choices, 1); %Generate the correct amount of image derivatives


%Now, to fill a tensor row with the image derivatives
%Build each r3 as a column vector
tmp = zeros(1, len);
for k = 1:length(tmp)
    bin_rep = dec2bin(k - 1);
    digits_bin = sscanf(strrep(num2str(bin_rep,8),'.',''),'%1d');
    num_y_derivs = sum(digits_bin);
    Iderivs_index = num_y_derivs + 1;
    tmp(k) = rand_choices(Iderivs_index);
end

r3 = tmp;


end
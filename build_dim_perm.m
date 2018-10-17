%This function takes in the effect of Mathematica's Transpose[array/list,
%transpose_vector] and makes the equivalent vector dim_perm so that
%Matlab's permute(array/list, dim_perm) has the same function.

%4/22/2015

function dim_perm = build_dim_perm(transpose_vector)
n = length(transpose_vector);
renamed_vector = zeros(n, 1);

for i = 1:n
    j = transpose_vector(i);
    renamed_vector(i) = n+1 - j;
end

[~, I] = sort(renamed_vector);
reverse = fliplr(1:n);
dim_perm = reverse(I);




end
% A possible check for this function: 6/3/2015
% Calculate DjN for j = 0:1:m for a given patch.
% Run this code for a single DjN and ensure that the output matches (up to numerical accuracy) a nT
% D_{j+1} N


%This function calculates the unfolded <DjN, n> from the sum of terms
%<DlN, DkN> where  l + k = j.  This constraint arises from differentiating
%the <n, n> = 1 constraint.

%We need to sum over "all" combinations at each differential level.  This
%involves a tensor product, then a contraction, then a permutation of
%dimensions.  See bottom of "normalization_investigation_1.nb".


function r1 = create_r1_from_normalization_constraints(DjN, k)

%Looking for the kth level, but note the index starts at 0.  So, k = 3
%means D2N.
sum = 0;
dirs = 1:k;

for i = 1:floor(k/2) %for each differential level
    %disp('Starting new differential level');
    %disp(i);
    possible_combs = nchoosek(dirs, i);
    if i == (k/2) %prevent double counting the middle terms
        combs = possible_combs(1:(length(possible_combs)/2), :);
    else
        combs = possible_combs;
    end
    
    base_term = calculate_baseterm(DjN, i, k);%going to be <DlN, DkN>
    diff_level_term = 0;
    for j = 1:length(combs) %for each dimension permutation
        transpose_vector = horzcat(combs(j, :), setdiff(1:k, combs(j, :)));
        dim_perm = build_dim_perm(transpose_vector);
        
        if isequal(dim_perm, dirs)
            new_perm = base_term;
        else
            %disp('Entering the permute thingy')
            new_perm = permute(base_term, dim_perm);
        end
        diff_level_term = diff_level_term + new_perm;
    end
    sum = sum + diff_level_term;
end

r1 = -reshape(sum, [1, 2^k]);

end
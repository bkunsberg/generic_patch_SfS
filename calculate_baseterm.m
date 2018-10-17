%This function calculates a base term <DiN, DlN> at the differential level
%i, where i + l = k.

%Assume each Djn is of the form 3x2x2 ... x2. (2^j).
%Then, you need to contract the first index and the j+2th index after the
%tensor product.

function base_term = calculate_baseterm(DjN, i, k)
%note that the starting index for DjN is 0.  That is DjN, j = 2, is actually the 3x2 DN.
%So if i = 2, k = 3, am I looking at the term <DN, N>?
%Best way to deal with it, I think, is to just build the DjN cell array to
%have first index as D1N.


l = k - i;
%disp('Base term is the tensor product of DjN{i} and DjN{l}:');
%disp([i, l]);
tensor_prod = ttt(tensor(DjN{l}), tensor(DjN{i}));
sz_prod = size(tensor_prod);


contract_indices = find(sz_prod == 3);
base_term = contract(tensor_prod, contract_indices(1), contract_indices(2)); %This dimension is scrambled compared to the Mathematica version

end
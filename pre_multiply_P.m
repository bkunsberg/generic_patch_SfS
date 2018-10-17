%This function create P^T DjN given {P, DjN}.  In shading terms, this
%places the solved for DjN from the normal, light source frame back into the {x, y, z} frame.

function adj_DjN = pre_multiply_P(P, DjN)

len = length(DjN);
adj_DjN = cell(len, 1);
%PT = transpose(P);
%I don't want the transpose, because I want to undo the effect of the
%transpose. The output is already in that basis, and I want to put it back
%into the {x, y, z} basis.

for i = 1:len
    %Will need to 1-mode unfold each DjN, left multiply by P^T and then
    %refold
    
    %Unfold
    curr_DjN = tensor(DjN{i});
    %unfolded_DjN = tenmat(curr_DjN, 1);
    
    unfolded_DjN = reshape(double(DjN{i}), [3, 2^i]);
    
    if sum(imag(unfolded_DjN(:))) ~= 0
       error('Unfolded DjN has imaginary parts!');
    end
    
    
    %Multiply
    P_DjN = P*unfolded_DjN;
    
    %Refold
    sz_DjN = size(curr_DjN);
    
    
    refolded_PDjN = reshape(P_DjN, sz_DjN); 
    adj_DjN{i} = tensor(refolded_PDjN);
    
end

end
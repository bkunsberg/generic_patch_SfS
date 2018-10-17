%This function needs as input:
%1. DjN, for all j < k, (cell array).
%2. It also needs DjI for j = k.
%3. It needs two parameters, which represent the linear combination of the two rows (r1 and r2)
%   which create r3.  Params will be a m x 2 matrix.

%The output is a modified DjN cell array, with the kth cell now defined.


function DjN = DjN_from_DjI(image, DjI, DjN, params, method_type)

if ~exist('method_type', 'var')
    method_type = 'generic'; %swich this back to 0 after done testing boundary
end

len = length(DjI);

for k = 2:len %Do not run the D1N case, but preprocess it instead
    r1 = double(create_r1_from_normalization_constraints(DjN, k)); %this is the hard part -- it contains the tensor unfoldings
    r2 = double(tensor_unfold_DjI(image, DjI, r1, k));
    
    if strcmp(method_type, 'boundary')
        
        tangent_dir_index = 1;
        scale = 0.1*rand(1); %let w = (1, 0), find 1-D solutions (param. by scale) such that e3 DjN(w^j) = 0
        if r2(tangent_dir_index) ~= 0
            ratio = r1(tangent_dir_index)/r2(tangent_dir_index);
            c1 = scale; c2 = - ratio*c1;
        else
            c1 = 0.1; c2 = 0;
        end
        r3 = c1 * r1 + c2 * r2;
        disp('Method type is boundary + generic');
        disp('Choosing r3 from boundary values: c1, c2');
        disp(c1); disp(c2);
        
        
    elseif strcmp(method_type, 'general')
        disp('Method type is general');
        general_norm = 1*mean(abs(r1));
        disp(['General_norm is: ', num2str(general_norm)]);
        r3 = choose_general_random_r3(r1, general_norm);
    elseif strcmp(method_type, 'generic')
        disp(k);
        disp('This differential level is chosen in generic fashion');
        r3 = double(params(k, 1) * r1 + params(k, 2) * r2);
    end
    
    DkN_as_matrix = [r1; r2; r3];
    
    if sum(imag(DkN_as_matrix(:))) ~= 0
        error('DkN_as_matrix has imaginary parts!');
    end
    
    %     DkN_as_matrix = zeros(3, length(r1));
    %     DkN_as_matrix(1, :) = r1;
    %     DkN_as_matrix(2, :) = r2;
    %     DkN_as_matrix(3, :) = r3;
    
    if k > 1
        previous_sz = size(DjN{k-1});
        new_sz = horzcat(previous_sz, 2);
    else
        new_sz = [3, 2];
    end
    
    DjN{k} = reshape(tensor(DkN_as_matrix), new_sz);
    %disp(DjN{k})
end


end
% computes Baker-Campbell-Hausdorff for n AlgebraElements
function [Z, Z_terms] = BCH_n_elements(algebra_elements, num_terms)
    % initialize storage of terms
    Z_terms = cell(1, num_terms);
    % enter first element just to bootstrap; this will be removed later
    Z_terms(1) = algebra_elements(1);
    % BCH is computed "nested" for ea. algebra element
    for e = 2:length(algebra_elements)
        terms_so_far = sum(~cellfun(@isempty, Z_terms));
        for order = 1:terms_so_far
            % BCH computed for ea. existing term to preserve knowledge of order
            [~,intermediate_terms] = BCH(Z_terms{order}, algebra_elements{e}, num_terms);
            % add resultant terms (of acceptable order) to existing terms
            for t = 1:length(intermediate_terms)
                if order+t-1 <= num_terms
                    % shift index, such that elements are correct order
                    Z_terms{order+t-1} = AlgebraElement(intermediate_terms{t}.vector);
                end
            end
        end
    end
    % create resultant element from terms
    % a cellfun would compact this, but I can't figure it out rn
    Z_vec = zeros(1,length(Z_terms{1}.vector))';
    for i = 1:length(Z_terms)
        Z_vec = Z_vec + Z_terms{i}.vector;
    end
    Z = AlgebraElement(Z_vec);
end
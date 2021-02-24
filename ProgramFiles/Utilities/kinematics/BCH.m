% computes Baker-Campbell-Hausdorff terms and result for two AlgebraElements
% pulled from https://en.wikipedia.org/wiki/Baker%E2%80%93Campbell%E2%80%93Hausdorff_formula#Explicit_forms
function [Z, Z_terms] = BCH(x, y, num_terms)
    % separate terms into cell array
    % this is so we don't waste compute time on unused terms
    z_series = {@(X,Y) X + Y,...
                @(X,Y) 1/2*(X*Y - Y*X),...
                @(X,Y) 1/12*(X^2*Y + X*Y^2 - 2*X*Y*X + Y^2*X + Y*X^2 - 2*Y*X*Y),...
                @(X,Y) 1/24*(X^2*Y^2 - 2*X*Y*X*Y - Y^2*X^2 + 2*Y*X*Y*X)};
    % compute each term, as well as aggregate
    Z_mat = 0;
    Z_terms = cell(1, num_terms);
    for t = 1:num_terms
        Z_term_func = z_series{t};
        current_term = Z_term_func(x.matrix, y.matrix);
        Z_terms{t} = AlgebraElement(current_term);
        Z_mat = Z_mat + current_term;
    end
    Z = AlgebraElement(Z_mat);
end
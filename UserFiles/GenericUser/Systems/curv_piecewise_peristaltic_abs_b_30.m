function output = curv_piecewise_peristaltic_abs_b_30(params,mode)

% 4/6/17 Jacquelin's attempt at a simple version of a curvature file, with
% eventual goal of adding elongation mode

% Turn params into a cell matrix
% params = num2cell(params);

switch mode
    
    case 'curvature'
        
        % kappas
        % get alphas for use in heaviside function
        a1 = params(1); 
        a2 = params(2);
        output = @(s) 0*(a1*heaviside(-s)+a2*heaviside(s));

    case 'angle'
        
        % thetas
        a1 = params(1); 
        a2 = params(2);
        output = @(s) 0*(a1*heaviside(-s).*s+a2*heaviside(s).*s);
        
    case 'dcurvature'
        
        % dk/da 
        output = @(s) 0*[heaviside(-s) heaviside(s)];
            
    case 'dcurvature_int'
        
        % int(dk/da)
        output = @(s) 0*[heaviside(-s).*s heaviside(s).*s];
        
    case 'stretch'
        
        % lambda term for elongation
        b = 0.3;       
        a1 = params(1);  % alphas
        a2 = params(2);
        e1 = b*abs(a1) + 1; % lambda, elongation factor
        e2 = b*abs(a2) + 1;
        output = @(s) e1*heaviside(-s)+e2*heaviside(s);
        
    case 'dstretch'
        
        % change in lambda
        b = 0.3;
        a1 = params(1);  % alphas
        a2 = params(2);
        e1 = b*sign(a1); % dlambda, elongation factor change
        e2 = b*sign(a2);
        output = @(s) [e1*heaviside(-s) e2*heaviside(s)];
                                   
end

end
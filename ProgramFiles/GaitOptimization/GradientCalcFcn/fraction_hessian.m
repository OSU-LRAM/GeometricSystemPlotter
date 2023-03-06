function gradfun = fraction_hessian(fun1,fun2,gradfun1,gradfun2,hessfun1,hessfun2)
    % fun1, fun2: a scalar value.
    % gradfun1, gradfun2: a nx1 vector.
    % Calculate d^2(fun1/fun2).
    
    gradfun = (hessfun1/fun2)-(2/fun2^2)*(gradfun1)*(gradfun2.')-...
    (fun1/fun2^2)*(hessfun2)-(fun1/fun2^3)*(gradfun2)*(gradfun2.');
end
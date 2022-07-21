function gradfun = fraction_gradient(fun1,fun2,gradfun1,gradfun2)
    % fun1, fun2: a scalar value.
    % gradfun1, gradfun2: a nx1 vector.
    % Calculate d(fun1/fun2).
    
    gradfun = gradfun1/fun2-fun1/fun2^2*gradfun2;
end
function [twoint,prev_integral,n,eval] = ctr(func, a, b, tol)
    n = 1; 
    prev_integral = 0;
    twoint = (b - a) * (func(a) + func(b)) / 2; 
    eval=0;
    while true
        h = (b - a) / n;
        sum_func = 0;
        for i = 1:n-1
            sum_func = sum_func + func(a + i * h);
        end
        integral_value_new = (b - a) * (func(a) + 2 * sum_func + func(b)) / (2 * n);
        if abs(integral_value_new - prev_integral) < tol
            twoint = integral_value_new;
            break;
        else
            prev_integral = integral_value_new;
            n = n * 2;
            eval=eval+1;
            twoint = integral_value_new;
        end
    end
end

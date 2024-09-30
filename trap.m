function [inttrap, trapiter,mtrap] = trap(a, b, fun, tol)
    inttrap = 0;
    trapiter = 0;
    mtrap=1;
    while abs(inttrap - (exp(b) - 1)) >= tol
        h = (b - a) / mtrap;
        x = a:h:b;
        dim = length(x);
        y = feval(fun, x);
        if size(y) == 1
            y = diag(ones(dim)) * y;
        end
        inttrap = h * (0.5 * y(1) + sum(y(2:end-1)) + 0.5 * y(end));
        trapiter = trapiter + 1;
        mtrap = mtrap +1; 
    end
end
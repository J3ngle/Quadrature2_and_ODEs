function [result, iterations] = gauss_legendre_2point(fun, a, b, tol)
    % Initialize result
    result = 0;
    iterations = 0;
    
    % Two-point Gauss-Legendre weights and nodes
    weights = [1, 1];
    nodes = [-sqrt(1/3), sqrt(1/3)];
    
    % Initialize previous result for error calculation
    prev_result = inf;
    
    % Iterate until the error tolerance is satisfied
    while abs(result - prev_result) >= tol
        % Update previous result
        prev_result = result;
        
        % Compute the quadrature approximation
        result = 0;
        for i = 1:numel(weights)
            x = (b - a) / 2 * nodes(i) + (b + a) / 2;
            result = result + weights(i) * fun(x);
        end
        result = result * (b - a) / 2;
        
        % Double the number of intervals and update result
        weights = repmat(weights, 1, 2);
        nodes = repmat(nodes, 1, 2);
        
        % Increment iteration counter
        iterations = iterations + 1;
    end
end


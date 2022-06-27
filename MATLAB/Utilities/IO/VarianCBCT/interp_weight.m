function [ii, weights] = interp_weight(x, X)
    % LOCATE Locate first node on grid below a given value.
    %
    %   [ii, weights] = locate(x, X) returns the first node in X that is below
    %   each element in x and the relative proximities to the two closest nodes.
    %
    %   X must be a monotonically increasing vector. x is a matrix (of any
    %   order).

    % Preallocate
    ii = ones(size(x));  % Indices of first node below (or 1 if no nodes below)
    weights = zeros([2, size(x)]);  % Relative proximity of the two closest nodes

    % Find indices
    for iX = 1:length(X) - 1
        ii(X(iX) <= x) = iX;
    end

    % Find weights
    below = x <= X(1);
    weights(1, below) = 1;  % All mass on the first node
    weights(2, below) = 0;

    above = x >= X(end);
    weights(1, above) = 0;
    weights(2, above) = 1;  % All mass on the last node

    interior = ~below & ~above;
    xInterior = x(interior)';
    iiInterior = ii(interior);
    XBelow = X(iiInterior)';
    XAbove = X(iiInterior + 1)';
    weights(:, interior) = ...
        [XAbove - xInterior; xInterior - XBelow] ./ (XAbove - XBelow);
end
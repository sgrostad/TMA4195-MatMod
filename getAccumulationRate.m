function q = getAccumulationRate(L, xs, xf, N, q0)
    % L = Length of glacier
    % xs = Starting point snow smelts in summer
    % xf = Toe of glacier
    % N = Number of points in discretization
    q = zeros(N,1);
    h = xf/N;
    % creating function describing the decreasing accumulation
    pointsBUntilXs = xs/h;
    meltingAtToe = -1.5*q0;
    gradient = (meltingAtToe - q0) / (xf - xs);
    a = q0 - gradient * xs;
    f = @(x) a + gradient * x;
    for i = 1:N
        if i <= pointsBUntilXs
            q(i) = q0;
        else
            q(i) = f(i*h);
        end
    end
end


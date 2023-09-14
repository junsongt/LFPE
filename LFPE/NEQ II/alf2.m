%=========================================================================
function result = alf2(E, n, xm1, xm2, xn1, xn2)
% WARNING: 
% s has 0 based index;
% xk has 0 based index;
% xj has 0 based index;
% A has 0 based index;
% alpha has 1 based index;


s = zeros(21,21);
xk = zeros(21,1);
xj = zeros(21,1);

nn = n + 5;
nnm1 = nn - 1;

% s = zeros(nn, nn);
% xk = zeros(nn, 1);

% Calculate the coefficients for the Sonine polynomials s(i,n)
s(1, 1) = 1;
for i = 1:nnm1+1
    xi = double(i-1); % index!!!
    s(i+1, 1) = (xi + 1.5) * s(i, 1) / (xi + 1);
end

for i = 2:nn+1
    xi = double(i-1);
    for k = 1:i
        xkk = double(k-1); % index!!!
        s(i, k+1) = -(xi - xkk) * s(i, k) / ((xkk + 1.5) * (xkk + 1));
    end
end

% Calculate the J(k) integrals
expe = exp(-E);
% xj = zeros(nn+1, 1);
xj(1) = expe / 2;
for i = 2:nn+2
    xi = double(i-1); % index!!!
    xj(i) = E^(i-1) * expe / 2 + xi * xj(i-1); % index!!!
    % xj(i) = E^(i) * expe / 2 + i * xj(i-1);
end

% Calculate the K(k) integrals
xk(1) = expe / 2;
for k = 2:nn+1
    xk(k) = xj(k+1) - E * xj(k);
    % xk(k) = xj(k) - E * xj(k-1);
end

% Calculate the A(j) integrals
% a = zeros(nn+1, 1);
for i = 1:nn+1
    sum = 0;
    for k = 1:i+1
        sum = sum + s(i, k) * xk(k);
    end
    A(i) = 4 * sum;
    
end



% Calculate the Alpha(i) quantities
alpha(1, 1) = (xm2 - xn1/(xn1 + xn2)) * A(2); % index!!!
alpha(2, 1) = (xm1 - xn2/(xn1 + xn2)) * A(2); % index!!!
% alpha(1, 1) = (xm2 - xm1*xn1/(xn1 + xn2)) * A(2); % index!!!
% alpha(2, 1) = (xm1 - xm2*xn2/(xn1 + xn2)) * A(2); % index!!!
for k = 2:nn
    alpha(1, k) = xm2^k * A(k+1); % index!!!
    alpha(2, k) = xm1^k * A(k+1); % index!!!
end


result.alpha = alpha;
result.A = A;
result.s = s;
end
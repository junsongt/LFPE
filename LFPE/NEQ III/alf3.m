%=========================================================================
function result = alf3(E, n, xm2)
% WARNING: 
% s has 0 based index;
% xk has 0 based index;
% xj has 0 based index;
% a has 0 based index;
% alpha has 1 based index;


s = zeros(21,21);
xk = zeros(21,1);
xj = zeros(21,1);

nn = n + 5;
nnm1 = nn - 1;


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
xj(1) = expe / 2;
for i = 2:nn+2
    xi = double(i-1); % index!!!
    xj(i) = E^(i-1) * expe / 2 + xi * xj(i-1); % index!!!
end

% Calculate the K(k) integrals

% For LOC cross section:
xk(1) = expe / 2;
for k = 2:nn+1
    xk(k) = xj(k+1) - E * xj(k);
end


% % For STEP cross section
% xk(1)=(1/2) *(E + 1) * expe;
% for k=2:nn+1
%     xk(k) = 0.5*E^k*expe + k*xk(k-1);
% end


% %  Calculate K(k)  integral
% for k = 1:nn+1
%     int_k = @(x) exp(-x.^2).*x.^(2*(k-1)+3).*(1 - E./x)./4;
%     xk(k) = 1./pi.*integral( @(x) exp(-x.^2).*x.^(2*(k-1)+3).*(1 - E./x)./4, E, Inf);
% end



% Calculate the A(j) integrals
% a = zeros(nn+1, 1);
for i = 1:nn+1
    sum = 0;
    for k = 1:i+1
        sum = sum + s(i, k) * xk(k);
    end
    % a(i) = 4 * sum;
    A(i) = 4 * (xm2^(i-1)) * sum;
end


% Calculate the Alpha(i) quantities
% alpha(1) = xm2 * a(2); % index!!!
alpha(1) = A(2); % index!!!
for k = 2:nn
    % alpha(k) = xm2^k * a(k+1); % index!!!
    alpha(k) = A(k+1);
end

result.alpha = alpha;
result.A = A;
result.s = s;
result.K = xk;
end
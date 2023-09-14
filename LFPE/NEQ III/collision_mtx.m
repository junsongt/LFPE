function M = collision_mtx(rm, N)
M1 = rm/(1+rm);
M2 = 1 - M1;
syms x y
f(x,y) = (-8) * (M1*M2) * ((x*y)/(1-x*y)^2) * (sqrt(1 - M2*(x+y) + (1-2*M1)*x*y) / (1 - (1-4*M1*M2)*x*y));


M = [];
for i = 1:N
    for j = 1:N
        if i <= j
            M(i,j) = taylorcoeff(f,i,j);
        else
            M(i,j) = M(j,i);
        end
    end
end

end


% Helper functions:

% partial derivative of any order wrt x, returns a symbolic expression
function res = partialx(f, order)
syms x y
F(x,y) = f(x,y);
while order ~= 0
    F = diff(F, x);
    order = order - 1;
end
res = F;
end

% partial derivative of any order wrt y, returns a symbolic expression
function res = partialy(f, order)
syms x y
F(x,y) = f(x,y);
while order ~= 0
    F = diff(F, y);
    order = order - 1;
end
res = F;
end


% calculate the coeff of term of power i,j in Taylor expansion of given function
function coeff = taylorcoeff(f,i,j)
syms x y
p(x,y) = partialy(partialx(f,i), j);
coeff = (1/factorial(i+j)) * nchoosek(i+j, i) * p(0,0);
end
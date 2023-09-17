% Author: Junsong Tang
% This function generates AA collision matrix at arbitrary size, based on
% Eq (36) in Shizgal, 1971 (Paper II)


function M = collision_mtx(rm, N)
% rm is mass ratio: m1/m2; N is matrix size
M1 = rm/(1+rm);
M2 = 1 - M1;
syms x y
% Refer to Eq (36) in Paper II
f(x,y) = (-8) * (M1*M2) * ((x*y)/(1-x*y)^2) * (sqrt(1 - M2*(x+y) + (1-2*M1)*x*y) / (1 - (1-4*M1*M2)*x*y));

% M is symmetric
% M(i,j) := coeff of (x^i y^j) term
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


%==========================================================================
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


% calculate the coeff of the term of power i,j in 
% Taylor expansion of given function f
% return a numerical value
function coeff = taylorcoeff(f,i,j)
syms x y
p(x,y) = partialy(partialx(f,i), j);
% For x^iy^j term, by Clairaut Thm, cross term partial derivative are all
% equal, there are (i+j) choose i such terms in total.
coeff = nchoosek(i+j, i) * (1/factorial(i+j)) * p(0,0);
end
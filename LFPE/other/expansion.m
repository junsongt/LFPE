rm = 1/2;
M1 = rm/(1+rm);
M2 = 1 - M1;
order = 4;
syms x y
f(x,y) = (-8) * (M1*M2) * ((x*y)/(1-x*y)^2) * (sqrt(1 - M2*(x+y) + (1-2*M1)*x*y) / (1 - (1-4*M1*M2)*x*y));
T = taylor(f, [x y], [0 0], "Order", 9);


% f(x,y) = exp(x+y);
% T = taylor(f, [x y], [0 0], "Order", 8);



% function res = partialx(f, order)
% syms x y
% F = f;
% while order ~= 0
%     F = diff(F, x);
%     order = order - 1;
% end
% res = F;
% end
% 
% function res = partialy(f, order)
% syms x y
% F = f;
% while order ~= 0
%     F = diff(F, y);
%     order = order - 1;
% end
% res = F;
% end
% 
% 
% 
% function coeff = taylorcoeff(f,i,j)
% syms x y
% p(x,y) = partialy(partialx(f,i), j);
% coeff = (1/factorial(i+j)) * nchoosek(i+j, i) * p(0,0);
% end







%mat11
%=============================================
function [a] = mat11(xm)
ym = 1 - xm;    
ym2 = ym^2;
ym3 = ym^3;
ym4 = ym^4;
xm2 = xm^2;
xm3 = xm^3;
xm4 = xm^4;
xm5 = xm^5;
xm6 = xm^6;

a(1, 1) = -8 * xm * ym;
a(1, 2) = 4 * xm * ym2;
a(1, 3) = xm * ym3;
% a(1, 4) = xm * ym * (xm3 - 3*xm2 + 3*xm - 1)/ 2; % from paper: a(1, 4) = xm * ym * (xm3 - 3*xm2 + 3*xm - 1)/ 2 = -xm * ym4/ 2
a(1, 4) = xm * ym4/ 2; % from code NEQ2
a(2, 2) = -2 * xm * ym * (15 * xm2 - 18 * xm + 13);
a(2, 3) = xm * ym2 * (35 * xm2 - 30 * xm + 23) / 2;
a(2, 4) = xm * ym * (21 * xm4 - 56 * xm3 + 60 * xm2 - 36 * xm + 11) / 4;
a(3, 3) = -xm * ym * (945 * xm4 - 2100 * xm3 + 2190 * xm2 - 1188 * xm + 433) / 8;
% a(3, 4) = xm * ym * (1155 * xm5 - 3255 * xm4 + 4130 * xm3 - 2970 * xm2 + 1299 * xm - 359) / 16; %from paper
a(3, 4) = -xm * ym * (1155 * xm5 - 3255 * xm4 + 4130 * xm3 - 2970 * xm2 + 1299 * xm - 359) / 16; % from code 
a(4, 4) = -xm * ym * (15015 * xm6 - 48510 * xm5 + 72135 * xm4 - 62300 * xm3 + 34485 * xm2 - 12102 * xm + 2957) / 32;

for i = 1:3
    for j = i+1:4
        a(j, i) = a(i, j);
    end
end

end
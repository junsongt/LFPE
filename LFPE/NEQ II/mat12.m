%mat11
%==============================================
function [a] = mat12(xm)
ym = 1 - xm;
ym2 = ym^2;
ym3 = ym^3;
ym4 = ym^4;
xm2 = xm^2;
xm3 = xm^3;
xm4 = xm^4;

a(1, 1) = 8 * xm * ym;
a(1, 2) = -4 * xm * ym2;
a(1, 3) = -xm * ym3;
a(1, 4) = -xm * ym4 / 2;
a(2, 2) = 30 * xm2 * ym2;
a(2, 3) = -35 * xm2 * ym3 / 2;
a(2, 4) = -21 * xm2 * ym4 / 4;
a(3, 3) = 945 * xm3 * ym3 / 8;
a(3, 4) = -1155 * xm3 * ym4 / 16;
a(4, 4) = 15015 * xm4 * ym4 / 32;

for i = 1:3
    for j = i+1:4
        a(j, i) = a(i, j);
    end
end
end

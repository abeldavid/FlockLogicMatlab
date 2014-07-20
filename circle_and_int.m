% Computes the intersection points of two circles with origins at o1, o2, and radii r.
% INPUT: origins(vectors 1x2) of circles; r = radius of circles (assume same radius)
% OUTPUT: intersection points
% 18/10/2013, Doris Voina.

function [int1, int2] = circle_and_int(o1, o2, r)

w1 = -o1(1) + o2(1);
w2 = -o1(2) + o2(2);

a = 1 + w1^2/w2^2;
b = -w1*(w1^2+w2^2)/w2^2;
c = ((w1^2+w2^2)/(2*w2))^2 - r^2;

if w2~=0
syms x;
int = solve(a*x^2 + b*x + c == 0);
int_x1 = int(1); int_x2 = int(2);
int_y1 = 1/w2*((w1^2+w2^2)/2 - int_x1*w1); int_y2 = 1/w2*((w1^2+w2^2)/2 - int_x2*w1);

int_x1 = int_x1 + o1(1);
int_x2 = int_x2 + o1(1);

int_y1 = int_y1 + o1(2);
int_y2 = int_y2 + o1(2);

int1 = [int_x1, int_y1];
int2 = [int_x2, int_y2];
else
if w1==0    
    int1=NaN;
    int2=NaN;
else
    int_x1 = o1(1)+w1/2;
    int_x2 = o1(1)+w1/2;
    
    int_y1 = sqrt(r^2-(w1/2)^2)+o1(2);
    int_y2 = -sqrt(r^2-(w1/2)^2)+o1(2);
    
    int1 = [int_x1, int_y1];
    int2 = [int_x2, int_y2];
end
end

end
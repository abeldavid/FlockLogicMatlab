function [x y theta] = transf_in_familiar_coord(x,y,theta)

xx(1) = x(3);
xx(2) = x(4);
xx(3) = x(2);
xx(4) = x(12);
xx(5) = x(1);
xx(6) = x(6);
xx(7) = x(8);
xx(8) = x(7);
xx(9) = x(10);
xx(10) = x(9);
xx(11) = x(11);
xx(12) = x(13);
xx(13) = x(5);

yy(1) = y(3);
yy(2) = y(4);
yy(3) = y(2);
yy(4) = y(12);
yy(5) = y(1);
yy(6) = y(6);
yy(7) = y(8);
yy(8) = y(7);
yy(9) = y(10);
yy(10) = y(9);
yy(11) = y(11);
yy(12) = y(13);
yy(13) = y(5);

thetat(1) = theta(3);
thetat(2) = theta(4);
thetat(3) = theta(2);
thetat(4) = theta(12);
thetat(5) = theta(1);
thetat(6) = theta(6);
thetat(7) = theta(8);
thetat(8) = theta(7);
thetat(9) = theta(10);
thetat(10) = theta(9);
thetat(11) = theta(11);
thetat(12) = theta(13);
thetat(13) = theta(5);

x = xx';
y = yy';
theta = thetat';

end


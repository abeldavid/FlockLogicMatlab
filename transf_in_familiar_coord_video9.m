function [x y theta] = transf_in_familiar_coord_video9(x,y,theta)

xx(1) = x(12);
xx(2) = x(13);
xx(3) = x(11);
xx(4) = x(1);
xx(5) = x(3);
xx(6) = x(4);
xx(7) = x(6);
xx(8) = x(9);
xx(9) = x(5);
xx(10) = x(8);
xx(11) = x(10);
xx(12) = x(7);
xx(13) = x(2);

yy(1) = y(12);
yy(2) = y(13);
yy(3) = y(11);
yy(4) = y(1);
yy(5) = y(3);
yy(6) = y(4);
yy(7) = y(6);
yy(8) = y(9);
yy(9) = y(5);
yy(10) = y(8);
yy(11) = y(10);
yy(12) = y(7);
yy(13) = y(2);

thetat(1) = theta(12);
thetat(2) = theta(13);
thetat(3) = theta(11);
thetat(4) = theta(1);
thetat(5) = theta(3);
thetat(6) = theta(4);
thetat(7) = theta(6);
thetat(8) = theta(9);
thetat(9) = theta(5);
thetat(10) = theta(8);
thetat(11) = theta(10);
thetat(12) = theta(7);
thetat(13) = theta(2);

x = xx';
y = yy';
theta = thetat';

end
function theta = anglerestrict(phi)
%ANGLERESTRICT function to restrict an angle to a value between -PI and PI
% ANGLERESTRICT(PHI) takes the angle PHI and computes the equivalent angle
% between -PI and PI

% George Young
% August 2011

theta = phi - 2*pi*ceil((phi-pi)/(2*pi));

end
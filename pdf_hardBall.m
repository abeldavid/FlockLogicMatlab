% applies MLE to find an optimum radius around some node. Inside the circle
% of the radius found no node is permitted and the circle marks a `private 
% space`

% Doris Voina, May 2014

function f = pdf_hardBall(x,positionMatrix, density, lb, ub, o, numberOfNodes)

r_0 = x; % x is the radius, denoted by r_0
f=0; %objective function

for i=lb:ub
    for s=[1:o-1, o+1:numberOfNodes]
        pos_dist = pdist([squeeze(positionMatrix(i,o,:))'; squeeze(positionMatrix(i,s,:))']); %distances from node o to every other node
        f = f + log(pdf_R(pos_dist, r_0, density, o)); %adding log(probability) for MLE to compute objective function
    end
end

f=-f;
end

function p = pdf_R(r,r_0,density,o)

%formulas from paper (see plotVisualizer.m for citation)
packFr = 4/3*pi*(r_0)^3*density(o);
A = 8*packFr*(1+4*packFr);
B = 18*packFr^2;
C = 24*(packFr)^3;

p = exp(-A*((r/(2*r_0))^3-1)+B*((r/(2*r_0))^2-1) - C*(r/(2*r_0)-1))*(-3*A*(r/(2*r_0))^2+2*B*(r/(2*r_0))-C/(2*r_0));

end
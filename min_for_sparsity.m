%MLE for causality algorithm inspired by Sznaier's paper
%Doris Voina, April 2014

function f = min_for_sparsity(x, weights_a, weights_u, q_x, q_y, Y, numberOfNodes, F, delay, lambda, delta)

%translate x into coordinates of interest: a and u

 for j=1:numberOfNodes
    a(j,:) = x((j-1)*(delay+1)+1:j*(delay+1));
 end

 Delta_u_x = x(numberOfNodes*(delay+1)+1:numberOfNodes*(delay+1) + F)';
 Delta_u_y = x(numberOfNodes*(delay+1)+F+1:numberOfNodes*(delay+1) + 2*F)';
 f=0;

  for j=1:numberOfNodes
       f = f + weights_a(j)*norm(a(j,:)); %computing objective function
  end
  
  f=f+lambda*norm(diag(weights_u)*[Delta_u_x; Delta_u_y], 1);
  
end


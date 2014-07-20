%computes betweenness centrality given an adjacency matrix (for one frame)
%Doris Voina, May 2014

function betweenness_c = betweenness(adj)

numberOfNodes = length(adj);
shortestLengths = zeros(numberOfNodes);

   for ii = 1:numberOfNodes
      for jj = [1:(ii-1) (ii+1):numberOfNodes]
          if adj(ii,jj) > 0
             shortestLengths(ii,jj) = adj(ii,jj);
          else
             shortestLengths(ii,jj) = Inf;
          end
      end
   end
   
   for kk = 1:numberOfNodes %interesting... computing the length of the shortest path
       for ii = 1:numberOfNodes
          for jj = 1:numberOfNodes
              shortestLengths(ii,jj) = min(shortestLengths(ii,jj), shortestLengths(ii,kk)+shortestLengths(kk,jj));
          end
       end
   end
   
   for kk = 1:numberOfNodes 
       for ii = 1:numberOfNodes
           if shortestLengths(kk,ii)==Inf
               shortestLengths(kk,ii)=0;
           end
       end
   end
   
   %computing how many shortest paths between two nodes

   possibleLengths = zeros(numberOfNodes);
   
   for j=1:numberOfNodes
       for k=1:numberOfNodes
           
         level=1;
              
         if shortestLengths(j,k) > 0
             
           neighbors = find(adj(j,:)==1);
           if shortestLengths(j,k)==1
             possibleLengths(j,k) = 1;  
           else    
             possibleLengths(j,k) = 0;
           end
           
        while level < shortestLengths(j,k)  
            
           neighbors2=[];
           
          for l=1:length(neighbors) 
              
           if shortestLengths(j,k) == level + shortestLengths(neighbors(l),k)
               if level == shortestLengths(j,k)-1
                 possibleLengths(j,k) = possibleLengths(j,k)+1;  
               else
                 neighbors2 = [neighbors2, find(adj(neighbors(l),:)==1)];
               end
           end
                     
          end
          
          neighbors = neighbors2;
          level=level+1;
          
        end
        
         end
       end
   end
 
   betweenness_c = zeros(1,numberOfNodes);
   
 for o=1:numberOfNodes
     
     for j=1:numberOfNodes-1
         for k=j+1:numberOfNodes
           
           if j~=o & k~=o
             if shortestLengths(j,k)==shortestLengths(j,o)+shortestLengths(o,k) & shortestLengths(j,k)~=0 & possibleLengths(j,k)~=0
                 betweenness_c(o) = betweenness_c(o) + possibleLengths(j,o)*possibleLengths(o,k)/possibleLengths(j,k);
             end
           end
         end
     end
 end
 
 betweenness_c
 [sorted_betw, sortedInd] = sort(betweenness_c)
 
end
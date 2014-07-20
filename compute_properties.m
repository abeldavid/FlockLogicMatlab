%computing properties given an adjacency matrix
% see properties computed displayed on the GUI (H2 norm, algebraic connectivity etc.)

%Doris Voina, May 2014

function compute_properties(adjacencyMatrix, numberOfNodes, min_frame, inputComputeChoice, inputWeightChoice, handles, hObject)

min_frame = size(adjacencyMatrix,1);
lb = 1;
ub = min_frame;

nodeStatus = zeros(ub-lb+1,numberOfNodes);
node_status_cum = zeros(1, numberOfNodes);

threshold=0.5;
cnt_ns = zeros(1, numberOfNodes);
numberOfNeighbors = 2;

for i = lb:ub
 
    file=strcat('Dance',num2str(i),'.dat');
    data_type = handles.data_type;

    if data_type==1
        data=handles.dancedata{i};
        set(handles.dataselection_editText,'String',i);
        x = data(:,1);
        y = data(:,2);
        theta = data(:,3);
        %[x y theta] = transf_in_familiar_coord(x,y,theta);        
        dim = 2;
        positionMatrix = horzcat(x,y);
    else
        disp('Cannot run algorithm on data that is not "Dance data"');
        return;
    end
    
    numberOfNodes = length(positionMatrix);
    distanceMatrix = ipdm(positionMatrix);
    
    relativeAngles = NaN(numberOfNodes);
    HeadingDiff = NaN(numberOfNodes);
    
    for s = 1:(numberOfNodes-1);
        for j = (s+1):numberOfNodes
            relativeAngles(s,j) = atan2(positionMatrix(j,2)-positionMatrix(s,2),positionMatrix(j,1)-positionMatrix(s,1));
            relativeAngles(j,s) = anglerestrict(relativeAngles(s,j) + pi);
        end
    end
        
    for s=1:numberOfNodes
        for j=1:numberOfNodes
            if distanceMatrix(s,j) == 0
                distanceMatrix(s,j) = Inf;
            end
            HeadingDiff(s,j) = abs(anglerestrict(theta(s)-theta(j)));
        end
    end 
   
    
    [sensor, sensed] = find(squeeze(adjacencyMatrix(i,:,:))==1);
    sortedNodePairs = sortrows([sensor,sensed],1);
%up-date weights    
    if inputWeightChoice == 1
        if i==lb
            weight = str2num(get(handles.fixedWeight,'String'));
        else weight = weight; %#ok<*ASGSL>
        end
    end
    
    if inputWeightChoice == 2
        if inputEdgeChoice == 1 | inputEdgeChoice==3
            weight = 1/numberOfNeighbors;
        elseif inputEdgeChoice == 2
            weight = zeros(numberOfNodes,1);
            for j=1:length(sortedNodePairs)
                weight(j) = 1/length(find(sortedNodePairs(:,1)==sortedNodePairs(j,1)));
            end
        end
    end
    
    if inputWeightChoice == 3
        if i==1
            bweight = str2num(cell2mat(inputdlg('Input fixed border edge weight.')));
            iweight = str2num(cell2mat(inputdlg('Input fixed interior edge weight.')));
        else
            bweight = bweight;
            iweight = iweight;
        end
    end
  %up-date adjacency matix according to weights      
    for s=1:numberOfNodes
        for j=1:numberOfNodes
            if adjacencyMatrix(i,s,j)==1
                if inputWeightChoice == 1
                    adjacencyMatrix(i,s,j) = weight;
                end
                if inputWeightChoice == 2
                    if inputEdgeChoice == 1 | inputEdgeChoice ==3
                        adjacencyMatrix(i,s,j) = weight;
                    elseif inputEdgeChoice == 2
                        for k=1:length(sortedNodePairs)
                            adjacencyMatrix(i,sensor(k),sensed(k)) = weight(k);
                        end
                    end
                end
                if inputWeightChoice == 3
                    if s <= divider
                        adjacencyMatrix(i,s,j) = bweight;
                    else adjacencyMatrix(i,s,j) = iweight;
                    end
                end
                if inputWeightChoice == 4
                    set(handles.edges_editText,'String','1/length');
                    adjacencyMatrix(i,s,j) = 1/distanceMatrix(s,j);
                end
                if inputWeightChoice == 5
                    weight(s,j) = 1/HeadingDiff(s,j);
                    adjacencyMatrix(i,s,j) = 1/HeadingDiff(s,j);
                end
            end
        end
    end
    
%     save(strcat('AdjacencyMatrix_t=',num2str(i)),'adjacencyMatrix');
   
    for t=1:length(sortedNodePairs)
        sensorNode = sensor(t);
        sensedNode = sensed(t);
        correspondingDistances(t,:) = distanceMatrix(sensorNode,sensedNode);
        correspondingAngles(t,:) = HeadingDiff(sensorNode,sensedNode);
    end   
    
    if inputComputeChoice ~= 3
        if dim==3
            axis([xmin xmax ymin ymax zmin zmax]);
            grid on;
            hold on;
            plot3(x,y,z,'g.');
        else
            axis([xmin xmax ymin ymax]);
            grid on;
            hold on;
            plot(x,y,'g.');
        end
        
        if inputGraphChoice == 1
            q=1;
            for j=1:length(sortedNodePairs)
                xVal1 = x(sortedNodePairs(j,1));
                xVal2 = x(sortedNodePairs(j,2));
                u(j) = xVal2 - xVal1;
           
                yVal1 = y(sortedNodePairs(j,1));
                yVal2 = y(sortedNodePairs(j,2));
                v(j) = yVal2 - yVal1;
                
                if dim==3
                    zVal1 = z(sortedNodePairs(j,1));
                    zVal2 = z(sortedNodePairs(j,2));
                    w(j) = zVal2 - zVal1;
                end
                
                if inputWeightChoice == 3
                    if j <= divider*numberOfNeighbors
                        ub(j) = u(j);
                        vb(j) = v(j);
                        if dim==3
                            wb(j) = w(j);
                        end
                        borderedges(j,:) = sortedNodePairs(j,:);
                    else
                        ui(q) = u(j);
                        vi(q) = v(j);
                        if dim==3
                            wi(q) = w(j);
                        end
                        q=q+1;
                        interioredges(j,:) = sortedNodePairs(j,:);
                    end
                end
            end
            
            xplot = x(sortedNodePairs(:,1));
            yplot = y(sortedNodePairs(:,1));
            u = transpose(u);
            v = transpose(v);
            if dim==3
                zplot = z(sortedNodePairs(:,1));
                w = transpose(w);
            end
            
            if inputWeightChoice == 1
                if dim ==3
                    quiver3(xplot,yplot,zplot,u,v,w,0,'AutoScale','off','Color','black','LineWidth',weight);
                else
                    quiver(xplot,yplot,u,v,0,'AutoScale','off','Color','black','LineWidth',weight);
                    %viscircles([x(11), y(11)], inputRadius);
                    %viscircles([x(12), y(12)], inputRadius);
                end
            end
            
            if inputWeightChoice == 2
                if dim ==3
                    quiver3(xplot,yplot,zplot,u,v,w,0,'AutoScale','off','Color','black','LineWidth',weight);
                else quiver(xplot,yplot,u,v,0,'AutoScale','off','Color','black','LineWidth',weight);
                end
            end
            
            if inputWeightChoice == 3
                [rw,cl]=find(borderedges==0);
                borderedges(rw,:)=[];
                
                [rw,cl]=find(interioredges==0);
                interioredges(rw,:)=[];
                
                xbplot = x(borderedges(:,1));
                xiplot = x(interioredges(:,1));
                ybplot = y(borderedges(:,1));
                yiplot = y(interioredges(:,1));
                ub = transpose(ub);
                vb = transpose(vb);
                ui = transpose(ui);
                vi = transpose(vi);
                if dim==3
                   zbplot = z(borderedges(:,1));
                    ziplot = z(interioredges(:,1));
                    wb = transpose(wb);
                    wi = transpose(wi);
                    quiver3(xbplot,ybplot,zbplot,ub,vb,wb,0,'AutoScale','off','LineWidth',bweight,'Color','blue');
                    quiver3(xiplot,yiplot,ziplot,ui,vi,wi,0,'AutoScale','off','LineWidth',iweight,'Color','magenta');
                else
                    quiver(xbplot,ybplot,ub,vb,0,'AutoScale','off','LineWidth',bweight,'Color','blue');
                    quiver(xiplot,yiplot,ui,vi,0,'AutoScale','off','LineWidth',iweight,'Color','magenta');
                end
            end
            
            if inputWeightChoice == 4
                for r=1:length(correspondingDistances)
                    if correspondingDistances(r) == Inf
                        correspondingDistances(r) = 0;
                    end
                end
                maxDistance = max(correspondingDistances(:));
                for k=1:length(correspondingDistances)
                    green(k) = correspondingDistances(k)/maxDistance;
                    if dim==3
                        obj = quiver3(xplot(k),yplot(k),zplot(k),u(k),v(k),w(k),0,'Color',[0,green(k),0],'AutoScale','off');
                        set(obj,'AutoScale','off');
                    else
                        obj = quiver(xplot(k),yplot(k),u(k),v(k),0,'Color',[0,green(k),0],'AutoScale','off');
                        set(obj,'AutoScale','off');
                    end
                end
            end
            
            if inputWeightChoice == 5
                for p = 1:length(correspondingAngles)
                lines(p) = 1/correspondingAngles(p);
                lines2(p) = lines(p)^2;
                end
                sqrt_sum_lines = sqrt(sum(lines2(:)));
            for p = 1:length(correspondingAngles)
                norm_lines(p) = lines(p)/sqrt_sum_lines;
            end
                for k=1:length(correspondingAngles)
                    if dim==3
                        obj = quiver3(xplot(k),yplot(k),zplot(k),u(k),v(k),w(k),0,'ShowArrowHead','off','LineWidth',15*norm_lines(k),'AutoScale','off');
                        set(obj,'AutoScale','off');
                    else
                        obj = quiver(xplot(k),yplot(k),u(k),v(k),0,'ShowArrowHead','off','LineWidth',15*norm_lines(k),'AutoScale','off');
                        
                        set(obj,'AutoScale','off');
                    end
                end
            end
            if data_type == 1
                title(strcat(avi_file_name,' at t = ',num2str(i/20),' sec'),'Interpreter','none');
            else
                title(strcat(avi_file_name,' at frame = ',num2str(i)),'Interpreter','none');
            end
            axis equal
%             saveas(gcf,strcat('NearestNeighborGraph_t=',num2str(i)),'fig');
            hold off;
        end
        if inputGraphChoice == 2
            reverseNodePairs = horzcat(sensed, sensor);
            check = ismember(sortedNodePairs, reverseNodePairs, 'rows');
            
            undirectedEdges = [];
            directedEdges = [];
            for k=1:length(reverseNodePairs)
                if check(k)==1
                    undirectedEdges(k,:) = sortedNodePairs(k,:);
                end
                if check(k)==0
                    directedEdges(k,:) = sortedNodePairs(k,:);
                end
            end
            
            [r,c]=find(undirectedEdges==0);
            undirectedEdges(r,:)=[];
            
            [r,c]=find(directedEdges==0);
            directedEdges(r,:)=[];
            
            if ~isempty(undirectedEdges)
                q=1;
                for j=1:size(undirectedEdges,1)
                    xVal1 = x(undirectedEdges(j,1));
                    xVal2 = x(undirectedEdges(j,2));
                    u(j) = xVal2 - xVal1;

                    yVal1 = y(undirectedEdges(j,1));
                    yVal2 = y(undirectedEdges(j,2));
                    v(j) = yVal2 - yVal1;

                    if dim==3
                        zVal1 = z(undirectedEdges(j,1));
                        zVal2 = z(undirectedEdges(j,2));
                        w(j) = zVal2 - zVal1;
                    end

                    if inputWeightChoice == 3
                        if j<=divider*numberOfNeighbors
                            ub(j) = u(j);
                            vb(j) = v(j);
                            if dim==3
                                wb(j) = w(j);
                            end
                            borderedges(j,:) = undirectedEdges(j,:);
                        else
                            ui(q) = u(j);
                            vi(q) = v(j);
                            if dim==3
                                wi(q) = w(j);
                            end
                            q=q+1;
                            interioredges(j,:) = undirectedEdges(j,:);
                        end
                    end
                end

                xplotu = x(undirectedEdges(:,1));
                yplotu = y(undirectedEdges(:,1));
                u = transpose(u);
                v = transpose(v);
                if dim==3
                    zplotu = z(undirectedEdges(:,1));
                    w = transpose(w);
                end

                if inputWeightChoice == 1
                    edgeWeight = weight;
                    if dim==3
                        quiver3(xplotu,yplotu,zplotu,u,v,w,'AutoScale','off','ShowArrowHead','off','Color','red','LineWidth',edgeWeight);
                    else quiver(xplotu,yplotu,u,v,'AutoScale','off','ShowArrowHead','off','Color','red','LineWidth',edgeWeight);
                        %viscircles([x(11), y(11)], inputRadius);
                        %viscircles([x(13), y(13)], inputRadius);
                        %viscircles([x(8),y(8)], inputRadius);
                        %viscircles([x(9), y(9)], inputRadius);
                    end
                end

                if inputWeightChoice == 2
                    edgeWeight = weight;
                    if dim==3
                        quiver3(xplotu,yplotu,zplotu,u,v,w,'AutoScale','off','ShowArrowHead','off','Color','red','LineWidth',edgeWeight);
                    else quiver(xplotu,yplotu,u,v,'AutoScale','off','ShowArrowHead','off','Color','red','LineWidth',edgeWeight);
                    end
                end

                if inputWeightChoice == 3
                    [rw,cl]=find(borderedges==0);
                    borderedges(rw,:)=[];

                    [rw,cl]=find(interioredges==0);
                    interioredges(rw,:)=[];

                    xbplot = x(borderedges(:,1));
                    xiplot = x(interioredges(:,1));
                    ybplot = y(borderedges(:,1));
                    yiplot = y(interioredges(:,1));
                    ub = transpose(ub);
                    vb = transpose(vb);
                    ui = transpose(ui);
                    vi = transpose(vi);

                    if dim==3
                        zbplot = z(borderedges(:,1));
                        ziplot = z(interioredges(:,1));
                        wb = transpose(wb);
                        wi = transpose(wi);
                    end

                    borderinput = bweight;
                    interiorinput = iweight;
                    if dim==3
                        quiver3(xbplot,ybplot,zbplot,ub,vb,wb,0,'ShowArrowHead','off','LineWidth',borderinput,'Color','red');
                        quiver3(xiplot,yiplot,ziplot,ui,vi,wi,0,'ShowArrowHead','off','LineWidth',interiorinput,'Color','red');
                    else
                        quiver(xbplot,ybplot,ub,vb,0,'ShowArrowHead','off','LineWidth',borderinput,'Color','red');
                        quiver3(xiplot,yiplot,ui,vi,0,'ShowArrowHead','off','LineWidth',interiorinput,'Color','red');
                    end
                end

                if inputWeightChoice == 4
                    sensoru = undirectedEdges(:,1);
                    sensedu = undirectedEdges(:,2);
                    correspondingUDistances = zeros(length(undirectedEdges),1);

                    for l=1:length(undirectedEdges)
                        sensorNode = sensoru(l);
                        sensedNode = sensedu(l);
                        correspondingUDistances(l,:) = distanceMatrix(sensorNode,sensedNode);
                    end

                    for l=1:length(correspondingUDistances)
                        if correspondingUDistances(l) == Inf
                            correspondingUDistances(l) = 0;
                        end
                    end
                    maxDistanceU = max(correspondingUDistances(:));
                    for k=1:length(correspondingUDistances)
                        red(k) = correspondingUDistances(k)/maxDistanceU;
                        if dim==3
                            obj = quiver3(xplotu(k),yplotu(k),zplotu(k),u(k),v(k),w(k),0,'Color',[red(k),0,0]);
                            set(obj,'ShowArrowHead','off');
                        else
                            obj = quiver(xplotu(k),yplotu(k),u(k),v(k),0,'Color',[red(k),0,0]);
                            set(obj,'ShowArrowHead','off');
                        end
                    end
                end
                if inputWeightChoice == 5
                   for p = 1:length(correspondingAngles)
                    lines(p) = 1/correspondingAngles(p);
                    lines2(p) = lines(p)^2;
                   end
                    sqrt_sum_lines = sqrt(sum(lines2(:)));
                   for p = 1:length(correspondingAngles)
                     norm_lines(p) = lines(p)/sqrt_sum_lines;
                   end
                    for k=1:length(correspondingAngles)
                        if dim==3
                            obj = quiver3(xplot(k),yplot(k),zplot(k),u(k),v(k),w(k),0,'ShowArrowHead','off','LineWidth',15*norm_lines(k),'AutoScale','off');
                            set(obj,'AutoScale','off');
                        else
                            obj = quiver(xplot(k),yplot(k),u(k),v(k),0,'ShowArrowHead','off','LineWidth',15*norm_lines(k),'AutoScale','off');
                        
                            set(obj,'AutoScale','off');
                        end
                    end
                end
            end
            
            if ~isempty(directedEdges)
                r=1;
                for j=1:length(directedEdges)
                    xVal1 = x(directedEdges(j,1));
                    xVal2 = x(directedEdges(j,2));
                    a(j) = xVal2 - xVal1;

                    yVal1 = y(directedEdges(j,1));
                    yVal2 = y(directedEdges(j,2));
                    b(j) = yVal2 - yVal1;

                    if dim==3
                        zVal1 = z(directedEdges(j,1));
                        zVal2 = z(directedEdges(j,2));
                        d(j) = zVal2 - zVal1;
                    end

                    if inputWeightChoice == 3
                        if j<=divider*numberOfNeighbors
                            ab(j) = a(j);
                            bb(j) = b(j);
                            if dim==3
                                db(j) = d(j);
                            end
                            borderedgesd(j,:) = directedEdges(j,:);
                        else
                            ai(r) = a(j);
                            bi(r) = b(j);
                            if dim==3
                                di(r) = d(j);
                            end
                            interioredgesd(j,:) = directedEdges(j,:);
                            r=r+1;
                        end
                    end
                end

                xplotd = x(directedEdges(:,1));
                yplotd = y(directedEdges(:,1));
                a = transpose(a);
                b = transpose(b);
                if dim==3
                    zplotd = z(directedEdges(:,1));
                    d = transpose(d);
                end

                if inputWeightChoice == 1
                    edgeWeight = weight;
                    if dim==3
                        quiver3(xplotd,yplotd,zplotd,a,b,d,'AutoScale','off','LineWidth',edgeWeight,'Color','blue');
                    else quiver(xplotd,yplotd,a,b,'AutoScale','off','LineWidth',edgeWeight,'Color','blue');
                    end
                end

                if inputWeightChoice == 2
                    edgeWeight = weight;
                    if dim==3
                        quiver3(xplotd,yplotd,zplotd,a,b,d,'AutoScale','off','LineWidth',edgeWeight,'Color','blue');
                    else quiver(xplotd,yplotd,a,b,'AutoScale','off','LineWidth',edgeWeight,'Color','blue');
                    end
                end

                if inputWeightChoice == 3
                    [rw,cl]=find(borderedgesd==0);
                    borderedgesd(rw,:)=[];

                    [rw,cl]=find(interioredgesd==0);
                    interioredgesd(rw,:)=[];

                    xbplotd = x(borderedgesd(:,1));
                    xiplotd = x(interioredgesd(:,1));
                    ybplotd = y(borderedgesd(:,1));
                    yiplotd = y(interioredgesd(:,1));
                    ab = transpose(ab);
                    bb = transpose(bb);
                    ai = transpose(ai);
                    bi = transpose(bi);
                    if dim==3
                        zbplotd = z(borderedgesd(:,1));
                        ziplotd = z(interioredgesd(:,1)); 
                        db = transpose(db);
                        di = transpose(di);
                    end

                    borderinput = bweight;
                    interiorinput = iweight;
                    if dim==3
                        quiver3(xbplotd,ybplotd,zbplotd,ab,bb,db,0,'AutoScale','off','LineWidth',borderinput,'Color','blue');
                        quiver3(xiplotd,yiplotd,ziplotd,ai,bi,di,0,'AutoScale','off','LineWidth',interiorinput,'Color','blue');
                    else
                        quiver(xbplotd,ybplotd,ab,bb,0,'AutoScale','off','LineWidth',borderinput,'Color','blue');
                        quiver(xiplotd,yiplotd,ai,bi,0,'AutoScale','off','LineWidth',interiorinput,'Color','blue');
                    end
                end

                if inputWeightChoice == 4
                    sensord = directedEdges(:,1);
                    sensedd = directedEdges(:,2);
                    correspondingDDistances = zeros(length(directedEdges),1);

                    for l=1:length(directedEdges)
                        sensorNode = sensord(l);
                        sensedNode = sensedd(l);
                        correspondingDDistances(l,:) = distanceMatrix(sensorNode,sensedNode);
                    end
                    for g=1:length(correspondingDDistances)
                        if correspondingDDistances(g) == Inf
                            correspondingDDistances(g) = 0;
                        end
                    end
                    maxDistance = max(correspondingDDistances(:));
                    for g=1:length(correspondingDDistances)
                        blue(g) = correspondingDDistances(g)/maxDistance;
                        if dim==3
                            obj = quiver3(xplotd(g),yplotd(g),zplotd(g),a(g),b(g),d(g),0,'Color',[0,0,blue(g)]);
                            set(obj,'AutoScale','off');
                        else
                            obj = quiver(xplotd(g),yplotd(g),a(g),b(g),0,'Color',[0,0,blue(g)]);
                            set(obj,'AutoScale','off');
                        end
                    end
                end
            end
            
            if data_type == 1
                title(strcat(avi_file_name,' at t = ',num2str(i/20),' sec'),'Interpreter','none');
            else
                title(strcat(avi_file_name,' at frame = ',num2str(i)),'Interpreter','none');
            end
            axis equal
            %             saveas(gcf,strcat('EdgeDistinctionGraph_t=',num2str(i)),'fig');
            hold off;
        end
        if inputGraphChoice == 3
            reverseNodePairs = horzcat(sensed, sensor);
            check = ismember(sortedNodePairs, reverseNodePairs, 'rows');
            
            for h=1:length(reverseNodePairs)
                if check(h)==0
                    directedEdges(h,:) = sortedNodePairs(h,:);
                end
            end
            
            [r,c]=find(directedEdges==0); %#ok<*NASGU>
            directedEdges(r,:)=[];
            
            r=1;
            for j=1:length(directedEdges)
                xVal1 = x(directedEdges(j,1));
                xVal2 = x(directedEdges(j,2));
                a(j) = xVal2 - xVal1;
                
                yVal1 = y(directedEdges(j,1));
                yVal2 = y(directedEdges(j,2));
                b(j) = yVal2 - yVal1;
                
                if dim==3
                    zVal1 = z(directedEdges(j,1));
                    zVal2 = z(directedEdges(j,2));
                    d(j) = zVal2 - zVal1;
                end
                
                if inputWeightChoice == 3
                    if j<=divider*numberOfNeighbors
                        ub(j) = a(j);
                        vb(j) = b(j);
                        if dim==3
                            wb(j) = d(j);
                        end
                        borderedges(j,:) = directedEdges(j,:);
                    else
                        ui(r) = a(j);
                        vi(r) = b(j);
                        if dim==3
                            wi(r) = d(j);
                        end
                        interioredges(j,:) = directedEdges(j,:);
                        r=r+1;
                    end
                end
            end
            
            xplotd = x(directedEdges(:,1));
            yplotd = y(directedEdges(:,1));
            a = transpose(a);
            b = transpose(b);
            if dim==3
                zplotd = z(directedEdges(:,1));
                d = transpose(d);
            end
            
            if inputWeightChoice == 1
                edgeWeight = weight;
                if dim==3
                    quiver3(xplotd,yplotd,zplotd,a,b,d,0,'AutoScale','off','Color','b','LineWidth',edgeWeight);
                else quiver(xplotd,yplotd,a,b,0,'AutoScale','off','Color','b','LineWidth',edgeWeight);
                end
            end
            
            if inputWeightChoice == 2
                edgeWeight = weight;
                if dim==3
                    quiver3(xplotd,yplotd,zplotd,a,b,d,0,'AutoScale','off','Color','b','LineWidth',edgeWeight);
                else quiver(xplotd,yplotd,a,b,0,'AutoScale','off','Color','b','LineWidth',edgeWeight);
                end
            end
            
            if inputWeightChoice == 3
                [rw,cl]=find(borderedges==0);
                borderedges(rw,:)=[];
                
                [rw,cl]=find(interioredges==0);
                interioredges(rw,:)=[];
                
                xbplot = x(borderedges(:,1));
                xiplot = x(interioredges(:,1));
                ybplot = y(borderedges(:,1));
                yiplot = y(interioredges(:,1));
                ub = transpose(ub);
                vb = transpose(vb);
                ui = transpose(ui);
                vi = transpose(vi);
                if dim==3
                    zbplot = z(borderedges(:,1));
                    ziplot = z(interioredges(:,1));
                    wb = transpose(wb);
                    wi = transpose(wi);
                end
                
                borderinput = bweight;
                interiorinput = iweight;
                if dim==3
                    quiver3(xbplot,ybplot,zbplot,ub,vb,wb,0,'AutoScale','off','LineWidth',borderinput,'Color','blue');
                    quiver3(xiplot,yiplot,ziplot,ui,vi,wi,0,'AutoScale','off','LineWidth',interiorinput,'Color','blue');
                else 
                    quiver(xbplot,ybplot,ub,vb,0,'AutoScale','off','LineWidth',borderinput,'Color','blue');
                    quiver(xiplot,yiplot,ui,vi,0,'AutoScale','off','LineWidth',interiorinput,'Color','blue');
                end
            
            end
            
            if inputWeightChoice == 4
                sensord = directedEdges(:,1);
                sensedd = directedEdges(:,2);
                correspondingDDistances = zeros(length(directedEdges),1);
                
                for g=1:length(directedEdges)
                    sensorNode = sensord(g);
                    sensedNode = sensedd(g);
                    correspondingDDistances(g,:) = distanceMatrix(sensorNode,sensedNode);
                end
                for g=1:length(correspondingDDistances)
                    if correspondingDDistances(g) == Inf
                        correspondingDDistances(g) = 0;
                    end
                end
                maxDistance = max(correspondingDDistances(:));
                for l=1:length(correspondingDDistances)
                    blue(l) = correspondingDDistances(l)/maxDistance;
                    red(l) = correspondingDDistances(l)/maxDistance;
                    if dim==3
                        obj = quiver3(xplotd(l),yplotd(l),zplotd(l),a(l),b(l),d(l),0,'Color',[red(l),0,blue(l)]);
                        set(obj,'AutoScale','off');
                    else
                        obj = quiver(xplotd(l),yplotd(l),a(l),b(l),0,'Color',[red(l),0,blue(l)]);
                        set(obj,'AutoScale','off');
                    end
                end
            end
            if data_type == 1
                title(strcat(avi_file_name,' at t = ',num2str(i/20),' sec'),'Interpreter','none');
            else
                title(strcat(avi_file_name,' at frame = ',num2str(i)),'Interpreter','none');
            end
            axis equal
%             saveas(gcf,strcat('DirectedGraph_t=',num2str(i)),'fig');
            hold off;
        end
        if inputGraphChoice == 4
            reverseNodePairs = horzcat(sensed, sensor);
            check = ismember(sortedNodePairs, reverseNodePairs, 'rows');
            
            for l=1:length(reverseNodePairs)
                if check(l)==1
                    undirectedEdges(l,:) = sortedNodePairs(l,:);
                end
            end
            
            [r,c]=find(undirectedEdges==0);
            undirectedEdges(r,:)=[];
            
            r=1;
            for j=1:length(undirectedEdges)
                xVal1 = x(undirectedEdges(j,1));
                xVal2 = x(undirectedEdges(j,2));
                u(j) = xVal2 - xVal1;
                
                yVal1 = y(undirectedEdges(j,1));
                yVal2 = y(undirectedEdges(j,2));
                v(j) = yVal2 - yVal1;
                
                if dim==3
                    zVal1 = z(undirectedEdges(j,1));
                    zVal2 = z(undirectedEdges(j,2));
                    w(j) = zVal2 - zVal1;
                end
                
                if inputWeightChoice == 3
                    if j <= divider*numberOfNeighbors
                        ub(j) = u(j);
                        vb(j) = v(j);
                        if dim==3
                            wb(j) = w(j);
                        end
                        borderedges(j,:) = undirectedEdges(j,:);
                    else
                        ui(r) = u(j);
                        vi(r) = v(j);
                        if dim==3
                            wi(r) = w(j);
                        end
                        interioredges(j,:) = undirectedEdges(j,:);
                        r = r+1;
                    end
                end
            end
            
            xplotu = x(undirectedEdges(:,1));
            yplotu = y(undirectedEdges(:,1));
            u = transpose(u);
            v = transpose(v);
            if dim==3
                zplotu = z(undirectedEdges(:,1));
                w = transpose(w);
            end
            
            if inputWeightChoice == 1
                edgeWeight = weight;
                if dim==3
                    quiver3(xplotu,yplotu,zplotu,u,v,w,0,'ShowArrowHead','off','Color','red','LineWidth',edgeWeight);
                else quiver(xplotu,yplotu,u,v,0,'ShowArrowHead','off','Color','red','LineWidth',edgeWeight);
                end
            end
            
            if inputWeightChoice == 2
                edgeWeight = weight;
                if dim==3
                    quiver3(xplotu,yplotu,zplotu,u,v,w,0,'ShowArrowHead','off','Color','red','LineWidth',edgeWeight);
                else quiver(xplotu,yplotu,u,v,0,'ShowArrowHead','off','Color','red','LineWidth',edgeWeight);
                end
            end
            
            if inputWeightChoice == 3
                [rw,cl]=find(borderedges==0);
                borderedges(rw,:)=[];
                
                [rw,cl]=find(interioredges==0);
                interioredges(rw,:)=[];
                
                xbplot = x(borderedges(:,1));
                xiplot = x(interioredges(:,1));
                ybplot = y(borderedges(:,1));
                yiplot = y(interioredges(:,1));
                ub = transpose(ub);
                vb = transpose(vb);
                ui = transpose(ui);
                vi = transpose(vi);
                if dim==3
                    zbplot = z(borderedges(:,1));
                    ziplot = z(interioredges(:,1));
                    wb = transpose(wb);
                    wi = transpose(wi);
                end
                
                borderinput = bweight;
                interiorinput = iweight;
                if dim==3
                    quiver3(xbplot,ybplot,zbplot,ub,vb,wb,0,'ShowArrowHead','off','LineWidth',borderinput,'Color','red');
                    quiver3(xiplot,yiplot,ziplot,ui,vi,wi,0,'ShowArrowHead','off','LineWidth',interiorinput,'Color','red');
                else
                    quiver(xbplot,ybplot,zbplot,ub,vb,0,'ShowArrowHead','off','LineWidth',borderinput,'Color','red');
                    quiver(xiplot,yiplot,ziplot,ui,vi,0,'ShowArrowHead','off','LineWidth',interiorinput,'Color','red');
                end
            end
            
            if inputWeightChoice == 4
                sensoru = undirectedEdges(:,1);
                sensedu = undirectedEdges(:,2);
                correspondingUDistances = zeros(length(undirectedEdges),1);
                for l=1:length(undirectedEdges)
                    sensorNode = sensoru(l);
                    sensedNode = sensedu(l);
                    correspondingUDistances(l,:) = distanceMatrix(sensorNode,sensedNode);
                end
                
                for l=1:length(correspondingUDistances)
                    if correspondingUDistances(l) == Inf
                        correspondingUDistances(l) = 0;
                    end
                end
                maxDistanceU = max(correspondingUDistances(:));
                for k=1:length(correspondingUDistances)
                    red(k) = correspondingUDistances(k)/maxDistanceU;
                    if dim==3
                        obj = quiver3(xplotu(k),yplotu(k),zplotu(k),u(k),v(k),w(k),0,'Color',[red(k),0,0]);
                        set(obj,'ShowArrowHead','off');
                    else
                        obj = quiver(xplotu(k),yplotu(k),u(k),v(k),0,'Color',[red(k),0,0]);
                        set(obj,'ShowArrowHead','off');
                    end
                end
            end
            if subject == 1
                title(strcat(avi_file_name,' at t = ',num2str(i/20),' sec'),'Interpreter','none');
            else
                title(strcat(avi_file_name,' at frame = ',num2str(i)),'Interpreter','none');
            end
            axis equal
%             as(gcf,strcat('UndirectedGraph_t=',num2str(i)),'fig');
            hold off;
        end
        if inputGraphChoice == 5
            for n=1:length(sortedNodePairs)
                j=1;
                nodePath(1,:) = sortedNodePairs(n,:);
                nodeConnection = nodePath(1,2);
                
                breakFlag = false;
                while nodeConnection ~= nodePath(1,1)
                    [nodeSensor,nodeSensed] = find(sortedNodePairs(:,1)==nodeConnection,1,'first');
                    j=j+1;
                    if nodeSensor > 0
                        nodePath(j,:) = sortedNodePairs(nodeSensor,:);
                        clear('nodeSensor');
                        nodeConnection = nodePath(j,2);
                        if nodePath(j,1) == nodePath(j,2)
                            break
                        end
                        for q=2:length(sortedNodePairs)
                            if j > q
                                if nodePath(j-q,1) == nodePath(j-1,2) && nodePath(j-1,2) == nodePath(j,1)
                                    breakFlag = true;
                                    break
                                end
                            end
                        end
                        if breakFlag == true
                            break
                        end
                        if size(nodePath,1) > 1
                            if nodePath(j-1,1) == nodePath(j,2) && nodePath(j-1,2) == nodePath(j,1)
                                nodeConnection = 0;
                            end
                        end
                    else
                        break
                    end
                end
                
                if nodeConnection == nodePath(1,1)
                    edgesCycle = nodePath;
                    for l=1:numberOfNeighbors
                        for k=1:length(edgesCycle)
                            xVal1 = x(edgesCycle(k,1));
                            xVal2 = x(edgesCycle(k,2));
                            f(k) = xVal2 - xVal1;
                        
                            yVal1 = y(edgesCycle(k,1));
                            yVal2 = y(edgesCycle(k,2));
                            g(k) = yVal2 - yVal1;
                            
                            if dim==3
                                zVal1 = z(edgesCycle(k,1));
                                zVal2 = z(edgesCycle(k,2));
                                h(k) = zVal2 - zVal1;
                            end
                        end
                    end
                    
                    xplotcycle = x(edgesCycle(:,1));
                    yplotcycle = y(edgesCycle(:,1));
                    if dim==3
                        zplotcycle = z(edgesCycle(:,1));
                    end
                    
                    if size(f,1) ~= size(xplotcycle,1)
                        f = transpose(f);
                    end
                    
                    if size(g,1) ~= size(yplotcycle,1)
                        g = transpose(g);
                    end
                    
                    if dim==3
                        if size(h,1) ~= size(zplotcycle,1)
                            h = transpose(h);
                        end
                    end
                    red = rand(1);
                    green = rand(1);
                    blue = rand(1);
                    if dim==3
                        quiver3(xplotcycle,yplotcycle,zplotcycle,f,g,h,0,'AutoScale','off','ShowArrowHead','on','Color',[red,green,blue]);
                    else quiver(xplotcycle,yplotcycle,f,g,0,'AutoScale','off','ShowArrowHead','on','Color',[red,green,blue]);
                    end
                    clear('nodePath','xplotcycle','yplotcycle','zplotcycle','f','g','h');
                end
                clear('nodePath');
            end
            disp('Done with an image.');
            if data_type == 1
                title(strcat(avi_file_name,' at t = ',num2str(i/20),' sec'),'Interpreter','none');
            else
                title(strcat(avi_file_name,' at frame = ',num2str(i)),'Interpreter','none');
            end
            axis equal
%             saveas(gcf,strcat('IndependentClustersGraph_t=',num2str(i)),'fig');
            hold off;
        end
        if inputGraphChoice == 6
            edges = sortrows([sensor,sensed],1);
            
            u = zeros(length(edges),1);
            v = zeros(length(edges),1);
            if dim==3
                w = zeros(length(edges),1);
            end
            
            for j=1:length(edges)
                xVal1 = x(edges(j,1));
                xVal2 = x(edges(j,2));
                u(j) = xVal2 - xVal1;
                
                yVal1 = y(edges(j,1));
                yVal2 = y(edges(j,2));
                v(j) = yVal2 - yVal1;
                
                if dim==3
                    zVal1 = z(edges(j,1));
                    zVal2 = z(edges(j,2));
                    w(j) = zVal2 - zVal1;
                end
                
                if inputWeightChoice == 3
                    xrow1 = find(x==xVal1);
                    if xrow1 < divider
                        ub(j) = u(j);
                        vb(j) = v(j);
                        if dim==3
                            wb(j) = w(j);
                        end
                        borderedges(j,:) = edges(j,:);
                    else
                        ui(j) = u(j);
                        vi(j) = v(j);
                        if dim==3
                            wi(j) = w(j);
                        end
                        interioredges(j,:) = edges(j,:);
                    end
                end
            end
            
            xplot = x(edges(:,1));
            yplot = y(edges(:,1));
            if dim==3
                zplot = z(edges(:,1));
            end
            
            if inputWeightChoice == 1
                edgeWeight = weight;
                if dim==3
                    quiver3(xplot,yplot,zplot,u,v,w,0,'ShowArrowHead','off','LineWidth',edgeWeight);
                else quiver(xplot,yplot,u,v,0,'ShowArrowHead','off','LineWidth',edgeWeight);
                end
            end
            
            if inputWeightChoice == 2
                edgeWeights = zeros(numberOfNodes,1);
                for l =1:length(edges)
                    edgeWeights(l) = 1/length(find(edges(:,1)==edges(l,1)));
                    if dim==3
                        quiver3(xplot(l),yplot(l),zplot(l),u(l),v(l),w(l),0,'ShowArrowHead','off','LineWidth',2*edgeWeights(l));
                    else quiver(xplot(l),yplot(l),u(l),v(l),0,'ShowArrowHead','off','LineWidth',2*edgeWeights(l));
                    end
                end
            end
            
            if inputWeightChoice == 3
                [rw,cl]=find(borderedges==0);
                borderedges(rw,:)=[];
                
                [rw,cl]=find(interioredges==0);
                interioredges(rw,:)=[];
                
                xbplot = x(borderedges(:,1));
                xiplot = x(interioredges(:,1));
                ybplot = y(borderedges(:,1));
                yiplot = y(interioredges(:,1));
                
                ub = transpose(ub);
                [rw,cl]=find(ub==0);
                ub(rw,:)=[];
                vb = transpose(vb);
                [rw,cl]=find(vb==0);
                vb(rw,:)=[];
                ui = transpose(ui);
                [rw,cl]=find(ui==0);
                ui(rw,:)=[];
                vi = transpose(vi);
                [rw,cl]=find(vi==0);
                vi(rw,:)=[];
                
                if dim==3
                    zbplot = z(borderedges(:,1));
                    ziplot = z(interioredges(:,1));                
                    wb = transpose(wb);
                    wi = transpose(wi);
                    [rw,cl]=find(wb==0);
                    
                    wb(rw,:)=[];
                    [rw,cl]=find(wi==0);
                    wi(rw,:)=[];
                end
                
                borderinput = bweight;
                interiorinput = iweight;
                if dim==3
                    quiver3(xbplot,ybplot,zbplot,ub,vb,wb,0,'ShowArrowHead','off','LineWidth',borderinput);
                    quiver3(xiplot,yiplot,ziplot,ui,vi,wi,0,'ShowArrowHead','off','LineWidth',interiorinput);
                else
                    quiver(xbplot,ybplot,ub,vb,0,'ShowArrowHead','off','LineWidth',borderinput);
                    quiver(xiplot,yiplot,ui,vi,0,'ShowArrowHead','off','LineWidth',interiorinput);
                end
            end

               if inputWeightChoice == 4
                correspondingDistances = zeros(length(edges),1);
                for l=1:length(edges)
                    correspondingDistances(l) = distanceMatrix(sensor(l),sensed(l));
                    edgeWeight(l) = 1/(correspondingDistances(l));
                    if dim==3
                        quiver3(xplot(l),yplot(l),zplot(l),u(l),v(l),w(l),0,'ShowArrowHead','off','LineWidth',2*edgeWeight(l));
                    else quiver(xplot(l),yplot(l),u(l),v(l),0,'ShowArrowHead','off','LineWidth',2*edgeWeight(l));
                    end
                end
            end
            if subject == 1
                title(strcat(avi_file_name,' at t = ',num2str(i/20),' sec'),'Interpreter','none');
            else
                title(strcat(avi_file_name,' at frame = ',num2str(i)),'Interpreter','none');
            end
            axis equal
%             saveas(gcf,strcat('DirectInputRadiusGraph_t=',num2str(i)),'fig');
        end
    end
   
    if inputComputeChoice == 1 || inputComputeChoice == 3

% CALCULATIONS QWERTY
            BG_Obj = biograph(squeeze(adjacencyMatrix(i,:,:)));
            [stronglyConnected, C2] = conncomp(BG_Obj,'Directed',true);
%             save(strcat('stronglyConnectedComponentsNN',num2str(i)),'stronglyConnected');
            graphProperties(i).stronglyConnected = stronglyConnected;
            
            [weaklyConnected, C3] = conncomp(BG_Obj,'Weak',true); %#ok<*ASGLU>
%             save(strcat('weaklyConnectedComponentsNN',num2str(i)),'weaklyConnected');
            graphProperties(i).weaklyConnected = weaklyConnected;
            
            set(handles.numberStronglyConnectedComponents_editText,'String',num2str(stronglyConnected));
            set(handles.numberWeaklyConnectedComponents_editText,'String',num2str(weaklyConnected));
            
        N = numberOfNodes;
        A = squeeze(adjacencyMatrix(i,:,:));
        D = diag(sum(A,2));
        L = D - A;
        Q = zeros(N-1,N);
        
        for k = 1:(N-1)
            for j = 1:k
                Q(k,j) = -1/sqrt(k*(k+1));
            end
            Q(k,k+1) = k/sqrt(k*(k+1));
        end
        
        % calculate algebraic connectivity
        Lbarsym = (1/2)*Q*(L+L')*Q';
        eigen = eig(Lbarsym);
        eigen = sort(eigen);
        algebraicConnectivity = roundn(eigen(1,1),-5);
        set(handles.algebraicConnectivity_editText,'String',num2str(algebraicConnectivity));
%         save(strcat('algebraicConnectivity',num2str(i)),'algebraicConnectivity');
        graphProperties(i).algebraicConnectivity = algebraicConnectivity;
        
        % calculate speed of convergence
        eigen = eig(L);
        eigen = sort(eigen);
        secondE = eigen(2,1);
        speedOfConvergence = roundn(real(secondE),-5);
        set(handles.speedOfConvergence_editText,'String',num2str(speedOfConvergence));
%         save(strcat('SpeedOfConvergence',num2str(i)),'speedOfConvergence');
        graphProperties(i).speedOfConvergence = speedOfConvergence;
        if speedOfConvergence > 0
            set(handles.connectedOrNot_editText,'String','Yes');
        else
            set(handles.connectedOrNot_editText,'String','No');
        end
        
        % calculate h2 norm
        evals = sort(real(eig(L)));
        if evals(2) > 1e-10
            H2 = H2_norm(L, Q);
        else
            H2 = Inf;
        end
        set(handles.h2Norm_editText,'String',num2str(H2));
%         save(strcat('H2Norm',num2str(i)),'H2');
        graphProperties(i).H2Norm = H2;
        
        % calculate l2 gain
        if algebraicConnectivity > 0
            L2 = 1/algebraicConnectivity;
            set(handles.l2gain_editText,'String',num2str(L2));
%             save(strcat('L2gain',num2str(i)),'L2');
            graphProperties(i).L2Gain = L2;
        else
            set(handles.l2gain_editText,'String','---');
            graphProperties(i).L2Gain = NaN;
        end
        
        % calculate h2 norm eigenvalue bounds
%         if algebraicConnectivity <= 0
%             msgbox('Eigenvalue bounds on H2 Norm cannot be computed because algebraic connectivity is not greater than zero.');
%         end
        
        if speedOfConvergence > 0
            eigenL = sort(eig(L));
            sumLB = 0;
            for k=2:N
                sumLB = sumLB + 1/(real(eigenL(k)));
            end
            lowerEBound = sqrt((1/2) * sumLB);
        else
            lowerEBound = Inf;
        end
        
        set(handles.lowerEigenBound_editText,'String',num2str(lowerEBound));
%         save(strcat('lowerEigenBound',num2str(i)),'lowerEBound');
        graphProperties(i).lowerEBound = lowerEBound;
        
        if algebraicConnectivity > 0
            Lbar = Q*L*Q';
            LBS = (1/2)*(Lbar+Lbar');
            eigenLBS = eig(LBS);
            sumLBS = 0;
            for k=1:N-1
                sumLBS = sumLBS + 1/(eigenLBS(k));
            end
            upperEBound = sqrt((1/2)*sumLBS);
        else
            upperEBound = Inf;

        end
        
        set(handles.upperEigenBound_editText,'String',num2str(upperEBound));
%         save(strcat('upperEigenBound',num2str(i)),'upperEBound');
        graphProperties(i).upperEBound = upperEBound;
        
        
        %Node Status Calculations
        
        shortestReversedLengths = zeros(numberOfNodes);
        for ii = 1:numberOfNodes
            for jj = [1:(ii-1) (ii+1):numberOfNodes]
                if adjacencyMatrix(i,jj,ii) > 0
                    shortestReversedLengths(ii,jj) = adjacencyMatrix(i,jj,ii);
                else
                    shortestReversedLengths(ii,jj) = Inf;
                end
            end
        end
        
        for kk = 1:numberOfNodes
            for ii = 1:numberOfNodes
                for jj = 1:numberOfNodes
                    shortestReversedLengths(ii,jj) = min(shortestReversedLengths(ii,jj), shortestReversedLengths(ii,kk)+shortestReversedLengths(kk,jj));
                end
            end
        end
        
        inverseDistances = 1./shortestReversedLengths;
        for ii = 1:numberOfNodes
            nodeStatus(i,ii) = 1/(numberOfNodes - 1)*sum(inverseDistances(ii,[1:(ii-1) (ii+1):numberOfNodes]));
            
            if nodeStatus(i,ii)>threshold
                if i>lb
                 if nodeStatus(i-1,ii)>threshold
                    nodeStatusMx(ii, 2, cnt_ns(ii)) = nodeStatusMx(ii, 2, cnt_ns(ii))+1;
                    nodeStatusMx(ii, 3, cnt_ns(ii)) = nodeStatusMx(ii, 3, cnt_ns(ii))+nodeStatus(i,ii);
                 else
                    cnt_ns(ii) = cnt_ns(ii)+1;
                    nodeStatusMx(ii, 1, cnt_ns(ii))=i;
                    nodeStatusMx(ii, 2, cnt_ns(ii))=1;
                    nodeStatusMx(ii, 3, cnt_ns(ii))=nodeStatus(i,ii);
                 end
                else
                    cnt_ns(ii) = cnt_ns(ii)+1;
                    nodeStatusMx(ii, 1, cnt_ns(ii))=i;
                    nodeStatusMx(ii, 2, cnt_ns(ii))=1;
                    nodeStatusMx(ii, 3, cnt_ns(ii))=nodeStatus(i,ii);
                end
            end
            %{
            allEdges = unique([find(squeeze(adjacencyMatrix(i,ii,:))) find(squeeze(adjacencyMatrix(i,:,ii)'))]);
            adjacentAngles = sort(relativeAngles(ii,allEdges));
            angleDifferences = [diff(adjacentAngles) (adjacentAngles(1) - adjacentAngles(end) + 2*pi)];
            [maxDiff locationOfMax] = max(angleDifferences);
            textDirection = anglerestrict(adjacentAngles(locationOfMax) + maxDiff/2); 
            %}
            if inputComputeChoice == 1 || inputComputeChoice == 2  
                text(x(ii) + 0.04, y(ii), num2str(roundn(nodeStatus(i,ii),-2)),'HorizontalAlignment','center');
                text(x(ii)-0.1, y(ii), num2str(ii),'HorizontalAlignment','left'); %code I added to display node number alongside node status
              %text(x(ii) + 0.04*lengthScale*cos(textDirection), y(ii) + 0.04*lengthScale*sin(textDirection), num2str(roundn(nodeStatus(ii),-2)),'HorizontalAlignment','center');             
              %text(x(ii) - 0.07*lengthScale*cos(textDirection), y(ii) + 0.04*lengthScale*sin(textDirection), num2str(ii),'HorizontalAlignment','center');
            end
        end
        
        node_status_cum = node_status_cum+ nodeStatus(i,:);
    end
    
    if inputComputeChoice == 1 || inputComputeChoice == 2  
       frame = getframe(gcf);
      aviobj = addframe(aviobj,frame);
    end
    cla;
end
hold off;

if inputComputeChoice == 1 || inputComputeChoice == 2
  aviobj = close(aviobj);
end

disp('mean node status for nodes')
node_status_cum = node_status_cum/(ub-lb+1)

for ii=1:numberOfNodes
    node_status_std(ii) = std(nodeStatus(:,ii));
end

figure(2)
plot(1:numberOfNodes, node_status_cum, 'x', 'Color', 'b')

disp('periods with higher node statuses')
nodeStatusMx
cnt_ns

disp('mean node status overall')
mean(node_status_cum)

disp('mean node status sorted')
[node_status_cum_sorted, ind_node_stat] = sort(node_status_cum);
node_status_cum_sorted
ind_node_stat

disp('std of node status sorted')
[node_status_std_sorted, ind_node_stat] = sort(node_status_std);
node_status_std_sorted
ind_node_stat

figure(1)
barwitherr(node_status_std, node_status_cum)

for i=lb:ub
Wconnect(i) = graphProperties(i).weaklyConnected;
end

states = 1;

Wconnect_mx(1,1)=lb;
Wconnect_mx(2,1)=1;
Wconnect_mx(3,1)=Wconnect(lb);

for i=lb+1:ub
    if Wconnect(i)==Wconnect(i-1)
      Wconnect_mx(2,states) = Wconnect_mx(2,states)+1;
    else
        states=states+1;
      Wconnect_mx(1,states) = i;
      Wconnect_mx(2,states) = 1;
      Wconnect_mx(3,states) = Wconnect(i);
    end
end

Wconnect_mx

size(find(Wconnect==1))
size(find(Wconnect==2))
size(find(Wconnect==3))

%{
disp('node Status Mx')
for ii=1:numberOfNodes
    ii
  squeeze(nodeStatusMx(ii,:,:))
end
%}

end
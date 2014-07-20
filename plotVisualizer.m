function varargout = plotVisualizer(varargin)
% PLOTVISUALIZER MATLAB code for plotVisualizer.fig
%      PLOTVISUALIZER, by itself, creates a new PLOTVISUALIZER or raises
%      the existing singleton*.
%
%      H = PLOTVISUALIZER returns the handle to a new PLOTVISUALIZER or the
%      handle to the existing singleton*.
%
%      PLOTVISUALIZER('CALLBACK',hObject,eventData,handles,...) calls the
%      local function named CALLBACK in PLOTVISUALIZER.M with the given
%      input arguments.
%
%      PLOTVISUALIZER('Property','Value',...) creates a new PLOTVISUALIZER
%      or raises the existing singleton*.  Starting from the left, property
%      value pairs are applied to the GUI before plotVisualizer_OpeningFcn
%      gets called.  An unrecognized property name or invalid value makes
%      property application stop.  All inputs are passed to
%      plotVisualizer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% 
% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @plotVisualizer_OpeningFcn, ...
    'gui_OutputFcn',  @plotVisualizer_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
    
end
end

function plotVisualizer_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;
set(hObject,'toolbar','figure');
set(handles.NeighOrRad2,'SelectionChangeFcn',@NeighOrRad2_SelectionChangeFcn);
guidata(hObject, handles);
end

function varargout = plotVisualizer_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
cla;
end

function configureWHAT_pushbutton_Callback(~, ~, ~)
msgbox...
    ('To use dancer data, click this pushbutton once to reformat the data to create plots and graphs and make computations.','"Configure Dance Data" Help','help');
end

function plotDatawithBorderwhat_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click this to plot position data. You will have the opportunity to select data to work with, and you can select a preference to plot data with or without distinction between border and interior individuals.','"Plot Data with Border Distinction" Help','help');
end

function saveImage_pushbutton_Callback(hObject, ~, handles)
% saves current figure

fileName = inputdlg('Please enter a name for your figure.');
directoryName = uigetdir('','Please select a folder to save to.');
if directoryName == 0      %User pressed the "Cancel" button
    directoryName = '';      %choose the empty string for folder
end
filePath = fullfile(directoryName,fileName{1});  %Create file path
extensions = {'fig'};
for k = 1:length(extensions)
    title(fileName,'Interpreter','none')
    saveas(gcf,filePath,extensions{k});  %Save the file
    set(gcf,'PaperPositionMode','auto');
end

guidata(hObject, handles);  %update the handles
end

function savewhat_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click this to save the current figure as a .fig file for further MATLAB manipulation. You will be asked to provide a name for the file, as well as a choice for the directory to save it to.','"Save Image" Help','help');
end

function autoPlot_pushbutton_Callback(hObject, ~, handles)
% plots position data in current directory -- making the video of the
% changing graph
data_type = handles.data_type;
lb = handles.lb;
ub = handles.ub;
if data_type ~= 1
     inputChoice = str2num(cell2mat(inputdlg('Select from the following two options:                                         1) Plot with distinction between border and interior individuals                                    2) Plot without distinction between border and interior individuals')));
end
avi_file_name = cell2mat(inputdlg('Please enter a name for the avi file to be created:'));
aviobj = avifile(avi_file_name,'compression','none','fps',10);
xmin=Inf;
ymin=Inf;
zmin=Inf;
xmax=-Inf;
ymax=-Inf;
zmax=-Inf;
for i=lb:ub
    if data_type==1
        data=handles.dancedata{i};
        xtempmin = min(data(:,1));
        if xtempmin < xmin
            xmin = xtempmin;
        end
        ytempmin = min(data(:,2));
        if ytempmin < ymin
            ymin = ytempmin;
        end
        xtempmax = max(data(:,1));
        if xtempmax > xmax
            xmax = xtempmax;
        end
        ytempmax = max(data(:,2));
        if ytempmax > ymax
            ymax = ytempmax;
        end
    else %%% This may need to be changed for other data types!
        borderFilename = strcat(num2str(i),'b','.dat');
        borderData = load(borderFilename);
        interiorFilename = strcat(num2str(i),'i','.dat');
        interiorData = load(interiorFilename);
        x = [borderData(:,1);interiorData(:,1)];
        y = [borderData(:,2);interiorData(:,2)];
        z = [borderData(:,3);interiorData(:,3)];
       
        xtempmin = min(x);
        if xtempmin < xmin
            xmin = xtempmin;
        end
        ytempmin = min(y);
        if ytempmin < ymin
            ymin = ytempmin;
        end
        ztempmin = min(z);
        if ztempmin < zmin
            zmin = ztempmin;
        end
        xtempmax = max(x);
        if xtempmax > xmax
            xmax = xtempmax;
        end
        ytempmax = max(y);
        if ytempmax > ymax
            ymax = ytempmax;
        end
        ztempmax = max(z);
        if ztempmax > zmax
            zmax = ztempmax;
        end
    end
end

for i = lb:ub
     if data_type == 1
        cla;
        data=handles.dancedata{i};
        x = data(:,1);
        y = data(:,2);
       
        axis([xmin xmax ymin ymax]);
        grid on;
        hold on;
        plot(x,y,'g.');
        F = getframe(gcf);
        aviobj = addframe(aviobj,F);
        saveas(gcf,fullfile(strcat('t=',num2str(i))),'fig');
     else %%This may need to be changed for other data types!
        borderFilename = strcat(num2str(i),'b','.dat');
        borderData = load(borderFilename);
        dim = size(borderData,2);
        axis([xmin xmax ymin ymax zmin zmax]);
        grid on;
        hold on;
        xb = borderData(:,1);
        yb = borderData(:,2);
       
        if dim == 3
            zb = borderData(:,3);
            if inputChoice == 1
                plot3(xb,yb,zb,'bx');
                drawnow;
            else
                plot3(xb,yb,zb,'g.');
                drawnow;
            end
        elseif inputChoice == 1
            plot(xb,yb,'bx');
            drawnow;
        else
            plot(xb,yb,'g.');
            drawnow;
        end
        
        interiorFilename = strcat(num2str(i),'i','.dat');
        interiorData = load(interiorFilename);
        
        xi = interiorData(:,1);
        yi = interiorData(:,2);
        if dim == 3
            zi = interiorData(:,3);
            plot3(xi,yi,zi,'g.');
        else plot(xi,y,'g.');
        end
        drawnow;
        F = getframe(gcf);
        aviobj = addframe(aviobj,F);
        saveas(gcf,fullfile(strcat('t=',num2str(i))),'fig');
        cla;
    end
end
hold off;

aviobj = close(aviobj);
guidata(hObject, handles);
end

function autoplotwhat_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click on this to plot multiple sets of position data from current directory.  All figures will be saved as "t= (i) .fig" files for further MATLAB manipulation. The function will then compile all figures and create a time series animation saved as PlotAnimation.avi.','"Automatic Plot/Movie" Help','help');
end

function nearestNeighborsWhat_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click on this to visualize edges on the graph dictated by your input number of nearest neighbors. The function plots the position data. Then, it computes the Euclidean distance between each node and every other node. You will be asked to input your desired number of nearest neighbors, and it will compile a list of the nearest neighbor node pairs, as well as their corresponding distances. You will also have the opportunity to choose from a list of options for assigning edge weight. An Adjacency Matrix will be saved in the current MATLAB directory, named "UnweightedAdjacencyMatrix.mat". Then, MATLAB will plot directed arrows connecting all nearest neighbor node pairs.','"See Nearest Neighbor Edges" Help','help');
end

function distinctionWhat_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click on this to visualize the distinction between node connections as being directed (if the sensing only goes in one direction) or undirected (if the sensing goes in both directions). You will be asked to select the data you wish to plot, as well as the number of nearest neighbors you wish to use, and how to weight the edges. The function computes your input number of minimum distances between each node and every other node, and then goes on to figure out whether the node pair appears again in the sensor/sensed pairs, but in the reverse order. If it does, it means that both nodes are sensing each other and the edge is therefore undirected. If the node pair sensing only goes in one direction, then the edge is directed (arrow going from the sensor node to the sensed node).','"Directed and Undirected Edge Distinction" Help','help');
end

function biographNetworkWhat_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click on this button to visualize the nearest neighbor edges in a different way. The function will still ask you to select the data sets you would like to use, as well as the number of nearest neighbors, but instead of plotting position data, this will open a new figure in the MATLAB biograph viewer. The position data is therefore not a part of this visualization, but you can more clearly see connections between nodes. Note that it can take a few minutes to produce the biograph object.','"See Network in Separate Figure" Help','help');
end

function eitherorWhat_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click on either of these to visualize ONLY the directed edges or ONLY the undirected edges between nodes. You will be asked to select the data sets you wish to plot, the number of nearest neighbors, the edge weight, and the function will go through the same process as in "Directed & Undirected Edge Distinction" to separate the node connections based on whether they are directed or undirected.','"Just Directed" and "Just Undirected" Help','help');
end

function independentGroupsWhat_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click on this button to visualize independent groups within the network of individual agents.  A group is independent if the nodes connect to each other in a cyclic fashion, meaning that if the group were to be removed from the network, no other nodes would be affected. You will be asked to select the number of nearest neighbors you wish to visualize, as well as the data you will be working with. The function will plot the position data and each independent cluster in a different color.','"See Independent Groups within Network" Help','help');
end

function sliderMIN_editText_Callback(~, ~, ~)
end
function sliderMIN_editText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function sliderMAX_editText_Callback(~, ~, ~)
end
function sliderMAX_editText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function directInputRadiusEdgeWhat_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click this button to create graphs with edges based on a user input radius. You will be asked to select data to plot, and you will also be given the opportunity to specify edge weights. The figure will be updated with position data, and then you must type in your desired input radius. Use the buttons in the "Input Radius Toolbox" below. The minimum and maximum edge lengths will appear next to "MIN" and "MAX". Watch the command window for detailed, step by step instructions for this particular visulaization. During the first input radius visualization, the figure will pause construction at various points to allow the user time to view different aspects of the graph. After viewing the graph, you can input new values into the text box, or slide the slider to get new graphs. These visualizations will be more rapid, but you will still see a sphere with the input radius, followed by the edges. Each visualization will feature differently colored edges.','"See Edges from Input Radius" Help','help');
end

function numberStronglyConnectedComponents_editText_Callback(~, ~, ~)
end
function numberStronglyConnectedComponents_editText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function numberWeaklyConnectedComponents_editText_Callback(~, ~, ~)
end
function numberWeaklyConnectedComponents_editText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function computeConnectedComponentsWHAT_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click this to compute the numbers of both strongly and weakly connected graph components.','"Compute Connected Components" Help','help');
end

function algebraicConnectivity_editText_Callback(~, ~, ~)
end
function algebraicConnectivity_editText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function algebraicConnectivityWHAT_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click this to compute the algebraic connectivity of the graph.','"Algebraic Connectivity" Help','help');
end

function speedOfConvergence_editText_Callback(~, ~, ~)
end
function speedOfConvergence_editText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function speedOfConvergenceWHAT_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click this to calculate the speed of convergence of the consensus of the graph.','"Speed of Convergence" Help','help');
end

function h2Norm_editText_Callback(~, ~, ~)
end
function h2Norm_editText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function h2NormWHAT_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click this to compute the H2 Norm of the graph.','"H2 Norm" Help','help');
end

function l2gain_editText_Callback(~, ~, ~)
end
function l2gain_editText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function l2GainWHAT_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click this to compute the L2 Gain of the graph.','"L2 Gain" Help','help');
end

function lowerEigenBound_editText_Callback(~, ~, ~)
end
function lowerEigenBound_editText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function upperEigenBound_editText_Callback(~, ~, ~)
end
function upperEigenBound_editText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function eigBoundsWHAT_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click this to calculate the lower and upper eigenvalue bounds on the H2 Norm.','"H2 Norm Eigenvalue Bounds" Help','help');
end


function autoPlotandCompute_pushbutton_Callback(hObject, ~, handles)
% automatically creates graph based on user specification and computes and
% prints all graph properties to command window and saves them as matlab
% variables

lb = handles.lb;
ub = handles.ub;
data_type = handles.data_type;
inputEdgeChoice = str2num(get(handles.NeighOrRad2,'UserData'));
inputNewModel = get(handles.NewModel, 'Value');
inputNewModel2 = get(handles.NewModel2, 'Value');

if inputNewModel == 1 
  inputEdgeChoice = 3;
end

if inputNewModel2 == 1 & inputNewModel==0
  inputEdgeChoice = 4;
end

if inputNewModel2 == 1 & inputNewModel == 1
    inputEdgeChoice = 5;
end

ProbCheck_Status = get(handles.ProbCheck,'Value');
if inputEdgeChoice == 1
    numberOfNeighbors = str2num(get(handles.nuneighbors_editText,'String'));
elseif inputEdgeChoice == 2
    inputRadius = str2num(get(handles.inputradius_editText,'String'));
elseif inputEdgeChoice == 3 %if NewModel, need both # of neighbors AND a radius
    numberOfNeighbors = str2num(get(handles.nuneighbors_editText,'String'));
    inputRadius = str2num(get(handles.inputradius_editText,'String'));
elseif inputEdgeChoice == 4
    numberOfNeighbors = str2num(get(handles.nuneighbors_editText,'String')); %these have now a different meaning: # of neighbors to consider, r = radius around nodes
    inputRadius = str2num(get(handles.inputradius_editText,'String'));
end

if data_type == 1
    viewAngle = str2num(get(handles.viewAngle_editText, 'String'));
    if viewAngle < 0 || viewAngle > 360
        msgbox('The angle you entered is invalid. Please enter a number between 0 and 360.');
        return;
    end  
else
    viewAngle = 360;
end

inputWeightChoice = (get(handles.edges_popup,'Value'));
inputComputeChoice = str2num(cell2mat(inputdlg...
    ('Please choose whether you would like to plot and/or compute:                  Select 1) to plot graphs and compute every graph property                   Select 2) to only plot graphs                                                                                                                             Select 3) to only compute graph properties')));
if inputComputeChoice < 1 || inputComputeChoice > 3
    msgbox('Your selection was invalid. Please select from options 1-3.')
end
avi_file_name = cell2mat(inputdlg('Please enter a name for the avi file to be created:'));

if inputComputeChoice ~= 3
    if inputEdgeChoice == 1 | inputEdgeChoice == 3 | inputEdgeChoice ==4 | inputEdgeChoice ==5
        inputGraphChoice = str2num(cell2mat(inputdlg...
            ('Please enter your choice for the type of graphs to be created:            Select 1) to plot nearest neighbor edges without directed/undirected distinction                                                                                                                                                      Select 2) to plot nearest neighbor edges with distinction between directed and undirected edges                                                                               Select 3) to only plot nearest neighbor directed edges                                  Select 4) to only plot nearest neighbor undirected edges                                        Select 5) to see independent groups within each network of individuals')));
        if inputGraphChoice < 1 || inputGraphChoice > 5
            msgbox('Your selection was invalid. Please select from options 1-6.')
        end
    elseif inputEdgeChoice == 2
        inputGraphChoice = str2num(cell2mat(inputdlg...
        ('Please enter your choice for the type of graphs to be created:            Select 1) to plot all edges         Select 2) to see independent groups')));
        if inputGraphChoice == 2
            inputGraphChoice = 5;
        end
    end
end

%CHECK:
numberPictures = ub-lb;
graphProperties = struct('UWadjacencyMatrix', cell(1,numberPictures),'adjacencyMatrix',cell(1,numberPictures),'stronglyConnected',cell(1,numberPictures),'weaklyConnected',cell(1,numberPictures),'algebraicConnectivity',cell(1,numberPictures),'speedOfConvergence',cell(1,numberPictures),'H2Norm',cell(1,numberPictures),'L2Gain',cell(1,numberPictures),'lowerEBound',cell(1,numberPictures),'upperEBound',cell(1,numberPictures),'status',cell(1,numberPictures),'probabilityMatrix',cell(1,numberPictures));

if inputComputeChoice == 1 || inputComputeChoice == 2
   aviobj = avifile(strcat(avi_file_name,'.avi'),'compression','none','fps',10); %fps = 10??!
    xmin=Inf;
    ymin=Inf;
    zmin=Inf;
    xmax=-Inf;
    ymax=-Inf;
    zmax=-Inf;
    for i=lb:ub
        if data_type==1
            file=strcat('Dance',num2str(i),'.dat');
            data=handles.dancedata{i};
            xtempmin = min(data(:,1));
            if xtempmin < xmin
                xmin = xtempmin;
            end
            ytempmin = min(data(:,2));
            if ytempmin < ymin
                ymin = ytempmin;
            end
            xtempmax = max(data(:,1));
            if xtempmax > xmax
                xmax = xtempmax;
            end
            ytempmax = max(data(:,2));
            if ytempmax > ymax
                ymax = ytempmax;
            end
        else
            borderFilename = strcat(num2str(i),'b','.dat');
            borderData = load(borderFilename);
            interiorFilename = strcat(num2str(i),'i','.dat');
            interiorData = load(interiorFilename);
            x = [borderData(:,1);interiorData(:,1)];
            y = [borderData(:,2);interiorData(:,2)];
            z = [borderData(:,3);interiorData(:,3)];
            xtempmin = min(x);
            if xtempmin < xmin
                xmin = xtempmin;
            end
            ytempmin = min(y);
            if ytempmin < ymin
                ymin = ytempmin;
            end
            ztempmin = min(z);
            if ztempmin < zmin
                zmin = ztempmin;
            end
            xtempmax = max(x);
            if xtempmax > xmax
                xmax = xtempmax;
            end
            ytempmax = max(y);
            if ytempmax > ymax
                ymax = ytempmax;
            end
            ztempmax = max(z);
            if ztempmax > zmax
                zmax = ztempmax;
            end
        end
    end
    
    if data_type == 1
        lengthScale = min([xmax - xmin, ymax - ymin]);
    else
        lengthScale = min([xmax - xmin, ymax - ymin, zmax - zmin]);
    end
end

%determine number of nodes outside for i=lb:ub loop, for practical reasons
data=handles.dancedata{handles.lb};
x = data(:,1);
numberOfNodes = length(x);
node_status_cum = zeros(1,numberOfNodes);
switch_nodes = zeros(1, numberOfNodes);
neighborNumber = zeros(1, numberOfNodes);

tau1=15; %3 tau's: long/short tau and limit at which neighbor is possible
tau2=50;
tau3=15;

neighbor_data = zeros(numberOfNodes,3,1);
cnt = ones(1, numberOfNodes);

if inputEdgeChoice==1 | inputEdgeChoice==2 | inputEdgeChoice==3
adjacencyMatrix = zeros(ub-lb+1,numberOfNodes, numberOfNodes);

for i = lb:ub
    
    i
    %open files, read etc.
    cla;
    lb = handles.lb;
    ub = handles.ub;
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
        
        borderFilename = strcat(num2str(i),'b','.dat');
        borderData = load(borderFilename);
        dim = size(borderData,2);
        interiorFilename = strcat(num2str(i),'i','.dat');
        interiorData = load(interiorFilename);
        set(handles.dataselection_editText,'String',strcat(borderFilename,' & ',interiorFilename));
        x = [borderData(:,1); interiorData(:,1)];
        y = [borderData(:,2); interiorData(:,2)];
        if dim == 3
            z = [borderData(:,3); interiorData(:,3)];
            positionMatrix = horzcat(x,y,z);
        else positionMatrix = horzcat(x,y);
        end
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
       
    if inputEdgeChoice == 1
%         organizedDistanceMatrix = sort(distanceMatrix, 2);
%         kminDistances = organizedDistanceMatrix(:,1:numberOfNeighbors);
        
         probabilityMatrix = zeros(numberOfNodes,numberOfNeighbors);
%        probabilityMatrix = zeros(numberOfNodes);
         br = zeros(numberOfNodes);
         
        for s=1:numberOfNodes;  %the beginning of a huge for loop
            if viewAngle == 360
                distanceRow = distanceMatrix(s,:);
                organizedDistanceRow = sort(distanceRow);
                [sortedDist, sortedNodes] = sort(distanceRow);
                visibleNodes = 1:numberOfNodes;
                keep = ismember(sortedNodes, visibleNodes);
                index = keep == 0;
                visible_sorted = sortedNodes;
                visible_sorted(index) = [];
                minimumRow = organizedDistanceRow(1:numberOfNeighbors);
                
                adjacencyRow = ismember(distanceRow,minimumRow);
                adjacencyMatrix(i,s,:) = adjacencyRow;
            else
                viewMin = anglerestrict(theta(s) - pi*viewAngle/360);
                viewMax = anglerestrict(theta(s) + pi*viewAngle/360);
                
                if viewMin < viewMax
                    visibleNodes = find((relativeAngles(s,:) <= viewMax) & (relativeAngles(s,:) >= viewMin));
                else
                    visibleNodes = find((relativeAngles(s,:) <= viewMax) | (relativeAngles(s,:) >= viewMin));
                end
                
                if length(visibleNodes) < numberOfNeighbors
                    angles = relativeAngles(s,[1:(s-1) (s+1):numberOfNodes]);
                    sortedAngles = sort(angles);
                    
                    viewAngleRequired = zeros(1,numberOfNodes-1);
                    for ii = 1:(numberOfNodes-1)
                        if ii <= numberOfNodes - numberOfNeighbors
                            viewAngleRequired(ii) = sortedAngles(ii+numberOfNeighbors-1) - sortedAngles(ii);
                        else
                            viewAngleRequired(ii) = anglerestrict(sortedAngles(ii+numberOfNeighbors-numberOfNodes) - sortedAngles(ii));
                        end
                    end
                    
                    %number of neighbors goes down until a possible viewing
                    %angle is found for tempNumberofNeighbors; heading
                    %could be ignored! -all that matters is diff in angle!
                    tempNumberOfNeighbors = numberOfNeighbors;
                    while min(viewAngleRequired) > pi*viewAngle/180
                        tempNumberOfNeighbors = tempNumberOfNeighbors - 1;
                        
                        viewAngleRequired = zeros(1,numberOfNodes-1);
                        for ii = 1:(numberOfNodes-1)
                            if ii <= numberOfNodes - tempNumberOfNeighbors
                                viewAngleRequired(ii) = sortedAngles(ii+tempNumberOfNeighbors-1) - sortedAngles(ii);
                            else
                                viewAngleRequired(ii) = anglerestrict(sortedAngles(ii+tempNumberOfNeighbors-numberOfNodes) - sortedAngles(ii));
                            end
                        end
                    end
                    
                    distanceRow = distanceMatrix(s,:);
                    [mindist closestnode] = min(distanceRow);
                    direction = 2*(relativeAngles(s,closestnode) > 0) - 1;
                    
                    if direction > 0
                        startloc = find(sortedAngles > viewMin,1);
                        if isempty(startloc)
                            startloc = 1;
                        end
                    else
                        startloc = find(sortedAngles > viewMax,1) - tempNumberOfNeighbors;
                        if isempty(startloc)
                            startloc = numberOfNodes - 1 - tempNumberOfNeighbors;
                        end
                        if startloc < 1
                            startloc = startloc + numberOfNodes - 1;
                        end
                    end
                    
                    neededViewAngle = pi*viewAngle/180;
                    currentViewAngle = neededViewAngle + 1;
                    stepnum = 0;
                    while currentViewAngle > neededViewAngle
                        if mod(stepnum,2) ~= 0
                            location = startloc + direction*(stepnum+1)/2;
                        else
                            location = startloc - direction*stepnum/2;
                        end
                        
                        if location < 1
                            location = location + numberOfNodes - 1;
                        end
                        if location >= numberOfNodes
                            location = location - numberOfNodes + 1;
                        end
                    
                        currentViewAngle = viewAngleRequired(location);
                        
                        stepnum = stepnum + 1;
                    end
                    
                    if location <= (numberOfNodes - tempNumberOfNeighbors)
                        anglesOfVisibleNodes = sortedAngles(location:(location + tempNumberOfNeighbors - 1));
                    else
                        anglesOfVisibleNodes = sortedAngles([location:(numberOfNodes-1) 1:(location + tempNumberOfNeighbors - numberOfNodes)]);
                    end
                    
                    visibleLogic = ismember(relativeAngles(s,:),anglesOfVisibleNodes);
                    visibleNodes = find(visibleLogic == 1);
                    
                end
                
                distanceRow = distanceMatrix(s,:);
                organizedDistances = sort(distanceRow(visibleNodes));
                [sortedDist, sortedNodes] = sort(distanceRow);
                keep = ismember(sortedNodes, visibleNodes);
                index = keep == 0;
                visible_sorted = sortedNodes;
                visible_sorted(index) = [];
                if length(visibleNodes) < numberOfNeighbors
                    minimumRow = organizedDistances;
                else
                    minimumRow = organizedDistances(1:numberOfNeighbors);
                end
            end
            h=1;
% PROBABILITY SWITCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
             if (ProbCheck_Status)
                m = str2num(get(handles.Increment,'String'));
                mm = str2num(get(handles.Increment2,'String'));
                L = str2num(get(handles.max_prob,'String'));
                
                if i == lb
                    adjacencyRow = ismember(distanceRow,minimumRow);
                    neighbors = find(adjacencyRow == 1);
                    probabilityMatrix(s,:)= zeros(1,numberOfNeighbors);                    
                else
                 
                 if s==h   
                    disp('frame and node')
                    i
                    s
                 end
                                 
                    probabilityMatrix = graphProperties(i-1).probabilityMatrix;
                    old_neighbors = graphProperties(i-1).neighbors(s,:); % node neighbors at previous frame
                    old_neighbors2 = old_neighbors;
                    closest = ismember(distanceRow,minimumRow);
                    num_closest = find(closest==1); % closest neighbors in view angle
                   
                    potential_new = num_closest;
                %{    
                    if s==h
                     disp('FIRST old neighbors')
                     old_neighbors
                     disp('potential_new2_first')
                     potential_new
                    end
                  %}  
                    for j = 1:numberOfNeighbors
                         if (ismember(old_neighbors(1,j),num_closest)) % if was neighbor and still closest in view
                                probabilityMatrix(s,j) = 0;
                         elseif ismember(old_neighbors(1,j),num_closest) == 0 % if was neighbor, still in view but not closest(2)
                            if probabilityMatrix(s,j) < (L - m) % increase pswitch towards saturation (possibly have this increase fasterthan case above)
                                probabilityMatrix(s,j) = probabilityMatrix(s,j) + m;
                            else
                                probabilityMatrix(s,j)= L;
                            end
                         end
                         if probabilityMatrix(s,j) > 0
                            change(s,j) = binornd(1,probabilityMatrix(s,j));
                         else 
                            change(s,j) = 0;
                         end
                    end
                    num_changes = sum(change(s,:));
                    
%                     if length(visibleNodes) == numberOfNeighbors
%                         neighbors = visibleNodes;
                    if num_changes == 0
                        neighbors = old_neighbors;
                    else
                        potential_new = visible_sorted;
                        for j = 1:numberOfNeighbors
                            if (ismember(old_neighbors(1,j),potential_new)) 
                                 potential_new(potential_new == old_neighbors(1,j)) = [];
                            end
                        end
                        if length(potential_new) <= num_changes
                              for j = 1:numberOfNeighbors
                                  if (change(s,j))
                                    neighbors(1,j) = potential_new(1);
                                    potential_new(potential_new == neighbors(1,j)) = [];
                                    probabilityMatrix(s,j) = 0; %reset probability switch
                                  else
                                    neighbors(1,j) = old_neighbors(1,j);
                                  end
                              end
                        elseif length(potential_new)> num_changes
                                for j = 1:numberOfNeighbors
                                    if (change(s,j))
                                        if numel(potential_new) == 0 %debug!
                                            stop = 1;
                                        end
                                        neighbors(1,j) = potential_new(1);
                                        potential_new(potential_new == neighbors(1,j)) = [];
                                        probabilityMatrix(s,j) = 0; %reset probability switch
                                    else
                                        neighbors(1,j) = old_neighbors(1,j);
                                    end
                                end
                        else
                            error = 1;
                        end
                    end
                    adjacencyRow = zeros(1,numberOfNodes);
                    adjacencyRow(1,neighbors) = 1;
                end
             else 
                 adjacencyRow = ismember(distanceRow,minimumRow);
                 neighbors = find(adjacencyRow == 1);
             end
             
                 adjacencyMatrix(i,s,:) = adjacencyRow;
                 graphProperties(i).probabilityMatrix(s,:) = probabilityMatrix(s,:);
                 graphProperties(i).neighbors(s,:) = neighbors;
              
                if i==lb
                    neighborNumber(s) = neighborNumber(s)+1;
                    Invneighbors(s, neighborNumber(s),:)=neighbors;
                    fr_neigh(s, neighborNumber(s))=1;
                else
                if sum(neighbors==graphProperties(i-1).neighbors(s,:))~=numberOfNeighbors %this remembers the neighbors and how many
                    neighborNumber(s) = neighborNumber(s)+1
                    Invneighbors(s, neighborNumber(s),:)=neighbors
                    fr_neigh(s, neighborNumber(s)) = 1
                else
                    fr_neigh(s, neighborNumber(s)) = fr_neigh(s, neighborNumber(s))+1;
                end
                end
                
             end

                graphProperties(i).UWadjacencyMatrix = squeeze(adjacencyMatrix(i,:,:)); %%unweighted!!!
                
             if ProbCheck_Status & i~=lb   
                clear neighbors
                clear change
             end
             
                if adjacencyRow(s) ~= 0
                    disp(' ')
                    disp('Problem: edge connecting a node to itself')
                    disp('-----------------------------------------')
                    disp(['Frame: ' num2str(i)])
                    disp(['Node: ' num2str(s)])
                    distanceRow
                    viewMin
                    viewMax
                    theta  
                    relativeAngles(s,:)
                    if viewMin < viewMax
                        find((relativeAngles(s,:) <= viewMax) & (relativeAngles(s,:) >= viewMin))
                    else
                        find((relativeAngles(s,:) <= viewMax) | (relativeAngles(s,:) >= viewMin))
                    end
                    visibleNodes
                    adjacencyRow
                    
                    pause
                end
                
                Ndist(i,s,:) = distanceMatrix(s, find(adjacencyMatrix(i,s,:)==1)); %store in N distances to neighbors. i=frame, o=node, :=dist to neighbors 
                Ndist = sort(Ndist,3);            
                    
                handles.adjacencyMatrix = adjacencyMatrix;
                
    elseif inputEdgeChoice == 2
        adjacencyMatrix = zeros(numberOfNodes);
        for o=1:numberOfNodes
            for j=1:numberOfNodes
                if viewAngle == 360
                    if distanceMatrix(o,j) < inputRadius
                        adjacencyMatrix(o,j) = 1;
                    end
                else
                    viewMin = anglerestrict(theta(o) - pi*viewAngle/360);
                    viewMax = anglerestrict(theta(o) + pi*viewAngle/360);
                    
                    if viewMin < viewMax
                        visible = (relativeAngles(o,j) <= viewMax) & (relativeAngles(o,j) >= viewMin);
                    else
                        visible = (relativeAngles(o,j) <= viewMax) | (relativeAngles(o,j) >= viewMin);
                    end
                    
                    if visible && distanceMatrix(o,j) < inputRadius
                        adjacencyMatrix(o,j) = 1;
                    end
                end
            end
            
            Ndist(i,o,:) = distanceMatrix(o, find(adjacencyMatrix(o,:)==1)); %store in N distances to neighbors. i=frame, o=node, :=dist to neighbors 
            Ndist = sort(Ndist,3);

        end
   
 elseif inputEdgeChoice==3 %NewModel stuff!!!
        
        for o=1:numberOfNodes
            
                if viewAngle == 360
                  distanceRow = distanceMatrix(o,:);
                  organizedDistanceRow = sort(distanceRow);
                  minimumRow = organizedDistanceRow(1:numberOfNeighbors);
                                                                         
                  adjacencyRow = ismember(distanceRow,minimumRow);
                  adjacencyMatrix(i,o,:) = adjacencyRow;
                  potential_neighbors = find(adjacencyRow==1);
                  
                else
                    
                    viewMin = anglerestrict(theta(o) - pi*viewAngle/360);
                    viewMax = anglerestrict(theta(o) + pi*viewAngle/360);
                    
                    if viewMin < viewMax
                        visible = find(relativeAngles(o,:) <= viewMax & relativeAngles(o,:) >= viewMin);
                    else
                        visible = find(relativeAngles(o,:) <= viewMax | relativeAngles(o,:) >= viewMin);
                    end
                    
                   for j=1:size(visible,2) 
                    if distanceMatrix(o,visible(j)) < inputRadius
                        adjacencyMatrix(i,o,visible(j)) = 1;
                    end
                   end
                end
          
            visible = find(adjacencyMatrix(i,o,:)==1);
            
            %{
            if o==8
                o
               find(adjacencyMatrix(o,:)==1)
            end
            %}     
                   
               if sum(adjacencyMatrix(i,o,:)) > numberOfNeighbors
                   ind_radii = find(adjacencyMatrix(i,o,:)==1);
                   distanceMatrix_withinangle = distanceMatrix(o,ind_radii);
                   distanceMatrix_withinangle_sorted = sort(distanceMatrix_withinangle);
                   Reduced_distMatrix_withinangle = ...
                   distanceMatrix_withinangle_sorted(1:numberOfNeighbors);
               
                   g=0;
                   neigh = [];
                   for f = 1:numberOfNodes
                      if adjacencyMatrix(i,o,f)==1 & ismember(distanceMatrix(o,f), Reduced_distMatrix_withinangle)  
                        neigh(g+1) = f;
                        g=g+1;  
                      end
                   end 
                   
                   potential_neighbors=neigh;
                   
              %{
                   if o==8
                distanceMatrix_withinangle
                distanceMatrix_withinangle_sorted
                Reduced_distMatrix_withinangle
                   end
              %}
               elseif sum(adjacencyMatrix(i,o,:)) < numberOfNeighbors
                     neighbors_left = numberOfNeighbors - sum(adjacencyMatrix(i,o,:));
                     ind2 = find(adjacencyMatrix(i,o,:)==0);
                     [distanceMatrix_sorted, ind_dM_sorted] = sort(distanceMatrix(o,:));
                     
                     j=0;
                     k=1; 
                     ind_dM_sorted2 = zeros(1, neighbors_left);
                     while j < neighbors_left
                       if ismember(ind_dM_sorted(k), ind2)  
                        ind_dM_sorted2(j+1) = ind_dM_sorted(k);
                        j=j+1;    
                       end
                     k=k+1;
                     end
                     
                     potential_neighbors = [visible, ind_dM_sorted2];
                     %{   
                     if o==8
                        neighbors_left
                        ind2
                        distanceMatrix_sorted
                        ind_dM_sorted
                        ind_dM_sorted2
                        adjacencyMatrix(o,:)
                     end
                      %}
                   else
                       potential_neighbors = visible;
                       
               end
              
               h=5;
               
               if (ProbCheck_Status)
                m = str2num(get(handles.Increment,'String'));      
                mm =  str2num(get(handles.Increment2,'String'));
                L = str2num(get(handles.max_prob,'String'));
                
                 if i == lb
                    neighbors = potential_neighbors;
                    probabilityMatrix(o,:)= zeros(1,numberOfNeighbors);                    
                    graphProperties(i).probabilityMatrix = probabilityMatrix;
                    adjacencyRow = zeros(1, numberOfNodes);
                    adjacencyRow(1,neighbors) = 1;
                else
                 
                 if o==h   
                    disp('frame and node')
                    i
                    o
                 end
                                 
                    probabilityMatrix = graphProperties(i-1).probabilityMatrix;
                    old_neighbors = graphProperties(i-1).neighbors(o,:); % node neighbors at previous frame
                    old_neighbors2 = old_neighbors;                   
                    potential_new = potential_neighbors;
                    
                    if o==h
                     disp('FIRST old neighbors')
                     old_neighbors
                     disp('potential_new2_first')
                     potential_new
                    end
                    
                    for j = 1:numberOfNeighbors
                         if (ismember(old_neighbors(1,j),potential_new)) % if was neighbor and still closest in view
                                probabilityMatrix(o,j) = 0;
                         elseif ismember(old_neighbors(1,j),potential_new) == 0 % if was neighbor, still in view but not closest(2)
                            if probabilityMatrix(o,j) < (L - m) % increase pswitch towards saturation (possibly have this increase fasterthan case above)
                                probabilityMatrix(o,j) = probabilityMatrix(o,j) + m;
                            else
                                probabilityMatrix(o,j)= L;
                            end
                         end
                         if probabilityMatrix(o,j) > 0
                            change(o,j) = binornd(1,probabilityMatrix(o,j));
                         else 
                            change(o,j) = 0;
                         end
                    end
                    num_changes = sum(change(o,:));
                    
%                     if length(visibleNodes) == numberOfNeighbors
%                         neighbors = visibleNodes;
                    if num_changes == 0
                        neighbors = old_neighbors;
                    else
                        for j = 1:numberOfNeighbors
                            if (ismember(old_neighbors(1,j),potential_new)) 
                                 potential_new(potential_new == old_neighbors(1,j)) = [];
                            end
                        end
                        if length(potential_new) <= num_changes
                              for j = 1:numberOfNeighbors
                                  if (change(o,j))
                                    neighbors(1,j) = potential_new(1);
                                    potential_new(potential_new == neighbors(1,j)) = [];
                                    probabilityMatrix(o,j) = 0; %reset probability switch
                                  else
                                    neighbors(1,j) = old_neighbors(1,j);
                                  end
                              end
                        elseif length(potential_new)> num_changes
                                for j = 1:numberOfNeighbors
                                    if (change(o,j))
                                        if numel(potential_new) == 0 %debug!
                                            stop = 1;
                                        end
                                        neighbors(1,j) = potential_new(1);
                                        potential_new(potential_new == neighbors(1,j)) = [];
                                        probabilityMatrix(o,j) = 0; %reset probability switch
                                    else
                                        neighbors(1,j) = old_neighbors(1,j);
                                    end
                                end
                        else
                            error = 1;
                        end
                    end
                    adjacencyRow = zeros(1,numberOfNodes);
                    adjacencyRow(1,neighbors) = 1;
                end
             else 
                 adjacencyRow = ismember(distanceRow,minimumRow);
                 neighbors = find(adjacencyRow == 1);
               end

                 adjacencyMatrix(i,o,:) = adjacencyRow;
                 graphProperties(i).probabilityMatrix(o,:) = probabilityMatrix(o,:);
                 graphProperties(i).neighbors(o,:) = neighbors;
              
                if i==lb
                    neighborNumber(o) = neighborNumber(o)+1;
                    Invneighbors(o, neighborNumber(o),:)=neighbors;
                    fr_neigh(o, neighborNumber(o))=1;
                else
                if sum(neighbors==graphProperties(i-1).neighbors(o,:))~=numberOfNeighbors %this remembers the neighbors and how many
                    neighborNumber(o) = neighborNumber(o)+1;
                    Invneighbors(o, neighborNumber(o),:)=neighbors;
                    fr_neigh(o, neighborNumber(o)) = 1;
                else
                    fr_neigh(o, neighborNumber(o)) = fr_neigh(o, neighborNumber(o))+1;
                end
                end
                
                if adjacencyMatrix(i,o,o) ~= 0
                    disp(' ')
                    disp('Problem: edge connecting a node to itself')
                    disp('-----------------------------------------')
                    disp(['Frame: ' num2str(i)])
                    disp(['Node: ' num2str(o)])
                    viewMin
                    viewMax
                    theta  
                    relativeAngles(o,:)
                    if viewMin < viewMax
                        find((relativeAngles(o,:) <= viewMax) & (relativeAngles(o,:) >= viewMin))
                    else
                        find((relativeAngles(o,:) <= viewMax) | (relativeAngles(o,:) >= viewMin))
                    end
                    o
                    pause
              end
                
            Ndist(i,o,:) = distanceMatrix(o, find(adjacencyMatrix(i,o,:)==1)); %store in N distances to neighbors. i=frame, o=node, :=dist to neighbors
            Ndist = sort(Ndist,3);
            
        end
        
      graphProperties(i).UWadjacencyMatrix = squeeze(adjacencyMatrix(i,:,:)); %%unweighted!
      
       if ProbCheck_Status & i~=lb   
           clear neighbors
           clear change
       end
            
  end
end

disp('average neighbor number/node and avg neighbor number')
neighborNumber
mean(neighborNumber)

for s=1:numberOfNodes
    s
   %squeeze(Invneighbors(s,:,:)  %display neighbors of each node
   mean_fr_neigh(s) = sum(fr_neigh(s,:))/neighborNumber(s);
end

%fr_neigh %displaying for how many frames a neighbor was kept
disp('mean frame number for which a neighbor was kept')
mean_fr_neigh
disp('mean frame number for this sensing graph')
mean(mean_fr_neigh)

handles.adjacencyMatrix = adjacencyMatrix;

end

 if inputEdgeChoice==4 %circle-intersection model 
     
  chosen=1; 
  neighborNumber = zeros(1, numberOfNodes);
  
   lb = handles.lb;
   ub = handles.ub;
   i=lb;

 %read data (before Read_coordinatesAndCenter was implemented)
  while i<= ub
    
    file=strcat('Dance',num2str(i),'.dat');
    data_type = handles.data_type;
    if data_type==1
        data=handles.dancedata{i};
        set(handles.dataselection_editText,'String',i);
        x = data(:,1);
        y = data(:,2);
        theta = data(:,3);

        %[x y theta] = transf_in_familiar_coord(x,y,theta);        
        pos_x(i,:)=x;
        pos_y(i,:)=y;
        dir_theta(i,:)=theta;
        positionMatrix(i-lb+1,:,:) = horzcat(pos_x(i,:)',pos_y(i,:)');
    end
    i=i+1;
  end
  
  [adj1, adj2, adj3] = simple_synchrony(handles, hObject, lb, ub);
  
  numberOfNodes = length(squeeze(positionMatrix(1,:,:)));
  
  for o = 1:1
      %numberOfNodes
            %for now, assume 360 viewing angle!
            %choose closest neighbors
   o         
   i=lb;
  
   while i <= ub-tau2   
       
     distanceMatrix = ipdm(squeeze(positionMatrix(i-lb+1,:,:)));
     distanceMatrix(o,o)=Inf;
     [ignore, posInd] = sort(distanceMatrix(o,:));
    
     for k=1:length(posInd)-1
             pos_h(k) = adj2{i}(o,posInd(k));
     end
              
%look through the closest 6 neighbors and first pick those that have closest direction to you, then average-closest direction, then no matching direction.         
         k=1; kk=1; stage=2;
         top_dist=6;
         
while k<=numberOfNeighbors
    
   if pos_h(kk)==stage
      closest_neighbors(k)=posInd(kk);
      k=k+1;
   end
   
   if kk<top_dist
    kk=kk+1;
   else
       kk=1;
       stage=stage-1;
   end
end

%closest_neighbors are the vector from which the neighbor nodes will be chosen

%look within a window of time tau
            sum_cnt = 1;
            sum_tau = [];
            %determine pair of neighbors with smallest integrated distance to
            %intersection pt
            
            for k = 0:tau1-1
                
                j = i+k; 
                
             %integrate after frame i for tau frames
                
          for l=1:numberOfNeighbors
             pos(l,:) =[pos_x(j,closest_neighbors(l)),pos_y(j,closest_neighbors(l))];   %position of neighbors
          end
         
             pos_node = [pos_x(j,o), pos_y(j,o)];            %position of node in question 
             
            % exclude the neighbors that are far apart at the beginning -
            % more than 3R away.
             if k==0
                 ind_pos_neigh = [];
                for f = 1:numberOfNeighbors-1
                  for g = f+1:numberOfNeighbors
                      
                      if pdist([pos_x(j,closest_neighbors(f)),pos_y(j,closest_neighbors(f)); pos_x(j,closest_neighbors(g)),pos_y(j,closest_neighbors(g))]) < 3*inputRadius  %if these are too far apart, we do not consider them
                        if closest_neighbors(f)<closest_neighbors(g)  
                         ind_pos_neigh(sum_cnt,:,:) = [closest_neighbors(f), closest_neighbors(g); f,g];
                        else
                         ind_pos_neigh(sum_cnt,:,:) = [closest_neighbors(g), closest_neighbors(f); g,f];  
                        end
                        
                         sum_cnt = sum_cnt+1;
                      end
                  end
                end
           
                %list of pairs of neighbors
                 listOfNeighbors = squeeze(ind_pos_neigh(:,1,:));
             if size(listOfNeighbors,2) == 1
                 listOfNeighbors=listOfNeighbors';
             end
             
               sum_cnt = size(squeeze(listOfNeighbors(:,1,:)),1); % how many pairs of neighbors are eventually considered
               
             end
             
           %if sum_cnt ==0, so no neighbors close to each other -- not done
           %code for this, since algoritm always predicts some neighbors,
           %at least for the first video
           
          %integrate the distances
          for h = 1:sum_cnt
              %pair of nodes considered
                   f = ind_pos_neigh(h,2,1);
                   g = ind_pos_neigh(h,2,2);
                   
                   %intersection of two circles
                   [int1, int2] = circle_and_int([pos_x(j,ind_pos_neigh(h,1,1)), pos_y(j, ind_pos_neigh(h,1,1))], [pos_x(j,ind_pos_neigh(h,1,2)), pos_y(j,ind_pos_neigh(h,1,2))], inputRadius); %intersection of circles around the 2 neighbors
                   int1 = double(int1);
                   int2 = double(int2);
                   
                   if isnan(int1)==0
                 if imag(int1(1))~=0 %if the circles don't intersect
                      
                    d1 = sqrt((pos(f,1) - pos_node(1))^2+(pos(f,2) - pos_node(2))^2); %distnces to circles
                    d2 = sqrt((pos(g,1) - pos_node(1))^2+(pos(g,2) - pos_node(2))^2);
                    
                    if d1 > inputRadius %outside or inside the circle?
                        d1 = d1 - inputRadius;
                    else
                        d1 = inputRadius - d1;    
                    end
                    
                    if d2>inputRadius
                        d2 = d2 - inputRadius;
                    else
                        d2 = inputRadius - d2;    
                    end
                    
                   if k==0
                    sum_tau(h) = d1+d2;
                   else
                    sum_tau(h) = sum_tau(h) + d1 + d2;  %adding up the distances
                   end
                   
                 else %if the circles intersect
                      
                   pos_mx_pt1 = [pos_node; int1];
                   pos_mx_pt2 = [pos_node; int2];
                   
                   if k==0
                       if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                         sum_tau(h) = pdist(pos_mx_pt1);
                       else
                         sum_tau(h) = pdist(pos_mx_pt2);   
                       end
                   else
                       if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                         sum_tau(h) = sum_tau(h) + pdist(pos_mx_pt1);
                       else
                        sum_tau(h) = sum_tau(h) + pdist(pos_mx_pt2);   
                       end
                   end
                 end
                   else
                       
                       if k==0
                           if  pdist([pos_node; pos_x(j,ind_pos_neigh(h,1,1)), pos_y(j, ind_pos_neigh(h,1,1))])>inputRadius
                             sum_cnt(h) = pdist([pos_node; pos_x(j,ind_pos_neigh(h,1,1)), pos_y(j, ind_pos_neigh(h,1,1))])-inputRadius;
                           else
                             sum_cnt(h) = inputRadius-pdist([pos_node; pos_x(j,ind_pos_neigh(h,1,1)), pos_y(j, ind_pos_neigh(h,1,1))]);  
                           end
                       else
                           if  pdist([pos_node; pos_x(j,ind_pos_neigh(h,1,1)), pos_y(j, ind_pos_neigh(h,1,1))])>inputRadius
                             sum_cnt(h) = sum_cnt(h) + pdist([pos_node; pos_x(j,ind_pos_neigh(h,1,1)), pos_y(j, ind_pos_neigh(h,1,1))])-inputRadius;
                           else
                             sum_cnt(h) = sum_cnt(h) + inputRadius-pdist([pos_node; pos_x(j,ind_pos_neigh(h,1,1)), pos_y(j, ind_pos_neigh(h,1,1))]);  
                           end
                       end
                   end
          end
            end
         
            sum_tau2=sum_tau;
            %find the minimums
            [min_dist, min_dist_ind] = min(sum_tau);
            potential_neighbors = ind_pos_neigh(min_dist_ind,1,1:2);
            potential_neighbors = sort(squeeze(potential_neighbors)'); %pair of neighbors that has the minimal distance after integration
            
            if i==lb
                
                %a simple approach at frame lb: known neighbors and then
                %proceed one more frame
                neighbors = squeeze(potential_neighbors);
                graphProperties(i).neighbors(o,:) = neighbors;
    
                neighbor_data(o,1, cnt(o)) = i;
                neighbor_data(o,2, cnt(o)) = 1;
                neighbor_data(o,3, cnt(o)) = 0;
                
                neighbor_data1(o,cnt(o),:) = neighbors;
                i = i+1;
            
           else
          
            old_neighbors = sort(graphProperties(i-1).neighbors(o,:)); %old_neighbors: what neighbors were before
            potential_neighbors = squeeze(potential_neighbors)';
            
            if size(potential_neighbors,1)~= size(old_neighbors,1)
                potential_neighbors = potential_neighbors';
            end
             
            if sum(eq(old_neighbors, potential_neighbors)) == size(old_neighbors,2) %if old neighbors are the same as new predicted ones, then change to old neighbor
                
                neighbors = potential_neighbors;
               for b=i:i+tau1-1 
                 graphProperties(b).neighbors(o,:) = neighbors;
               end
               
                neighbor_data(o, 2, cnt(o)) = neighbor_data(o, 2, cnt(o))+tau1;
                i = i+tau1;  
                
            else
                
               for k = tau1:tau2 %if not, enlarge window and integrate again the distances; no more comments since exactly the same procedure for integration as before
                   
                j = i+k;
                
           pos = [];     
          for l=1:numberOfNeighbors
             pos(l,:) =[pos_x(j,closest_neighbors(l)),pos_y(j,closest_neighbors(l))];   %position of neighbors
          end
          
             pos_node = [pos_x(j,o), pos_y(j,o)];            %position of node in question
              
           %if sum_cnt ==0, so no neighbors close to each other -- have not
           %analyzed this case
           
          % for each frame k, consider each pair of neighbors   
          
          %integrate the distances
          for h = 1:sum_cnt
                   
              f = ind_pos_neigh(h,2,1);
              g = ind_pos_neigh(h,2,2);
              
                   [int1, int2] = circle_and_int([pos_x(j,ind_pos_neigh(h,1,1)), pos_y(j, ind_pos_neigh(h,1,1))],[pos_x(j,ind_pos_neigh(h,1,2)), pos_y(j, ind_pos_neigh(h,1,2))], inputRadius); %intersection of circles around the 2 neighbors
                   int1 = double(int1);
                   int2 = double(int2);
                   
                 if isnan(int1)==0  
                 if imag(int1(1))~=0
                      
                    d1 = sqrt((pos(f,1) - pos_node(1))^2+(pos(f,2) - pos_node(2))^2);
                    d2 = sqrt((pos(g,1) - pos_node(1))^2+(pos(g,2) - pos_node(2))^2);
                    
                    if d1>inputRadius
                        d1 = d1 - inputRadius;
                    else
                        d1 = inputRadius - d1;    
                    end
                    
                    if d2>inputRadius
                        d2 = d2 - inputRadius;
                    else
                        d2 = inputRadius - d2;    
                    end
                    
                    sum_tau(h) = sum_tau(h) + d1 + d2;  
                   
                  else
                      
                   pos_mx_pt1 = [pos_node; int1];
                   pos_mx_pt2 = [pos_node; int2];
                 
                   if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                    sum_tau(h) = sum_tau(h) + pdist(pos_mx_pt1);
                   else
                    sum_tau(h) = sum_tau(h) + pdist(pos_mx_pt2);   
                   end
                   
                 end
                 else
                     
                     if k==0
                           if  pdist([pos_node; pos_x(j,ind_pos_neigh(h,1,1)), pos_y(j, ind_pos_neigh(h,1,1))])>inputRadius
                             sum_cnt(h) = pdist([pos_node; pos_x(j,ind_pos_neigh(h,1,1)), pos_y(j, ind_pos_neigh(h,1,1))])-inputRadius;
                           else
                             sum_cnt(h) = inputRadius-pdist([pos_node; pos_x(j,ind_pos_neigh(h,1,1)), pos_y(j, ind_pos_neigh(h,1,1))]);  
                           end
                     else
                           if  pdist([pos_node; pos_x(j,ind_pos_neigh(h,1,1)), pos_y(j, ind_pos_neigh(h,1,1))])>inputRadius
                             sum_cnt(h) = sum_cnt(h) + pdist([pos_node; pos_x(j,ind_pos_neigh(h,1,1)), pos_y(j, ind_pos_neigh(h,1,1))])-inputRadius;
                           else
                             sum_cnt(h) = sum_cnt(h) + inputRadius-pdist([pos_node; pos_x(j,ind_pos_neigh(h,1,1)), pos_y(j, ind_pos_neigh(h,1,1))]);  
                           end
                     end
                       
                 end
          end
               end
             
               [min_dist, min_dist_ind] = min(sum_tau);
               potential_neighbors2 = ind_pos_neigh(min_dist_ind,1,1:2);
               potential_neighbors2 = sort(squeeze(potential_neighbors2));
             
               if size(potential_neighbors,1)~=size(potential_neighbors2,1)
                   potential_neighbors2 = potential_neighbors2';
               end
               
               %if neighbors after larger integration are the same as the neighbors after the first integration
               if sum(eq(potential_neighbors, potential_neighbors2)) == size(potential_neighbors2,2) 
                  %when tau1==tau2 => switch!
                   
                   neighbors = potential_neighbors2;
                  if sum(old_neighbors==0)==length(old_neighbors)
                   switch_pt = i;
                   
                  else
                  %determine the switch point  
                   switch_pt = determine_switch(handles, i, squeeze(neighbor_data(o,2,cnt(o))), old_neighbors, potential_neighbors, o, tau1, inputRadius);
                                      
                  end
                   
                  for b=switch_pt:i+tau1-1 
                     graphProperties(b).neighbors(o,:) = neighbors; %save this pair of neighbors in graphProperties
                  end
                   
                   cnt(o) = cnt(o)+1; %increase number of neighbors
                   
                   %neighbor_data is the matrix that remebers when the
                   %switches happen and whether the pair of neighbors was
                   %conclusive or inconclusive
                   neighbor_data(o,1,cnt(o))=switch_pt;
                   neighbor_data(o,2,cnt(o))=i+tau1 - switch_pt;
                   neighbor_data(o,3,cnt(o))=0;
                   
                   %neighbor_data1 remembers the neighbors; neighbor_data
                   %and neighbor_data1 will be stored in a file
                   neighbor_data1(o,cnt(o),:) = neighbors;
                   
                   %all cases must be treated, even when switch point is
                   %before previous neighbors (old_neighbors) came about.
                   if i>=switch_pt
                       neighbor_data(o,2,cnt(o)-1)=neighbor_data(o,2,cnt(o)-1)-(i-switch_pt);
                       
                       if neighbor_data(o,2,cnt(o)-1)<=0
                          
                           neighbor_data(o,:,cnt(o)-1) = neighbor_data(o,:, cnt(o));
                           neighbor_data(o,:, cnt(o)) = zeros(1,1,3);
                           cnt(o) = cnt(o)-1;
                           neighbor_data1(o,cnt(o),:) = neighbors;
                           neighbor_data1(o,cnt(o)+1,:) = [0, 0];
                          
                        if cnt(o)-1>=1   
                           if sum(eq(neighbor_data1(o,cnt(o),:), neighbor_data1(o,cnt(o)-1,:))) == length(neighbor_data1(o,cnt(o),:))
                                neighbor_data(o,2,cnt(o)-1) = neighbor_data(o,2,cnt(o)-1) + neighbor_data(o,2, cnt(o));
                                neighbor_data(o,:, cnt(o)) = zeros(1,1,3);
                                cnt(o) = cnt(o)-1;
                                neighbor_data1(o,cnt(o),:) = neighbors;
                                neighbor_data1(o,cnt(o)+1,:) = [0, 0];        
                           end
                        end
                        
                       end
                       
                   else
                       neighbor_data(o,2,cnt(o)-1)=neighbor_data(o,2,cnt(o)-1)+(switch_pt-i);
                       
                       for b=i:switch_pt-1
                           graphProperties(b).neighbors(o,:) = graphProperties(i-1).neighbors(o,:);
                       end
                       
                       if neighbor_data(o, 2,cnt(o))<=0
                           neighbor_data(o,:,cnt(o)) = zeros(1,3)';
                           cnt(o)=cnt(o)-1;                           
                       end
                   end
 
                  i = i+tau1;        %skip tau1 frames 
                  
               else  %when in the long run not the main pair (tau1 ~= tau2)
                   
                if sum(eq(potential_neighbors2, old_neighbors)) == size(old_neighbors,2) %... but old_neighbors are the same as postential_neighbors2
                       
                       ok=0;
                       %check is in the list of neighbors that are
                       %potential the old_neighbors
                       for v = 1: size(listOfNeighbors,1)
                           if sum(eq(old_neighbors, listOfNeighbors(v,:)))==2
                               ok=1;
                               ind_neighbor = v;
                           end
                       end

                       %if answer is yes to the ? posed in the previous
                       %comment and if neighbors sufficiently close, then
                       %neighbors are the same another frame
                       if ok==1 & sum_tau2(ind_neighbor)/tau1 < 0.4
                           
                          neighbors = old_neighbors;
                          graphProperties(i).neighbors(o,:) = neighbors;

                          neighbor_data(o,2,cnt(o)) = neighbor_data(o,2,cnt(o))+1; %another frame
                          i=i+1;
                          
                       else
                           
                          neighbors = [0 0]; %inconclusive, marked by '1' in neighbor_data
                          graphProperties(i).neighbors(o,:) = neighbors;
                         
                          if neighbor_data(o,3,cnt(o))==1
                           neighbor_data(o,2,cnt(o)) = neighbor_data(o,2,cnt(o))+1;
                          else
                              
                           cnt(o)=cnt(o)+1;
                           neighbor_data(o,1,cnt(o))=i; %where it begins
                           neighbor_data(o,2,cnt(o))=1; %how long it lasts
                           neighbor_data(o,3,cnt(o))=1; %conclusive/inconclusive
                           
                           neighbor_data1(o,cnt(o),:) = neighbors;
                           
                          end
                               
                          i=i+1; %when inconclusive, increment frame number by only 1

                       end
                     
                else
                       neighbors = [0 0]; % again inconclusive
                       graphProperties(i).neighbors(o,:) = neighbors;
                       
                       if neighbor_data(o,3,cnt(o))==1
                           neighbor_data(o,2,cnt(o)) = neighbor_data(o,2,cnt(o))+1;
                       else
                           cnt(o)=cnt(o)+1;
                           neighbor_data(o,1,cnt(o))=i; %where it begins
                           neighbor_data(o,2,cnt(o))=1; %how long it lasts
                           neighbor_data(o,3,cnt(o))=1; %conclusive/inconclusive
                           
                           neighbor_data1(o,cnt(o),:) = neighbors;
                       end
                     
                       i=i+1; %again increment by 1
                       
                 end
                   end
               end
              
            end
   end
   
  end
  
%go back... and modify wherever needed; if inconclusive for a small number
%of frames, then we can suppose neighbors are either the previous or the
%next

adjacencyMatrix = zeros(ub-tau2-lb+1, numberOfNodes, numberOfNodes);
neighborNumber = zeros(1, numberOfNodes);
neighborNumber2 = zeros(1, numberOfNodes);

jj = ones(1, numberOfNodes);
neighbor_data4 = zeros(numberOfNodes,1,2);

for o=1:numberOfNodes
    %display neighbor switches and neighbor pairs before modification
    o
    squeeze(neighbor_data(o,:,:))
    squeeze(neighbor_data1(o,:,:))
    
    final = find(neighbor_data(o,1,:)==0,1);
    final
    
    if size(final,1)==0
        final = size(squeeze(neighbor_data1(o,:,:)),1)
    else
        final = final-1;
    end
    
    %final represents how many neighbor pairs there are
    
    neighborNumber(o) = final;
    numberOfZeros(o) = length(find(squeeze(neighbor_data(o,3,1:final))==1));
    neighborNumber2(o) = final-numberOfZeros(o);
    mean_fr_neigh(o) = mean(neighbor_data(o,2,1:final));
    
    %neighbor_data3 and neighbor_data4 store the new neighbor switches and
    %the new neighbor pairs
    j=1;
    while j<=final
      
        if sum(eq(squeeze(neighbor_data1(o,j,:))',[0 0]))~=2  % when neighbors are not [0 0], then it's ok, nothing should be changed
            
        if jj(o)~=1 | j~=1
          if   sum(eq(neighbor_data1(o,j,:), neighbor_data4(o,jj(o),:)))~=length(squeeze(neighbor_data1(o,j,:)))
              jj(o)=jj(o)+1;
          end
        end
        
            neighbor_data3(o,:,jj(o))=neighbor_data(o,1:3,j);
            neighbor_data4(o,jj(o),:) = neighbor_data1(o,j,:);
            j=j+1;
            
            adjacencyMatrix(neighbor_data3(o,1,jj(o)):neighbor_data3(o,1,jj(o))+neighbor_data3(o,2,jj(o))-1,o,neighbor_data4(o,jj(o),1))=1;
            adjacencyMatrix(neighbor_data3(o,1,jj(o)):neighbor_data3(o,1,jj(o))+neighbor_data3(o,2,jj(o))-1,o,neighbor_data4(o,jj(o),2))=1;
            
        else
            
           if j<final %when neighbors are [0 0], but are not at the end
               
            if sum(eq(squeeze(neighbor_data1(o,j-1,:)), squeeze(neighbor_data1(o,j+1,:))))==2 & neighbor_data(o,2,j)<=tau2 %if neighbors before and after are the same
                          %then divide the inconclusive period into 2:
                          %a period when the neighbors are old_neighbors, and another period when the neighbors are the future ones 
       
                neighbor_data3(o,2,jj(o)) = neighbor_data3(o,2,jj(o)) + neighbor_data(o,2,j) + neighbor_data(o,2,j+1);
                
                j=j+2;
                adjacencyMatrix(neighbor_data3(o,1,jj(o)):neighbor_data3(o,1,jj(o))+neighbor_data3(o,2,jj(o))-1, o, neighbor_data4(o,jj(o),1))=1;
                adjacencyMatrix(neighbor_data3(o,1,jj(o)):neighbor_data3(o,1,jj(o))+neighbor_data3(o,2,jj(o))-1, o, neighbor_data4(o,jj(o),2))=1;
                
            elseif neighbor_data(o,2,j)<=tau2 % if neighbors before and after are different, but the time of inconclusiveness is short
                
                r1 = floor(neighbor_data(o,2,j)/2);
                r2 = floor(neighbor_data(o,2,j)/2);
                
                if r1*2 ~= neighbor_data(o,2,j)
                    r1=r1+1;
                end
                
                neighbor_data3(o,2,jj(o)) = neighbor_data3(o,2,jj(o)) + r1;
                jj(o) = jj(o)+1;
                neighbor_data3(o,1,jj(o)) = neighbor_data(o,1,j+1)-r2;
                neighbor_data3(o,2,jj(o)) = neighbor_data(o,2,j+1)+r2;
                neighbor_data4(o,jj(o),:) = neighbor_data1(o, j+1,:);
                j = j+2;
                adjacencyMatrix(neighbor_data3(o,1,jj(o)-1):neighbor_data3(o,1,jj(o)-1)+neighbor_data3(o,2,jj(o)-1)-1,o,neighbor_data4(o,jj(o)-1,1))=1;
                adjacencyMatrix(neighbor_data3(o,1,jj(o)-1):neighbor_data3(o,1,jj(o)-1)+neighbor_data3(o,2,jj(o)-1)-1,o,neighbor_data4(o,jj(o)-1,2))=1;
                
            else %else it remains inconclusive
                
               jj(o)=jj(o)+1;
               neighbor_data(o,1:3,j)
               neighbor_data3(o,1:3,jj(o))=neighbor_data(o,1:3,j);
               neighbor_data4(o,jj(o),:) = neighbor_data1(o,j,:);
               j=j+1;
               
            end
            
           else %last frame, so there are only neighbors before 
      
               if neighbor_data(o,2,j)<=tau2
                 neighbor_data3(o,2,jj(o)) = neighbor_data3(o,2,jj(o)) + neighbor_data(o,2,j);
                 adjacencyMatrix(neighbor_data3(o,1,jj(o)):neighbor_data3(o,1,jj(o))+neighbor_data3(o,2,jj(o))-1,o,neighbor_data4(o,jj(o),1))=1;
                 adjacencyMatrix(neighbor_data3(o,1,jj(o)):neighbor_data3(o,1,jj(o))+neighbor_data3(o,2,jj(o))-1,o,neighbor_data4(o,jj(o),2))=1;
               else
                   jj(o) = jj(o)+1;  
                   neighbor_data3(o,:,jj(o)) = neighbor_data(o,1:3,j);
                   neighbor_data4(o,jj(o),:) = neighbor_data1(o,j,:);
               end
               j=j+1;
           end
     
        end

             for i = neighbor_data3(o,1,jj(o)):neighbor_data3(o,1,jj(o)) + neighbor_data3(o,2,jj(o))-1
               distanceMatrix = ipdm(squeeze(positionMatrix(i-lb+1,:,:)));
               if isempty(find(squeeze(adjacencyMatrix(neighbor_data3(o,1,jj(o)),o,:))==1))==0
                  Ndist(i,o,:) = distanceMatrix(o, find(squeeze(adjacencyMatrix(neighbor_data3(o,1,jj(o)),o,:))==1)); %store in N distances to neighbors. i=frame, o=node, :=dist to neighbors 
               else
                  Ndist(i,o,:) = [NaN, NaN];
               end
               
             end
             Ndist = sort(Ndist,3);
    end
    
    mean_fr_neigh(o) = mean(neighbor_data3(o,2,1:jj(o)));
end

switch_nodes = jj;

%neighbor_data3
%neighbor_data4

%{
disp('number of neighbors of each:')
neighborNumber
disp('mean number of neighbors for this sensing graph')
mean(neighborNumber)
disp('number of neighbors without [0,0]s')
neighborNumber2
disp('mean number of neighbors for this sensing graph without [0 0]s')
mean(neighborNumber2)

disp('mean frame number for every node')
mean_fr_neigh
disp('mean frame number for sensing graph')
mean(mean_fr_neigh)
%}
 end
%saving all this information, of neighbor pair switching and pair of neighbors

%save('senseGraph_r03_N4.mat', 'neighbor_data', 'neighbor_data1', 'neighbor_data3', 'neighbor_data4', 'adjacencyMatrix');
%save('NewModel2_r03_N4_Vreala.mat','adjacencyMatrix', 'neighbor_data', 'neighbor_data1', 'neighbor_data3', 'neighbor_data4')
%save('adjacencyMatrix_NewModel_180_0.3_video2.mat', 'adjacencyMatrix')
%save('adjacencyMatrix_NewModel_180_0.4_video2.mat','adjacencyMatrix')
%save('adjacencyMatrix_NewModel_180_0.5_video2.mat','adjacencyMatrix')
%save('180_r_infty_video2.mat','adjacencyMatrix')
%save('adjacencyMatrix_360_video2s.mat','adjacencyMatrix')

 if inputEdgeChoice>3
     ub = ub-tau2;
 end
 
nodeStatus = zeros(ub-lb+1,numberOfNodes);
handles.adjacencyMatrix = adjacencyMatrix; 

threshold=0.5;
cnt_ns = zeros(1, numberOfNodes);

%evaluate properties of sensing graph
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
        graphProperties(i).status = nodeStatus(i,:);
    end
    
    graphProperties(i).Theta = theta;
    graphProperties(i).ViewAngle = viewAngle;
    if (ProbCheck_Status)
    graphProperties(i).ProbSlope = m;
    end
    graphProperties(i).PositionMatrix = positionMatrix;
    clearvars -except *weight avi_file_name inputEdgeChoice inputGraphChoice inputWeightChoice adjacencyMatrix lb ub tau2 time_len numberOfNeighbors numberOfNodes Ndist tau switch_nodes node_status_cum inputRadius viewAngle numberPictures inputComputeChoice aviobj subject xmin xmax ymin ymax zmin zmax handles graphProperties lengthScale ProbCheck_Status m nodeStatus hObject nodeStatusMx cnt_ns threshold
    if inputComputeChoice == 1 || inputComputeChoice == 2  
       frame = getframe(gcf);
      aviobj = addframe(aviobj,frame);
    end
    cla;
end
hold off;
'here2'
if inputComputeChoice == 1 || inputComputeChoice == 2
  aviobj = close(aviobj);
end

for j=1:numberOfNodes
    disp(j);
    for k=1:2
      mean_dist(j,k) = mean(Ndist(:,j,k));
    end  
      %mean(m(j,:))
      disp(mean_dist(j,:));
end

for j=1:numberOfNodes
    disp(j);
    mean_dist_avg(j) = mean(mean_dist(j,:));
    disp(mean_dist_avg(j));
end

[mean_dist_sorted, ind_avg] = sort(mean_dist_avg);

ind_avg
mean_dist_sorted

%how often ppl switch
switch_nodes
[switch_nodes_sorted, ind_switch] = sort(switch_nodes);
switch_nodes_sorted
ind_switch

%}

disp('mean node status for nodes')
node_status_cum = node_status_cum/(handles.ub-handles.lb+1)

figure(2)
plot(1:numberOfNodes, node_status_cum, 'x', 'Color', 'b')

disp('periods with higher node statuses')
nodeStatusMx
cnt_ns

disp('mean node status overall')
mean(node_status_cum)
[node_status_cum_sorted, ind_node_stat] = sort(node_status_cum);
node_status_cum_sorted
ind_node_stat

for i=lb:ub
Wconnect(i) = graphProperties(i).weaklyConnected;
end

states = 1;

Wconnect_mx(1,1)=lb;
Wconnect_mx(2,1)=1;
Wconnect_mx(3,1)=Wconnect(i);

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

%{
figure(1)
xmin = lb; xmax = ub;
ymin = 0; ymax = 3;
axis([xmin, xmax, ymin, ymax]);
plot(lb:ub, graphProperties(lb:ub).stronglyConnected);
%}

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
%save(strcat(avi_file_name,'.mat'),'graphProperties', 'adjacencyMatrix');
%}
guidata(hObject,handles)

end

%without creating a video, computing properties that are shown in Flock
%Grapher: H2 norm, connectedness (strong/weak), L2 gain and so on.

function compute_properties_withoutVideo(handles, hObject)

adjacencyMatrix = handles.adjacencyMatrix;
numberOfNodes = handles.numberOfNodes;
min_frame = handles.min_frame;

inputWeightChoice = (get(handles.edges_popup,'Value'));
inputComputeChoice = 3;

compute_properties(adjacencyMatrix, numberOfNodes, min_frame, inputComputeChoice, inputWeightChoice, handles, hObject)
end

%from an adjacency matrix with undecided neighbors (denoted by [0 0])
%assign neighbors so that half of the frames when neighbors are undecided
%are assigned to the neighbors before the undecided period, and the other
%half to the neighbors after.

function [adjacencyMatrix min_frame] = build_adjMatrix(handles, hObject)

    file = ['NewModel2_r03_N4_Vreala_v2', '.mat'];
    load(file);
    
    numberOfNodes = size(neighbor_data5,1);
    min_frame =Inf;
    
    for o=1:numberOfNodes
        final2 = find(neighbor_data5(o,1,:)==0,1);
    
        if size(final2,1)==0
          final(o) = length(squeeze(neighbor_data5(o,1,:)));
        else
          final(o) = final2-1;
        end
    
        if min_frame > neighbor_data5(o,1,final(o))+neighbor_data5(o,2,final(o))-1;
            min_frame = neighbor_data5(o,1,final(o))+neighbor_data5(o,2,final(o))-1;
        end
        
    end
 
    disp('min_frame')
    min_frame
    adjacencyMatrix = zeros(min_frame, numberOfNodes, numberOfNodes);
    
    for o=1:numberOfNodes
        for i=1:final(o)
            adjacencyMatrix(neighbor_data5(o,1,i):neighbor_data5(o,1,i)+neighbor_data5(o,2,i)-1, o, neighbor_data6(o,i, 1)) = 1;
            adjacencyMatrix(neighbor_data5(o,1,i):neighbor_data5(o,1,i)+neighbor_data5(o,2,i)-1, o, neighbor_data6(o,i, 2)) = 1;
        end
    end
    
    handles.adjacencyMatrix = adjacencyMatrix;
    handles.min_frame = min_frame;
    save('NewModel2_r03_N4_Vreala_v2', 'adjacencyMatrix', '-append'); 
    guidata(hObject, handles);
 
end

%integrate appears in the implementation of circle-intersection model as
%well
function [sum, winner] = integrate(handles, start_frame, duration, pairOfNodes, node, inputRadius)

  for i = start_frame : start_frame + duration-1
%read data
      file=strcat('Dance',num2str(i),'.dat');
      data_type = handles.data_type;
      if data_type==1
         data=handles.dancedata{i};
         set(handles.dataselection_editText,'String',i);
         pos_x = data(:,1);
         pos_y = data(:,2);
         theta = data(:,3);
         %[pos_x pos_y] = transf_in_familiar_coord(pos_x,pos_y,theta); 
      end
      
      pos_node = [pos_x(node), pos_y(node)];
      for j = 1:size(pairOfNodes,1)
          pos_pairOfNodes(j,1,:) = [pos_x(pairOfNodes(j,1)), pos_y(pairOfNodes(j,1))];
          pos_pairOfNodes(j,2,:) = [pos_x(pairOfNodes(j,2)), pos_y(pairOfNodes(j,2))];
          
          [int1, int2] = circle_and_int(pos_pairOfNodes(j,1,:), pos_pairOfNodes(j,2,:), inputRadius); %intersection of circles around the 2 neighbors
          int1 = double(int1);
          int2 = double(int2);
                   
           if imag(int1(1))~=0 %if the circles do not intersect
                      
              d1 = sqrt((pos_pairOfNodes(j,1,1) - pos_node(1))^2+(pos_pairOfNodes(j,1,2) - pos_node(2))^2);
              d2 = sqrt((pos_pairOfNodes(j,2,1) - pos_node(1))^2+(pos_pairOfNodes(j,2,2) - pos_node(2))^2);
                    %if inside or outside the circle
                    if d1>inputRadius
                        d1 = d1 - inputRadius;
                    else
                        d1 = inputRadius - d1;    
                    end
                    
                    if d2>inputRadius
                        d2 = d2 - inputRadius;
                    else
                        d2 = inputRadius - d2;    
                    end
               %integrate the sum     
                if i==start_frame
                    sum_cum(j) = d1+d2;
                else
                    sum_cum(j) = sum_cum(j) + d1 + d2;  
                end 
                
           else  % if intersection exists
                      
                   pos_mx_pt1 = [pos_node; int1];
                   pos_mx_pt2 = [pos_node; int2];
                   
                if i==start_frame
                    if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                    sum_cum(j) = pdist(pos_mx_pt1);
                   else
                    sum_cum(j) = pdist(pos_mx_pt2);   
                   end
                else
                   if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                    sum_cum(j) = sum_cum(j) + pdist(pos_mx_pt1);
                   else
                    sum_cum(j) = sum_cum(j) + pdist(pos_mx_pt2);   
                   end
                end   
            end
        
      end 
      
      %the pair of nodes chosen as neighbors is the one for which the
      %distance to the intersection point is minimal
      [sum, ind] = min(sum_cum);
      winner = pairOfNodes(ind,:);
      
  end
end

%since many times, I increment i with tau1 (which was set to 15 frames), switching pairs of neighbors would lack precision, 
%and would not give the exact frame at which neighbors are switched. This
%algorithm is rather complicated, and could most likely be simplified.

%there were also a few mistakes, corrected after the algorithm was
%implemented with the results (these mistakes however are corrected through commented code; the original code was left as it is). 
% This affects ns SG, but intuition is that
%the algorithm change would not modify results VERY much since all that
%determine_switch does is determine the exact frame of the switch, which is 
%around frame pt.

%it seems like a long function, but actually there is a lot of repetition,
%for which separate function could have been considered

function [switch_pt] = determine_switch(handles, pt, duration, neighbors1, neighbors2, node, tau1, inputRadius)
      
        %look for neighbors within middle time frame and see where to look
        %first; search only among neighbors1,2
        
        sum_tau =[0 0];
        duration2 = min(duration, floor(tau1/2)); %the interval on which we consider a switch could happen, always <tau1/2
        
        %pt is the frame we're at
           for i = pt-duration2:pt+floor(tau1/2) %looking before and after pt
        %at first looking in the whole interval and establishing to which pair of nodes/neighbors the node is closest to               
                    data=handles.dancedata{i}; %read data
                    pos_x = data(:,1);
                    pos_y = data(:,2);
                    theta = data(:,3);
                    %[pos_x pos_y] = transf_in_familiar_coord(pos_x,pos_y,theta); 
      
                pos_node = [pos_x(node), pos_y(node)];
                
                f = neighbors1(1);
                g = neighbors1(2);
                
                pos1(1,:) = [pos_x(f), pos_y(f)];  %position of nodes that are neighbors before frame pt
                pos1(2,:) = [pos_x(g) pos_y(g)];
                
                f = neighbors2(1);
                g = neighbors2(2);
                 
                pos2(1,:) = [pos_x(f), pos_y(f)];   %position of nodes that are neighbors after frame pt
                pos2(2,:) = [pos_x(g), pos_y(g)];   
                
                %for first pair of neighbors, from before
                   [int1, int2] = circle_and_int(pos1(1,:), pos1(2,:), inputRadius); %intersection of circles around the 2 neighbors (from before pt)
                   int1 = double(int1);
                   int2 = double(int2);
                   
                 if imag(int1(1))~=0 %when circles don't intersect
                      
                    d1 = sqrt((pos1(1,1) - pos_node(1))^2+(pos1(1,2) - pos_node(2))^2);
                    d2 = sqrt((pos1(2,1) - pos_node(1))^2+(pos1(2,2) - pos_node(2))^2);
                    
                    if d1 > inputRadius
                        d1 = d1 - inputRadius;
                    else
                        d1 = inputRadius - d1;    
                    end
                    
                    if d2>inputRadius
                        d2 = d2 - inputRadius;
                    else
                        d2 = inputRadius - d2;    
                    end
                    
                   if i==pt-duration2  %integration over time frame [pt-duration2, pt+duration2]
                    sum_tau(1) = d1+d2;
                   else
                    sum_tau(1) = sum_tau(1) + d1 + d2;  
                   end
                   
                 else  %when circles intersect
                      
                   pos_mx_pt1 = [pos_node; int1];
                   pos_mx_pt2 = [pos_node; int2];
                   
                   if i==pt-duration2
                       if pdist(pos_mx_pt1) < pdist(pos_mx_pt2) %integrating with distance from node to intersection which is closest
                         sum_tau(1) = pdist(pos_mx_pt1);
                       else
                         sum_tau(1) = pdist(pos_mx_pt2);   
                       end
                   else
                       if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                         sum_tau(1) = sum_tau(1) + pdist(pos_mx_pt1);
                       else
                        sum_tau(1) = sum_tau(1) + pdist(pos_mx_pt2);   
                       end
                   end
                 end
                
                 %for second pair of neighbors, from after
                   [int1, int2] = circle_and_int(pos2(1,:), pos2(2,:), inputRadius); %intersection of circles around the 2 neighbors (from after frame pt)
                   int1 = double(int1);
                   int2 = double(int2);
                   
                 if imag(int1(1))~=0 %same thing: if circles don't intersect
                      
                    d1 = sqrt((pos2(1,1) - pos_node(1))^2+(pos2(1,2) - pos_node(2))^2);
                    d2 = sqrt((pos2(2,1) - pos_node(1))^2+(pos2(2,2) - pos_node(2))^2);
                    
                    if d1 > inputRadius
                        d1 = d1 - inputRadius;
                    else
                        d1 = inputRadius - d1;    
                    end
                    
                    if d2>inputRadius
                        d2 = d2 - inputRadius;
                    else
                        d2 = inputRadius - d2;    
                    end
                    
                   if i==pt-duration2
                    sum_tau(2) = d1+d2;
                   else
                    sum_tau(2) = sum_tau(2) + d1 + d2;  
                   end
                   
                  else
                      
                   pos_mx_pt1 = [pos_node; int1];
                   pos_mx_pt2 = [pos_node; int2];
                   
                   if i==pt-duration2
                       
                       if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                         sum_tau(2) = pdist(pos_mx_pt1);
                       else
                         sum_tau(2) = pdist(pos_mx_pt2);   
                       end
                       
                   else
                       if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                         sum_tau(2) = sum_tau(2) + pdist(pos_mx_pt1);
                       else
                        sum_tau(2) = sum_tau(2) + pdist(pos_mx_pt2);   
                       end
                   end
                 end 
             
           end
          
           %determining which pair of nodes is closer
           if sum_tau(1)>=sum_tau(2)
               look = -1;  %pair of neighbors from after are closer
           else
               look = 1;  %pair of neighbors from before are closer
           end
          
           %switch can't be lower than pt_end1 and higher than pt_end2
           pt_end1 = pt-min(tau1, duration);
           pt_end2 = pt+tau1;

           if look == -1 %since pair of neighbors from after frame pt are closer, switch should be before pt, so we look in the interval [pt-duration2, pt]
               
              pt2 = pt;  %ends of the first interval
              pt1 = pt - duration2;
              
              End = 0;
              iteration=1; change=0;
              
              sum_tau = [0 0]; %sum of integration for both pairs of neighbors, at first initialized to 0
              
              while pt2-pt1 >=1 & End == 0 % look for switch while the ends of the interval are different; I am not looking for the switch at the frames where neighbors from before/after are (before neighbors1, after neighbors2)
                
                  for i = pt1:pt2  %integrating distances between pt1 and pt2
                     
                        data=handles.dancedata{i};
                        pos_x = data(:,1);
                        pos_y = data(:,2);
                        theta = data(:,3);
                        %[pos_x pos_y] = transf_in_familiar_coord(pos_x,pos_y,theta); 
      
                pos_node = [pos_x(node), pos_y(node)];
                
                f = neighbors1(1);
                g = neighbors1(2);
                
                pos1(1,:) = [pos_x(f), pos_y(f)];
                pos1(2,:) = [pos_x(g) pos_y(g)];
                
                f = neighbors2(1);
                g = neighbors2(2);
                
                pos2(1,:) = [pos_x(f), pos_y(f)];
                pos2(2,:) = [pos_x(g), pos_y(g)];   
                
                   [int1, int2] = circle_and_int(pos1(1,:), pos1(2,:), inputRadius); %intersection of circles around the 2 neighbors
                   int1 = double(int1);
                   int2 = double(int2);
                   
                 if imag(int1(1))~=0
                      
                    d1 = sqrt((pos1(1,1) - pos_node(1))^2+(pos1(1,2) - pos_node(2))^2);
                    d2 = sqrt((pos1(2,1) - pos_node(1))^2+(pos1(2,2) - pos_node(2))^2);
                    
                    if d1 > inputRadius
                        d1 = d1 - inputRadius;
                    else
                        d1 = inputRadius - d1;    
                    end
                    
                    if d2>inputRadius
                        d2 = d2 - inputRadius;
                    else
                        d2 = inputRadius - d2;    
                    end
                    
                    %integrating for first neighbor pair, neighbors1, in
                    %sum_tau(1)
                   if i==pt1
                    sum_tau(1) = d1+d2;
                   else
                    sum_tau(1) = sum_tau(1) + d1 + d2;  
                   end
                   
                  else
                      
                   pos_mx_pt1 = [pos_node; int1];
                   pos_mx_pt2 = [pos_node; int2];
                   
                   if i==pt1
                       if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                         sum_tau(1) = pdist(pos_mx_pt1);
                       else
                         sum_tau(1) = pdist(pos_mx_pt2);   
                       end
                       
                   else
                       if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                         sum_tau(1) = sum_tau(1) + pdist(pos_mx_pt1);
                       else
                        sum_tau(1) = sum_tau(1) + pdist(pos_mx_pt2);   
                       end
                   end
                 end
                 
                   [int1, int2] = circle_and_int(pos2(1,:), pos2(2,:), inputRadius); %intersection of circles around the 2 neighbors
                   int1 = double(int1);
                   int2 = double(int2);
                   
                 if imag(int1(1))~=0
                      
                    d1 = sqrt((pos2(1,1) - pos_node(1))^2+(pos2(1,2) - pos_node(2))^2);
                    d2 = sqrt((pos2(2,1) - pos_node(1))^2+(pos2(2,2) - pos_node(2))^2);
                    
                    if d1 > inputRadius
                        d1 = d1 - inputRadius;
                    else
                        d1 = inputRadius - d1;    
                    end
                    
                    if d2>inputRadius
                        d2 = d2 - inputRadius;
                    else
                        d2 = inputRadius - d2;    
                    end
                    
                    %integrating for second neighbor pair, neighbors2, in
                    %sum_tau(2)
                   if i==pt1
                    sum_tau(2) = d1+d2;
                   else
                    sum_tau(2) = sum_tau(2) + d1 + d2;  
                   end
                   
                  else
                      
                   pos_mx_pt1 = [pos_node; int1];
                   pos_mx_pt2 = [pos_node; int2];
                   
                   if i==pt1
                       if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                         sum_tau(2) = pdist(pos_mx_pt1);
                       else
                         sum_tau(2) = pdist(pos_mx_pt2);   
                       end
                   else
                       if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                         sum_tau(2) = sum_tau(2) + pdist(pos_mx_pt1);
                       else
                        sum_tau(2) = sum_tau(2) + pdist(pos_mx_pt2);   
                       end
                   end
                 end
                  end  %end of integration from pt1 to pt2 of distances to pairs of neighbors, neighbor1 and neighbor2
                  
                  if iteration == 1  %first time
                      
                    if sum_tau(1) >= sum_tau(2) %neighbors2 are closer
                        status = 2; %status=1 when neighbors1 are closer; status=2 when neighbors2 are closer
                        m = floor((pt2-pt1)/2);
                        pt2 = pt1+m;  %moving more towards pt1, since neighbors2 are closer, makes sense
                        pt1 = pt1-m;
                    else 
                        status = 1;
                        m = floor((pt2-pt1)/2);
                        pt2 = pt2+m;
                        pt1 = pt2-m;
                    end
                    
                    if pt1 <= pt_end1 %if reached the limit for pt1, or pt2
                       pt1 = pt_end1; 
                    end
                              
                    if pt2 >=pt_end2
                       pt2 = pt_end2;
                    end
                         
                    iteration=iteration+1;
                    
                  else
                      
                      if sum_tau(1) >= sum_tau(2) %neighbors2 are closer
                          switch status
                              case 1 %if this has changed, change=1; if still the same, change=0.
                                  change = 1;
                                  status = 2;
                                  pt1 = pt2 - floor((pt2-pt1)/2);  %coming closer to pt2, since neighbors2 are closer, looking for the switch
                                  %pt2 = pt1 + floor((pt2-pt1)/2);                             
      %coming closer to pt1, since neighbors2 are closer, looking for the switch more towards pt1... this is correct, but was not implemented
      
                              case 2
                                  change=0;  %if the smae neighbors2 are closer move pt1 and pt2 more to the left
                                  m = floor((pt2-pt1)/2);
                                  pt1 = pt1-m;
                                  pt2 = pt1+m;
                           end
                          
                          if (pt1 <= pt_end1 | pt2 >= pt_end2) & change == 0 %if there is no change and have reached the end of the interval, stop loop
                              End = 1;
                          end
                          
                           if pt1 <= pt_end1
                                 pt1 = pt_end1; 
                           end
                              
                           if pt2 >=pt_end2
                                  pt2 = pt_end2;
                           end
                          
                      else %when neighbors1 are closer
                          
                         switch status
                             case 1
                               change = 0;
                               m = floor((pt2-pt1)/2);
                               pt1 = pt2-m;
                               pt2 = pt2+m;
                             case 2
                               change = 1;
                               status = 1;
                               pt2 = pt1 + floor((pt2-pt1)/2);  
                               %pt1 = pt2 - floor((pt2-pt1)/2);  %when
                               %neighbors1 are closer, look more in the
                               %right, towards pt2
                               %correct version, but not the one
                               %implemented
                           end
                         
                         if (pt1 <= pt_end1 | pt2 >= pt_end2) & change == 0
                              End = 1;
                         end
                         
                             if pt1 <= pt_end1
                                 pt1 = pt_end1; 
                             end
                              
                             if pt2 >=pt_end2
                                 pt2 = pt_end2;
                             end
                         
                      end
                  
                  end
                 
              
              end
              
           else  %if first pair of neighbors is closer, look in interval [pt, pt+duration2] for a switch
               %same thing as before
              pt1 = pt;
              pt2 = pt + floor(tau1/2);
              
              End = 0;
              iteration=1; change=0;
              
              sum_tau = [0 0];
              
              while pt2-pt1 >=1 & End == 0
               
                  for i = pt1:pt2
                     
                        data=handles.dancedata{i};
                        set(handles.dataselection_editText,'String',i);
                        pos_x = data(:,1);
                        pos_y = data(:,2);
                        theta = data(:,3);
                        %[pos_x pos_y] = transf_in_familiar_coord(pos_x,pos_y,theta); 
      
                pos_node = [pos_x(node), pos_y(node)];
                
                f = neighbors1(1);
                g = neighbors1(2);
                
                pos1(1,:) = [pos_x(f), pos_y(f)];
                pos1(2,:) = [pos_x(g) pos_y(g)];
                
                f = neighbors2(1);
                g = neighbors2(2);
                
                pos2(1,:) = [pos_x(f), pos_y(f)];
                pos2(2,:) = [pos_x(g), pos_y(g)];   
                
                   [int1, int2] = circle_and_int(pos1(1,:), pos1(2,:), inputRadius); %intersection of circles around the 2 neighbors
                   int1 = double(int1);
                   int2 = double(int2);
                   
                 if imag(int1(1))~=0
                      
                    d1 = sqrt((pos1(1,1) - pos_node(1))^2+(pos1(1,2) - pos_node(2))^2);
                    d2 = sqrt((pos1(2,1) - pos_node(1))^2+(pos1(2,2) - pos_node(2))^2);
                    
                    if d1 > inputRadius
                        d1 = d1 - inputRadius;
                    else
                        d1 = inputRadius - d1;    
                    end
                    
                    if d2>inputRadius
                        d2 = d2 - inputRadius;
                    else
                        d2 = inputRadius - d2;    
                    end
                    
                   if i==pt-duration
                    sum_tau(1) = d1+d2;
                   else
                    sum_tau(1) = sum_tau(1) + d1 + d2;  
                   end
                   
                  else
                      
                   pos_mx_pt1 = [pos_node; int1];
                   pos_mx_pt2 = [pos_node; int2];
                   
                   if i==pt-duration
                       if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                         sum_tau(1) = pdist(pos_mx_pt1);
                       else
                         sum_tau(1) = pdist(pos_mx_pt2);   
                       end
                   else
                       if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                         sum_tau(1) = sum_tau(1) + pdist(pos_mx_pt1);
                       else
                        sum_tau(1) = sum_tau(1) + pdist(pos_mx_pt2);   
                       end
                   end
                 end
               
                   [int1, int2] = circle_and_int(pos2(1,:), pos2(2,:), inputRadius); %intersection of circles around the 2 neighbors
                   int1 = double(int1);
                   int2 = double(int2);
                   
                 if imag(int1(1))~=0
                      
                    d1 = sqrt((pos2(1,1) - pos_node(1))^2+(pos2(1,2) - pos_node(2))^2);
                    d2 = sqrt((pos2(2,1) - pos_node(1))^2+(pos2(2,2) - pos_node(2))^2);
                    
                    if d1 > inputRadius
                        d1 = d1 - inputRadius;
                    else
                        d1 = inputRadius - d1;    
                    end
                    
                    if d2>inputRadius
                        d2 = d2 - inputRadius;
                    else
                        d2 = inputRadius - d2;    
                    end
                    
                   if i==pt-floor(tau1/2)
                    sum_tau(2) = d1+d2;
                   else
                    sum_tau(2) = sum_tau(2) + d1 + d2;  
                   end
                   
                 else
                      
                   pos_mx_pt1 = [pos_node; int1];
                   pos_mx_pt2 = [pos_node; int2];
                   
                   if i==pt1
                       if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                         sum_tau(2) = pdist(pos_mx_pt1);
                       else
                         sum_tau(2) = pdist(pos_mx_pt2);   
                       end
                   else
                       if pdist(pos_mx_pt1) < pdist(pos_mx_pt2)
                         sum_tau(2) = sum_tau(2) + pdist(pos_mx_pt1);
                       else
                        sum_tau(2) = sum_tau(2) + pdist(pos_mx_pt2);   
                       end
                   end
                 end
                  end
          
                  %after evaluating the distances to the pairs of
                  %neighbors, compare them and move pt1 and pt2 accordingly
                  
                  if iteration == 1
                      
                    if sum_tau(1) >= sum_tau(2) 
                        status = 2;
                        m = floor((pt2-pt1)/2);
                        pt2 = pt1+m;  %if closer to neighbors2, then move towards pt1, to the left
                        pt1 = pt1-m;
                    else 
                        status = 1;
                        m = floor((pt2-pt1)/2);
                        pt2 = pt2+m;  %if closer to neighbors1, then move towards pt1, to the right
                        pt1 = pt2-m;
                    end
                    
                         if pt1 <= pt_end1
                                 pt1 = pt_end1; 
                         end
                              
                        if pt2 >=pt_end2
                                 pt2 = pt_end2;
                        end
                         
                    iteration = iteration+1;
                  else
                      
                      if sum_tau(1) >= sum_tau(2) 
                          switch status
                              case 1
                                  change = 1;
                                  status = 2;
                                  pt1 = pt2 - floor((pt2-pt1)/2);
                                  % pt2 = pt1 + floor((pt2-pt1)/2);
                                  % correct version, not implemented
                              case 2
                                  change=0;
                                  m = floor((pt2-pt1)/2);
                                  pt1 = pt1-m;
                                  pt2 = pt1+m;
                                  
                           end
                          
                          if (pt1 <= pt_end1 | pt2 >= pt_end2) & change == 0
                              End = 1;
                          end
                          
                              if pt1 <= pt_end1
                                 pt1 = pt_end1; 
                              end
                              
                             if pt2 >=pt_end2
                                 pt2 = pt_end2;
                             end
                          
                      else
                          
                         switch status
                             case 1
                               change = 0;
                               m = floor((pt2-pt1)/2);
                               pt1 = pt2-m;
                               pt2 = pt2+m;
                             case 2
                               change = 1;
                               status = 1;
                               pt2 = pt1 + floor((pt2-pt1)/2);
                               %pt1 = pt2 - floor((pt2-pt1)/2);  
                               %correct version
                           end
                         
                         if (pt1 <= pt_end1 | pt2 >= pt_end2) & change == 0
                              End = 1;
                         end
                         
                              if pt1 <= pt_end1
                                 pt1 = pt_end1; 
                              end
                              
                             if pt2 >=pt_end2
                                 pt2 = pt_end2;
                             end
                         
                  end
                  end
   end
 end
  
     if pt1 == pt_end1 %make all neighbor2
       switch_pt = pt_end1;
     else
      if pt2 == pt_end2 %make all neighbor1
       switch_pt = pt_end2;
      else   %else there is a switch point
        switch_pt = floor((pt1+pt2)/2); 
      end
     end
     
end

% read coordinates of nodes and computes the mean of these coordinates which
% represents the center of the group

%read coordinates and store them in positionMatrix; also look at
%coordinates of the center
function [positionMatrix, center, numberOfNodes] = Read_coordinatesAndCenter(start_frame, end_frame, handles, hObject)
    
lb = handles.lb;
ub = handles.ub;

if start_frame<lb | end_frame>ub
    error('Provide different start and end frames');
    return 
end

    for i=start_frame:end_frame
      
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
        positionMatrix(i-start_frame+1,:,:) = horzcat(x,y);
        center(i-start_frame+1,:) = [mean(x), mean(y)];
        
        if i==start_frame
            numberOfNodes = length(x);
            handles.numberOfNodes = numberOfNodes;
        end
        
      end
      
    end
    
end

function choose_graph_pushbutton_Callback(hObject, ~, handles)
choose_graph(hObject, handles,0);
end

function A = choose_graph(hObject, handles, graph)

%in case graph is not specified at input ("graph"), choose sensing graph
if graph ==0
    
prompt = {'Choose graph to use:                                                         Select 1). 360 degrees modes                                                                                        Select 2). 180 degrees, r=0.3                                                                                                                               Select 3). 180 degrees, r=0.4                                                                                                                           Select 4). 180 degrees, r=0.5                                                                                           Select 5). 180 degrees, R=infty                                                                                                                     Select 6). New Model'};
graph = str2num(cell2mat(inputdlg(prompt)));

switch graph
    case 1
        load('360.mat');
    case 2
        load('180_r03.mat');
    case 3
        load('180_r04.mat');
    case 4
        load('180_r05.mat');
    case 5
        load('180_r_infty.mat');
    case 6
        load('circleInt_N4_r03.mat');
end

%when what sensing graph to use is specified
else
    
        %1- 360 degress model; 2- 180 degrees model, r=0.3; 3- 180 degrees
        %model, r=0.4; 4- 180 degrees model, r=0.4; 5- 180 degrees model, r=0.infty;
        %6 - circle intersection model
    
    switch graph
    case 1
        load('360.mat');
    case 2
        load('180_r03.mat');
    case 3
        load('180_r04.mat');
    case 4
        load('180_r05.mat');
    case 5
        load('180_r_infty.mat');
    case 6
        load('circleInt_N4_r03.mat');
    end

end

handles.adjacencyMatrix = adjacencyMatrix;
handles.numberOfNodes = length(squeeze(adjacencyMatrix(1,:,:)));
A = adjacencyMatrix;

guidata(hObject,handles)
end

%to compare different sensing graphs and node statuses of these sensing
%graphs
function compare_graph_properties(hObject, handles, start_frame, end_frame)

lb = handles.lb;
ub = handles.ub;

%further properties of the graph, besides node status, can be compared
prompt = {'Select what property you want to compare:                                                                                        Select 1). Node status'};
property =str2num(cell2mat(inputdlg(prompt)))

S = [{'1'},{'2'},{'3'},{'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'9'}, {'10'}, {'11'}, {'12'}, {'13'}];
nodes = listdlg('ListString', S, 'Name', 'Node status & correlations', 'PromptString', 'Select the node(s)');

switch property
    case 1
        
        %1- 360 degress model; 2- 180 degrees model, r=0.3; 3- 180 degrees
        %model, r=0.4; 4- 180 degrees model, r=0.4; 5- 180 degrees model, r=0.infty;
        %6 - circle intersection model
        
        A1 = choose_graph(hObject, handles, 1);
        A2 = choose_graph(hObject, handles, 2);
        A3 = choose_graph(hObject, handles, 3);
        A4 = choose_graph(hObject, handles, 4);
        A5 = choose_graph(hObject, handles, 5);
        A6 = choose_graph(hObject, handles, 6);
        
        %compute node status for each of these sensing graphs
        [nodeStatus1, nodeStatus1_avg] = compute_nodeStatus(hObject, handles, A1, start_frame, end_frame);
        [nodeStatus2, nodeStatus2_avg] = compute_nodeStatus(hObject, handles, A2, start_frame, end_frame);
        [nodeStatus3, nodeStatus3_avg] = compute_nodeStatus(hObject, handles, A3, start_frame, end_frame);
        [nodeStatus4, nodeStatus4_avg] = compute_nodeStatus(hObject, handles, A4, start_frame, end_frame);
        [nodeStatus5, nodeStatus5_avg] = compute_nodeStatus(hObject, handles, A5, start_frame, end_frame);
        [nodeStatus6, nodeStatus6_avg] = compute_nodeStatus(hObject, handles, A6, start_frame, end_frame);
        
        colors = jet(6);
        
        %average so that irregularities get smoothed out
        nodeStatus1_smooth = averaging(nodeStatus1(:, nodes(1)), 30, 1, 1);
        nodeStatus2_smooth = averaging(nodeStatus2(:, nodes(1)), 30, 1, 1);
        nodeStatus3_smooth = averaging(nodeStatus3(:, nodes(1)), 30, 1, 1);
        nodeStatus4_smooth = averaging(nodeStatus4(:, nodes(1)), 30, 1, 1);
        nodeStatus5_smooth = averaging(nodeStatus5(:, nodes(1)), 30, 1, 1);
        nodeStatus6_smooth = averaging(nodeStatus6(:, nodes(1)), 30, 1, 1);
        
        %plot node status without averaging
        figure(1)
        plot(start_frame:end_frame, nodeStatus1(:, nodes(1)), 'Color', colors(1,:), 'linewidth', 4)
        hold on
        plot(start_frame:end_frame, nodeStatus2(:, nodes(1)), 'Color', colors(2,:), 'linewidth', 4)
        hold on
        plot(start_frame:end_frame, nodeStatus3(:, nodes(1)), 'Color', colors(3,:), 'linewidth', 4)
        hold on
        plot(start_frame:end_frame, nodeStatus4(:, nodes(1)), 'Color', colors(4,:), 'linewidth', 4)
        hold on
        plot(start_frame:end_frame, nodeStatus5(:, nodes(1)), 'Color', colors(5,:), 'linewidth', 4)
        hold on
        plot(start_frame:end_frame, nodeStatus6(:, nodes(1)), 'Color', colors(6,:), 'linewidth', 4)
        
        legend('360', '180 deg, r03', '180 deg, r04', '180 deg, r05', '180 deg infty', 'New Model')
        
        %plot node status with averaging
        figure(2)
        plot(start_frame:end_frame, nodeStatus1_smooth, 'Color', colors(1,:), 'linewidth', 4)
        hold on
        plot(start_frame:end_frame, nodeStatus2_smooth, 'Color', colors(2,:), 'linewidth', 4)
        hold on
        plot(start_frame:end_frame, nodeStatus3_smooth, 'Color', colors(3,:), 'linewidth', 4)
        hold on
        plot(start_frame:end_frame, nodeStatus4_smooth, 'Color', colors(4,:), 'linewidth', 4)
        hold on
        plot(start_frame:end_frame, nodeStatus5_smooth, 'Color', colors(5,:), 'linewidth', 4)
        hold on
        plot(start_frame:end_frame, nodeStatus6_smooth, 'Color', colors(6,:), 'linewidth', 4)
        
        legend('360', '180 deg, r03', '180 deg, r04', '180 deg, r05', '180 deg infty', 'New Model')
        
        %correlation between node statuses for different sensing graphs
        disp('node status 360 vs. node status 180 0.3')
        [xcf, lags, bounds] = crosscorr(nodeStatus1_smooth(1:end_frame-start_frame+1), nodeStatus2_smooth(1:end_frame-start_frame+1));
        xcf
        bounds
        [c p] = corrcoef(nodeStatus1_smooth(1:end_frame-start_frame+1), nodeStatus2_smooth(1:end_frame-start_frame+1));
         c
         p
         
         disp('node status 360 vs. node status 180 0.5')
        [xcf, lags, bounds] = crosscorr(nodeStatus1_smooth(1:end_frame-start_frame+1), nodeStatus4_smooth(1:end_frame-start_frame+1));
        xcf
        bounds
        [c p] = corrcoef(nodeStatus1_smooth(1:end_frame-start_frame+1), nodeStatus4_smooth(1:end_frame-start_frame+1));
         c
         p
         
         disp('node status 180 0.5 vs. node status 180 infty')
        [xcf, lags, bounds] = crosscorr(nodeStatus4_smooth(1:end_frame-start_frame+1), nodeStatus5_smooth(1:end_frame-start_frame+1));
        xcf
        bounds
        [c p] = corrcoef(nodeStatus4_smooth(1:end_frame-start_frame+1), nodeStatus5_smooth(1:end_frame-start_frame+1));
         c
         p
         
         disp('node status 360 vs. node status 180 infty')
        [xcf, lags, bounds] = crosscorr(nodeStatus1_smooth(1:end_frame-start_frame+1), nodeStatus5_smooth(1:end_frame-start_frame+1));
        xcf
        bounds
        [c p] = corrcoef(nodeStatus1_smooth(1:end_frame-start_frame+1), nodeStatus5_smooth(1:end_frame-start_frame+1));
         c
         p
         
         disp('node status 360 vs. node status New Model')
        [xcf, lags, bounds] = crosscorr(nodeStatus1_smooth(1:end_frame-start_frame+1), nodeStatus6_smooth(1:end_frame-start_frame+1));
        xcf
        bounds
        [c p] = corrcoef(nodeStatus1_smooth(1:end_frame-start_frame+1), nodeStatus6_smooth(1:end_frame-start_frame+1));
         c
         p
         
         disp('node status 180 infty vs. node status New Model')
        [xcf, lags, bounds] = crosscorr(nodeStatus5_smooth(1:end_frame-start_frame+1), nodeStatus6_smooth(1:end_frame-start_frame+1));
        xcf
        bounds
        [c p] = corrcoef(nodeStatus5_smooth(1:end_frame-start_frame+1), nodeStatus6_smooth(1:end_frame-start_frame+1));
         c
         p
         
end

end

function autoPlotandComputeWHAT_pushbutton_Callback(~, ~, ~)
msgbox...
    ('Click this to automatically process data from multiple sets in the current directory. You will have the option to plot graphs for every data set and/or compute graph properties (which will be saved in the current directory as MATLAB variables). Graphs can be based on a number of nearest neighbors or a direct input radius, and you will have the opportunity to specify how to weight the edges. Be aware that processing many data sets and calculating graph properties and creating complex graphs can be time consuming.','"Auto Plot & Compute Graphs" Help','help');
end

% text boxes where selections are updated
function dataselection_editText_Callback(hObject, ~, handles)
guidata(hObject, handles);
end

function dataselection_editText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function nuneighbors_editText_Callback(~, ~, ~)
end
function nuneighbors_editText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function inputradius_editText_Callback(~, ~, ~)
end
function inputradius_editText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function connectedOrNot_editText_Callback(~, ~, ~)
end
function connectedOrNot_editText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function viewAngle_editText_Callback(~, ~, ~)
end
function viewAngle_editText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function LoadData_Callback(hObject, ~, handles)
str = {'Dance Data'};
[data_type,ok] = listdlg('PromptString','Select Data Type:', 'SelectionMode',...
    'single','ListString',str,'ListSize',[160 160]);
handles.data_type = data_type;
h = msgbox('Please indicate the range of configured data you would like to load');
uiwait(h);

%Gets the lower boundary 
lb = str2num(cell2mat(inputdlg('Please type the number of the first file to load')));
set(handles.lb_frame,'String',lb);
ub = str2num(cell2mat(inputdlg('Please type the number of the last file to load')));
set(handles.ub_frame,'String',ub);
if data_type == 1
    old_dir = cd;
    
    location = uigetdir('','Please browse to Dance folder');
    addpath(location);
    cd(location);
    for n = lb:ub
    file =[num2str(n) '.dat'];
    handles.dancedata{n} = importdata(file);
    end
    cd(old_dir);
    data = handles.dancedata{n};
    x=data(:,1);
    numberOfNodes = length(x);
end

dataInventory = [ub-lb+1; lb; ub];
handles.lb = lb;
handles.ub = ub;
handles.numberOfNodes = numberOfNodes;
handles.dataInventory = dataInventory;

guidata(hObject,handles);

end

%loads data in addition to the files loaded at one point
function LoadData2_pushbutton_Callback(hObject, ~, handles)

lb = handles.lb;
ub = handles.ub;

str = {'Dance Data'};
[data_type,ok] = listdlg('PromptString','Select Data Type:', 'SelectionMode',...
    'single','ListString',str,'ListSize',[160 160]);
handles.data_type = data_type;
h = msgbox('Please indicate the range of configured data you would like to load');
uiwait(h);

%Gets the lower boundary 
lb2 = str2num(cell2mat(inputdlg('Please type the number of the first file to load')));
ub2 = str2num(cell2mat(inputdlg('Please type the number of the last file to load')));
set(handles.ub_frame,'String',ub+ub2-lb2+1);

if data_type == 1
    old_dir = cd;
    
    location = uigetdir('','Please browse to Dance folder');
    addpath(location);
    cd(location);
    for n = lb:ub
    file =[num2str(n) '.dat'];
    handles.dancedata{ub+n} = importdata(file);
    end
    cd(old_dir);
end

%store this information so that we know data files could come from
%different videos
handles.ub = ub+ub2-lb2+1;
handles.dataInventory = [handles.dataInventory, [ub2-lb2+1; lb2; ub2]];
guidata(hObject,handles);

end

function ConfigureData_Callback(hObject, ~, handles)

% Configures Data, currently only set up for dance data but other types of
% data sets can be added
str = {'Dance Data'};
[data_type,ok] = listdlg('PromptString','Select Data Type:', 'SelectionMode',...
    'single','ListString',str,'ListSize',[160 160]);
handles.data_type = data_type;
if data_type == 1
    curr_dir = pwd;
    %cd('../Dance')
    location = uigetdir('','Please browse to Dance folder');
    addpath(location);
    cd(location);
    ReadFlock1;
    %cd(curr_dir)
    for i=1:length(X)
        xdata = transpose(X(i,:)); %#ok<*NODEF>
        ydata = transpose(Y(i,:));
        headingdata = transpose(Theta(i,:));
        position = horzcat(xdata,ydata,headingdata);
        save(strcat(num2str(i),'.dat'),'position','-ascii');
    end
% else %%% Here is where you add scripts to configure other data types, 
end    
guidata(hObject,handles);
end

function update_Callback(hObject, ~, handles)

%%%%% Plots Selection and Updates Calculations-----------------------------
cla;
data_type = handles.data_type;
Near_Neigh_Status = get(handles.nearNeigh,'Value');
plot_position_Status = get(handles.plot_position,'Value');
direct_edges_Status = get(handles.direct_edges,'Value');
undir_edge_Status = get(handles.undir_edge,'Value');
network_fig_Status = get(handles.network_fig,'Value');
indep_grp_Status = get(handles.indep_grp,'Value');
disp_status_Status = get(handles.disp_status,'Value');

inputChoice = (get(handles.edges_popup,'Value'));
inputEdgeChoice = str2num(get(handles.NeighOrRad2,'UserData'));
numberOfNeighbors = str2num(get(handles.nuneighbors_editText,'String'));

if inputEdgeChoice == 1
    numberOfNeighbors = str2num(get(handles.nuneighbors_editText,'String'));
elseif inputEdgeChoice == 2
    inputRadius = str2num(get(handles.inputradius_editText,'String'));
end

if data_type==1
    i = str2num(get(handles.dataselection_editText, 'String'));
    data = handles.dancedata{i};
    x=data(:,1);
    y=data(:,2);

    if plot_position_Status == 1
    plot(x,y,'g.');
    hold on;
    grid on;
    end
    positionMatrix=horzcat(x,y);
    dim = 2;
else
    bfilename = uigetfile('*.dat');
    border = load(bfilename);
    ifilename = uigetfile('*.dat');
    interior = load(ifilename);
    set(handles.dataselection_editText,'String',strcat(bfilename,' & ',ifilename));
    x = [border(:,1); interior(:,1)];
    y = [border(:,2); interior(:,2)];
    dim = size(border,2);
    divider = length(border);
    
    if dim==3
        z = [border(:,3); interior(:,3)];
        if plot_position_Status == 1
        plot3(x,y,z,'g.');
        end
        positionMatrix = horzcat(x,y,z);
    else
        if plot_position_Status == 1
        plot(x,y,'g.');
        end
        positionMatrix = horzcat(x,y);
    end
    hold on;
    grid on;
end
handles.dim = dim;
handles.x = x;
handles.y = y;
if dim ==3
    handles.z = z;
end
numberOfNodes = length(positionMatrix);
handles.numberOfNodes = numberOfNodes;

%position matrix => distance matrix => organized dist matrix
%=>kminDistances => adjancency matrix
distanceMatrix = ipdm(positionMatrix);
for i=1:numberOfNodes
    for j=1:numberOfNodes
        if distanceMatrix(i,j) == 0
            distanceMatrix(i,j) = Inf;
        end
    end
end

if inputEdgeChoice == 1
    organizedDistanceMatrix = sort(distanceMatrix, 2);
    kminDistances = organizedDistanceMatrix(:,1:numberOfNeighbors);
    
    adjacencyMatrix = zeros(numberOfNodes);
    for s=1:numberOfNodes
        distanceRow = distanceMatrix(s,:);
        minimumRow = kminDistances(s,:);
        adjacencyRow = ismember(distanceRow,minimumRow); 
        adjacencyMatrix(s,:) = adjacencyRow; 
    end
    save('UnweightedAdjacencyMatrix.mat','adjacencyMatrix');
elseif inputEdgeChoice == 2
    distanceMatrix = ipdm(positionMatrix);
%identify min and max overall radius 
    for i=1:numberOfNodes
        for j=1:numberOfNodes
            if distanceMatrix(i,j) == 0
              distanceMatrix(i,j) = Inf;
            end
        end
    end

    minimumRadius = min(distanceMatrix(:));
    set(handles.sliderMIN_editText,'String',num2str(minimumRadius));

    for i=1:numberOfNodes
        for j=1:numberOfNodes
            if distanceMatrix(i,j) == Inf
                distanceMatrix(i,j) = 0;
            end
        end
    end

    maximumRadius = max(distanceMatrix(:));
    set(handles.sliderMAX_editText,'String',num2str(maximumRadius));
%making sure the right radius
    inputRadius = str2num(get(handles.inputradius_editText,'String'));
    if (inputRadius == 0) || (inputRadius > maximumRadius) || (inputRadius < minimumRadius)
        inputRadius = str2num(cell2mat(inputdlg('Please input your desired radius.')));
        set(handles.inputradius_editText,'String',num2str(inputRadius));
    end
    guidata(hObject,handles);

    for i=1:numberOfNodes
        for j=1:numberOfNodes
            if distanceMatrix(i,j) == 0
                distanceMatrix(i,j) = Inf;
            end
        end
    end
%creating adjacency matrix
   adjacencyMatrix = zeros(numberOfNodes);
    for i=1:numberOfNodes
        for j=1:numberOfNodes
            if distanceMatrix(i,j) < inputRadius
                adjacencyMatrix(i,j) = 1;
            end
        end
    end
    save('UnweightedAdjacencyMatrix.mat','adjacencyMatrix');
    %elseif inputEdgeChoice == 3
end

sparseMatrix = zeros(numberOfNodes);
for i=1:numberOfNodes
    for j=1:numberOfNodes
        if adjacencyMatrix(i,j)==1
            sparseMatrix(i,j) = distanceMatrix(i,j);
        end
    end
end

[sensor,sensed] = find(adjacencyMatrix==1);
sortedNodePairs = sortrows([sensor,sensed],1);

%up-date weights
if inputChoice == 1
    weight = str2num(get(handles.fixedWeight,'String'));
end
if inputChoice == 2
    if inputEdgeChoice == 1
        weight = 1/numberOfNeighbors;
    elseif inputEdgeChoice == 2
        weight = zeros(numberOfNodes,1);
        for i =1:length(sortedNodePairs)
            weight(i) = 1/length(find(sortedNodePairs(:,1)==sortedNodePairs(i,1)));
        end
    end
end
if inputChoice == 3  %look for where it is established whether the node is border/interior
    bweight = str2num(cell2mat(inputdlg('Input fixed border edge weight.')));
    iweight = str2num(cell2mat(inputdlg('Input fixed interior edge weight.')));
end  % where are the other options?

%%%%% Display Node Status -------------------------------------------------
threshold = 0.5;
cnt_ns = zeros(1, numberOfNodes);
if disp_status_Status == 1
            shortestReversedLengths = zeros(numberOfNodes);
        for ii = 1:numberOfNodes
            for jj = [1:(ii-1) (ii+1):numberOfNodes]
                if adjacencyMatrix(jj,ii) > 0
                    shortestReversedLengths(ii,jj) = adjacencyMatrix(jj,ii);
                else
                    shortestReversedLengths(ii,jj) = Inf;
                end
            end
        end
        
        for kk = 1:numberOfNodes %interesting...
            for ii = 1:numberOfNodes
                for jj = 1:numberOfNodes
                    shortestReversedLengths(ii,jj) = min(shortestReversedLengths(ii,jj), shortestReversedLengths(ii,kk)+shortestReversedLengths(kk,jj));
                end
            end
        end
        
        inverseDistances = 1./shortestReversedLengths;
        nodeStatus = zeros(ub-lb+1,numberOfNodes);
        for ii = 1:numberOfNodes
            nodeStatus(i,ii) = 1/(numberOfNodes - 1)*sum(inverseDistances(ii,[1:(ii-1) (ii+1):numberOfNodes]));
            text(x(ii) + 0.04, y(ii), num2str(roundn(nodeStatus(i,ii),-2)),'HorizontalAlignment','center');
            text(x(ii)-0.1, y(ii), num2str(ii),'HorizontalAlignment','left'); %code I added to display node number alongside node status
%             allEdges = unique([find(adjacencyMatrix(ii,:)) find(adjacencyMatrix(:,ii)')]);
%             adjacentAngles = sort(relativeAngles(ii,allEdges));
%             angleDifferences = [diff(adjacentAngles) (adjacentAngles(1) - adjacentAngles(end) + 2*pi)];
%             [maxDiff locationOfMax] = max(angleDifferences);
%             textDirection = anglerestrict(adjacentAngles(locationOfMax) + maxDiff/2);
%             
%             text(x(ii) + 0.04*lengthScale*cos(textDirection), y(ii) + 0.04*lengthScale*sin(textDirection), num2str(roundn(nodeStatus(ii),-2)),'HorizontalAlignment','center');
% %             text(x(ii) + 0.04*lengthScale*cos(textDirection), y(ii) + 0.04*lengthScale*sin(textDirection), num2str(ii),'HorizontalAlignment','center');
        end
end

%%%%% SEE Nearest Neighbor Edges ------------------------------------------
if Near_Neigh_Status == 1
if inputEdgeChoice == 1
correspondingDistances = zeros((numberOfNeighbors*numberOfNodes),1);
for i=1:length(sortedNodePairs)
    sensorNode = sensor(i);
    sensedNode = sensed(i);
    correspondingDistances(i,:) = distanceMatrix(sensorNode,sensedNode);
end

%create u, v, and w matrices to use with quiver3 plot function where u is
%a sortedNodePairs x 1 matrix with each entry being the x vector between
%the ith and k nearest neighbor, and the same for y in the v matrix and z
%in the w matrix.
q=1;
for j=1:length(sortedNodePairs)
    xVal1 = x(sortedNodePairs(j,1));
    xVal2 = x(sortedNodePairs(j,2));
    u(j) = xVal2 - xVal1; %#ok<*AGROW>
    
    yVal1 = y(sortedNodePairs(j,1));
    yVal2 = y(sortedNodePairs(j,2));
    v(j) = yVal2 - yVal1;
    
    if dim == 3
        zVal1 = z(sortedNodePairs(j,1));
        zVal2 = z(sortedNodePairs(j,2));
        w(j) = zVal2 - zVal1;
    end
    
    if inputChoice == 3
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

%need sortedNodePairs x 1 matrix for each x, y, and z to copy position
%coordinates k number of times to use in quiver3 plot function.

xplot = x(sortedNodePairs(:,1));
u = transpose(u);
yplot = y(sortedNodePairs(:,1));
v = transpose(v);
if dim==3
    zplot = z(sortedNodePairs(:,1));
    w = transpose(w);
end

% direct input fixed weight option
if inputChoice == 1
    if dim==3
        quiver3(xplot,yplot,zplot,u,v,w,0,'AutoScale','off','Color','black','LineWidth',weight);
    else quiver(xplot,yplot,u,v,0,'AutoScale','off','Color','black','LineWidth',weight);
    end
end

if inputChoice == 2
    if dim==3
        quiver3(xplot,yplot,zplot,u,v,w,0,'AutoScale','off','Color','black','LineWidth',weight);
    else quiver(xplot,yplot,u,v,0,'AutoScale','off','Color','black','LineWidth',weight);
    end
end

% border distinction - 2 direct input weights option
if inputChoice == 3
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

    if dim==3
        quiver3(xbplot,ybplot,zbplot,ub,vb,wb,0,'AutoScale','off','LineWidth',bweight,'Color','blue');
        quiver3(xiplot,yiplot,ziplot,ui,vi,wi,0,'AutoScale','off','LineWidth',iweight,'Color','magenta');
    else
        quiver(xbplot,ybplot,ub,vb,0,'AutoScale','off','LineWidth',bweight,'Color','blue');
        quiver(xiplot,yiplot,ui,vi,0,'AutoScale','off','LineWidth',iweight,'Color','magenta');
    end
end

if inputChoice == 4
    for i=1:length(correspondingDistances)
        if correspondingDistances(i) == Inf
            correspondingDistances(i) = 0;
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
elseif inputEdgeChoice ==2
% have to do things all over again because we are in the case when neighbors are
% determined by the radius
[sensor,sensed] = find(adjacencyMatrix == 1);
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
    
    if inputChoice == 3
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

% direct input fixed weight option
if inputChoice == 1
    edgeWeight = str2num(get(handles.fixedWeight,'String'));
    if dim==3
        quiver3(xplot,yplot,zplot,u,v,w,0,'ShowArrowHead','off','LineWidth',edgeWeight);
    else quiver(xplot,yplot,u,v,0,'ShowArrowHead','off','LineWidth',edgeWeight);
    end
end
% inverse out degree weight option
if inputChoice == 2
    % compute weights for each node based on how many connections they
    % make, and store a matrix with these weights to use in quiver function
    edgeWeights = zeros(numberOfNodes,1);
    for i =1:length(edges)
        edgeWeights(i) = 1/length(find(edges(:,1)==edges(i,1)));
        if dim==3
            quiver3(xplot(i),yplot(i),zplot(i),u(i),v(i),w(i),0,'ShowArrowHead','off','LineWidth',2*edgeWeights(i));
        else quiver(xplot(i),yplot(i),u(i),v(i),0,'ShowArrowHead','off','LineWidth',2*edgeWeights(i));
        end
    end
end

% border distinction - 2 direct input weights option
if inputChoice == 3
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
    
    [rw,cl]=find(ub==0);
    ub(rw,:)=[];
    [rw,cl]=find(vb==0);
    vb(rw,:)=[];
    [rw,cl]=find(ui==0);
    ui(rw,:)=[];
    [rw,cl]=find(vi==0);
    vi(rw,:)=[];
    if dim==3
        [rw,cl]=find(wb==0);
        wb(rw,:)=[];
        [rw,cl]=find(wi==0);
        wi(rw,:)=[];
    end

    borderinput = str2num(cell2mat(inputdlg('Input fixed weight for border edges.'))); %#ok<*ST2NM>
    interiorinput = str2num(cell2mat(inputdlg('Input fixed weight for interior edges.')));
    if dim==3
        quiver3(xbplot,ybplot,zbplot,ub,vb,wb,0,'ShowArrowHead','off','LineWidth',borderinput);
        quiver3(xiplot,yiplot,ziplot,ui,vi,wi,0,'ShowArrowHead','off','LineWidth',interiorinput);
    else 
        quiver(xbplot,ybplot,ub,vb,0,'ShowArrowHead','off','LineWidth',borderinput);
        quiver(xiplot,yiplot,ui,vi,0,'ShowArrowHead','off','LineWidth',interiorinput);
    end
end

% function input weight option
if inputChoice == 4
    correspondingDistances = zeros(length(edges),1);
    for i=1:length(edges)
        correspondingDistances(i) = distanceMatrix(sensor(i),sensed(i));
        edgeWeight(i) = 1/(correspondingDistances(i));
        if dim==3
            quiver3(xplot(i),yplot(i),zplot(i),u(i),v(i),w(i),0,'ShowArrowHead','off','LineWidth',2*edgeWeight(i));
        else quiver(xplot(i),yplot(i),u(i),v(i),0,'ShowArrowHead','off','LineWidth',2*edgeWeight(i));
        end
    end
end
end
end
clear u; clear v; clear w; clear ub; clear vb; clear wb;
clear ui; clear vi; clear wi;
%%%%% SEE Directed Edges --------------------------------------------------
reverseNodePairs = horzcat(sensed, sensor);
check = ismember(sortedNodePairs, reverseNodePairs, 'rows');
if direct_edges_Status == 1
for i=1:length(reverseNodePairs)
    if check(i)==0
        directedEdges(i,:) = sortedNodePairs(i,:);
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
    
    if inputChoice == 3
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

if inputChoice == 1
    edgeWeight = str2num(get(handles.df_weight, 'String'));
    if dim==3
        quiver3(xplotd,yplotd,zplotd,a,b,d,0,'AutoScale','off','Color','b','LineWidth',edgeWeight);
    else quiver(xplotd,yplotd,a,b,0,'AutoScale','off','Color','b','LineWidth',edgeWeight);
    end
end

if inputChoice == 2
    edgeWeight = 1/numberOfNeighbors;
    if dim==3
        quiver3(xplotd,yplotd,zplotd,a,b,d,0,'AutoScale','off','Color','b','LineWidth',edgeWeight);
    else quiver(xplotd,yplotd,a,b,0,'AutoScale','off','Color','b','LineWidth',edgeWeight);
    end
end

if inputChoice == 3
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
    
    borderinput = str2num(cell2mat(inputdlg('Input fixed weight for border individuals.'))); %#ok<*ST2NM>
    interiorinput = str2num(cell2mat(inputdlg('Input fixed weight for interior individuals.')));
    if dim==3
        quiver3(xbplot,ybplot,zbplot,ub,vb,wb,0,'AutoScale','off','LineWidth',borderinput,'Color','blue');
        quiver3(xiplot,yiplot,ziplot,ui,vi,wi,0,'AutoScale','off','LineWidth',interiorinput,'Color','blue');
    else
        quiver(xbplot,ybplot,ub,vb,0,'AutoScale','off','LineWidth',borderinput,'Color','blue');
        quiver(xiplot,yiplot,ui,vi,0,'AutoScale','off','LineWidth',interiorinput,'Color','blue');
    end
end

if inputChoice == 4
    sensord = directedEdges(:,1);
    sensedd = directedEdges(:,2);
    correspondingDDistances = zeros(length(directedEdges),1);
    
    for i=1:length(directedEdges)
        sensorNode = sensord(i);
        sensedNode = sensedd(i);
        correspondingDDistances(i,:) = distanceMatrix(sensorNode,sensedNode);
    end
    for g=1:length(correspondingDDistances)
        if correspondingDDistances(g) == Inf
            correspondingDDistances(g) = 0;
        end
    end
    maxDistance = max(correspondingDDistances(:));
    for i=1:length(correspondingDDistances)
        blue(i) = correspondingDDistances(i)/maxDistance;
        red(i) = correspondingDDistances(i)/maxDistance;
        if dim==3
            obj = quiver3(xplotd(i),yplotd(i),zplotd(i),a(i),b(i),d(i),0,'Color',[red(i),0,blue(i)]);
            set(obj,'AutoScale','off');
        else
            obj = quiver(xplotd(i),yplotd(i),a(i),b(i),0,'Color',[red(i),0,blue(i)]);
            set(obj,'AutoScale','off');
        end
    end
end
end
clear a; clear b; clear d; clear ub; clear vb; clear wb;
clear ui; clear vi; clear wi;
%%%%% SEE Undirected Edges-------------------------------------------------
if undir_edge_Status == 1
    for i=1:length(reverseNodePairs)
    if check(i)==1
        undirectedEdges(i,:) = sortedNodePairs(i,:);
    end
    end

[r,c]=find(undirectedEdges==0);
undirectedEdges(r,:)=[];

r=1;
clear u; clear v;
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
    
    if inputChoice == 3
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

if inputChoice == 1
    edgeWeight = str2num(get(handles.uf_weight,'String'));
    if dim==3
        quiver3(xplotu,yplotu,zplotu,u,v,w,0,'ShowArrowHead','off','Color','red','LineWidth',edgeWeight);
    else quiver(xplotu,yplotu,u,v,0,'ShowArrowHead','off','Color','red','LineWidth',edgeWeight);
    end
end

if inputChoice == 2
    edgeWeight = 1/numberOfNeighbors;
    if dim==3
        quiver3(xplotu,yplotu,zplotu,u,v,w,0,'ShowArrowHead','off','Color','red','LineWidth',edgeWeight);
    else quiver(xplotu,yplotu,u,v,0,'ShowArrowHead','off','Color','red','LineWidth',edgeWeight);
    end
end

if inputChoice == 3
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
    
    borderinput = str2num(cell2mat(inputdlg('Input fixed weight for undirected border edges.'))); %#ok<*ST2NM>
    interiorinput = str2num(cell2mat(inputdlg('Input fixed weight for undirected interior edges.')));
    if dim==3
        quiver3(xbplot,ybplot,zbplot,ub,vb,wb,0,'ShowArrowHead','off','LineWidth',borderinput,'Color','red');
        quiver3(xiplot,yiplot,ziplot,ui,vi,wi,0,'ShowArrowHead','off','LineWidth',interiorinput,'Color','red');
    else 
        quiver(xbplot,ybplot,ub,vb,0,'ShowArrowHead','off','LineWidth',borderinput,'Color','red');
        quiver(xiplot,yiplot,ui,vi,0,'ShowArrowHead','off','LineWidth',interiorinput,'Color','red');
    end
end

if inputChoice == 4
    sensoru = undirectedEdges(:,1);
    sensedu = undirectedEdges(:,2);
    correspondingUDistances = zeros(length(undirectedEdges),1);
    for i=1:length(undirectedEdges)
        sensorNode = sensoru(i);
        sensedNode = sensedu(i);
        correspondingUDistances(i,:) = distanceMatrix(sensorNode,sensedNode);
    end
    
    for i=1:length(correspondingUDistances)
        if correspondingUDistances(i) == Inf
            correspondingUDistances(i) = 0;
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
end
clear u; clear v; clear w; clear ub; clear vb; clear wb;
clear ui; clear vi; clear wi;    
%%%%% SEE Network in separate Figure --------------------------------------
%CREATES SEPARATE FIGURE WITH BOXES CORRESPONDING TO EACH NODE IN THE DATA
%SET (CHOSEN BY THE USER- BORDER DATA AND INTERIOR DATA), AND EDGES
%CONNECTING (USER INPUT NUMBOR OF) NEAREST NEIGHBORS. ALSO SAVES THE
%ADJACENCY MATRIX AS AdjacencyMatrix.mat TO CURRENT DIRECTORY (MATLAB).
%DOES NOT REPRESENT POSITION DATA, SO NOT AN EFFECTIVE VISUALIZATION TOOL,
%BUT CREATES AN INTERESTING IMAGE THAT SHOWS CONNECTED GROUPS WITHIN THE
%LARGER NETWORK. NOTE: VERY TIME INTENSIVE- ONLY USE FOR SMALL NUMBERS OF
%NEIGHBORS. ALSO, BE PATIENT WHILE THE FIGURE LOADS- MAY TAKE A WHILE.
if (network_fig_Status)
    nodeNames = (1:numberOfNodes);
    nodeNames = transpose(nodeNames);

    sparseMatrix = zeros(numberOfNodes);
    for i=1:numberOfNodes
        for j=1:numberOfNodes
            if adjacencyMatrix(i,j)==1
                sparseMatrix(i,j) = distanceMatrix(i,j);
            end
        end
    end
    gObj = biograph(sparseMatrix,nodeNames);
    view(gObj);
end

%%%%% SEE Independent Groups ----------------------------------------------
if indep_grp_Status == 1
for i=1:length(sortedNodePairs)
    j=1;
    nodePath(1,:) = sortedNodePairs(i,:);
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
            quiver3(xplotcycle,yplotcycle,zplotcycle,f,g,h,0,'MaxHeadSize',0.01,'ShowArrowHead','on','Color',[red,green,blue]);
        else quiver(xplotcycle,yplotcycle,f,g,0,'MaxHeadSize',0.01,'ShowArrowHead','on','Color',[red,green,blue]);
        end
        clear('nodePath','xplotcycle','yplotcycle','zplotcycle','f','g','h');
    end
    clear('nodePath');
end
end

%create matrices for directed/undirected connections and corresponding adjacency matrices
if inputEdgeChoice == 1
    reverseNodePairs = horzcat(sensed, sensor);
    check = ismember(sortedNodePairs, reverseNodePairs, 'rows');
    
    for i=1:length(reverseNodePairs)
        if check(i)==1
            undirectedEdges(i,:) = sortedNodePairs(i,:);
        end
        if check(i)==0
            directedEdges(i,:) = sortedNodePairs(i,:);
        end
    end
    
    [r,c]=find(undirectedEdges==0);
    undirectedEdges(r,:)=[];
    
    [r,c]=find(directedEdges==0); %#ok<*NASGU>
    directedEdges(r,:)=[];
    
    adjacencyMatrixU = adjacencyMatrix;
    adjacencyMatrixU(directedEdges) = 0;
    adjacencyMatrixD = adjacencyMatrix;
    adjacencyMatrixD(undirectedEdges) = 0;
    
    sparseMatrixU = zeros(numberOfNodes);
    sparseMatrixD = zeros(numberOfNodes);
    for i=1:numberOfNodes
        for j=1:numberOfNodes
            if adjacencyMatrixU(i,j)==1
                sparseMatrixU(i,j) = distanceMatrix(i,j);
            end
            if adjacencyMatrixD(i,j)==1
                sparseMatrixD(i,j) = distanceMatrix(i,j);
            end
        end
    end
    
%%%%% calculates number of strongly connected components in the graph.-----    
    BGUObj = biograph(sparseMatrixU);
    [stronglyConnectedU, C] = conncomp(BGUObj,'Directed',false);
    BGDObj = biograph(sparseMatrixD);
    [stronglyConnectedD, C2] = conncomp(BGDObj,'Directed',true);
    set(handles.numberStronglyConnectedComponents_editText,'String',num2str(stronglyConnectedU+stronglyConnectedD));
    
    BGWObj = biograph(sparseMatrix);
    [weaklyConnected, C3] = conncomp(BGWObj,'Weak',true);
    set(handles.numberWeaklyConnectedComponents_editText,'String',num2str(weaklyConnected));
elseif inputEdgeChoice == 2
    BGSObj = biograph(sparseMatrix);
    [stronglyConnected, C4] = conncomp(BGSObj,'Directed',false);
    set(handles.numberStronglyConnectedComponents_editText,'String',num2str(stronglyConnected));
    set(handles.numberWeaklyConnectedComponents_editText,'String',num2str(stronglyConnected));
end
%%%%%%% Computes algebraic connectivity of graph--------------------------- 
for i=1:numberOfNodes
    for j=1:numberOfNodes
        if adjacencyMatrix(i,j)==1
            % direct input for edge weights
            if inputChoice == 1
                adjacencyMatrix(i,j) = weight;
            end
            % inverse of number of neighbors/out degree
            if inputChoice == 2
                if inputEdgeChoice == 1
                    adjacencyMatrix(i,j) = weight;
                elseif inputEdgeChoice == 2
                    for k=1:length(sortedNodePairs)
                        adjacencyMatrix(sensor(k),sensed(k)) = weight(k);
                    end
                end
            end
            % border & interior distinction
            if inputChoice == 3
                if i <= divider
                    adjacencyMatrix(i,j) = bweight;
                else adjacencyMatrix(i,j) = iweight;
                end
            end
            % 1/distance decay function
            if inputChoice == 4
                adjacencyMatrix(i,j) = 1/distanceMatrix(i,j);
            end
        end
    end
end
save('AdjacencyMatrix.mat','adjacencyMatrix');
 
N = numberOfNodes;
A = adjacencyMatrix;
D = diag(sum(A,2));
L = D - A;
Q = zeros(N-1,N);
for i = 1:(N-1)
    for j = 1:i
        Q(i,j) = -1/sqrt(i*(i+1));
    end
    Q(i,i+1) = i/sqrt(i*(i+1));
end
Lbarsym = (1/2)*Q*(L+L')*Q';
eigen = eig(Lbarsym);
eigen = sort(eigen);
algebraicConnectivity = roundn(eigen(1,1),-3);
set(handles.algebraicConnectivity_editText,'String',num2str(algebraicConnectivity));

%%%%% Computes speed of convergence of consensus---------------------------
eigenl = eig(L);
eigenl = sort(eigenl);
secondE = eigenl(2,1);
speedOfConvergence = roundn(real(secondE),-3);
        if speedOfConvergence > 0
            set(handles.connectedOrNot_editText,'String','Yes');
        else
            set(handles.connectedOrNot_editText,'String','No');
        end
set(handles.speedOfConvergence_editText,'String',num2str(speedOfConvergence));

%%%%% Computes H2 Norm-----------------------------------------------------
evals = sort(real(eig(L)));
if evals(2) > 1e-10
    H2 = H2_norm(L, Q);
else
    H2 = Inf;
end
set(handles.h2Norm_editText,'String',num2str(H2));

%%%%% Computes L2 Gain-----------------------------------------------------
algebraicConnectivity_L = eigen(2,1);
if algebraicConnectivity_L > 0
    L2 = roundn(1/algebraicConnectivity_L,-3);
    set(handles.l2gain_editText,'String',num2str(L2));
else msgbox('L2 gain could not be computed because algebraic connectivity is not positive.');
end

%%%%% Calculates upper and lower eigenvalue bounds for H2 Norm-------------
algebraicConnectivity_b = eigen(1,1);
if algebraicConnectivity_b <= 0
    msgbox('Eigenvalue bounds on H2 Norm cannot be computed because algebraic connectivity is not greater than zero.','ERROR','error');
    return;
end

eigenL = sort(eig(L));
sumLB = 0;
for i=2:N
    sumLB = sumLB + 1/(real(eigenL(i)));
end
lowerEBound = roundn(sqrt((1/2) * sumLB),-3);
set(handles.lowerEigenBound_editText,'String',num2str(lowerEBound));

Lbar = Q*L*Q';
LBS = (1/2)*(Lbar+Lbar');
eigenLBS = eig(LBS);
sumLBS = 0;
for i=1:N-1
    sumLBS = sumLBS + 1/(eigenLBS(i));
end
upperEBound = roundn(sqrt((1/2)*sumLBS),-3);
set(handles.upperEigenBound_editText,'String',num2str(upperEBound));

guidata(hObject, handles);
end

function edges_popup_Callback(hObject, ~, handles)
    if get(handles.edges_popup,'Value') == 1
     set(handles.fixedWeight,'Visible','on');  
     set(handles.df_weight,'Visible','on'); 
     set(handles.uf_weight,'Visible','on'); 
     set(handles.fixedtext,'Visible','on'); 
     set(handles.dirfixedtext,'Visible','on'); 
     set(handles.undirfixedtext,'Visible','on');   
    else
     set(handles.fixedWeight,'Visible','off'); 
     set(handles.df_weight,'Visible','off'); 
     set(handles.uf_weight,'Visible','off'); 
     set(handles.fixedtext,'Visible','off'); 
     set(handles.dirfixedtext,'Visible','off'); 
     set(handles.undirfixedtext,'Visible','off'); 
    end
guidata(hObject,handles);
end

function edges_popup_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function lb_frame_Callback(~, ~, ~)
end
function lb_frame_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function ub_frame_Callback(~, ~, ~)
end
function ub_frame_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function fixedWeight_Callback(hObject, ~, handles)
guidata(hObject,handles);
end
function fixedWeight_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function NeighOrRad2_CreateFcn(hObject, ~, handles)
set(hObject,'UserData','1')
guidata(hObject,handles);
end

function NeighOrRad2_Callback(hObject, ~, handles)
guidata(hObject,handles);
end

function NeighOrRad2_SelectionChangeFcn(hObject, eventdata, ~)
handles = guidata(hObject);
switch get(eventdata.NewValue,'Tag')
    case 'numNeigh'
      set(handles.NeighOrRad2,'UserData','1');
      set(handles.nearNeigh,'String','Nearest Neighbors')
      set(handles.direct_edges,'Visible','on')
      set(handles.undir_edge,'Visible','on')
      set(handles.dirfixedtext,'Visible','on')
      set(handles.df_weight,'Visible','on')
      set(handles.undirfixedtext,'Visible','on')
      set(handles.uf_weight,'Visible','on')
%      set(handles.eitherorWhat_pushbutton,'Visible','on')
    case 'radiusInput'
      set(handles.NeighOrRad2,'UserData','2')
      set(handles.nearNeigh,'String','Edges')
      set(handles.direct_edges,'Visible','off')
      set(handles.undir_edge,'Visible','off')
      set(handles.dirfixedtext,'Visible','off')
      set(handles.df_weight,'Visible','off')
      set(handles.undirfixedtext,'Visible','off')
      set(handles.uf_weight,'Visible','off')
%      set(handles.eitherorWhat_pushbutton,'Visible','off')
end
guidata(hObject, handles);
end

function NewModel_Callback(hObject, eventdata, handles)
end

function NewModel2_Callback(hObject, eventdata, handles)
end

function plot_position_Callback(~, ~, ~)
end

function nearNeigh_Callback(~, ~, ~)
end

function direct_edges_Callback(~, ~, ~)
end

function undir_edge_Callback(~, ~, ~)
end

function network_fig_Callback(~, ~, ~)
end

function indep_grp_Callback(~, ~, ~)
end

function df_weight_Callback(~, ~, ~)
end
function df_weight_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function uf_weight_Callback(~, ~, ~)
end
function uf_weight_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function disp_status_Callback(~, ~, ~)
end

%plotting the trajectory of nodes of choice then choosing what time
%frames to display

function plot_trace_pushbutton_Callback(hObject, ~, handles)

lb = handles.lb;
ub = handles.ub;

S = [{'1'},{'2'},{'3'},{'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'9'}, {'10'}, {'11'}, {'12'}, {'13'}, {'center'}];
nodes = listdlg('ListString', S, 'Name', 'Plot trace', 'PromptString', 'Select the node(s)');

prompt = {'Enter starting frame:','Enter end frame:', 'How frequently to display frame:'};
answer = inputdlg(prompt, 'Select start/end frame');
lb1 = str2num(answer{1});
ub1 = str2num(answer{2});
incr = str2num(answer{3});

if lb1<lb | ub1>ub
    msgbox('Wrong start or end frames', 'Error');
    return;
else
    lb = lb1;
    ub = ub1;
end

if incr<0 | incr>ub
    msgbox('Incorrect increment', 'Error');
    return;
end

numberOfNodes = handles.numberOfNodes;
color = jet(numberOfNodes+1);

for j=1:length(nodes)
node=nodes(j);

for i=lb:ub
    
    data=handles.dancedata{i};
    x = data(:,1);
    y = data(:,2);
    theta = data(:,3);
    
%[x y theta] = transf_in_familiar_coord(x,y,theta);        

if node<=numberOfNodes
    xx1(i-lb+1)=x(node);
    yy1(i-lb+1)=y(node);
else
    xx1(i-lb+1)=mean(x);
    yy1(i-lb+1)=mean(y);
end

end

fig = figure(1);
plot(xx1, yy1, 'Color', color(j,:), 'LineStyle', '--');
hold on;
plot(xx1(1), yy1(1), 'Marker', '*', 'MarkerSize', 20, 'MarkerFaceColor', 'r');
hold on;
text(xx1(1)+0.1, yy1(1), num2str(nodes(j)), 'HorizontalAlignment', 'center')

if incr>0
for i=lb:incr:ub
  plot(xx1(i-lb+1), yy1(i-lb+1), 'Marker', 'x', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
  hold on;
  str = num2str(i);
  text(xx1(i-lb+1)+0.02, yy1(i-lb+1), str, 'HorizontalAlignment', 'center');
  hold on;
end
end

title('Trajectories (traces) of the nodes selected');
hold on;
end

%close(fig);
%imwrite(fig, 'trace11.jpg');  %how to save the figure?

guidata(hObject, handles)
end

%computes coordinates of the center, direction of the center and the angle
%formed by the line connecting the center and some node of choice

function [center, center_avg, center2, dir_center, dir_node, dir_node2] = center_dir(hObject, handles, delay, start_frame, end_frame)
        
numberOfNodes = handles.numberOfNodes;
lb = handles.lb;
ub = handles.ub;

if end_frame+delay>ub
   u=ub;
else
   u = end_frame+delay;
end

    for i=start_frame:u
        data=handles.dancedata{i};
        x(:,i-start_frame+1) = data(:,1);
        y(:,i-start_frame+1) = data(:,2);
        theta(:,i-start_frame+1) = data(:,3);
        
        %coordinates of center
        center(i-start_frame+1,:) = [mean(x(:,i-start_frame+1)), mean(y(:,i-start_frame+1))];
    end
    
    %smoothing the coordinates of the center out
    [center_avg(:,1), ignore, ignore, ignore, ignore] = averaging(center(:,1),delay,1,1);
    [center_avg(:,2), ignore, ignore, ignore, ignore] = averaging(center(:,2),delay,1,1);
    
    for i=start_frame:end_frame
        
        if i+delay>ub
            u =end_frame;
        else
            u = i+delay;
        end
        
        %direction of movement of center
            center2(i-start_frame+1,:) = center_avg(u-start_frame+1,:)-center_avg(i-start_frame+1,:); %vector
            dir_center(i-start_frame+1) = atan2(center2(i-start_frame+1,2), center2(i-start_frame+1,1)); %direction of movement
          for j=1:numberOfNodes  
            angle(j, i-start_frame+1) = atan2(y(j, i-start_frame+1)-center_avg(i-start_frame+1,2), x(j, i-start_frame+1)-center_avg(i-start_frame+1,1));
            dir_node(i-start_frame+1, j) = anglerestrict(dir_center(i-start_frame+1)-angle(j,i-start_frame+1));
          end
    end
   
    for j=1:numberOfNodes
       dir_node2(:,j) = averaging(dir_node(:,j),delay,1,1);
    end
        
    
 figure(1)   

    plot(start_frame:end_frame-delay, dir_center(1:end_frame-start_frame+1-delay)) 
    colors = jet(numberOfNodes);
    
 figure(2)
    for j=1:numberOfNodes
        plot(start_frame:end_frame-delay, dir_node2(1:end_frame-start_frame+1-delay,j), 'Color', colors(j,:))
        hold on
    end
      %}  
    %[Distance2, Distance3, distance, leader, backer, N_leader, N_backer] = LeaderBacker(hObject, handles, start_frame, end_frame,1);
 %{   
    frame = 300;
    figure(100)
    plot(x(:, frame-start_frame+1), y(:, frame-start_frame+1), 'x')
    hold on
    quiver(center_avg(frame-start_frame+1,1), center_avg(frame-start_frame+1,2), center2(frame-start_frame+1,1), center2(frame-start_frame+1,2), 'color', 'g')
    quiver(center_avg(frame-start_frame+1,1), center_avg(frame-start_frame+1,2), x(7,frame-start_frame+1)-center_avg(frame-start_frame+1,1), y(7,frame-start_frame+1)-center_avg(frame-start_frame+1,2), 'color', 'b')
    
    disp('momentu adevarului')
    dir_node2(frame-start_frame+1,7)*180/pi
    dir_node(frame-start_frame+1,7)*180/pi
    dir_center(frame-start_frame+1)*180/pi
    angle(7, frame-start_frame+1)*180/pi
    Distance3(frame-start_frame+1,7)
    distance(frame-start_frame+1,7)
   %} 
    
end

%smoothing  out coordinates, and other measures like node status SG/IG
function [mat2 , mat3, xax, mat4, xx] = averaging(mat, interval,interval2, space)

   n = length(mat);
   i=1;
   cnt=1;
   
   while i<=n
       xax(cnt)=i;
       if i+interval>n
           mat2(i:n) = mean(mat(i:n));
           mat3(cnt) = mean(mat(i:n));
       else 
           
           mat2(i:i+interval2-1) = mean(mat(i:i+interval-1));
           mat3(cnt) = mean(mat(i:i+interval-1));
       end
       
       i=i+interval2;
       cnt=cnt+1;
       
   end
   
   xx = xax(1):space:n;
   mat4 = spline(xax, mat3, xx);

   %there are different ways of smoothing functions; I used mat2 most often
   %interval determines how many frames to average and interval2 how many
   %frames to move forward after averaging
end

%tranform into coordinates of the center
function [positionMatrix2, center] = transform_body_frame(positionMatrix, handles, hObject)

end_frame = size(positionMatrix,1);
delay=20;
numberOfNodes = handles.numberOfNodes;

for i = 1:end_frame
    
    center(i,:) = [mean(positionMatrix(i,:,1)), mean(positionMatrix(i,:,2))];
    
  if i<=end_frame-delay  
    center2(i,:) = [mean(positionMatrix(i+delay,:,1)), mean(positionMatrix(i+delay,:,2))];
  else
    center2(i,:) = [mean(positionMatrix(end_frame,:,1)), mean(positionMatrix(end_frame,:,2))];
  end
  
    dir_center(i,:)= center2(i,:) - center(i,:);
    angle_center(i) = anglerestrict(atan2(dir_center(i,2), dir_center(i,1)));
    
    for j=1:numberOfNodes
        positionMatrix2(i,j,:) = [positionMatrix(i,j,1)*cos(angle_center(i))+ positionMatrix(i,j,2)*sin(angle_center(i)), -positionMatrix(i,j,1)*sin(angle_center(i))+ positionMatrix(i,j,2)*cos(angle_center(i))];
    end
    
end
end

%compute inertia matrix, find first moment of inertia (for normalization purposes)
%see wikipedia.com for formulas on inertia matrix used here

function  [I, V, D] = inertia_matrix(positionMatrix, handles, hObject)

    numberOfNodes = handles.numberOfNodes;
    end_frame = size(positionMatrix,1);
    delay = 20;
    
    I = zeros(end_frame, 3, 3);
    
    [positionMatrix, center] = transform_body_frame(positionMatrix, handles, hObject);
    
    for i=1:end_frame
      for j=1:numberOfNodes  
       x = center(i,1) - positionMatrix(i,j,1);
       y = center(i,2) - positionMatrix(i,j,2);
       I(i,:,:) = squeeze(I(i,:,:)) + 1/numberOfNodes*[-y*y, x*y, 0; x*y, -x*x, 0; 0, 0, -x*x-y*y];
      end
      
      [V(i,:,:),D(i,:,:)] = eig(squeeze(I(i,:,:)));
      
    end
    
end

% average distance to other nodes in time and over all the video, average distance to closest nodes,
% average angle to other nodes and to closest nodes

function [positionMatrix, distanceMatrix, avg_distance, avg_distance_inTime, avg_distance_inTime_closest1, avg_distance_inTime_closest2, avg_distanceAll, avg_dist_index, closest_neighbors, relativeAngles, HeadingDiff, angleToclosest, Avg_angleToclosest, Avg_angleOverall] = computeAll_pos_dist(handles, start, ending, option)

  for i=start:ending
   
        data=handles.dancedata{i};
        x = data(:,1);
        y = data(:,2);
        theta = data(:,3);
        
        %[x y theta] = transf_in_familiar_coord(x,y,theta);        

        dim = 2;
        positionMatrix(i,:,:) = horzcat(x,y);
      
      numberOfNodes = length(squeeze(positionMatrix(i,:,:)));
      distanceMatrix(i,:,:) = ipdm(squeeze(positionMatrix(i,:,:)));
    
    if i==start
        avg_distance = zeros(numberOfNodes, numberOfNodes);
        avg_angle = zeros(1,numberOfNodes);
        for j=1:numberOfNodes
         avg_distance(j,j) = Inf;
        end
    end
    
    %avg_distance computes the average distance to all other nodes
    avg_distance = avg_distance + squeeze(distanceMatrix(i,:,:));
    
   for s=1:numberOfNodes 
       %average distance of some node to other nodes at a certain frame i
       avg_distance_inTime(i,s) = mean(distanceMatrix(i,s,:));
       %determine closest nodes
       [ignore, Ind_dist_s_sorted] = sort(squeeze(distanceMatrix(i,s,:)));
       %mean distance to closest 2 nodes
       avg_distance_inTime_closest2(i,s) = mean([distanceMatrix(i,s,Ind_dist_s_sorted(1)), distanceMatrix(i,s,Ind_dist_s_sorted(2))]);
       %mean distance to closest node
       avg_distance_inTime_closest1(i,s) = distanceMatrix(i,s,Ind_dist_s_sorted(1));
   end 
   
   %relative angles from one node to the others
   % Warning: considering the average of these angles is a mistake since a
   % we need to consider the BORDER EFFECTS. 
   % This code needs to be modified to consider border effects!
   % see Cavagna A. et al., The STARFLAG handbook on collective animal behaviour:
   % Part II, three-dimensional analysisThe STARFLAG handbook on collective
   % animal behaviour.
   
    for s = 1:(numberOfNodes-1)
        for j = (s+1):numberOfNodes
            relativeAngles(i,s,j) = atan2(squeeze(positionMatrix(i,j,2))-squeeze(positionMatrix(i,s,2)),squeeze(positionMatrix(i,j,1))-squeeze(positionMatrix(i,s,1)));
            relativeAngles(i,j,s) = anglerestrict(relativeAngles(i,s,j) + pi);
        end
    end
        
    for s=1:numberOfNodes
        for j=1:numberOfNodes
            
            if distanceMatrix(i,s,j) == 0
                distanceMatrix(i,s,j) = Inf;
            end
            %heading differences of nodes
            HeadingDiff(i,s,j) = abs(anglerestrict(theta(s)-theta(j)));
        end
    end
    
    for j=1:numberOfNodes
         [ignore close_neighbors] = sort(distanceMatrix(i,j,:));
         closest_neighbors(i,j,:) = close_neighbors(1:2);
    end
   
  if option==1  
    for j=1:numberOfNodes
         angleToclosest(i,j,1) = anglerestrict(theta(j) - relativeAngles(i,j,closest_neighbors(i,j,1)))*180/pi; %angle to closest node
         angleToclosest(i,j,2) = anglerestrict(theta(j) - relativeAngles(i,j,closest_neighbors(i,j,2)))*180/pi; % angle to second closest node
         angleToclosest(i,j,3) = (angleToclosest(i,j,1)+angleToclosest(i,j,2))/2; % mean between the two angles
         angleToclosest(i,j,4) = (abs(angleToclosest(i,j,1))+abs(angleToclosest(i,j,2)))/2; %average of absolute values between 2; could be of use!
         
         if angleToclosest(i,j,1)<angleToclosest(i,j,2) %smaller and larger angle
             angleToclosest(i,j,5) = angleToclosest(i,j,1); 
             angleToclosest(i,j,6) = angleToclosest(i,j,2);
         else
             angleToclosest(i,j,5) = angleToclosest(i,j,2);
             angleToclosest(i,j,6) = angleToclosest(i,j,1);
         end
    end
  end
  
  end
  
  avg_distance = avg_distance/(ending-start+1);
  
   for j=1:numberOfNodes
        avg_distance2 = squeeze(avg_distance(j,:));
        avg_distance2(j)=[];
        %the mean over all frames
        avg_distanceAll(j) = mean(avg_distance2);
    end
    
    [ignore, avg_dist_index] = sort(avg_distanceAll);
  
  for j=1:numberOfNodes
      %over all frames
    Avg_angleToclosest(j,1) = mean(angleToclosest(:,j,1));
    Avg_angleToclosest(j,2) = mean(angleToclosest(:,j,2));
    Avg_angleToclosest(j,3) = mean((angleToclosest(:,j,1)+angleToclosest(:,j,2))/2);
    Avg_angleToclosest(j,4) = mean(angleToclosest(:,j,4));
    Avg_angleToclosest(j,5) = mean(angleToclosest(:,j,5));
    Avg_angleToclosest(j,6) = mean(angleToclosest(:,j,6));
  end
  
  Avg_angleOverall(1) = mean(Avg_angleToclosest(:,1));
  Avg_angleOverall(2) = mean(Avg_angleToclosest(:,2));
  Avg_angleOverall(3) = mean(Avg_angleToclosest(:,3));
  Avg_angleOverall(4) = mean(Avg_angleToclosest(:,4));
  Avg_angleOverall(5) = mean(Avg_angleToclosest(:,5));
  Avg_angleOverall(6) = mean(Avg_angleToclosest(:,6));
  
end
 
function position_in_group_pushbutton_Callback(hObject, ~, handles)

S = [{'Plot trace of center'},{'Position of nodes'},{'Summary of front/back positions'}, {'Average relative angles'}, {'Average position matrix'}];
option = listdlg('ListString',S,'Name','Options','PromptString','Select a task');

if option==1
    
    %compute position of center
 prompt = {'Enter starting frame:','Enter end frame:'};   
 answer = inputdlg(prompt, 'Select start/end frame');
 lb1 = str2num(answer{1});
 end_frame = str2num(answer{2});
 
 [positionMatrix, c, numberOfNodes] = Read_coordinatesAndCenter(lb1, end_frame, handles, hObject); 

 %plot trajectory of center
fig = figure(1);
plot(c(:,1), c(:,2), 'Color', 'b', 'LineStyle', '--');
hold on;
plot(c(1,1), c(1,2), 'Marker', '*', 'MarkerSize', 20, 'MarkerFaceColor', 'r');

elseif option == 2
    %computing distance to center and the average of it over all frames
 S = [{'1'},{'2'},{'3'},{'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'9'}, {'10'}, {'11'}, {'12'}, {'13'}];
 nodes = listdlg('ListString', S, 'Name', 'Plot trace', 'PromptString', 'Select the node(s)');

 prompt = {'Enter starting frame:','Enter end frame:'};   
 answer = inputdlg(prompt, 'Select start/end frame');
 start_frame = str2num(answer{1});
 end_frame = str2num(answer{2});
 
[positionMatrix, center, numberOfNodes] = Read_coordinatesAndCenter(start_frame, end_frame, handles, hObject);
[dist_toCenter, avg_toCenter] = distToCenter(handles, hObject, start_frame, end_frame, nodes,1);
 
 elseif option == 3
     
 lb = handles.lb;
 ub = handles.ub;
 prompt = {'Enter starting frame:','Enter end frame:'};   
 answer = inputdlg(prompt, 'Select start/end frame');
 start_frame = str2num(answer{1});
 end_frame = str2num(answer{2});

 %distance to center-axis and determining whether node is in front or in the back
 [Distance2, Distance3, distance, leader, backer, N_leader, N_backer] = LeaderBacker(hObject, handles, start_frame, end_frame,1);
     
elseif option==4
    
    prompt = {'Enter starting frame:','Enter end frame:'};   
    answer = inputdlg(prompt, 'Select start/end frame');
    lb1 = str2num(answer{1});
    end_frame = str2num(answer{2});
     [positionMatrix, distanceMatrix, avg_distance, avg_distance_inTime, avg_distance_inTime_closest1, avg_distance_inTime_closest2,  avg_distanceAll, avg_dist_index, closest_neighbors, relativeAngles, HeadingDiff, angleToclosest, Avg_angleToclosest, Avg_angleOverall]= computeAll_pos_dist(handles, lb1, end_frame, 1);
    
    prompt = {'Enter node:'};
    answer = inputdlg(prompt, 'Select node');
    node = str2num(answer{1});
    
    %print results, plot figures
    figure(1)
    plot(lb1:end_frame,angleToclosest(lb1:end_frame,node,1), 'b')
    hold on;
    plot(lb1:end_frame,angleToclosest(lb1:end_frame,node,2), 'k')
    title('Angle to closest and second closest of chosen node')
    xlabel('frames')
    ylabel('angle')
    figure(2)
    plot(lb1:end_frame,angleToclosest(lb1:end_frame,node,3), 'g')
    title('Average of angles to closest and second closest node')
    xlabel('frames')
    ylabel('angle')
    figure(3)
    plot(lb1:end_frame,angleToclosest(lb1:end_frame,node,4), 'c')
    title('Average of absolute values of angles to closest and second closest nodes')
    xlabel('frames')
    ylabel('angle')
    figure(4)
    plot(lb1:end_frame,angleToclosest(lb1:end_frame,node,5), 'b')
    hold on;
    plot(lb1:end_frame,angleToclosest(lb1:end_frame,node,6), 'r')
    title('Average smallest and largest angle to closest 2 nodes')
    xlabel('frames')
    ylabel('angle')
    legend('smallest angle', 'largest angle')
    
    figure(5)
    bar(Avg_angleToclosest(:,1))
    title('Average angle of nodes to closest node')
    figure(6)
    bar(Avg_angleToclosest(:,2))
    title('Average angle of nodes to second closest node')
    figure(7)
    bar(Avg_angleToclosest(:,3))
    title('Average angle of nodes closest nodes (both)')
    figure(8)
    bar(Avg_angleToclosest(:,4))
    title('abs value of Average angle of nodes')
    figure(9)
    bar(Avg_angleToclosest(:,5))
    title('Smallest average angle of nodes to closest nodes')
    figure(10)
    bar(Avg_angleToclosest(:,6))
    title('Largest average angle of nodes to closest nodes')
    
    disp('Overall node Average angle is:')
    Avg_angleOverall(1)
    Avg_angleOverall(2)
    Avg_angleOverall(3)
    Avg_angleOverall(4)
    Avg_angleOverall(5)
    Avg_angleOverall(6)
    
elseif option ==5
    
    prompt = {'Enter starting frame:','Enter end frame:'};   
    answer = inputdlg(prompt, 'Select start/end frame');
    start_frame = str2num(answer{1});
    end_frame = str2num(answer{2});
    [positionMatrix, distanceMatrix, avg_distance, avg_distance_inTime, avg_distance_inTime_closest1, avg_distance_inTime_closest2,  avg_distanceAll, avg_dist_index, closest_neighbors, relativeAngles, HeadingDiff, angleToclosest, Avg_angleToclosest, Avg_angleOverall] = computeAll_pos_dist(handles, start_frame, end_frame, 1);
    
    prompt = {'Enter node:'};
    answer = inputdlg(prompt, 'Select node');
    node = str2num(answer{1});
   
    numberOfNodes = length(avg_distance);
    disp('Average distance matrix for chosen node is:')
    squeeze(avg_distance(node,:))
    [ignore,avg_distance_sorted_ind] = sort(squeeze(avg_distance(node,:)));
    disp('Closest nodes to chosen nodes are:')
    avg_distance_sorted_ind(1:2)
   
    for j=1:numberOfNodes
        avg_distance2 = squeeze(avg_distance(j,:));
        avg_distance2(j)=[];
        avg_distanceAll(j) = mean(avg_distance2); %average over all nodes
    end
    
    disp('Overall how far away are nodes from each other')
    avg_distanceAll
    figure(1)
    bar(avg_distanceAll)
    
    [ignore, avg_dist_index] = sort(avg_distanceAll);
    disp('Sorted nodes according to average distance to all other nodes  -- a measure of isolation/centrality')
    avg_dist_index
end
    
guidata(hObject, handles)
end

% distance to center (first moment) on each frame and average distance to center over all frames
function [dist_to_center, avg_toCenter] = distToCenter(handles, hObject, start_frame, end_frame, nodes, n)

[positionMatrix, center, numberOfNodes] = Read_coordinatesAndCenter(start_frame, end_frame, handles, hObject);
colorPlot = jet(numberOfNodes);

if nodes==numberOfNodes
    nodes = 1:numberOfNodes;
end

avg_toCenter = zeros(1, length(nodes));

[I, V, D] = inertia_matrix(positionMatrix, handles, hObject);

for i=1:end_frame-start_frame+1
    sorted_eig(i,:) = sort([D(i,1,1), D(i,2,2), D(i,3,3)]);
    ratio_eig(i) = sorted_eig(i,1)/sorted_eig(i,3);
end

for i = start_frame:end_frame
    
   for j= 1:length(nodes)
        pos_mx = [center(i-start_frame+1,:); squeeze(positionMatrix(i-start_frame+1, nodes(j), :))'];
        dist_to_center(i-start_frame+1,j) = pdist(pos_mx);  %distance to center of group of node j at frame i
   end
       
   max_dist = max(dist_to_center(i-start_frame+1,:));
   
   for j=1:length(nodes)
    if n==1       
        %computing distances and normalizing
      dist_to_center(i-start_frame+1,j) = dist_to_center(i-start_frame+1,j)/max_dist;
      avg_toCenter(j) = avg_toCenter(j) + dist_to_center(i-start_frame+1,j)/max_dist;
  
      %dist_to_center(i-start_frame+1,j) = dist_to_center(i-start_frame+1,j)/ratio_eig(i-start_frame+1);
      %avg_toCenter(j) = avg_toCenter(j) + dist_to_center(i-start_frame+1,j)/ratio_eig(i-start_frame+1);
    else
        
      avg_toCenter(j) = avg_toCenter(j) + dist_to_center(i-start_frame+1,j);
    end
   end
      
end   

for j=1:length(nodes)
    std_dist_to_Center(nodes(j)) = std(dist_to_center(:,j));
end

 %{
 f = figure(1);
 title('Average position tot center of several nodes')
 xlabel('frames')
 ylabel('distance from center')
 
 for j=1:size(nodes,2)
     plot(start_frame:end_frame, dist_to_center(:,j), 'Color', colorPlot(j,:));
   hold on;
 end
 %}

 avg_toCenter = avg_toCenter/(end_frame-start_frame+1);

disp('sorted distances')
 %}
 disp('hierarchy of distance to center(point)')
[distances, sorted_dist] = sort(avg_toCenter);
sorted_dist
distances

 disp('hierarchy of STD distance to center(point)')
[std_distances, sorted_std_dist] = sort(std_dist_to_Center);
sorted_std_dist
std_distances

f = figure(2);
 barwitherr(std_dist_to_Center, avg_toCenter);
 title('Average distance from center at particular frame')

end

%distance from center-axis, denoted by D, i.e. the projection of the node
%position from center on the direction of movement  of the whole group;
%also, absolute vaue of D, and node most up-front and most in the back
%how many times node is in the back/front
function [Distance2, Distance3, distance, leader, backer, N_leader, N_backer] = LeaderBacker(hObject, handles, start_frame, end_frame, n)

c=[]; c2 = [];
leader = [];
backer = [];
 
 data=handles.dancedata{start_frame};
 x = data(:,1);
 y = data(:,2);
 N = size(y,1);
 
 [dist_to_center, avg_toCenter] = distToCenter(handles, hObject, start_frame, end_frame, 13, 0);
 [positionMatrix, center, numberOfNodes] = Read_coordinatesAndCenter(start_frame, end_frame, handles, hObject); 
 [I, V, D] = inertia_matrix(positionMatrix, handles, hObject);

for i=1:end_frame-start_frame+1
    sorted_eig(i,:) = sort([D(i,1,1), D(i,2,2), D(i,3,3)]);
    ratio_eig(i) = sorted_eig(i,1)/sorted_eig(i,3);
end

%read data; must end 20 frames before end_frame because direction at frame
%i is taken as position at frame i+20 - position at frame i. 
 for i = start_frame:end_frame-20
    
    distance1 = zeros(1,N);
    distance2 = zeros(1,N);
     
    data=handles.dancedata{i};
    x = data(:,1);
    y = data(:,2);
    theta = data(:,3);
    %[x y theta] = transf_in_familiar_coord(x,y,theta);
    c = [mean(x), mean(y)];  %compute center
     
    data=handles.dancedata{i+20};
    x2 = data(:,1);
    y2 = data(:,2);
    theta2 = data(:,3);
    %[x2 y2 theta2] = transf_in_familiar_coord(x2, y2, theta2);
    c2 = [mean(x2), mean(y2)];
    
    %two ways of computing D: by analytical geometry formulas of distances
    %to axes; OR easier, by computing projection on the direction of movement
    
    m_old = (c2(2)-c(2))/(c2(1)-c(1));
    m_new = -1/m_old;
    n_new = c(2) - m_new*c(1);
    
    a = 1; b= -m_new; d = -n_new;
    
    %which side of the center-axis the node is
    if a*c2(2)+b*c2(1)+d >=0 
     for j = 1:N 
         
       Distance2(i-start_frame+1,j) = (a*y(j)+b*x(j)+d)/(norm([a,b]));
       
       angle1(j) = atan2(y(j)-c(2), x(j)-c(1));
       angle2(j) = atan2(c2(2)-c(2), c2(1)-c(1));
       angle3(j) = anglerestrict(angle1(j)-angle2(j)); 
       Distance3(i-start_frame+1,j) = dist_to_center(i-start_frame+1,j)*cos(angle3(j)); %second way of computing D
       
      %distance(i-start_frame+1,j) = abs(a*y(j)Distance3+b*x(j)+d)/(norm([a,b]));
      distance(i-start_frame+1,j) = abs(Distance3(i-start_frame+1,j));
       
      if a*y(j)+b*x(j)+d >=0   
          distance1(j) = abs(Distance3(i-start_frame+1,j)); %|D|, absolute value of D
          %distance1(j) = (a*y(j)+b*x(j)+d)/(norm([a,b])); %same side
      else
          %distance1(j) = abs(Distance3(i-start_frame+1,j));
          distance2(j) = (a*y(j)+b*x(j)+d)/(norm([a,b])); %opposite side
      end
     end
     
     %maximum of D, d, in order to normalize per frame
     max_Distance2 = max(Distance2(i-start_frame+1,:));
     max_distance = max(distance(i-start_frame+1,:));
     max_Distance3 = max(Distance3(i-start_frame+1,:));
     
   if n==1   
     Distance2(i-start_frame+1,:) = Distance2(i-start_frame+1,:)/max_Distance2;
     distance(i-start_frame+1, :) = distance(i-start_frame+1, :)/max_distance;
     Distance3(i-start_frame+1,:) = Distance3(i-start_frame+1,:)/max_distance;
     
     %normalizing by first moment of inertia
     %Distance2(i-start_frame+1,:) = Distance2(i-start_frame+1,:)/ratio_eig(i-start_frame+1);
     %distance(i-start_frame+1, :) = distance(i-start_frame+1, :)/ratio_eig(i-start_frame+1);
     %Distance3(i-start_frame+1,:) = Distance3(i-start_frame+1,:)/ratio_eig(i-start_frame+1);
     
   end  
     %the opposite side of the center-axis
    else
        
     for j = 1:N 
         
        Distance2(i-start_frame+1,j) = -(a*y(j)+b*x(j)+d)/(norm([a,b]));
        
        angle1(j) = atan2(y(j)-c(2), x(j)-c(1));
        angle2(j) = atan2(c2(2)-c(2), c2(1)-c(1));
        angle3(j) = anglerestrict(angle1(j)-angle2(j)); 
        Distance3(i-start_frame+1,j) = dist_to_center(i-start_frame+1,j)*cos(angle3(j));
        
        %distance(i-start_frame+1,j) = abs(a*y(j)+b*x(j)+d)/(norm([a,b]));
        distance(i-start_frame+1,j) = abs(Distance3(i-start_frame+1,j));
        
      if a*y(j)+b*x(j)+d <0   
        distance1(j) = -(a*y(j)+b*x(j)+d)/(norm([a,b])); %same side
      else
        distance2(j) = -(a*y(j)+b*x(j)+d)/(norm([a,b])); %opposite side
      end
      
     end
     
     max_Distance2 = max(Distance2(i-start_frame+1,:));
     max_distance = max(distance(i-start_frame+1,:));
     max_Distance3 = max(Distance3(i-start_frame+1,:));
     
   if n==1  
       
     Distance2(i-start_frame+1,:) = Distance2(i-start_frame+1,:)/max_Distance2;
     Distance3(i-start_frame+1,:) = Distance3(i-start_frame+1,:)/max_distance;
     distance(i-start_frame+1, :) = distance(i-start_frame+1, :)/max_distance;
     
     %Distance2(i-start_frame+1,:) = Distance2(i-start_frame+1,:)/ratio_eig(i-start_frame+1);
     %distance(i-start_frame+1, :) = distance(i-start_frame+1, :)/ratio_eig(i-start_frame+1);
     %Distance3(i-start_frame+1,:) = Distance3(i-start_frame+1,:)/ratio_eig(i-start_frame+1);
     
   end   
     
    end
    
    [ignore, leader(i-start_frame+1)] = max(distance1);
    [ignore, backer(i-start_frame+1)] = min(distance2);
    end

 disp('front and back means')
 for j=1:N
     
     Distance3_mean(j) = mean(Distance3(:,j));
     Distance3_std(j) = std(Distance3(:,j));
     Distance3_mean(j)
     mean(Distance3(j,Distance3(j,:)>=0))
     mean(Distance3(j,Distance3(j,:)<=0))
 end
 
 figure(1)
 barwitherr(Distance3_std, Distance3_mean)
 [ignore, ind_Distance3_mean] = sort(Distance3_mean);
 disp('hierarchy of nodes to center as axis, NOT absolute value')
 ind_Distance3_mean
 ignore
 [ignore, ind_Distance3_std] = sort(Distance3_std);
 disp('hierarchy of nodes STD to center as axis, NOT absolute value')
 ind_Distance3_std
 ignore
 
 disp('average distance to center of each node')
 %}
for j=1:N
    j
    disp(mean(distance(:,j)));
    mean_distToCenter(j) = mean(distance(:,j));
    std_distToCenter(j) = std(distance(:,j));
end

[ignore, ind_mean_distToCenter] = sort(mean_distToCenter);
disp('Hierarchy of distances to axis of center')
ind_mean_distToCenter
ignore

[ignore, ind_std_distToCenter] = sort(std_distToCenter);
disp('Hierarchy of distances to axis of center')
ind_std_distToCenter
ignore

disp('max and min average distance to center')
[max_distToCenter, max_ind]=max(mean_distToCenter)
[min_distToCenter, min_ind]=min(mean_distToCenter)

% how many times a node leads
N_leader = zeros(1, N); N_backer = zeros(1, N);

     for i=start_frame:end_frame-20
       N_leader(leader(i-start_frame+1)) = N_leader(leader(i-start_frame+1))+1;
       N_backer(backer(i-start_frame+1)) = N_backer(backer(i-start_frame+1))+1;
     end
  
     N_leader = N_leader/(end_frame-start_frame-19);
     N_backer = N_backer/(end_frame-start_frame-19);
     
     %showing how much % of the time the respective node is leader/backer
     figure(2)
     bar(N_leader);
     title('How much % of the time node is in front')
     figure(3);
     bar(N_backer);
     title('How much % of the time node is in the back');
  
end

%for each frame largest distance between nodes
function d = diameter(handles, hObject, start_frame, end_frame)

    for i =start_frame:end_frame 
        data=handles.dancedata{i};
        positions = ipdm(data(:,1:2));
        d(i-start_frame+1) = max(max(positions));
    end
end

function distTo_other(handles, hObject, start_frame, end_frame)
end

function [directions PosAndAngles] = directionTable(hObject, handles) % a table with directions of each node

%there are 2 alternatives for this vector of directions: either in regular,
%Euclidean coordinates of the room OR coordinates determined by how the
%center of the group moves. The later seems more in handy.

    lb = handles.lb;
    ub = handles.ub;
    delay2 = 1;
    
      file=strcat('Dance',num2str(lb),'.dat');
      data_type = handles.data_type;
      if data_type==1
        data=handles.dancedata{lb};
        x = data(:,1);
        y = data(:,2);
        theta = data(:,3);
        %[x y theta] = transf_in_familiar_coord(x,y,theta);        
        dim = 2;
        PosMatrix1_lb = horzcat(x,y);
        centerGroup1 = [mean(x), mean(y)];
      end
      
      numberOfNodes = length(x);
      
      file=strcat('Dance',num2str(lb+delay2),'.dat');
      data_type = handles.data_type;
      if data_type==1
        data=handles.dancedata{lb+delay2};
        x = data(:,1);
        y = data(:,2);
        theta = data(:,3);
        %[x y theta] = transf_in_familiar_coord(x,y,theta);        
        dim = 2;
        PosMatrix2_lb = horzcat(x,y);
        centerGroup2 = [mean(x), mean(y)];
      end
      
      PosMatrix_lb = PosMatrix2_lb - PosMatrix1_lb;
      centerDir_lb = centerGroup2 - centerGroup1;
   
      C_min = 0.9;
      
  for j=1:numberOfNodes
    
   index = 1;   
      
      theta1 = atan2(PosMatrix_lb(j,2), PosMatrix_lb(j,1));
      theta2 = atan2(centerDir_lb(2), centerDir_lb(1));
      thetaFinal = anglerestrict(theta2 - theta1);
      cosTheta = cos(thetaFinal);
      sinTheta = sin(thetaFinal);
      
      directions(j,1,1)=lb;
      directions(j,2,1)=1;
      %directions(j,3:4,1)= [PosMatrix(j,1),PosMatrix(j,2)]/norm([PosMatrix(j,1),PosMatrix(j,2)]);
      directions(j,3:4,1) = [cosTheta, sinTheta];
      directions(j,5,1) = thetaFinal;

      %group coordinates
      %PosAndAngles(lb, j, 1:2) = [cosTheta, sinTheta]';
      %room coordinates
      PosAngles(lb,j,1:2) = [PosMatrix_lb(j,1),PosMatrix_lb(j,2)]/norm([PosMatrix_lb(j,1),PosMatrix_lb(j,2)]);
      
   for i=lb+1:ub-1
       
      file=strcat('Dance',num2str(i),'.dat');
      data_type = handles.data_type;
      if data_type==1
        data=handles.dancedata{i};
        x = data(:,1);
        y = data(:,2);
        theta = data(:,3);
        %[x y theta] = transf_in_familiar_coord(x,y,theta);        
        dim = 2;
        PosMatrix1 = horzcat(x,y);    
        center_group1 = [mean(x), mean(y)];
      end
    
    if i+delay2>ub
        u=ub;
    else
        u=i+delay2;
    end
    
      file=strcat('Dance',num2str(u),'.dat');
      data_type = handles.data_type;
      if data_type==1
        data=handles.dancedata{u};
        x = data(:,1);
        y = data(:,2);
        theta = data(:,3);
        %[x y theta] = transf_in_familiar_coord(x,y,theta);        
        dim = 2;
        PosMatrix2 = horzcat(x,y);    
        center_group2 = [mean(x), mean(y)];
      end
      
      PosMatrix = PosMatrix2 - PosMatrix1;
      centerDir = centerGroup2 - centerGroup1;
      
      theta1 = atan2(PosMatrix(j,2), PosMatrix(j,1));
      theta2 = atan2(centerDir(2), centerDir(1));
      thetaFinal = anglerestrict(theta2 - theta1);
      cosTheta = cos(thetaFinal);
      sinTheta = sin(thetaFinal);
      
      %group coordinates
      PosAndAngles(i, j, 1:2) = [cosTheta, sinTheta]';
      %room coordinates
      %PosAndAngles(i,j,1:2) = [PosMatrix(j,1), PosMatrix(j,2)]/norm([PosMatrix(j,1), PosMatrix(j,2)]);
      
      if  dot(squeeze(PosAndAngles(i,j,1:2)), squeeze(directions(j,3:4,index))) >= C_min 
          directions(j, 2, index) = directions(j, 2, index)+1;
      else
          index=index+1;
          directions(j, 1, index) = i;
          directions(j, 2, index) = 1;
          directions(j, 3, index) = PosAndAngles(i,j,1);
          directions(j, 4, index) = PosAndAngles(i,j,2);
          directions(j, 5, index) = thetaFinal;
      end
      
  end
  end
  
end

function [LeadAndLag, frame, frame1, frame2] = lead_and_lag(hObject, handles) 

%directionTable used below when implementing lead_and_lag for one frame; uncertain if it
%is necessary

   %[directions, PosAndAngles] = directionTable(hObject, handles)
  
   numberOfNodes = handles.numberOfNodes;
  
   frame = 0; frame1 = 0; frame2 = 0;
      
      lb = handles.lb;
      ub = handles.ub;
   
  %two options: (a). run the algorithm for one frame; (b). run the algorithm for multiple frames and then average   
    prompt = {'Run lead & lag algorithm for:                                             Select 1). One frame                                                                   Select 2). Average of more frames'};
    answer_initial = cell2mat(inputdlg(prompt));
    answer = str2num(answer_initial);
    
    %one frame: there are two versions
   if answer==1
    
    prompt = {'Enter  frame to analyze lead and lag network for:'};   
    answer = inputdlg(prompt, 'Select frame');
    frame = str2num(answer{1});
    
      file=strcat('Dance',num2str(frame),'.dat');
      data_type = handles.data_type;
      if data_type==1
        data=handles.dancedata{frame};
        set(handles.dataselection_editText,'String',frame);
        x = data(:,1);
        y = data(:,2);
        theta = data(:,3);
        %[x y theta] = transf_in_familiar_coord(x,y,theta);        
      end
      
    lb = handles.lb;
    ub = handles.ub;
    %{
    %using directionTable (version #1)
  %parameters of method: C_min and delay  
    C_min=0.99;
    delay=25;
    
    numberOfNodes = size(directions,1);
    LeadAndLag = zeros(numberOfNodes, numberOfNodes);  %matrix of interactions: LeadAndLag(j,k)=1 means that j follows k

    %explore pairs (j,k)
    for j=1:numberOfNodes
          
            ind = find(squeeze(directions(j,1,:))>frame,1)-1;  %find frame in directionTable
            pos_j = [directions(j,3,ind),directions(j,4,ind)];
            %pos_j = [PosAndAngles(frame, j, 1), PosAndAngles(frame, j, 2)];
            k=2; 
            
           while k<=numberOfNodes
            if k~=j  
         
             ind2 = find(squeeze(directions(k,1,:))>frame,1)-1;
             if directions(k,1,ind2)-delay<=0
                start=lb;
             else
                start = directions(k,1,ind2)-delay;
             end
           
             ind3 = find(squeeze(directions(k,1,:))>start,1)-1;
             fr = directions(k,1,ind3);
    
         %go back and see where the trajectories mathced most
             while fr <= directions(j,1,ind)   
                pos_k = [directions(k,3,ind3),directions(k,4,ind3)];
                
                %if dot product is higher than C_min (threshold)
                
                if dot(pos_j, pos_k)>=C_min
                    LeadAndLag(j,k)=1;
                end
             
                ind3 = ind3+1;
                fr = directions(k,1,ind3);
                
             end
            end

             k=k+1;
            end  
    end
   
    for j=1:numberOfNodes
        for k=j+1:numberOfNodes
    %if the following is mutual, we are not interested
            if LeadAndLag(j,k)==1 & LeadAndLag(k,j)==1
                LeadAndLag(j,k)=0;
                LeadAndLag(k,j)=0;
            end
        end
    end
   LeadAndLag
  %}  
 
    %version #2, when I do not use directionTable
  delay=25;
  C_min=0.99;
  
  [theta0, theta, theta2, theta3, theta4, theta_avg, theta_avg2] = direction(handles, hObject);
  [center, center_avg, center2, dir_center, dir_node, dir_node2] = center_dir(hObject, handles, delay, frame-delay, frame);
  
    for i=max(frame-delay,lb):frame
     
        %in the coordinate frame of the room
      cosTheta = cos(theta4(i-lb+1, :))';
      sinTheta = sin(theta4(i-lb+1, :))';

      %projections in the direction the group moves for each frame,
      %actually not necessary
    
      projections(i-max(frame-delay, lb)+1,:,:) = [cosTheta, sinTheta];
      
      %in the coordinate frame of the center of the group
      %{
      for node = 1:numberOfNodes
        projections(i-max(frame-delay, lb)+1,node,:) = [cos(anglerestrict(theta4(i-lb+1,node)-dir_center(i-max(frame-delay, lb)+1))), sin(anglerestrict(theta4(i-lb+1,node)-dir_center(i-max(frame-delay, lb)+1)))];
      end
     %}
    end
    
      numberOfNodes = length(x);
      LeadAndLag = zeros(numberOfNodes,numberOfNodes);
      
      for j=1:numberOfNodes
          for k = 1:numberOfNodes
       
        if j~=k    
            C_tau=[];  
            
            %evaluating how similar directions (and implicitly,
            %trajectories) are when comparing two nodes and looking at
            %history in the last tau frames
            
            for tau=max(-delay,lb-frame):0
              C_tau(tau-max(-delay,lb-frame)+1) = dot(projections(min(delay, frame-lb)+1,j,:), projections(tau-max(-delay,lb-frame)+1,k,:));
            end
            
            C_tau
              [C(j,k), delay_tau(j,k)] = max(C_tau);
              delay_tau(j,k) = delay_tau(j,k)+max(-delay,lb-frame)-1;
              
              if max(C_tau)>=C_min & delay_tau(j,k)>-delay  & delay_tau(j,k)<delay
                  
                    LeadAndLag(j,k) = 1;
                  
                    if LeadAndLag(k,j)==1 & LeadAndLag(j,k)==1 
                        LeadAndLag(j,k) = 0;
                        LeadAndLag(k,j) = 0;
                    end
                   
              end
        end
        
          end
        
      end
       
    LeadAndLag
    delay_tau
    C
    gObj = biograph(LeadAndLag);
    view(gObj);
    
    %}
  
   else
      
       %this version of lead and lag consideres the average trajectory
       %match over more frames
       
    prompt = {'Enter start frame to analyze lead and lag network for:', 'Enter end frame'};   
    answer = inputdlg(prompt, 'Select frames');
    frame1 = str2num(answer{1});
    frame2 = str2num(answer{2});
    
    lb = handles.lb;
    ub = handles.ub;
   
    C_min=0.99;
    delay=25;
    delay2 = 5;
    
    [theta0, theta, theta2, theta3, theta4, theta_avg, theta_avg2] = direction(handles, hObject);
    [center, center_avg, center2, dir_center, dir_node, dir_node2] = center_dir(hObject, handles, delay2, max(frame1-delay,lb), frame2);
    
    for i = max(frame1-delay,lb):frame2
      
       %using coordinates of the room 
       
      cosTheta = cos(theta4(i-lb+1,:))';
      sinTheta = sin(theta4(i-lb+1,:))';
      
      projections(i-frame1+delay+1,:,:) = [cosTheta, sinTheta];
       
      %using coordinates of the center of the group
      %{
      for node = 1:numberOfNodes
        projections(i-max(frame1-delay,lb)+1,node,:) = [cos(anglerestrict(theta4(i-lb+1,node)-dir_center(i-max(frame1-delay,lb)+1))), sin(anglerestrict(theta4(i-lb+1,node)-dir_center(i-max(frame1-delay,lb)+1)))];
      end
      %}
    end
   
    numberOfNodes = length(squeeze(projections(1,:,:)));
    LeadAndLag = zeros(numberOfNodes,numberOfNodes);
    sum_C = zeros(numberOfNodes,numberOfNodes);
    sum_delay = zeros(numberOfNodes,numberOfNodes);
    
    for i=frame1:frame2
      
        for j = 1:numberOfNodes
          for k = 1:numberOfNodes
            
              %compare trajectories, over multiple frames, one more
              %dimension added to C_tau
            C_tau=[];  
            for tau=max(-delay,lb-i):0  
              C_tau(i-frame1+1,tau-max(-delay,lb-i)+1) = dot(projections(i-frame1+delay+1,j,:), projections(i-frame1+delay+1+tau,k,:));
            end
            
              [C(i-frame1+1,j,k), delay_tau(i-frame1+1,j,k)] = max(C_tau(i-frame1+1,:));
              delay_tau(i-frame1+1,j,k) = delay_tau(i-frame1+1,j,k)+max(-delay,lb-i)-1;
              
              sum_C(j,k) = sum_C(j,k)+C(i-frame1+1,j,k);
              sum_delay(j,k) = sum_delay(j,k) + delay_tau(i-frame1+1,j,k);
              
          end
      end
        
    end
      
    %average over many frames
    sum_C = sum_C/(frame2-frame1+1);
    sum_delay = sum_delay/(frame2-frame1+1);
     
    for j=1:numberOfNodes
        for k=1:numberOfNodes
            
            if sum_C(j,k)>C_min 
                %& sum_delay(j,k)>-delay  & sum_delay(j,k)<delay
                LeadAndLag(j,k) = 1;
                if LeadAndLag(k,j)==1 & LeadAndLag(j,k)==1
                    LeadAndLag(k,j) = 0;
                    LeadAndLag(j,k) = 0;
                end
            end
            
        end
    end
    
    LeadAndLag
    sum_C
    sum_delay
   
end
  %}  

end

function [following, following2] = correlation_withCenter(handles, hObject, start_frame, end_frame)

%correlation_withCenter compares all trajectories of nodes with that of the center;
%it applies similar ideas

lb = handles.lb;
ub = handles.ub;

%select the two parameters, a delay and a threshold
delay = 40;
threshold = 0.8;

[positionMatrix, center, numberOfNodes] = Read_coordinatesAndCenter(lb, ub, handles, hObject);
[theta0, theta, theta2, theta3, theta4, theta_avg, theta_avg2] = direction(handles, hObject);
[center, center_avg, center2, dir_center, dir_node, dir_node2] = center_dir(hObject, handles, delay, lb, ub);

%following2 keeps track of the leading node at every frame, while following
%makes an average
following = zeros(1, numberOfNodes);
following2 = zeros(end_frame-start_frame+1, numberOfNodes);

for o =1:numberOfNodes

    for i=start_frame:end_frame
        
        %before frame i
        max_correl1(o,i-start_frame+1) = -1;
        for j=i:-1:max(i-delay, lb)
            %angle wrt direction of center; cos(angle) = dot product of
            %direction vector of node and that of center
            angle = anglerestrict(theta4(i-lb+1, o)-dir_center(j-lb+1));  
            %maximum dot product over a period of delay frames
            if cos(angle) > max_correl1(o,i-start_frame+1)   %when center is before which means center leads and node follows
                max_correl1(o,i-start_frame+1) = cos(angle);
                delay_correl1(o,i-start_frame+1) = i-j;
            end
        end
        
        %after frame i, same as before, only frames considered AFTER frame
        %i
        max_correl2(o,i-start_frame+1) = -1;
        for j=i:min(i+delay, ub)
            angle = anglerestrict(theta4(i-lb+1, o)-dir_center(j-lb+1));
            if cos(angle) > max_correl2(o,i-start_frame+1) %when center is after so center follows so node leads
                max_correl2(o,i-start_frame+1) = cos(angle);
                delay_correl2(o,i-start_frame+1) = i-j;
            end
        end
    
        if max_correl2(o, i-start_frame+1)>threshold
            following2(i-start_frame+1, o)=1; %leader
        elseif max_correl1(o, i-start_frame+1)>threshold
            following2(i-start_frame+1, o)=2; %follower
        end
    
    end
    
    %taking the mean of the dot products
     max_correl_avg_past(o) = mean(max_correl1(o,:));
     max_correl_avg_future(o) = mean(max_correl2(o,:));
     
     %mean of delays, which is the average time the node was correlated to the
     %direction of center before/after frame i
     avg_delay_past(o) = mean(delay_correl1(o,:));
     avg_delay_future(o) = mean(delay_correl2(o,:));
     
     %is the mean of correlations higher than some threshold?
    if max_correl_avg_future(o)>threshold
            following(o)=1;
            avg_delay(o) = -avg_delay_future(o);
    elseif max_correl_avg_past(o)>threshold
            following(o)=2;     
            avg_delay(o) = -avg_delay_past(o);
    else
        avg_delay(o) = -mean([avg_delay_past(o), avg_delay_future(o)]);
    end
    
end
     
%if delays are high, then node was ahead of the center by a high value, which
%makes it be in a so-called 'leading' position
[ignore, nodes_ordered]=sort(avg_delay);
ignore
nodes_ordered

end

function [following, following2, LeadAndLag2] = lead_and_lag2(hObject, handles, start_frame, end_frame)

lb = handles.lb;
ub = handles.ub;
delay =25;
threshold = 0.7;

[positionMatrix, center, numberOfNodes] = Read_coordinatesAndCenter(lb, ub, handles, hObject);
[theta0, theta, theta2, theta3, theta4, theta_avg, theta_avg2] = direction(handles, hObject);
following = zeros(numberOfNodes, numberOfNodes);
following2 = zeros(end_frame-start_frame+1, numberOfNodes, numberOfNodes);

for o =1:numberOfNodes
 for s = 1:numberOfNodes
     
   if o~=s  
    for i=start_frame:end_frame
        
        %before
        %finding the maximum dot product between the trajectory of node s
        %and node o
     
        max_correl1(o,s,i-start_frame+1) = -1;
        for j=i:-1:max(i-delay, lb)
            angle = anglerestrict(theta4(j-lb+1,s)-theta4(i-lb+1,o));
            if cos(angle) > max_correl1(o,s,i-start_frame+1)  %when s is leading and o is following
                max_correl1(o,s,i-start_frame+1) = cos(angle);
                delay_correl1(o,s,i-start_frame+1) = i-j;
            end
        end
        
        %after
        max_correl2(o,s,i-start_frame+1) = -1;
        for j=i:min(i+delay, ub)
            angle = anglerestrict(theta4(j-lb+1,s)-theta4(i-lb+1,o));
            if  cos(angle) > max_correl2(o,s,i-start_frame+1)  %when s is following and o is leading
                max_correl2(o,s,i-start_frame+1) = cos(angle);
                delay_correl2(o,s,i-start_frame+1) = i-j;
            end
        end
 
        % leader/follower on each frame
        if max_correl1(o, s, i-start_frame+1)>threshold
            following2(i-start_frame+1, o, s)=1;  % o is following s
        elseif max_correl2(o, i-start_frame+1)>threshold
            following2(i-start_frame+1, o, s)=2;  % s is following o
        end
    
    end
    
    %compute averages over all frames
     max_correl_avg_past(o,s) = mean(max_correl1(o,:));
     max_correl_avg_future(o,s) = mean(max_correl2(o,:));
     
     avg_delay_past(o) = mean(delay_correl1(o,:));
     avg_delay_future(o) = mean(delay_correl2(o,:));
     
     %leader/follower relation if above some threshold
    if max_correl_avg_past(o,s)>threshold
            following(o,s)=1;  %o is following s
            avg_delay(o) = avg_delay_past(o);
    elseif max_correl_avg_future(o,s)>threshold
            following(o,s)=2;  %s is following o
            avg_delay(o) = avg_delay_future(o);
    else
            avg_delay(o) = mean([avg_delay_past(o), avg_delay_future(o)]);
    end
 
    end
    
   end
end

max_correl_avg_past
max_correl_avg_future
following

LeadAndLag2 = zeros(numberOfNodes, numberOfNodes);

for o=1:numberOfNodes
    for s=o+1:numberOfNodes

        if following(o,s)~=following(s,o)
            if following(o,s)==1
                LeadAndLag2(o,s)=1; %o is following s
            elseif following(o,s)==2
                LeadAndLag2(s,o) = 1; %s is following o
            end
        end
    end
end

LeadAndLag2
gObj = biograph(LeadAndLag2);
view(gObj);        

end

%establishes how many main changes of direction (i.e. when node has shifted
%direction more than pi/6) there are

function direction = significant_difference(handles, hObject, start_frame, end_frame)

[positionMatrix, center, numberOfNodes] = Read_coordinatesAndCenter(lb, ub, handles, hObject);
[theta0, theta, theta2, theta3, theta4, theta_avg, theta_avg2] = direction(handles, hObject);
delay = 10;

for j=1:numberOfNodes

i = start_frame;
k = start_frame+1;
cnt =1;

 while k <= end_frame 
    
        if theta4(i-start_frame+1, j)-theta4(k-start_frame+1, j) > pi/6
            direction(j, cnt) = k; %at what frame the node changes direction
            cnt = cnt+1;
            i=k;
        else
             k=k+1;
        end
        
 end
end

end

%computes distance from center of group
function [r, numberOfNodes] = distToCenter2(handles, hObject, start_frame, end_frame)

      file=strcat('Dance',num2str(start_frame),'.dat');
      data_type = handles.data_type
      if data_type==1
        data=handles.dancedata{start_frame};
        set(handles.dataselection_editText,'String',start_frame);
        x = data(:,1);
        y = data(:,2);
        theta = data(:,3);
        %[x y theta] = transf_in_familiar_coord(x,y,theta);        
        dim = 2;
        positionMatrix = horzcat(x,y);    
        center_group = [mean(x), mean(y)];
      end
      
      numberOfNodes = length(x);

  for i=start_frame:end_frame
      
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
        center_group = [mean(x), mean(y)];
      end
 
      %after reading all data, distance from center is computed
      for j= 1:numberOfNodes
         r(i-start_frame+1,j) = norm(positionMatrix(j,:) - center_group);
      end
  end
  
end

%moments of inertia of the group
function moments = MomentsOfInertia(r, start_frame, end_frame, numberOfNodes)

   for i=1:end_frame-start_frame+1
       for j=1:5
           
           switch j
               case 1
                   moments(i,j)=mean(r(i,:)); % mean of distance from center of the group is first moment
               case 2
                   moments(i,j)=mean((r(i,:) - mean(r(i,:))).^2); %mean of distance squared
               case 3
                   moments(i,j)=mean((r(i,:) - mean(r(i,:))).^3); 
               case 4
                   moments(i,j)=mean((r(i,:) - mean(r(i,:))).^4);
               case 5
                   moments(i,j)=mean((r(i,:) - mean(r(i,:))).^5);
           end
       end
   end
   
end

%the pushbutton of moments of inertia
function momentsOfInertia_pushbutton_Callback(hObject, ~, handles)

    prompt = {'Enter start frame to analyze lead and lag network for:', 'Enter end frame'};   
    answer = inputdlg(prompt, 'Select frames');
    start_frame = str2num(answer{1});
    end_frame = str2num(answer{2});
 
    [r, numberOfNodes] = distToCenter2(handles, hObject, start_frame, end_frame);
    moments = MomentsOfInertia(r, start_frame, end_frame, numberOfNodes);
    
    figure(1);
    plot(start_frame:end_frame, moments(:,1), 'r', start_frame:end_frame, moments(:,2), 'b', start_frame:end_frame, moments(:,3), 'y', start_frame:end_frame, moments(:,4), 'c', start_frame:end_frame, moments(:,5), 'k')
    legend('Mean', 'Variance', 'moment 3', 'moment 4', 'moment 5');
    
end

%compute polarization of each node 
function [Pol, Polreal, Polimg] = polariz_moment(m, theta)

numberOfNodes= size(theta,2);

for n=1:numberOfNodes
   Polreal(n) = cos(theta(n));
   Polimg(n) = sin(theta(n));
end

Pol = (sum(Polreal)+i*sum(Polimg))/numberOfNodes;
Polreal=Polreal';
Polimg=Polimg';

end

%formula for Laplacian
function L = Laplacian(A)
   
 D = diag(sum(A,2));
 L = D - A;

end

%computing polarization, laplacian phase potential, angular momentum
function [W, Pol, angular_m] = polariz_momentum_shape_pushbutton_Callback(hObject, ~, handles)

m=1;
lb = handles.lb;
ub = handles.ub;

choose_graph(hObject, handles,0);
adj = handles.adjacencyMatrix;

ub = size(adj,1);

    for i=lb:ub %read the data
    
      file=strcat('Dance',num2str(i),'.dat');
      data_type = handles.data_type;
      if data_type==1
        data=handles.dancedata{i};
        set(handles.dataselection_editText,'String',i);
        x(i,:) = data(:,1);
        y(i,:) = data(:,2);
        theta(i,:) = data(:,3);
        %[x(i,:), y(i,:), theta(i,:)] = transf_in_familiar_coord(x(i,:),y(i,:),theta(i,:));        
        dim = 2;
        positionMatrix(i,:,:) = horzcat(x(i,:)',y(i,:)');    
        centerMatrix(i,:) = [mean(x(i,:)), mean(y(i,:))];
      end
    
      numberOfNodes = size(theta(i,:),2);
      [Pol(i), Polreal, Polimg] = polariz_moment(m, theta(i,:)); 
      Evector = [Polreal,Polimg]; %polarization
      %adj = handles.adjacencyMatrix;
      A = squeeze(adj(i,:,:));
      L = Laplacian(A);
      Evector2(:,1) = L*Evector(:,1); %Laplacian * polarization
      Evector2(:,2) = L*Evector(:,2);
      
      dotProd = dot(Evector(:,1),Evector2(:,1)) + dot(Evector(:,2),Evector2(:,2));
      W(i) = 1/(2*m*numberOfNodes)*sum(dotProd); % Laplacian phase potential (known formula)
 
      Pol_norm(i) = norm(Pol(i));
    
      r_ic = squeeze(positionMatrix(i,:,:)) - [ones(numberOfNodes,1)*centerMatrix(i,1), ones(numberOfNodes,1)*centerMatrix(i,2)]; %corrdinates given that center of the group would be origin
      r_ic = [zeros(numberOfNodes,1), r_ic];
      v = [zeros(numberOfNodes,1),Evector];
      
      sum_cross=[0,0,0]; 
      for j=1:numberOfNodes
        CrossPr(j,:,:,:) = cross(r_ic(j,:,:), v(j,:,:));
        sum_cross = sum_cross + squeeze(CrossPr(j,:,:,:));
      end
      
      angular_m(i) = norm(sum_cross)/numberOfNodes; %angular momentum
      
    end
    
    [W_avg, ignore, ignore, ignore, ignore] = averaging(W, 20, 1, 1);
    [Pol_norm_avg, ignore, ignore, ignore, ignore] = averaging(Pol_norm, 20, 1, 1);
    [angular_m_avg, ignore, ignore, ignore, ignore] = averaging(angular_m, 20, 1, 1);
 
 save('D:\WorkStuff\FlockLogic\version 11\W4_SG1.mat', 'W_avg')
 save('D:\WorkStuff\FlockLogic\version 11\Pol_norm_SG1.mat', 'Pol_norm_avg')
 save('D:\WorkStuff\FlockLogic\version 11\angular_m_SG1.mat', 'angular_m_avg')
 
    figure(1)
    plot(lb:ub, W, 'r', lb:ub, W_avg, 'b')
    title('Laplacian phase potential')
    xlabel('frame')
    ylabel('W')
    legend('Laplacian phase potential', 'average LPP')
    
    disp('mean W')
    mean(W)
    
    %Pol_norm and angular_m are not dependent on the sensing graph!, but W is
    
    figure(2)
    plot(lb:ub, Pol_norm, 'r', lb:ub, Pol_norm_avg, 'b')
    title('phase order parameter p, measuring synchrony')
    xlabel('frames')
    ylabel('p')
    legend('Phase order', 'average phase order')
    
    disp('mean Pol')
    mean(Pol_norm)
    
    figure(3)
    plot(lb:ub, angular_m, 'r', lb:ub, angular_m_avg, 'b')
    title('angular momementum')
    xlabel('frames')
    ylabel('m')
    legend('angular momentum', 'average angular momentum')
    
end

%speed, accelration
function speed_acc_or_pushbutton_Callback(hObject, ~, handles)

lb=handles.lb;
ub = handles.ub;

S = [{'1'},{'2'},{'3'},{'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'9'}, {'10'}, {'11'}, {'12'}, {'13'}];
node = listdlg('ListString', S, 'Name', 'Speed-acc-orientation graph', 'PromptString', 'Select the node(s)');

if size(node,2)>1
    msgbox('Please select only one node!')
    return
end

FPS=20;
i=lb+1;
cnt=1;
for i=lb+1:ub
    
    file=strcat('Dance',num2str(i-1),'.dat');
      data_type = handles.data_type;
      if data_type==1
        data=handles.dancedata{i-1};
        set(handles.dataselection_editText,'String',i-1);
        x = data(:,1);
        y = data(:,2);
        theta1 = data(:,3);
        %[x y theta1] = transf_in_familiar_coord(x,y,theta1);        
        x = x(node); y=y(node); theta1=theta1(node);
        dim = 2;
        positionMatrix1(i-1,:) = horzcat(x,y);    
      end
      
      file=strcat('Dance',num2str(i),'.dat');
      data_type = handles.data_type;
      if data_type==1
        data=handles.dancedata{i};
        set(handles.dataselection_editText,'String',i);
        x = data(:,1);
        y = data(:,2);
        theta2 = data(:,3);
        %[x y theta2] = transf_in_familiar_coord(x,y,theta2);        
        x = x(node); y=y(node); theta2=theta2(node);
        dim = 2;
        positionMatrix(i,:) = horzcat(x,y);    
      end
     
      %direction of nodes normalized
      %velocity and speed
       dr(i-lb,:) = positionMatrix(i,:) - positionMatrix(i-1,:);
       dr_norm(i-lb) = norm(dr(i-lb,:));
       dtheta(i-lb) = anglerestrict(theta1-theta2);
       
       if i>lb+1
           %acceleration
           dv(i-lb-1,:) = dr(i-lb,:)-dr(i-lb-1,:);
           dv_norm(i-lb-1) = norm(dv(i-lb-1,:));
       end
   
end

choose_graph(hObject, handles,0);
adjacencyMatrix = handles.adjacencyMatrix;

neighbors = find(adjacencyMatrix(lb,node,:)==1);

%changeMx is a matrix describing how many frames and what neighbors it has
neighborNumber=1;
changeMx(1,neighborNumber)=lb;
changeMx(2,neighborNumber)=1;
changeMx(3,neighborNumber)=neighbors(1);
changeMx(4,neighborNumber)=neighbors(2);

ub = size(adjacencyMatrix,1)
for i=lb+1:ub
    neighbors2 = find(adjacencyMatrix(i,node,:)==1);
    if sum(neighbors == neighbors2)== 2
        changeMx(2, neighborNumber)=changeMx(2, neighborNumber)+1;
    else
        neighborNumber = neighborNumber+1;
        changeMx(1, neighborNumber) = i;
        changeMx(2, neighborNumber)=1;
        changeMx(3, neighborNumber)=neighbors(1);
        changeMx(4, neighborNumber)=neighbors(2);
    end
    neighbors=neighbors2;
end

 figure(1)
 plot(lb+9:ub, dr_norm(9:ub-lb), 'b');
 hold on;
 plot(lb+10:ub, dv_norm(9:ub-lb-1), 'r')
 hold on;
 %ymin = min(min(dr_norm), min(dv_norm));
 %ymax = max(max(dr_norm), max(dv_norm));
ymin=-0.1; ymax=0.2;
 
 %legend('speed','acceleration')
 
 for i=1:neighborNumber
  xcoord(i) = changeMx(1,i);
  s = size(ymin:0.01:ymax,2);
  plot(xcoord(i)*ones(1,s), ymin:0.01:ymax, 'g')
  hold on;
 end
 
 title('Speed, Acceleration')
 xlabel('frames')
 ylabel('speed, acc')
 legend('Speed', 'Acceleration', 'Neighbor change')
 
 figure(2)
 plot(lb+1:ub, dtheta(1:ub-lb), 'k');
 hold on;
 ymin = min(dtheta);
 ymax = max(dtheta);
 
 for i=1:neighborNumber
  xcoord(i) = changeMx(1,i);
  s = size(ymin:ymax,2);
  plot(xcoord(i)*ones(1,s), ymin:ymax, 'g')
  hold on;
 end
 
 title('Changing theta')
 xlabel('frames')
 ylabel('dtheta')
 
end

%distances of node to center, to other nodes, to neighbors, to closest nodes
function Individual_DistanceAndMoments_pushbutton_Callback(hObject, ~, handles)

S = [{'1'},{'2'},{'3'},{'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'9'}, {'10'}, {'11'}, {'12'}, {'13'}];
nodes = listdlg('ListString', S, 'Name', 'Individual s distances', 'PromptString', 'Select the node(s)');

prompt = {'Enter start frame to analyze Distances of Individual Nodes for:', 'Enter end frame'};
answer1 = inputdlg(prompt, 'Select frames');
start_frame = str2num(answer1{1});
end_frame = str2num(answer1{2});

prompt = {'Select options of what to do:                                             Select 1). Distance of node to center                                                                   Select 2). Distance of node to other nodes                                        Select 3). Distance of node to its to closest node(s)                                            Select 4). Distance of node to its neighbors' };
answer2 = cell2mat(inputdlg(prompt));
answer2 = str2num(answer2);

switch answer2
    
    case 1
        %distance of node to center is represented by first moment of
        %inertia
        
        [r, numberOfNodes] = distToCenter2(handles, hObject, start_frame, end_frame);
        moments = MomentsOfInertia(r, start_frame, end_frame, numberOfNodes);

        %comparing distance from center of node as compared to other nodes
        
        figure(1)
        plot(start_frame:end_frame, moments(:,1), 'b', 'LineWidth', 0.8);
        hold on;
        for i=1:size(nodes,2)
           plot(start_frame:end_frame, r(1:end_frame-start_frame+1, nodes(i)), 'r', 'LineWidth', 0.7);
           hold on;
        end
        
        legend('Avg distance to center', 'Single node dist to center')
        title('Comparison of first moment of inertia w/ single node dist to center')
        xlabel('frames')
        ylabel('Dist to center')
        
      disp('Average distances to center are:');  
      for i=1:size(nodes,2)
        disp(nodes(i))  
        disp(sum(r(1:end_frame-start_frame+1, nodes(i)))/(end_frame-start_frame+1))
        distToCenter_avg(i) = mean(r(:, nodes(i)));
      end
      
      [sorted_distToCenter, ind_distToCenter] = sort(distToCenter_avg);
      disp('hierarchy of distances to center')
      ind_distToCenter
      
    case 2
    
    for i=start_frame:end_frame
      
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
        positionMatrix(i,:,:) = horzcat(x,y);
    
      end
      
    numberOfNodes = length(squeeze(positionMatrix(i,:,:)));
    distanceMatrix(i,:,:) = ipdm(squeeze(positionMatrix(i,:,:)));
    
    for j=1:size(nodes,2)
      distToOthers(j,i) = (sum(distanceMatrix(i,nodes(j),:)))/(numberOfNodes-1);    %average dist to others at each frame
    end
    
    end
    
    disp('Average distance to all other nodes is:')
    for j=1:size(nodes,2)
        disp(nodes(j));
        avg_distToOthers(j) = mean(distToOthers(j,:));
        disp(avg_distToOthers(j));
    end
    
    [r, numberOfNodes] = distToCenter2(handles, hObject, start_frame, end_frame);
   
    figure(1)
   for j=1:size(nodes,2) 
    plot(start_frame:end_frame, distToOthers(j,1:end_frame-start_frame+1), 'r', start_frame:end_frame, r(1:end_frame-start_frame+1,j), 'b');
    hold on;
   end
   
   title('Comparing mean distance to others to distance to center')
   xlabel('frames')
   ylabel('distance')
   legend('Mean distance to others', 'Distance to center')
   
    case 3
        
        for i=start_frame:end_frame
      
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
                positionMatrix(i,:,:) = horzcat(x,y);
    
            end
      
       numberOfNodes = length(squeeze(positionMatrix(i,:,:)));
       distanceMatrix(i,:,:) = ipdm(squeeze(positionMatrix(i,:,:)));

       for j=1:size(nodes,2)
         dist(i,j,:) = min(distanceMatrix(i, nodes(j),:),2); %identify closest nodes
         dist2(i,j) = mean(dist(i,j,:)); %mean distance to closest nodes
       end
       
        end
      
     disp('Average distance to nearest nodes');
      for j=1:size(nodes,2)
          nodes(j)
          dist_avg(j) = mean(dist(:,j,1))
          disp(dist_avg(j))
      end
      
      figure(1);
      for j=1:size(nodes,2)
          plot(start_frame:end_frame, dist2(start_frame:end_frame,j));
          hold on;
      end
      
    case 4
   
        choose_graph(hObject, handles,0);
        adjacencyMatrix = handles.adjacencyMatrix;

        for i=start_frame:end_frame
      
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
                positionMatrix(i,:,:) = horzcat(x,y);
    
            end
      
        numberOfNodes = length(squeeze(positionMatrix(i,:,:)));
        distanceMatrix(i,:,:) = ipdm(squeeze(positionMatrix(i,:,:)));
        
        end
        
        for j=1:size(nodes,2)
            for i=start_frame:end_frame
                nodes2 = find(adjacencyMatrix(i,nodes(j),:)==1,2); %2 neighbors
                distNeigh(j,i) = mean(distanceMatrix(i,j, nodes2)); %distance to neighbors
                distNeigh2(j,i,:) = distanceMatrix(i,j, nodes2);
            end
        end
        
        figure(1)
        for j=1:size(nodes,2)
            plot(start_frame:end_frame, distNeigh(j,start_frame:end_frame))
            hold on;
        end
end

end

%distance to neighbors: to each and average
function [distNeigh, distNeigh_avg] = dist_to_neighb(hObject, handles)

    prompt = {'Enter start frame to analyze Distances of Individual Nodes for:', 'Enter end frame'};
    answer1 = inputdlg(prompt, 'Select frames');
    start_frame = str2num(answer1{1});
    end_frame = str2num(answer1{2});

    [positionMatrix, center, numberOfNodes] = Read_coordinatesAndCenter(start_frame, end_frame, handles, hObject);
    numberOfNodes = handles.numberOfNodes;
    
   for i=start_frame:end_frame 
    
       %matrix of distances from one node to the other
       distanceMatrix(i,:,:) = ipdm(squeeze(positionMatrix(i,:,:)));
       
       for j=1:numberOfNodes
           nodes = find(adjacencyMatrix(i,j,:)==1,2);
           %average distance to neighbors
           distNeigh_avg(j,i) = mean(distanceMatrix(i,j,nodes));
           %distance to neighbors
           distNeigh(j,i,:) = distanceMatrix(i,j,nodes);
       end
   end
    
end

%dot product with center of group OR simply form group of synchronized
%nodes
function synchrony_pushbutton_Callback(hObject, ~, handles)

lb = handles.lb;
ub = handles.ub;

prompt = {'Select 1). Synchrony w.r.t. center of nodes                                                                       Select 2). Synchrony groups'};
answer = str2num(cell2mat(inputdlg(prompt, 'Select analysis')));

if answer == 1
S = [{'1'},{'2'},{'3'},{'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'9'}, {'10'}, {'11'}, {'12'}, {'13'}];
nodes = listdlg('ListString', S, 'Name', 'Individual s distances', 'PromptString', 'Select the node(s)');

prompt = {'Select start frame:','Select end frame:'}
answer = inputdlg(prompt, 'Select range of frames');
start_frame = str2num(answer{1})
end_frame = str2num(answer{2})

%read data and store
for i=start_frame+1:end_frame
    
          file=strcat('Dance',num2str(i-1),'.dat');
          data_type = handles.data_type;
            if data_type==1
                data=handles.dancedata{i-1};
                set(handles.dataselection_editText,'String',i-1);
                x = data(:,1);
                y = data(:,2);
                theta = data(:,3);
        
                %[x y theta] = transf_in_familiar_coord(x,y,theta);        
                dim = 2;
                positionMatrix(i-1,:,:) = horzcat(x,y);
                centerMatrix(i-1,:) = [mean(x), mean(y)];
            end
            
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
                positionMatrix(i,:,:) = horzcat(x,y);
                centerMatrix(i,:) = [mean(x), mean(y)];
            end
    
            %direction of the center normalized
              vector_c = centerMatrix(i,:)-centerMatrix(i-1,:);
              vector_c = vector_c/norm(vector_c);
              
              %for each node, dot product of direction vector of each node with
              %direction of center
            for j=1:size(nodes,2)
              vector = squeeze(positionMatrix(i,nodes(j),:)) - squeeze(positionMatrix(i-1,nodes(j),:));
              vector = vector/norm(vector);
              dotProd(i-start_frame, j) = dot(vector', vector_c);
            end

end

            figure(1);
            for j=1:size(nodes,2)
                plot(start_frame+1:end_frame, dotProd(1:end_frame-start_frame, j));
                hold on;
            end
else
    
    prompt = {'Select start frame:','Select end frame:'};
    answer = inputdlg(prompt, 'Select range of frames');
    start_frame = str2num(answer{1});
    end_frame = str2num(answer{2});

    group = synchronized_groups(handles, hObject, start_frame, end_frame);
    s = size(group,2);
    
    %form groups and display
    for i=start_frame:end_frame  
        if s>=2
        if sum(squeeze(group(i-start_frame+1,2,:)))~=0
            squeeze(group(i-start_frame+1,:,:))
        end
        end
    end
    
    end
end

%include nodes in groups according to whether they move in roughly the same
%direction
function group = synchronized_groups(handles, hObject, start_frame, end_frame)

lb = handles.lb;
ub = handles.ub;
numberOfNodes = handles.numberOfNodes;

%parameter, how far away in the past to look
tau = 15;
%parameter defining 'almost the same' direction
c_min = 0.99;

%read data and store
for i=max(lb,start_frame-tau):min(ub,end_frame+tau)
    
            file=strcat('Dance',num2str(i),'.dat');
            data_type = handles.data_type;
            data=handles.dancedata{i};
            x = data(:,1);
            y = data(:,2);
            theta = data(:,3);
        
                %[x y theta] = transf_in_familiar_coord(x,y,theta);        

            dim = 2;
            positionMatrix(i,:,:) = horzcat(x,y);
end
start_frame
for i=start_frame:end_frame
    for o=1:numberOfNodes-1
        for s=o+1:numberOfNodes
            
            for x = max(i-tau, lb):min(i+tau, ub)
                pos_o = squeeze(positionMatrix(x,o,:))/norm(squeeze(positionMatrix(x,o,:)));
                pos_s = squeeze(positionMatrix(x,s,:))/norm(squeeze(positionMatrix(x,s,:)));
                dotProd(x-max(i-tau,lb)+1) = dot(pos_o, pos_s);
            end
  
            if max(dotProd)>=c_min
                adj(o,s)=1;
                adj(s,o)=1;
            end
        end       
    end

      queue = [1];
      cnt=0; %numberOfNodes
      cnt2=1; %group number
      
      verify = zeros(1,numberOfNodes); 
      verify(1)=1;
      change = 1;
      index=1;
      
      while cnt <= numberOfNodes %while not all nodes have been included in a group
       
        if ~isempty(queue) %are there nodes left that could be part of the group cnt2?
           
          node = queue(1);
          group(i-start_frame+1,cnt2,index) = node; %include node in group, i.e. include the first element of queue in the group
          
          cnt=cnt+1;
          queue(1) = []; 
          neighbors1 = find(adj(node,:)==1); neighbors2 = find(adj(:,node)==1); %nodes that node follows or is followed by
          int = ismember(neighbors1, neighbors2); neighbors1(int) = []; %not to take nodes twice if followed and is followed by, then delete from one vector
          
          neighbors = [neighbors1, neighbors2]; % paste neighbors 1,2 together
          neighbors = neighbors(find(verify(neighbors)==0)); %nodes not yet put together in a group
          
          queue = [queue, neighbors'];
          verify(neighbors)=1; %nodes already in the group
          index = index+1;
          
        else
            
          queue = [find(verify==0,1)];  
          cnt2 = cnt2+1; %increase number of groups since no other node is added to the subgroup
          index = 1; %the new subgroup contains 0 members and the next member will be the first member of the group
          
          if cnt==numberOfNodes %to exit while loop
              cnt = cnt+1;
          end
          
        end
        
      end
      
end
end

%selecting nodes most probably to be neighbors; used in the
%circle-intersecttion model
%adj1: remembers all the max dot products between directions of vectors
%adj2: remembers for each pair of nodes 1, if nodes are headed in
%roughly the same direction; 2, if the nodes head in almost the same
%direction
%adj3: 1, if the direction is roughly the same and the distance is close (close top_neighbor=4)
%2, if the direction is almost the same and the distance is close (close top_neighbor=4)
function [adj1, adj2, adj3] = simple_synchrony(handles, hObject, start_frame, end_frame)

lb = handles.lb;
ub = handles.ub;
numberOfNodes=handles.numberOfNodes;

c_min1=0.9; c_min2 = 0.5;
tau=20;

[theta0, theta, theta2, theta3, theta4, theta_avg, theta_avg2] = direction(handles, hObject);

%read data and store it in a matrix
    for i=max(lb,start_frame-tau):min(ub,end_frame+tau)
    
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
                positionMatrix(i,:,:) = horzcat(x,y);
                distanceMatrix(i,:,:) = ipdm(squeeze(positionMatrix(i,:,:)));
            end
    end

top_neighbor=4;

for i=start_frame:end_frame
    
    adj1{i} = zeros(numberOfNodes);
    adj2{i} = zeros(numberOfNodes);
    adj3{i} = zeros(numberOfNodes);
    
    for o=1:numberOfNodes-1
        
          adj1{i}(o,o)=Inf;
          adj2{i}(o,o)=Inf;
          
          %sort distanceMatrix (matrix of distances between nodes) because
          %closer nodes have more priority than nodes farther away
          positionMatrix(i,o,o)=Inf;
          [ignore, posIndex] = sort(distanceMatrix(i,o,:));
          
        for s=o+1:numberOfNodes
            
            %compute dot product between the directions of the nodes
            %normalized
            for x = max(i-tau, lb):min(i+tau, ub)
                dotProd(x-max(i-tau,lb)+1) = cos(anglerestrict(theta4(i-start_frame+1,o)- theta4(i-start_frame+1,s)));
            end
  
            adj1{i}(o,s) = max(dotProd);
            adj1{i}(s,o) = max(dotProd);
            
            if max(dotProd)>=c_min1
                adj2{i}(o,s)=2;
                adj2{i}(s,o)=2;
                
              if find(posIndex==s)<=top_neighbor
                  adj3{i}(o,s)=2;
              end
              
            elseif max(dotProd)>=c_min2
                adj2{i}(o,s)=1;
                adj2{i}(s,o)=1; 
              if find(posIndex==s)<=top_neighbor
                  adj3{i}(o,s)=1;
              end
              
            end
        end       
    end

end

end

% direction of moveent, averaging in different ways
function [theta0, theta, theta2, theta3, theta4, theta_avg, theta_avg2] = direction(handles, hObject)

    start_frame = handles.lb;
    end_frame = handles.ub;
    delay = 5;
    delay2 = 5;
    
    for i=start_frame:end_frame
           data=handles.dancedata{i};
           x(i,:) = data(:,1);
           y(i,:) = data(:,2);
           theta0(i-start_frame+1,:) = data(:,3);
    end

    numberOfNodes = length(theta0(1,:));
    
    for i=start_frame:end_frame-1
        %from frame to frame difference in positions
        xx = x(i+1,:)-x(i,:);
        yy = y(i+1,:)-y(i,:);
        v(i,:,:) = [xx', yy'];
        
     for j=1:numberOfNodes   
         %angle given by the direction of movement of each node frame by
         %frame
           theta(i-start_frame+1,j) = atan2(v(i,j,2), v(i,j,1));
     end   
     
     %adding some delay to obtain a difference in positions, so that frame
     %by frame irregularity in direction is avoided. Not sure it is useful
     if i<=end_frame-delay
        xx = (x(i+delay,:)-x(i,:));
        yy = (y(i+delay,:)-y(i,:));
        v2(i,:,:) = [xx', yy']; 
     
        %angle
     for j=1:numberOfNodes   
           theta2(i-start_frame+1,j) = atan2(v2(i,j,2), v2(i,j,1));
     end
     
     end
     
   end
    
    theta(end_frame-start_frame+1,:) = theta(end_frame-start_frame,:);
    
 for j=1:numberOfNodes   
    theta2(end_frame-delay+1-start_frame+1:end_frame-start_frame+1,j) = theta2(end_frame-delay-start_frame+1,j);
 end
 
   for j=1:numberOfNodes 
    [theta3(:,j), ignore, ignore, theta_avg(:,j), ignore] = averaging(theta(:,j), delay2, 1, 1);
    [theta4(:,j), ignore, ignore, theta_avg2(:,j), ignore] = averaging(theta2(:,j), delay2, 1, 1);
   end 

   S = [{'1'},{'2'},{'3'},{'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'9'}, {'10'}, {'11'}, {'12'}, {'13'}];
   nodes = listdlg('ListString', S, 'Name', 'Node status & correlations', 'PromptString', 'Select the node(s)');

   %plot all 4 options of angles
  figure(1) 
    plot(start_frame:end_frame, theta0(:,nodes(1)), 'y-*')
    hold on
    plot(start_frame:end_frame, theta(:,nodes(1)), 'b-*')
    hold on
    plot(start_frame:end_frame, theta2(:,nodes(1)), 'r-*')
    hold on
    plot(start_frame:end_frame, theta3(:,nodes(1)), 'g-*')
    hold on
    plot(start_frame:end_frame, theta4(:,nodes(1)), 'k-*')
    
end

% want for each frame a measure of nodes' likelihood to be compatible:
% given by both distance to each other AND difference in direction of movement
function [F, Ftotal] = flockingIndex(handles, hObject, start_frame, end_frame)

[positionMatrix, center, numberOfNodes] = Read_coordinatesAndCenter(start_frame, end_frame, handles, hObject);
[theta0, theta, theta2, theta3, theta4, theta_avg, theta_avg2] = direction(handles, hObject);

for i=1:end_frame-start_frame+1
    for o=1:numberOfNodes
        for s=1:numberOfNodes
            delta = anglerestrict(theta_avg2(i+start_frame-1,o) - theta_avg2(i+start_frame-1,s)); 
            dist_norm = norm(positionMatrix(i,o,:)-positionMatrix(i,s,:));
            F(i,o,s) = (1-delta/pi)*1/(1+exp(-2*(dist_norm-1/2))); % a function of both distance and orientation
            Ftotal(i) = sum(sum(squeeze(F(i,:,:))))/(numberOfNodes*(numberOfNodes-1)); %normalization
        end
    end
end

end

%compute node status of an influence graph, given by matrix adjacency matrix adj
function [nodeStatus, nodeStatus_avg] = compute_nodeStatus(hObject, handles, adj, start_frame, end_frame)

lb = handles.lb;
ub = handles.ub;

if start_frame<lb | end_frame>ub
    return
end

%adjacencyMatrix = handles.adjacencyMatrix; %despending on whether
%adjMatrix is stored in the handle (if program was run) OR imported from a
%.mat file

adjacencyMatrix = adj;
numberOfNodes = handles.numberOfNodes;

for i=start_frame:end_frame
    
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
        
        %algorithm for computing shortest distances
       for kk = 1:numberOfNodes
            for ii = 1:numberOfNodes
                for jj = 1:numberOfNodes  
                    shortestReversedLengths(ii,jj) = min(shortestReversedLengths(ii,jj), shortestReversedLengths(ii,kk)+shortestReversedLengths(kk,jj));
                end
            end
       end
       
        inverseDistances = 1./shortestReversedLengths;
    
        %formula for node status
        for ii = 1:numberOfNodes
            nodeStatus(i-start_frame+1,ii) = 1/(numberOfNodes - 1)*sum(inverseDistances(ii,[1:(ii-1) (ii+1):numberOfNodes]));
        end
        
end

        for ii = 1:numberOfNodes
            nodeStatus_avg(ii) = mean(nodeStatus(:,ii));
        end
       
end

%putting together all the graphs after computing ns IG, each is 25 frames
%long, at pasting everything together
function [nodeStatus_IG, nodeStatus_avg2_IG, degree_in, degree_out] = nodeStatusIG(handles, hObject, start_frame, end_frame)

name = 'D:\WorkStuff\FlockLogic\version 8\IG\.mat files\IG_25-40.mat';
load(name);

%Interaction2 is matrix of who is following who
for i=25:40
Interaction2(i,:,:) = Interaction;
end

%compute node status of this matrix
[nodeStatus_IG(25:40,:), nodeStatus_avg_IG(1,:)] = compute_nodeStatus(hObject, handles, Interaction2, 25, 40);
[degree_in(25:40,:), degree_out(25:40,:)] = degreeInOut(Interaction2, 25, 40);
cnt = 2;  

%doing the same thing for all other frames
for i=40:20:920
  name = ['D:\WorkStuff\FlockLogic\version 8\IG\.mat files\IG_', num2str(i), '-', num2str(i+20)];
  load(name);
  
  for j=i:i+20
   Interaction2(j,:,:) = Interaction;
  end

  [nodeStatus_IG(i:i+20,:), nodeStatus_avg_IG(cnt,:)] = compute_nodeStatus(hObject, handles, Interaction2, i, i+20);
  cnt = cnt+1;
  [degree_in(i:i+20,:), degree_out(i:i+20,:)] = degreeInOut(Interaction2, i, i+20);
  
end

name = 'D:\WorkStuff\FlockLogic\version 8\IG\.mat files\IG_940-953.mat';
load(name);

%compute everyting for the last frames
for i=940:953
Interaction2(i,:,:) = Interaction;
end

[nodeStatus_IG(940:953,:), nodeStatus_avg_IG(cnt,:)] = compute_nodeStatus(hObject, handles, Interaction2, 940, 953);
[degree_in(940:953,:), degree_out(940:953,:)] = degreeInOut(Interaction2, 940, 953);
%disp('nodeStatus_avg_IG')
%nodeStatus_avg_IG

%compute mean and std
disp('average node status IG')
nodeStatus_avg2_IG = mean(nodeStatus_avg_IG,1);
disp('std of ns IG')
nodeStatus_std2_IG = std(nodeStatus_avg_IG,1);

[center, center_avg, center2, dir_center, dir_node, dir_node2] = center_dir(hObject, handles, 20, 25, 903);

[sorted_ns_IG, sorted_nodes] = sort(nodeStatus_avg2_IG);
disp('sorted nodes by mean')
sorted_ns_IG
sorted_nodes

[Distance2, Distance3, distance, leader, backer, N_leader, N_backer] = LeaderBacker(hObject, handles, 1, 953,1);

figure(3)
barwitherr(nodeStatus_std2_IG, nodeStatus_avg2_IG)

disp('average node status IG')
degree_in_avg2_IG = mean(degree_in,1)
disp('std of ns IG')
degree_in_std2_IG = std(degree_in,1)

[sorted_dgIn_IG, sorted_nodes] = sort(degree_in_avg2_IG)
disp('sorted nodes by mean')
sorted_nodes
disp('sorted nodes by std')
[sorted_std_dgIn_IG, sorted_nodes] = sort(degree_in_std2_IG)
sorted_nodes

figure(4)
barwitherr(degree_in_std2_IG, degree_in_avg2_IG)

end

function L = leadershipIndex(hObject, handles, adj, start_frame, end_frame)

lb = handles.lb;
ub = handles.ub;

if start_frame<lb | end_frame>ub
    return
end

%A = handles.adjacencyMatrix;
A = adj;
numberOfNodes = handles.numberOfNodes;
U = ones(1, numberOfNodes)';

L = zeros(numberOfNodes, end_frame-start_frame+1);

for i=start_frame:end_frame
    
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
                positionMatrix(i,:,:) = horzcat(x,y);
            end
            
    %compute F
    for j=1:numberOfNodes
        for k=1:numberOfNodes
            %using formula from paper - see LeadershipMeasures for name of paper
            F(j,k) = (1-min(abs(theta(j)-theta(k)), 2*pi-abs(theta(j)-theta(k)))/pi)*(1 - 1/(1+exp(-2*norm([x(j)-x(k), y(j)-y(k)])+1)));
        end
    end
        
    a=0.5;
   
     for j=1:1000
       L(:, i-start_frame+1) =L(:, i-start_frame+1)+(a*(squeeze(A(i,:,:)).*F)')^j*U;
     end  
       
    end
    
end

function [xx, exitflag, start_frame, end_frame, delay, nodes] = causality(hObject, handles)

  prompt = {'Select delay:'};
  answer2 = inputdlg(prompt, 'Delay');
  delay = str2num(answer2{1});

  prompt = {'Select start frame:','Select end frame:'};
  answer = inputdlg(prompt, 'Select range of frames');
  start_frame = str2num(answer{1});
  end_frame = str2num(answer{2});

  S = [{'1'},{'2'},{'3'},{'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'9'}, {'10'}, {'11'}, {'12'}, {'13'}];
  nodes = listdlg('ListString', S, 'Name', 'Individual s distances', 'PromptString', 'Select the node(s)');
  
  F = end_frame-start_frame+1;
  
  lb = handles.lb;
  ub = handles.ub;
  
  if start_frame<lb+delay-1 | end_frame>ub
      disp('start frame not chosen properly, mult be >delay')
      return
  end
    
  numberOfNodes = length(nodes);
  %numberOfNodes = handles.numberOfNodes;
            tic;
   %read data of trajectories
   for i=start_frame-delay:end_frame
       
       file=strcat('Dance',num2str(i),'.dat');
       data_type = handles.data_type;
       if data_type==1
          data=handles.dancedata{i};
          set(handles.dataselection_editText,'String', i);
          %traces
          x(:,i-start_frame+delay+1) = data(:,1);
          y(:,i-start_frame+delay+1) = data(:,2);
       end
       
   end
   
   q_final=[];
   B_final = [];
   
   for j=1:length(nodes)
      
       %constructing q
       %position of nodes at F frames
        q_x(1:F,j) = x(nodes(j), F+delay:-1:delay+1)'; %what needs to be used in the problem as q_x and q_y
        q_y(1:F,j) = y(nodes(j), F+delay:-1:delay+1)';
        
       %position of nodes at F frames, with delay
        q_x2(1:F+delay, j) = x(nodes(j), F+delay:-1:1)'; %as before, but with first frames that go only in the delay
        q_y2(1:F+delay, j) = y(nodes(j), F+delay:-1:1)';
        
        q(:,j)= [q_x(:,j); q_y(:,j)];
        q_final = [q_final; q(:,j)];
        
        %constructing B
      for i = F+delay:-1:delay+1
        hankel_x(j, F+delay+1-i, :) = q_x2(F+delay+1-i:F+delay+1-i+delay,j)';
        hankel_y(j, F+delay+1-i, :) = q_y2(F+delay+1-i:F+delay+1-i+delay,j)';
      end
        
      B(j, 1:2*F, 1:(delay+1)) = [squeeze(hankel_x(j,:,:)); squeeze(hankel_y(j,:,:))];
      B_final = [B_final, squeeze(B(j,:,:))]; %for each node B is the same, composed of hankel submatrices
     
   end 
   
   S=3/2;
   lambda = 1;
   eps1 = 0.00001; %for convergence when using fmincon
   eps2 = 0.0001; % as constraint, |eta|_inf <= eps2
   delta = 0.1; 
   
   options = optimoptions('fmincon','Algorithm','interior-point', 'MaxFunEvals', 20000);
   verify = zeros(1, numberOfNodes);
   
  for j = 1:length(nodes)
      
   I = eye(2*F);   %untity matrix of dimensions end_frame-start_frame+1
   Y = [B_final, I];
   
  %equality and inequality constraints 
   A = [Y; -Y]; b = [[q_x(:,j); q_y(:,j)]; [-q_x(:,j); -q_y(:,j)]] + eps2*ones(1, 4*F)';
   Aeq1 = ones(1,numberOfNodes*(delay+1)); Aeq2 = zeros(1,2*F); Aeq = [Aeq1, Aeq2]; 
   beq = 1;
   
   %weights
   weights_a = ones(1, numberOfNodes);
   weights_u = ones(1, 2*F);
   s = ones(1, numberOfNodes); 
   s(j) = S; 
  
   %initial point, where optimization algorithm begins looking
   
   %initial = [zeros(1,delay+1), 1/((numberOfNodes-1)*(delay+1))*ones(1, (numberOfNodes-1)*(delay+1)), zeros(1, 2*F)];
   %initial = [100*ones(1, numberOfNodes*(delay+1)), zeros(1, 2*F)];
   
   initial = [(1/((numberOfNodes-1)*(delay+1)))*ones(1, numberOfNodes*(delay+1)), zeros(1, 2*F)];
   %initial = [zeros(1, numberOfNodes*(delay+1)), zeros(1, 2*F)];
   
   for jj=1:numberOfNodes
      if jj==j 
        initial((jj-1)*(delay+1)+1:jj*(delay+1)) = 0;
      end
   end
   
   %initial = [zeros(1, numberOfNodes*(delay+1)), zeros(1,2*F)];
   
   old_x = Inf;
   xx(j,:) = initial;  
   cnt =1;
   
   %while the solution does not converge...
 while norm(old_x - xx(j,:), Inf) > eps1
     
          old_x = xx(j,:);
          
          %finding x...
          [xx(j,:), fval(nodes(j)), exitflag(j)] = fmincon(@(x) min_for_sparsity(x, weights_a, weights_u, q_x, q_y, Y, numberOfNodes, F, delay, lambda, delta), [old_x], A, b, Aeq, beq, zeros(1,numberOfNodes*(delay+1)+2*F), ones(1,numberOfNodes*(delay+1)+2*F), [], options);          
          
       for jj=1:numberOfNodes   
           weights_a(jj) = 1/(norm(xx(j,(jj-1)*(delay+1)+1:jj*(delay+1)))+delta);
       end
      
       %change weights
       weights_a = weights_a.*s;
       weights_u = 1./(abs(xx(numberOfNodes*delay+1:numberOfNodes*delay+2*F)) +delta);
       
       norm(old_x - xx(j,:), Inf)
       exitflag(j)
       
       cnt = cnt+1;
       
       %if it goes for too long, quit
       if cnt==30 
           %& exitflag(j) ~= 1
           verify(j)=-1;
          break;  
       end
     
      end  
  end
  
  toc
end

%betweenness centrality
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
   
   for kk = 1:numberOfNodes %computing the length of the shortest path
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

   %the algorithm goes from each node j to each node k on each of the
   %shortest paths counting them
   possibleLengths = zeros(numberOfNodes);
   
   for j=1:numberOfNodes
       for k=1:numberOfNodes
           
         level=1; %level 1 analyzes nodes that are neighbors; level k analyzes nodes that are at a distance k.
              
         if shortestLengths(j,k) > 0
             
           neighbors = find(adj(j,:)==1);
           if shortestLengths(j,k)==1
             possibleLengths(j,k) = 1;  %when the nodes are neighbors, the # of possible shortest paths is one
           else    
             possibleLengths(j,k) = 0;
           end
              
        while level < shortestLengths(j,k)  %if nodes are not neighbors
            
           neighbors2=[];
           
          for l=1:length(neighbors) %looking at all the nodes in the vector neighbors
              
           if shortestLengths(j,k) == level + shortestLengths(neighbors(l),k) %check if on the path of shortest length
               if level == shortestLengths(j,k)-1 %if on the level almost of the shortest path (right before reaching node k)
                 possibleLengths(j,k) = possibleLengths(j,k)+1;  
               else
                 neighbors2 = [neighbors2, find(adj(neighbors(l),:)==1)]; %else get all the neighbors of the neighbors of l; these are possibly on the shortest path
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
   
   %computing betweenness centrality, given by computing, for each node o
   %and pair (j,k), the number of shortest paths from j to k passing
   %through o.
   
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
  
end

%various methods of quantifying leadership
function LeadershipMeasures_pushbutton_Callback(hObject, ~, handles)

prompt = {'Assesing leadership by one of the 7 options:                                                                        Select1). Node status                                                   Select2). Leadership index                                                             Select 3). Lead and Lag                                                               Select 4). Lead and Lag v.2                                                                                                             Select 5). Correlation with center                                                                                                                          Select 6). Sparsity connections (Sznaier)                                                                                                      Select 7). Other centrality measures'};
option = str2num(cell2mat(inputdlg(prompt)));

numberOfNodes = handles.numberOfNodes;
colors = jet(numberOfNodes);

choose_graph(hObject, handles,0);
adj = handles.adjacencyMatrix;

switch option    
    case 1  %node status SG
        
        prompt = {'Select start frame:','Select end frame:'};
        answer = inputdlg(prompt, 'Select range of frames');
        start_frame = str2num(answer{1});
        end_frame = str2num(answer{2});

        prompt = {'For a few nodes or for multiple?                                                                                                       Select 1) if for a few nodes                                                                                                                              Select 2) if for all nodes'};
        answer = str2num(cell2mat(inputdlg(prompt)));
        
        [nodeStatus, nodeStatus_avg] = compute_nodeStatus(hObject, handles, adj, start_frame, end_frame);
        
      switch answer
          
          case 1 
              
        S = [{'1'},{'2'},{'3'},{'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'9'}, {'10'}, {'11'}, {'12'}, {'13'}];
        nodes = listdlg('ListString', S, 'Name', 'Node status & correlations', 'PromptString', 'Select the node(s)');
      
        leader = zeros(1, 1:end_frame-start_frame+1);
        colors = jet(length(nodes));
        
        figure(1)
        for j=1:length(nodes)
            %plot ns SG
            plot(start_frame:end_frame, nodeStatus(1:end_frame-start_frame+1, nodes(:,j)), 'Color', colors(j,:))
            hold on;
            j
            disp('average node status & std of node status')
            %display mean and std
            nodeStatus_mean(j) = mean(nodeStatus(1:end_frame-start_frame+1, nodes(:,j)))
            nodeStatus_std(j) = std(nodeStatus(1:end_frame-start_frame+1, nodes(:,j)))
        end
        
          case 2
              
              %for all nodes
          for j=1:numberOfNodes
            nodeStatus_mean(j) = mean(nodeStatus(1:end_frame-start_frame+1, j));
            nodeStatus_std(j) = std(nodeStatus(1:end_frame-start_frame+1, j));
          end    
          
          %leader is the one with the maxium node status
              for i=start_frame:end_frame
            [maxNodeStatus, indexNode] = max(nodeStatus(i-start_frame+1,:));
            leader(i) = indexNode;
              end
        
        disp('leader')
        leader
        
        [nodeStatus_mean_sorted, Index_nodeStatus_sorted] = sort(nodeStatus_mean);
        
        disp('sort as a function of node status')
        nodeStatus_mean_sorted
        Index_nodeStatus_sorted
        
        [nodeStatus_std_sorted, Index_nodeStatus_sorted] = sort(nodeStatus_std);
        
        disp('sort as a function of std of node status')
        nodeStatus_std_sorted
        Index_nodeStatus_sorted

        figure(1)
        barwitherr(nodeStatus_std, nodeStatus_mean)
        
      end
        
    case 2 %leadership index, according to formula in Quera V., Beltran F.S., Dolado R., Flocking behaviour: Agent-based simulation and hierarchical leadership.

        
        prompt = {'Select start frame:','Select end frame:'};
        answer = inputdlg(prompt, 'Select range of frames');
        start_frame = str2num(answer{1});
        end_frame = str2num(answer{2});

        L = leadershipIndex(hObject, handles, adj, start_frame, end_frame);
        leader = zeros(1, 1:end_frame-start_frame+1);
        
        S = [{'1'},{'2'},{'3'},{'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'9'}, {'10'}, {'11'}, {'12'}, {'13'}];
        nodes = listdlg('ListString', S, 'Name', 'Node status & correlations', 'PromptString', 'Select the node(s)');

        figure(2)
        for j=1:length(nodes)
            plot(start_frame:end_frame, L(nodes(j), :), 'Color', colors(j,:));
            hold on;
            j
            disp('Average Leadership Index')
            mean_L(j) = mean(L(nodes(j),:))
            std_L(j) = std(L(nodes(j),:))
        end
        
        %leader at frame i is the one with highest LI at frame i
      for i=start_frame:end_frame
          [maxLI, indexNode] = max(L(:,i-start_frame+1));
          leader(i) = indexNode;
      end
      
      leader
      
      [L_mean_sorted, Index_L_sorted] = sort(mean_L);
      L_mean_sorted
      Index_L_sorted  
      
      [L_std_sorted, Index_L_sorted] = sort(std_L);
      L_std_sorted
      Index_L_sorted  
      
    case 3
        
      [LeadAndLag, frame, frame1, frame2] = lead_and_lag(hObject, handles) 
      
        gObj = biograph(LeadAndLag);
        view(gObj);        
      
    case 4
        
        prompt = {'Select start frame:','Select end frame:'};
        answer = inputdlg(prompt, 'Select range of frames');
        start_frame = str2num(answer{1});
        end_frame = str2num(answer{2});

        [following, following2, LeadAndLag] = lead_and_lag2(hObject, handles, start_frame, end_frame);
        following
        following2
        
    case 5

        prompt = {'Select start frame:','Select end frame:'};
        answer = inputdlg(prompt, 'Select range of frames');
        start_frame = str2num(answer{1});
        end_frame = str2num(answer{2});

        [following, following2] = correlation_withCenter(handles, hObject, start_frame, end_frame)
        
    case 6 %compute ns IG
        
        [x, exitflag, start_frame, end_frame, delay, nodes] = causality(hObject, handles);
        disp('x')
        x
        disp('exitflag')
        exitflag
        
        F = end_frame-start_frame+1;
        threshold = 0.005;
        %threshold = 0.00248;
        numberOfNodes = length(nodes);
        Interaction = zeros(numberOfNodes);
        
        for i=1:numberOfNodes
            for j=1:numberOfNodes
                    
                if max(x(i,(j-1)*(delay+1)+1:j*(delay+1))) >= threshold & exitflag(i)==1 & exitflag(j)==1
                    %sum(x(i,(j-1)*(delay+1)+1:j*(delay+1)))>=(delay+1)*0.002 
                    Interaction(i,j) = 1;
                    
                    if Interaction(j,i)==1   %if both follow each other, we are not interested...
                        Interaction(i,j)=0;
                        Interaction(j,i)=0;
                    end
                end
                
            end
        end
      
        save('IG_940-953.mat', 'Interaction')
      
         for j=1:length(nodes)
            NodeID{j,:} = num2str(nodes(j));
         end
        
         figure(1)
         gObj = biograph(Interaction, NodeID);
         view(gObj);    

        Interaction2(handles.lb,:,:) = Interaction;   
        nodeStatus = compute_nodeStatus(hObject, handles, Interaction2, handles.lb, handles.lb);
        
        disp('node statuses')
        nodeStatus(handles.lb,:)
        
        [nodeStatus_sorted, Index_nodeStatus_sorted] = sort(nodeStatus(handles.lb,:));
        
        nodeStatus_sorted
        Index_nodeStatus_sorted
        
        plot(1:numberOfNodes, nodeStatus(handles.lb,:), 'x', 'Color', 'b', 'MarkerSize', 14)
        
    case 7
        
        prompt = {'Choose some other centrality measures:                               Select 1). Betweenness (how often node is bridge within shortest path)                                                                                    Select 2). Eigenvector centrality                                                        Select 3). Katz centrality version 1                                                                   Select 4). Katz centrality version 2                                                         Select 5). Page rank                                                                     Select 6). Degree'};
        option2 = str2num(cell2mat(inputdlg(prompt)));
        
        prompt = {'Select start frame:','Select end frame:'};
        answer = inputdlg(prompt, 'Select range of frames');
        start_frame = str2num(answer{1});
        end_frame = str2num(answer{2});

        switch option2
            
            case 1 %betweenness centrality, see wikipedia article for how to compute it
                
                Betweenness = [];
                for i = start_frame:end_frame
                    
                    betweenness_c = betweenness(squeeze(adj(i,:,:)));
                    Betweenness = [Betweenness; betweenness_c];
                    
                end
                
                b_mean = mean(Betweenness)
                b_std = std(Betweenness)
                [sorted_b, sortedInd_b] = sort(b_mean);
                disp('Sorted mean of betweenness')
                sorted_b
                sortedInd_b
                disp('Sorted std of betweenness')
                sorted_b
                [sorted_b, sortedInd_b] = sort(b_std);
                sortedInd_b
                
                figure(1)
                barwitherr(b_std, b_mean)
                
            case 2 %eigenvector centrality
            
                for i=start_frame:end_frame
                    A = squeeze(adj(i,:,:));
                    lambda = eig(A);
                    [V,D] = eig(A);
                    
                    [ignore, maxInd_lambda] = max(lambda);
                    eig_vector(i-start_frame+1,:) = V(:,maxInd_lambda)';
                end
                
                eig_mean = mean(eig_vector);
                eig_std = std(eig_vector);
                [ignore, Ind_sorted] = sort(eig_mean);
                'Mean sorted'
                Ind_sorted
                [ignore, Ind_sorted] = sort(eig_std);
                'Std sorted'
                Ind_sorted
                
                figure(1)
                %bar(eig_mean)
                barwitherr(eig_std, eig_mean)
                
            case 3 %Katz centrality
                
                m = 100;
                for i=start_frame:end_frame
                    
                    A(1,:,:) = squeeze(adj(i,:,:));
          
                    for k=2:m
                        A(k,:,:) = squeeze(A(k-1,:,:))*squeeze(A(1,:,:));
                    end
                    
                    lambda=eig(squeeze(A(1,:,:)));
                    max_lambda = max(lambda);
                    alpha = 1/max_lambda-1/(max_lambda*100);
               
                    for j=1:numberOfNodes
                        x(i-start_frame+1,j) = 0;
                       for k=1:m 
                          x(i-start_frame+1,j) = x(i-start_frame+1,j) + alpha^k*sum(squeeze(A(k,:,j)));
                       end
                    end
                   
                end
                
                x
                x_mean = mean(x)
                [ignore, Ind_sorted]=sort(x_mean);
                'Mean sorted'
                Ind_sorted
                x_std = std(x);
                [ignore, Ind_sorted]=sort(x_std);
                'Std sorted'
                Ind_sorted
                
                figure(1)
                barwitherr(x_std, x_mean)
                
            case 4 %Katz centrality v. 2, does it work?
                
                one_v  = ones(1, numberOfNodes)';
                
                for i=start_frame:end_frame
                    
                    lambda=eig(squeeze(adj(i,:,:)));
                    max_lambda = max(lambda);
                    alpha = 1/max_lambda-1/(max_lambda*100);
                    
                    A = eye(numberOfNodes) - alpha*squeeze(adj(i,:,:));
                    B = alpha*squeeze(adj(i,:,:))*one_v;
                    x(i-start_frame+1,:) = linsolve(A,B);
                end
                
                x
                x_mean = mean(x);
                [ignore, Ind_sorted]=sort(x_mean);
                'Mean sorted'
                Ind_sorted
                
                x_std = std(x);
                [ignore, Ind_sorted]=sort(x_std);
                'Std sorted'
                Ind_sorted
                
                'x_mean'
                x_mean
                'x_std'
                x_std
                
                figure(2)
                barwitherr(x_std, x_mean)
                
            case 5 %Page rank; does it work?
                
                for i=start_frame:end_frame
                    
                    a = squeeze(adj(i,:,:));
                    lambda=eig(a);
                    max_lambda = max(lambda);
                    alpha = 1/max_lambda-1/(max_lambda*10);
                    
                    for k=1:numberOfNodes
                        degree(k) = length(find(a(:,k)==1));
                        a(:,k) = a(:,k)/degree(k);
                    end
                    
                    A = eye(numberOfNodes) - alpha*a';
                    B = (1-alpha)/(numberOfNodes) * ones(1, numberOfNodes)';
                    
                    x(i-start_frame+1,:) = linsolve(A,B);
                    %x(i-start_frame+1,:) = A\B;
                end
                
                x
                x_mean = mean(x);
                x_std = std(x);
                
                'x_mean'
                x_mean
                'x_std'
                x_std
                
                [ignore, Ind_sorted]=sort(x_mean);
                'Mean sorted'
                Ind_sorted
                [ignore, Ind_sorted]=sort(x_std);
                'Std sorted'
                Ind_sorted
                
                figure(3)
                barwitherr(x_std, x_mean);
                
            case 6 %degree in/out
                                
                [degree_in, degree_out] = degreeInOut(adj, start_frame, end_frame)
                squeeze(adj(start_frame,:,:))
                
              if size(degree_in,1)~=1 & size(degree_in,2)~=1 
                degree_in_mean = mean(degree_in);
                degree_in_std = std(degree_in);
              else
                   degree_in_mean = degree_in;
                   degree_in_std = zeros(1, length(degree_in)); 
              end
              
                degree_in_mean
                degree_in_std
                
                'Sort in-degree mean'
                [ignore, Ind_sorted] = sort(degree_in_mean);
                Ind_sorted
                ignore
                'Sort degree std'
                [ignore, Ind_sorted] = sort(degree_in_std);
                Ind_sorted
                ignore
                
                figure(2)
                barwitherr(degree_in_std, degree_in_mean);
        end
end
end

%degree in and degree out given an adjacency matrix
function [degree_in, degree_out] = degreeInOut(adj, start_frame, end_frame)

numberOfNodes = size(adj,2)

    for i=start_frame:end_frame
                    
        a = squeeze(adj(i,:,:));
        for k=1:numberOfNodes
            degree_in(i-start_frame+1,k) = length(find(a(:,k)==1));
            degree_out(i-start_frame+1,k) = length(find(a(k,:)==1));
        end
                  
    end
    
end

%determining the shape of the group
function [cavity, density, boundary_nodes] = alpha_shape(handles, hObject)

lb = handles.lb;
ub = handles.ub;

prompt = {'Select start frame:','Select end frame:'};
answer = inputdlg(prompt, 'Select range of frames');
start_frame = str2num(answer{1});
end_frame = str2num(answer{2});

for i=start_frame:end_frame
    
  file=strcat('Dance',num2str(i),'.dat');
  data_type = handles.data_type;
  data=handles.dancedata{i};
  x = data(:,1);
  y = data(:,2);
  theta = data(:,3);
        
     %[x y theta] = transf_in_familiar_coord(x,y,theta);        
  dim = 2;
  positionMatrix(i,:,:) = horzcat(x,y);
  
  numberOfNodes = length(x);
  
 cnt=0;  
 array = 0.6:-0.1:0.1; 
 
 %use algorithm to determine alpha shape
 for R = array
  cnt = cnt+1;
  [V_cum(cnt),S_cum(cnt)] = alphavol(squeeze(positionMatrix(i,:,:)),R, 1);
 end

 %select proper parameter by where the difference is greatest
 [ignore, ind_dV]=sort(abs(diff(V_cum)));
 cavity(i-start_frame+1) = array(ind_dV(cnt-1));
 
 V(i-start_frame+1) = V_cum(ind_dV(cnt-1));
 S(i-start_frame+1) = S_cum(ind_dV(cnt-1));
 
 %V(i-start_frame+1)
 %compute density
 density(i-start_frame+1) = numberOfNodes/V(i-start_frame+1);
 
 %which nodes are on the boundary
  boundary_nodes{i-start_frame+1}=[];
for j=1:numberOfNodes 
    
 if ismember(j, S(i-start_frame+1).bnd)
  boundary_nodes{i-start_frame+1}=[boundary_nodes{i-start_frame+1}, j];
  boundary_nodes2(i-start_frame+1,j)=1;
 else
  boundary_nodes2(i-start_frame+1,j)=0;   
 end
 
end
end

disp('boundary nodes')
for i=start_frame:end_frame
boundary_nodes{i-start_frame+1}
end

end

%find radius inside which no neighbors are (the 'private space')
function r_0 = fit_R(handles, hObject)

lb = handles.lb;
ub = handles.ub;

for i=lb:ub
    
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
     positionMatrix(i,:,:) = horzcat(x,y);
  end
  
end

numberOfNodes = length(x);
[cavity, density, boundary_nodes] = alpha_shape(handles, hObject);

%use 'hard ball' idea from the paper Cavagna A. et al., The STARFLAG handbook on collective animal behaviour:
% Part II, three-dimensional analysisThe STARFLAG handbook on collective animal behaviour.

%I am not sure if pdf_hardBall does not need to be modified; the results outputted by this algorithm seem too small.

for o=1:numberOfNodes
    %search for radius
    r_0(o) = fminsearch(@(x) pdf_hardBall(x, positionMatrix, density, lb, ub, o, numberOfNodes), [0.1]);
end

r_0
end

function Miscellaneous_pushbutton_Callback(hObject, ~, handles)

prompt = {'Select what to do:                                                                               Select 1). fit R                                                                                                                                    Select 2) compare node status among graphs                                                                                                                                                          Select 3). Create huge matrix for ML                                                                                                                                                      Select 4). Plot direction of center                                                                                                                                               Select 5). compute inertia matrix and principal moments of inertia                                                                                                                                                          Select 6). Compute density                                                                                                                  Select 7). Plot direction of a single node of choice                                                                                                        Select 8). Build adjacenct matrix from New Model data                                                                                                                                                   Select 9). Compute properties and make video'};    
Misc_choice = str2num(cell2mat(inputdlg(prompt)))

  delay=20;
  
switch Misc_choice
    
    case 1
        
r_0 = fit_R(handles, hObject) % fit radius of each node where no neighbor can be

    case 2
  
  prompt = {'Select start frame:','Select end frame:'};
  answer = inputdlg(prompt, 'Select range of frames');
  start_frame = str2num(answer{1});
  end_frame = str2num(answer{2});

compare_graph_properties(hObject, handles, start_frame, end_frame) %compare node statuses from different sensing  graphs

    case 3
        
group_properties(hObject,handles)        

    case 4
        
prompt = {'Select start frame:','Select end frame:'};
answer = inputdlg(prompt, 'Select range of frames');
start_frame = str2num(answer{1});
end_frame = str2num(answer{2});
        
center_dir(hObject, handles, 30, start_frame, end_frame) %direction of center, angle of each node from the center

    case 5
        
        %computing first moment of inertia
  prompt = {'Select start frame:','Select end frame:'};
  answer = inputdlg(prompt, 'Select range of frames');
  start_frame = str2num(answer{1});
  end_frame = str2num(answer{2});
  
[positionMatrix, center, numberOfNodes] = Read_coordinatesAndCenter(start_frame, end_frame, handles, hObject);
[I, V, D] = inertia_matrix(positionMatrix, handles, hObject);

for i=1:end_frame-start_frame+1
    i
    [num2str(D(i,1,1)), '  ', num2str(D(i,2,2)), '  ', num2str(D(i,3,3))]
    sorted_eig(i,:) = sort([D(i,1,1), D(i,2,2), D(i,3,3)]);
    ratio_eig(i) = sorted_eig(i,1)/sorted_eig(i,3);
end
ratio_eig

    case 6
        %shape of group
[cavity, density, boundary_nodes] = alpha_shape(handles, hObject);

    case 7
        %compare direction of nodes smoothed out or not to see which method
        %of averaging might be best
[theta, theta2, theta3, theta4, theta_avg, theta_avg2] = direction(handles, hObject);

    case 8
[adjacencyMatrix min_frame]= build_adjMatrix(handles, hObject);

    case 9
[adjacencyMatrix min_frame]= build_adjMatrix(handles, hObject);
compute_properties_withoutVideo(handles, hObject)

end

end

%important function to draw all the conclusions about position and leadership in the group...
function correlations_pushbutton_Callback(hObject, ~, handles)
  
% choose a graph which gives the adj matrix
choose_graph(hObject, handles,0);
adj = handles.adjacencyMatrix;
numberOfNodes = handles.numberOfNodes;

  prompt = {'Select start frame:','Select end frame:'};
  answer = inputdlg(prompt, 'Select range of frames');
  start_frame = str2num(answer{1});
  end_frame = str2num(answer{2});
 
  %whether to analyze averages of node status and position in group (distances) for all nodes (option 2) OR
  %node status and position within bins for ONE node (option 2)
  prompt = {'Node status & positions. Choose what to analyze:                                                                               Select 1). Analyze correlations for one individual node                                                                                                         Select 2). Analyze correlations averaging over time for all nodes'};
  answer = str2num(cell2mat((inputdlg(prompt))));
  
  if answer==1 
    
      S = [{'1'},{'2'},{'3'},{'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'9'}, {'10'}, {'11'}, {'12'}, {'13'}];
      nodes = listdlg('ListString', S, 'Name', 'Node status & correlations', 'PromptString', 'Select the node(s)');
      
  if length(nodes)>1
      error('In this scenario, please choose only one node');
      return
  end
  
  delay = 20;

  node = 1:numberOfNodes;
%{
  %this loop selects a sensing graph (from options 1,2,3) and computes
  %correlations between node status SG or IG and: distance to center (d),
  %distance to center-axis (D) and absolute value of distance to
  %center-axis (|D|); also I compute the correlation between ns SG and ns
  %IG.
  
for adj_graph=1:3  
    
    switch adj_graph
        case 1
            disp('360 degrees')
            adj = choose_graph(hObject, handles,1);

        case 2
            disp('180 degrees, R= infty')
            adj = choose_graph(hObject, handles,5);
            
        case 3
            disp('New Model')
            adj = choose_graph(hObject, handles,6);
            end_frame = 903;      
    end
    
    %compute node status SG, IG, distance to center, distance to
    %center-axis
  [nodeStatus, nodeStatus_avg] = compute_nodeStatus(hObject, handles, adj, start_frame, end_frame);
  [dist_to_center, avg_toCenter] = distToCenter(handles, hObject, start_frame, end_frame, numberOfNodes,0);
  [Distance2, Distance3, distance, leader, backer, N_leader, N_backer] = LeaderBacker(hObject, handles, start_frame, end_frame,1);
  
  [nodeStatus_IG, nodeStatus_avg2_IG] = nodeStatusIG(handles, hObject, start_frame, end_frame);
  [center, center_avg, center2, dir_center, dir_node, dir_node2] = center_dir(hObject, handles, delay, start_frame, end_frame);
  
  for i=1:numberOfNodes  
     
      %smooth out these functions
  [ignore, ignore, ignore, nodeStatus_smooth2, xax] = averaging(nodeStatus(:, node(i)), 20, 1, 1);
  [ignore, ignore, ignore, nodeStatus_IG_smooth2, xax2] = averaging(nodeStatus_IG(:,node(i)), 20, 1, 1);
  
  %do all these correlation analyses with corrcoef
  disp('node status and distance to center correlated')
  [c p] = corrcoef(nodeStatus_smooth2(1:end_frame-start_frame+1), dist_to_center(1:end_frame-start_frame+1, node(i)));
  c(2)
  p(2)
  
  disp('node status and distance to center-axis')
  [c p] = corrcoef(nodeStatus_smooth2(1:end_frame-start_frame-19), Distance3(1:end_frame-start_frame-19, node(i)));
  c(2)
  p(2)
  
  disp('node status and abs value of distance to center-axis correlated')
  [c p] = corrcoef(nodeStatus_smooth2(1:end_frame-start_frame-19), distance(1:end_frame-start_frame-19, node(i)));
  c(2)
  p(2)
 
  
  disp('node status SG and node status IG correlated')
  [c p] = corrcoef(nodeStatus_smooth2(1:end_frame-start_frame), nodeStatus_IG_smooth2(1:end_frame-start_frame));
  c(2)
  p(2)
  
  disp('node status IG and distance to center correlated')
  [c p] = corrcoef(dist_to_center(1:end_frame-start_frame+1, node(i)), nodeStatus_IG_smooth2(start_frame:end_frame));
  c(2)
  p(2)
  
  disp('node status IG and distance from center-axis correlated')
  [c p] = corrcoef(nodeStatus_IG_smooth2(1:end_frame-start_frame-19), Distance3(1:end_frame-start_frame-19, node(i)));
  c(2)
  p(2)
  
  disp('node status IG and absolute value of distance to center-axis correlated')
  [c p] = corrcoef(nodeStatus_IG_smooth2(1:end_frame-start_frame-19), distance(1:end_frame-start_frame-19, node(i)));
  c(2)
  p(2)
 
  end
end
  %}
 
  %compute node status SG, IG, distance to center, distance to
  %center-axis and smoothing these functions out...
  [nodeStatus, nodeStatus_avg] = compute_nodeStatus(hObject, handles, adj, start_frame, end_frame);
  
  for o=1:numberOfNodes
      [nodeStatus_all(:,o),ignore,ignore,ignore,ignore] = averaging(nodeStatus(:,o),20,1,1);
  end
  
  [Distance2, Distance3, distance, leader, backer, N_leader, N_backer] = LeaderBacker(hObject, handles, start_frame, end_frame,1);
  [dist_to_center, avg_toCenter] = distToCenter(handles, hObject, start_frame, end_frame, numberOfNodes, 1);
  
  [nodeStatus_IG, nodeStatus_avg2_IG] = nodeStatusIG(handles, hObject, start_frame, end_frame);
  [center, center_avg, center2, dir_center, dir_node, dir_node2] = center_dir(hObject, handles, delay, start_frame, end_frame);
  
  for o=1:numberOfNodes
      [nodeStatus_IG_all(:,o), ignore, ignore, ignore, ignore] = averaging(nodeStatus_IG(:,o),20,1,1);
  end
  
  node = nodes(1);
  
  [nodeStatus_smooth2, ignore, ignore, ignore, xax] = averaging(nodeStatus(:, node), 20, 1, 1);
  [nodeStatus_IG_smooth2, ignore, ignore, ignore, xax2] = averaging(nodeStatus_IG(:,node), 20, 1, 1);
  
  nodeStatus = nodeStatus_smooth2;
  nodeStatus_IG = nodeStatus_IG_smooth2;
  
  %plot node status SG. IG versus distances (d, D, |D|; d = distToCenter; D
  %= Distance3; |D| = distance) representing position
  
  figure(1)
  plot(start_frame:end_frame, nodeStatus_smooth2(1:end_frame-start_frame+1), 'b-', 'linewidth', 2)
  hold on
  plot(start_frame:end_frame, dist_to_center(1:end_frame-start_frame+1, node), 'r-', 'linewidth', 2)
  title('Node status and distance to center as a function of time')
  xlabel('frames')
  legend('node status', 'distance to center')
  
  disp('node status and distance to center correlated')
  
  [xcf, lags] = xcorr(nodeStatus(1:end_frame-start_frame+1)-mean(nodeStatus(1:end_frame-start_frame+1)), dist_to_center(1:end_frame-start_frame+1, node)-mean(dist_to_center(1:end_frame-start_frame+1, node)), 20, 'coeff');
  xcf
  [c p] = corrcoef(nodeStatus_smooth2(1:end_frame-start_frame+1), dist_to_center(1:end_frame-start_frame+1, node));
  c
  p
  
  figure(2)
  plot(start_frame:end_frame-20, nodeStatus_smooth2(1:end_frame-start_frame-19), 'b-', 'linewidth', 2)
  hold on
  plot(start_frame:end_frame-20, Distance3(1:end_frame-start_frame-19), 'r-', 'linewidth', 2)
  title('Node status and distance from center-axis as a function of time')
  
  xlabel('frames')
  legend('node status', 'distance from center-axis')
  
  disp('node status and distance to center-axis')
  [xcf, lags] = xcorr(nodeStatus(1:end_frame-start_frame-19)-mean(nodeStatus(1:end_frame-start_frame-19)), Distance3(1:end_frame-start_frame-19, node)-mean(Distance3(1:end_frame-start_frame-19, node)), 20, 'coeff');
  xcf
  [c p] = corrcoef(nodeStatus_smooth2(1:end_frame-start_frame-19), Distance3(1:end_frame-start_frame-19, node));
  c
  p
  
  figure(3)
  plot(start_frame:end_frame-20, nodeStatus_smooth2(1:end_frame-start_frame-19), 'b-', 'linewidth', 2)
  hold on
  plot(start_frame:end_frame-20, distance(1:end_frame-start_frame-19, node), 'r-', 'linewidth', 2)
  title('node status and absolute value of distance to center-axis as a function of time')
  xlabel('frames')
  legend('node status', '|distance to center-axis|')
  
  disp('node status and abs value of distance to center-axis correlated')
  [xcf, lags] = xcorr(nodeStatus(1:end_frame-start_frame-19)-mean(nodeStatus(1:end_frame-start_frame-19)), distance(1:end_frame-start_frame-19, node)-mean(distance(1:end_frame-start_frame-19, node)), 20, 'coeff');
  xcf
  [c p] = corrcoef(nodeStatus_smooth2(1:end_frame-start_frame-19), distance(1:end_frame-start_frame-19, node));
  c
  p
  
  figure(4)
  
 % plot(start_frame:end_frame, nodeStatus_smooth, 'b-', 'linewidth', 2)
 % hold on
  %{
  for i=start_frame:end_frame-20
      if leader(i-start_frame+1)==node
          %plot(i, nodeStatus(i-start_frame+1,node), 'r-*', 'linewidth', 4)
          plot(i, 0.5, 'r-*', 'linewidth', 4)
          hold on
      end
  end
  
  title('Node status and periods of being at front end, as a function of time')
  xlabel('frames')
  %}
  
  %xx = xax(1):1:xax(length(xax));
  %yy = spline(xax, nodeStatus_smooth2, xx);
  
  plot(start_frame:end_frame, nodeStatus_smooth2, 'b-', 'linewidth', 2);
  hold on
  %plot(xx, yy, 'g-', 'linewidth', 2)
  
   plot(start_frame:end_frame, nodeStatus_IG_smooth2(start_frame:end_frame), 'r-', 'linewidth', 2)
   hold on
   title('Node status and periods of being at front end, as a function of time')
   xlabel('frames')
   legend('node status SG', 'node status IG')
  
  for i=start_frame:end_frame-20
      if leader(i-start_frame+1)==node
          %plot(i, nodeStatus(i-start_frame+1,node), 'r-*', 'linewidth', 4)
          plot(i, 0.5, 'g-*', 'linewidth', 2)
          hold on
      end
  end
  
   figure(5)
 plot(start_frame:end_frame, nodeStatus_smooth2, 'b-', 'linewidth', 2) 
 hold on
 plot(start_frame:end_frame, nodeStatus_IG_smooth2(start_frame:end_frame), 'r-', 'linewidth', 2) 
 title('node status SG vs node status IG')
 legend('node status SG', 'node status IG')
 
  disp('node status SG and node status IG correlated')
  [xcf, lags] = xcorr(nodeStatus(1:end_frame-start_frame)-mean(nodeStatus(1:end_frame-start_frame)), nodeStatus_IG(1:end_frame-start_frame)-mean(nodeStatus_IG(1:end_frame-start_frame)), 20, 'coeff');
  xcf
  [c p] = corrcoef(nodeStatus_smooth2(1:end_frame-start_frame), nodeStatus_IG_smooth2(1:end_frame-start_frame));
  c
  p
  
 figure(6)
 plot(start_frame:end_frame, nodeStatus_IG_smooth2(start_frame:end_frame), 'b-', 'linewidth', 2) 
 hold on
 plot(start_frame:end_frame, dist_to_center(1:end_frame-start_frame+1, node), 'r-', 'linewidth', 2) 
 title('node status IG vs. distance to center (point)')
 legend('node status IG', 'distance to center')
 
 disp('node status IG and distance to center correlated')
 [xcf, lags] = xcorr(nodeStatus_IG(start_frame:end_frame)-mean(nodeStatus_IG(start_frame:end_frame)), dist_to_center(1:end_frame-start_frame+1, node)-mean(dist_to_center(1:end_frame-start_frame+1, node)), 20, 'coeff');
  xcf
  [c p] = corrcoef(dist_to_center(1:end_frame-start_frame+1, node), nodeStatus_IG_smooth2(start_frame:end_frame));
  c
  p
 
  figure(7)
 plot(start_frame:end_frame-20, nodeStatus_IG_smooth2(start_frame:end_frame-20), 'b-', 'linewidth', 2) 
 hold on
 plot(start_frame:end_frame-20, Distance3(1:end_frame-start_frame-19, node), 'r-', 'linewidth', 2) 
 title('Node status IG vs. distance to center-axis')
 legend('node status IG', 'Distance2')
 
 disp('node status IG and distance from center-axis correlated')
  [xcf, lags] = xcorr(nodeStatus_IG(start_frame:end_frame-20)-mean(nodeStatus_IG(start_frame:end_frame-20)), Distance3(1:end_frame-start_frame-19, node)-mean(Distance3(1:end_frame-start_frame-19, node)), 20, 'coeff');
  xcf
  [c p] = corrcoef(Distance3(1:end_frame-start_frame-19, node), nodeStatus_IG_smooth2(start_frame:end_frame-20));
  c
  p
  
  figure(8)
 plot(start_frame:end_frame-20, nodeStatus_IG_smooth2(start_frame:end_frame-20), 'b-', 'linewidth', 2) 
 hold on
 plot(start_frame:end_frame-20, distance(1:end_frame-start_frame-19, node), 'r-', 'linewidth', 2) 
 title('Node status IG vs. absolute value of distance from center-axis')
 legend('node status IG', 'abs distance')
 
  disp('node status IG and absolute value of distance to center-axis correlated')
  [xcf, lags] = xcorr(nodeStatus_IG(start_frame:end_frame-20)-mean(nodeStatus_IG(start_frame:end_frame-20)), distance(1:end_frame-start_frame-19, node)-mean(distance(1:end_frame-start_frame-19, node)), 20, 'coeff');
  xcf
  [c p] = corrcoef(distance(1:end_frame-start_frame-19, node), nodeStatus_IG_smooth2(start_frame:end_frame-20));
  c
  p
 
  
  dist_to_center2 = dist_to_center(1:end_frame-delay-start_frame+1,:);
 
  %PART I: for ALL nodes, bins contain data for any node
  
  max_D3 = max(max(Distance3));
  min_D3 = min(min(Distance3));
  max_d  = max(max(distance));
  min_d = min(min(distance));
  min_dist_to_center = min(min(dist_to_center2));
  max_dist_to_center = max(max(dist_to_center2));
  
  D3_interval = (max_D3-min_D3)/10;
  d_interval = (max_d-min_d)/10;
  dcenter_interval = (max_dist_to_center - min_dist_to_center)/10;
  
  ss_D3_SG = NaN*ones(10,10);
  ss_d_SG = NaN*ones(10,10);
  ss_dCenter_SG = NaN*ones(10,10);
  ss_dCenter_front_SG = NaN*ones(10,10);
  ss_dCenter_middle_SG = NaN*ones(10,10);
  ss_dCenter_back_SG = NaN*ones(10,10);
  
  ss_D3_IG = NaN*ones(10,10);
  ss_d_IG = NaN*ones(10,10);
  ss_dCenter_IG = NaN*ones(10,10);
  ss_dCenter_front_IG = NaN*ones(10,10);
  ss_dCenter_middle_IG = NaN*ones(10,10);
  ss_dCenter_back_IG = NaN*ones(10,10);
  
  dist_to_center_longV =[];
  distance_longV = [];
  Distance3_longV = [];
  corr_row = [];
  corr_col = [];
  
  dim1 = size(distance,1);
  dim2 = size(distance,2);
  
  %put vectors of distances together for all nodes and then divide into
  %bins
  for s = 1:size(dist_to_center2,1)
      dist_to_center_longV = [dist_to_center_longV, dist_to_center2(s,:)];
      distance_longV = [distance_longV, distance(s,:)];
      Distance3_longV = [Distance3_longV, Distance3(s,:)];
  
  %vector of nodes and frames
      corr_row = [corr_row, s*ones(1,dim2)]; %rows = nodes
      corr_col = [corr_col, 1:dim2]; %columns = frames
  
  end
  
  %sort these vectors that were put together from smallest distances to
  %largest
  [dist_to_center_longV ind_dist_to_center_longV] = sort(dist_to_center_longV);
  [distance_longV ind_distance_longV] = sort(distance_longV);
  [Distance3_longV ind_Distance3_longV] = sort(Distance3_longV);
  
  corr_row_dist_to_center = corr_row(ind_dist_to_center_longV);
  corr_row_distance = corr_row(ind_distance_longV);
  corr_row_Distance3 = corr_row(ind_Distance3_longV);
  
  corr_col_dist_to_center = corr_col(ind_dist_to_center_longV);
  corr_col_distance = corr_col(ind_distance_longV);
  corr_col_Distance3 = corr_col(ind_Distance3_longV);
  
  m = length(distance_longV);
  
  %splitting data into bins, k is the bin
 for k=1:10    %quantifying overall, for all nodes, how the distance from center influences node status
     %rows are frames, cols are nodes
  %ONE idea is to split data so that each bin spans the same distance, this
  %is not detailed in the report
  [rows_d{k}, cols_d{k}] = find(distance>=min_d+(k-1)*d_interval & distance<=min_d+k*d_interval);
  [rows_D3{k}, cols_D3{k}] = find(Distance3>=min_D3+(k-1)*D3_interval & Distance3<=min_D3+k*D3_interval);
  [rows_dcenter{k}, cols_dcenter{k}] = find(dist_to_center2 >= min_dist_to_center + (k-1)*dcenter_interval & dist_to_center2 <= min_dist_to_center+k*dcenter_interval);

  %second idea, and the one actually implemented for the study is to split
  %data so that each bin contains equal number of trials
  
 if k~=10 
  rows_d2{k} = corr_row_distance((k-1)*floor(m/10)+1:k*floor(m/10));    
  rows_D32{k} = corr_row_Distance3((k-1)*floor(m/10)+1:k*floor(m/10));
  rows_dcenter2{k} = corr_row_dist_to_center((k-1)*floor(m/10)+1:k*floor(m/10));
  
  cols_d2{k} = corr_col_distance((k-1)*floor(m/10)+1:k*floor(m/10));    
  cols_D32{k} = corr_col_Distance3((k-1)*floor(m/10)+1:k*floor(m/10));
  cols_dcenter2{k} = corr_col_dist_to_center((k-1)*floor(m/10)+1:k*floor(m/10));
  
 else
  
   %the last bin can contain more trials
  rows_d2{k} = corr_row_distance((k-1)*floor(m/10)+1:m);    
  rows_D32{k} = corr_row_Distance3((k-1)*floor(m/10)+1:m);
  rows_dcenter2{k} = corr_row_dist_to_center((k-1)*floor(m/10)+1:m);   
  
  cols_d2{k} = corr_col_distance((k-1)*floor(m/10)+1:m);    
  cols_D32{k} = corr_col_Distance3((k-1)*floor(m/10)+1:m);
  cols_dcenter2{k} = corr_col_dist_to_center((k-1)*floor(m/10)+1:m);   
  
 end
 
  nodeStatus_SG_cluster_d{k}=[];
  nodeStatus_SG_cluster_D3{k}=[];
  nodeStatus_SG_cluster_dcenter{k} = [];
  nodeStatus_IG_cluster_d{k}=[];
  nodeStatus_IG_cluster_D3{k}=[];
  nodeStatus_IG_cluster_dcenter{k} = [];
  
  nodeStatus_SG_cluster_d2{k}=[];
  nodeStatus_SG_cluster_D32{k}=[];
  nodeStatus_SG_cluster_dcenter2{k} = [];
  nodeStatus_IG_cluster_d2{k}=[];
  nodeStatus_IG_cluster_D32{k}=[];
  nodeStatus_IG_cluster_dcenter2{k} = [];
  
  %node status in bin k for d,D,|D|
  for j=1:length(rows_d{k})
      nodeStatus_SG_cluster_d{k} = [nodeStatus_SG_cluster_d{k}, nodeStatus_all(rows_d{k}(j), cols_d{k}(j))];
      nodeStatus_IG_cluster_d{k} = [nodeStatus_IG_cluster_d{k}, nodeStatus_IG_all(rows_d{k}(j)+start_frame-1, cols_d{k}(j))];
  end
 
  %node status is 0 if there are no trials in the bin
  if isempty(rows_d{k})
      nodeStatus_SG_cluster_d{k}=[0];
      nodeStatus_IG_cluster_d{k}=[0];
  end
 
  for j=1:length(rows_d2{k})
      nodeStatus_SG_cluster_d2{k} = [nodeStatus_SG_cluster_d2{k}, nodeStatus_all(rows_d2{k}(j), cols_d2{k}(j))];
      nodeStatus_IG_cluster_d2{k} = [nodeStatus_IG_cluster_d2{k}, nodeStatus_IG_all(rows_d2{k}(j)+start_frame-1, cols_d2{k}(j))];
  end
  
  if isempty(rows_d2{k})
      nodeStatus_SG_cluster_d2{k}=[0];
      nodeStatus_IG_cluster_d2{k}=[0];
  end
  
  for j=1:length(rows_D3{k})
      nodeStatus_SG_cluster_D3{k} = [nodeStatus_SG_cluster_D3{k}, nodeStatus_all(rows_D3{k}(j), cols_D3{k}(j))];
      nodeStatus_IG_cluster_D3{k} = [nodeStatus_IG_cluster_D3{k}, nodeStatus_IG_all(rows_D3{k}(j)+start_frame-1, cols_D3{k}(j))];
  end
  
  if isempty(rows_D3{k})
      nodeStatus_SG_cluster_D3{k}=[0];
      nodeStatus_IG_cluster_D3{k}=[0];
  end
  
  for j=1:length(rows_D32{k})    
      nodeStatus_SG_cluster_D32{k} = [nodeStatus_SG_cluster_D32{k}, nodeStatus_all(rows_D32{k}(j), cols_D32{k}(j))];
      nodeStatus_IG_cluster_D32{k} = [nodeStatus_IG_cluster_D32{k}, nodeStatus_IG_all(rows_D32{k}(j)+start_frame-1, cols_D32{k}(j))];
  end
  
   if isempty(rows_D32{k})
      nodeStatus_SG_cluster_D32{k}=[0];
      nodeStatus_IG_cluster_D32{k}=[0];
   end
  
  nodeStatus_SG_cluster_dcenter_Front{k} = [];
  nodeStatus_SG_cluster_dcenter_Middle{k} = [];
  nodeStatus_SG_cluster_dcenter_Back{k} = [];
  nodeStatus_IG_cluster_dcenter_Front{k} = [];
  nodeStatus_IG_cluster_dcenter_Middle{k} = [];
  nodeStatus_IG_cluster_dcenter_Back{k} = [];
  
  nodeStatus_SG_cluster_dcenter_Front2{k} = [];
  nodeStatus_SG_cluster_dcenter_Middle2{k} = [];
  nodeStatus_SG_cluster_dcenter_Back2{k} = [];
  nodeStatus_IG_cluster_dcenter_Front2{k} = [];
  nodeStatus_IG_cluster_dcenter_Middle2{k} = [];
  nodeStatus_IG_cluster_dcenter_Back2{k} = [];
  
  for j=1:length(rows_dcenter{k})
      
     nodeStatus_SG_cluster_dcenter{k} = [nodeStatus_SG_cluster_dcenter{k}, nodeStatus_all(rows_dcenter{k}(j), cols_dcenter{k}(j))];
     nodeStatus_IG_cluster_dcenter{k} = [nodeStatus_IG_cluster_dcenter{k}, nodeStatus_IG_all(rows_dcenter{k}(j)+start_frame-1, cols_dcenter{k}(j))];
         %divide into front, middle, back, trials of d
     if dir_node2(rows_dcenter{k}(j), cols_dcenter{k}(j))>=-pi/3 & dir_node2(rows_dcenter{k}(j), cols_dcenter{k}(j))<=pi/3
         nodeStatus_SG_cluster_dcenter_Front{k} = [nodeStatus_SG_cluster_dcenter_Front{k}, nodeStatus_all(rows_dcenter{k}(j), cols_dcenter{k}(j))];
         nodeStatus_IG_cluster_dcenter_Front{k} = [nodeStatus_IG_cluster_dcenter_Front{k}, nodeStatus_IG_all(rows_dcenter{k}(j)+start_frame-1, cols_dcenter{k}(j))];
            
     elseif (dir_node2(rows_dcenter{k}(j), cols_dcenter{k}(j))>=pi/3 & dir_node2(rows_dcenter{k}(j), cols_dcenter{k}(j))<=2*pi/3) | (dir_node2(rows_dcenter{k}(j), cols_dcenter{k}(j))<=-pi/3 & dir_node2(rows_dcenter{k}(j), cols_dcenter{k}(j))>=-2*pi/3)
         nodeStatus_SG_cluster_dcenter_Middle{k} = [nodeStatus_SG_cluster_dcenter_Middle{k}, nodeStatus_all(rows_dcenter{k}(j), cols_dcenter{k}(j))];
         nodeStatus_IG_cluster_dcenter_Middle{k} = [nodeStatus_IG_cluster_dcenter_Middle{k}, nodeStatus_IG_all(rows_dcenter{k}(j)+start_frame-1, cols_dcenter{k}(j))];
         
     elseif dir_node2(rows_dcenter{k}(j), cols_dcenter{k}(j))>=2*pi/3 | dir_node2(rows_dcenter{k}(j), cols_dcenter{k}(j))<=-2*pi/3
         
         nodeStatus_SG_cluster_dcenter_Back{k} = [nodeStatus_SG_cluster_dcenter_Back{k}, nodeStatus_all(rows_dcenter{k}(j), cols_dcenter{k}(j))];
         nodeStatus_IG_cluster_dcenter_Back{k} = [nodeStatus_IG_cluster_dcenter_Back{k}, nodeStatus_IG_all(rows_dcenter{k}(j)+start_frame-1, cols_dcenter{k}(j))];
         
     end 
  end
    
  if isempty(rows_dcenter{k})
      nodeStatus_SG_cluster_dcenter{k} = [0];
      nodeStatus_IG_cluster_dcenter{k} = [0];
  end
  
  %compute the mean of these node statuses (ns SG and then ns IG)
   nodeStatus_SG_clusterAVG_d(k) = mean(nodeStatus_SG_cluster_d{k});
   nodeStatus_SG_clusterAVG_D3(k) = mean(nodeStatus_SG_cluster_D3{k});
   nodeStatus_SG_clusterAVG_dcenter(k) = mean(nodeStatus_SG_cluster_dcenter{k});
   nodeStatus_SG_clusterAVG_dcenter_Front(k) = mean(nodeStatus_SG_cluster_dcenter_Front{k});
   nodeStatus_SG_clusterAVG_dcenter_Middle(k) = mean(nodeStatus_SG_cluster_dcenter_Middle{k});
   nodeStatus_SG_clusterAVG_dcenter_Back(k) = mean(nodeStatus_SG_cluster_dcenter_Back{k});
  
   nodeStatus_IG_clusterAVG_d(k) = mean(nodeStatus_IG_cluster_d{k});
   nodeStatus_IG_clusterAVG_D3(k) = mean(nodeStatus_IG_cluster_D3{k});
   nodeStatus_IG_clusterAVG_dcenter(k) = mean(nodeStatus_IG_cluster_dcenter{k});
   nodeStatus_IG_clusterAVG_dcenter_Front(k) = mean(nodeStatus_IG_cluster_dcenter_Front{k});
   nodeStatus_IG_clusterAVG_dcenter_Middle(k) = mean(nodeStatus_IG_cluster_dcenter_Middle{k});
   nodeStatus_IG_clusterAVG_dcenter_Back(k) = mean(nodeStatus_IG_cluster_dcenter_Back{k});
  
  %how many trials for each bin
   nodeStatus_SG_clusterSize_d(k) = length(nodeStatus_SG_cluster_d{k});
   nodeStatus_SG_clusterSize_D3(k) = length(nodeStatus_SG_cluster_D3{k});
   nodeStatus_SG_clusterSize_dcenter(k) = length(nodeStatus_SG_cluster_dcenter{k});
   nodeStatus_SG_clusterSize_dcenter_Front(k) = length(nodeStatus_SG_cluster_dcenter_Front{k});
   nodeStatus_SG_clusterSize_dcenter_Middle(k) = length(nodeStatus_SG_cluster_dcenter_Middle{k});
   nodeStatus_SG_clusterSize_dcenter_Back(k) = length(nodeStatus_SG_cluster_dcenter_Back{k});
   
  %same thing again for bins of equal number of trials
 for j=1:length(rows_dcenter2{k})
     
     nodeStatus_SG_cluster_dcenter2{k} = [nodeStatus_SG_cluster_dcenter2{k}, nodeStatus_all(rows_dcenter2{k}(j), cols_dcenter2{k}(j))];
     nodeStatus_IG_cluster_dcenter2{k} = [nodeStatus_IG_cluster_dcenter2{k}, nodeStatus_IG_all(rows_dcenter2{k}(j)+start_frame-1, cols_dcenter2{k}(j))];
 
     if dir_node2(rows_dcenter2{k}(j), cols_dcenter2{k}(j))>=-pi/3 & dir_node2(rows_dcenter2{k}(j), cols_dcenter2{k}(j))<=pi/3  
      nodeStatus_SG_cluster_dcenter_Front2{k} = [nodeStatus_SG_cluster_dcenter_Front2{k}, nodeStatus_all(rows_dcenter2{k}(j), cols_dcenter2{k}(j))];
      nodeStatus_IG_cluster_dcenter_Front2{k} = [nodeStatus_IG_cluster_dcenter_Front2{k}, nodeStatus_IG_all(rows_dcenter2{k}(j)+start_frame-1, cols_dcenter2{k}(j))];
   elseif (dir_node2(rows_dcenter2{k}(j), cols_dcenter2{k}(j))>=pi/3 & dir_node2(rows_dcenter2{k}(j), cols_dcenter2{k}(j))<=2*pi/3) | (dir_node2(rows_dcenter2{k}(j), cols_dcenter2{k}(j))<=-pi/3 & dir_node2(rows_dcenter2{k}(j), cols_dcenter2{k}(j))>=-2*pi/3)    
      nodeStatus_SG_cluster_dcenter_Middle2{k} = [nodeStatus_SG_cluster_dcenter_Middle2{k}, nodeStatus_all(rows_dcenter2{k}(j), cols_dcenter2{k}(j))];
      nodeStatus_IG_cluster_dcenter_Middle2{k} = [nodeStatus_IG_cluster_dcenter_Middle2{k}, nodeStatus_IG_all(rows_dcenter2{k}(j)+start_frame-1, cols_dcenter2{k}(j))];
   elseif dir_node2(rows_dcenter2{k}(j), cols_dcenter2{k}(j))>=2*pi/3 | dir_node2(rows_dcenter2{k}(j), cols_dcenter2{k}(j))<=-2*pi/3   
      nodeStatus_SG_cluster_dcenter_Back2{k} = [nodeStatus_SG_cluster_dcenter_Back2{k}, nodeStatus_all(rows_dcenter2{k}(j), cols_dcenter2{k}(j))];
      nodeStatus_IG_cluster_dcenter_Back2{k} = [nodeStatus_IG_cluster_dcenter_Back2{k}, nodeStatus_IG_all(rows_dcenter2{k}(j)+start_frame-1, cols_dcenter2{k}(j))];
     end
 end
 
 if isempty(rows_dcenter2{k})

     nodeStatus_SG_cluster_dcenter2{k} = [0];
      nodeStatus_IG_cluster_dcenter2{k} = [0];
 end
  
   nodeStatus_SG_clusterAVG_d2(k) = mean(nodeStatus_SG_cluster_d2{k});
   nodeStatus_SG_clusterAVG_D32(k) = mean(nodeStatus_SG_cluster_D32{k});
   nodeStatus_SG_clusterAVG_dcenter2(k) = mean(nodeStatus_SG_cluster_dcenter2{k});
   nodeStatus_SG_clusterAVG_dcenter_Front2(k) = mean(nodeStatus_SG_cluster_dcenter_Front2{k});
   nodeStatus_SG_clusterAVG_dcenter_Middle2(k) = mean(nodeStatus_SG_cluster_dcenter_Middle2{k});
   nodeStatus_SG_clusterAVG_dcenter_Back2(k) = mean(nodeStatus_SG_cluster_dcenter_Back2{k});
   
   nodeStatus_IG_clusterAVG_d2(k) = mean(nodeStatus_IG_cluster_d2{k});
   nodeStatus_IG_clusterAVG_D32(k) = mean(nodeStatus_IG_cluster_D32{k});
   nodeStatus_IG_clusterAVG_dcenter2(k) = mean(nodeStatus_IG_cluster_dcenter2{k});
   nodeStatus_IG_clusterAVG_dcenter_Front2(k) = mean(nodeStatus_IG_cluster_dcenter_Front2{k});
   nodeStatus_IG_clusterAVG_dcenter_Middle2(k) = mean(nodeStatus_IG_cluster_dcenter_Middle2{k});
   nodeStatus_IG_clusterAVG_dcenter_Back2(k) = mean(nodeStatus_IG_cluster_dcenter_Back2{k});
   
   nodeStatus_SG_clusterSize_d2(k) = length(nodeStatus_SG_cluster_d2{k});
   nodeStatus_SG_clusterSize_D32(k) = length(nodeStatus_SG_cluster_D32{k});
   nodeStatus_SG_clusterSize_dcenter2(k) = length(nodeStatus_SG_cluster_dcenter2{k});
   nodeStatus_SG_clusterSize_dcenter_Front2(k) = length(nodeStatus_SG_cluster_dcenter_Front2{k});
   nodeStatus_SG_clusterSize_dcenter_Middle2(k) = length(nodeStatus_SG_cluster_dcenter_Middle2{k});
   nodeStatus_SG_clusterSize_dcenter_Back2(k) = length(nodeStatus_SG_cluster_dcenter_Back2{k});
  
 end
 
  %how statistically significant these results are
 for k=1:10
     for j=1:10
       ss_D3_SG(k,j) = kstest2(nodeStatus_SG_cluster_D3{k}, nodeStatus_SG_cluster_D3{j});  
       ss_d_SG(k,j) = kstest2(nodeStatus_SG_cluster_d{k}, nodeStatus_SG_cluster_d{j});
       ss_dCenter_SG(k,j) = kstest2(nodeStatus_SG_cluster_dcenter{k}, nodeStatus_SG_cluster_dcenter{j});
       
       ss_D3_SG2(k,j) = kstest2(nodeStatus_SG_cluster_D32{k}, nodeStatus_SG_cluster_D32{j});  
       ss_d_SG2(k,j) = kstest2(nodeStatus_SG_cluster_d2{k}, nodeStatus_SG_cluster_d2{j});
       ss_dCenter_SG2(k,j) = kstest2(nodeStatus_SG_cluster_dcenter2{k}, nodeStatus_SG_cluster_dcenter2{j});
       
     if ~isempty(nodeStatus_SG_cluster_dcenter_Front{k})  & ~isempty(nodeStatus_SG_cluster_dcenter_Front{j})
       ss_dCenter_front_SG(k,j) = kstest2(nodeStatus_SG_cluster_dcenter_Front{k}, nodeStatus_SG_cluster_dcenter_Front{j});
     end
     
     if ~isempty(nodeStatus_SG_cluster_dcenter_Front2{k})  & ~isempty(nodeStatus_SG_cluster_dcenter_Front2{j})
         ss_dCenter_front_SG2(k,j) = kstest2(nodeStatus_SG_cluster_dcenter_Front2{k}, nodeStatus_SG_cluster_dcenter_Front2{j});
     end
     
     if ~isempty(nodeStatus_SG_cluster_dcenter_Middle{k})  & ~isempty(nodeStatus_SG_cluster_dcenter_Middle{j})
       ss_dCenter_middle_SG(k,j) = kstest2(nodeStatus_SG_cluster_dcenter_Middle{k}, nodeStatus_SG_cluster_dcenter_Middle{j});
     end
     
     if ~isempty(nodeStatus_SG_cluster_dcenter_Middle2{k})  & ~isempty(nodeStatus_SG_cluster_dcenter_Middle2{j})
       ss_dCenter_middle_SG2(k,j) = kstest2(nodeStatus_SG_cluster_dcenter_Middle2{k}, nodeStatus_SG_cluster_dcenter_Middle2{j});
     end
     
     if ~isempty(nodeStatus_SG_cluster_dcenter_Back{k}) & ~isempty(nodeStatus_SG_cluster_dcenter_Back{j})
       ss_dCenter_back_SG(k,j) = kstest2(nodeStatus_SG_cluster_dcenter_Back{k}, nodeStatus_SG_cluster_dcenter_Back{j});
     end  
       
     if ~isempty(nodeStatus_SG_cluster_dcenter_Back2{k}) & ~isempty(nodeStatus_SG_cluster_dcenter_Back2{j})
       ss_dCenter_back_SG2(k,j) = kstest2(nodeStatus_SG_cluster_dcenter_Back2{k}, nodeStatus_SG_cluster_dcenter_Back2{j});
     end
     
       ss_D3_IG(k,j) = kstest2(nodeStatus_IG_cluster_D3{k}, nodeStatus_IG_cluster_D3{j});  
       ss_d_IG(k,j) = kstest2(nodeStatus_IG_cluster_d{k}, nodeStatus_IG_cluster_d{j});
       ss_dCenter_IG(k,j) = kstest2(nodeStatus_IG_cluster_dcenter{k}, nodeStatus_IG_cluster_dcenter{j});
       
       ss_D3_IG2(k,j) = kstest2(nodeStatus_IG_cluster_D32{k}, nodeStatus_IG_cluster_D32{j});  
       ss_d_IG2(k,j) = kstest2(nodeStatus_IG_cluster_d2{k}, nodeStatus_IG_cluster_d2{j});
       ss_dCenter_IG2(k,j) = kstest2(nodeStatus_IG_cluster_dcenter2{k}, nodeStatus_IG_cluster_dcenter2{j});
       
      if ~isempty(nodeStatus_SG_cluster_dcenter_Front{k})  & ~isempty(nodeStatus_SG_cluster_dcenter_Front{j})
       ss_dCenter_front_IG(k,j) = kstest2(nodeStatus_IG_cluster_dcenter_Front{k}, nodeStatus_IG_cluster_dcenter_Front{j});
      end
      
      if ~isempty(nodeStatus_SG_cluster_dcenter_Front2{k})  & ~isempty(nodeStatus_SG_cluster_dcenter_Front2{j})
       ss_dCenter_front_IG2(k,j) = kstest2(nodeStatus_IG_cluster_dcenter_Front2{k}, nodeStatus_IG_cluster_dcenter_Front2{j});
      end
      
      if ~isempty(nodeStatus_IG_cluster_dcenter_Middle{k})  & ~isempty(nodeStatus_IG_cluster_dcenter_Middle{j})  
       ss_dCenter_middle_IG(k,j) = kstest2(nodeStatus_IG_cluster_dcenter_Middle{k}, nodeStatus_IG_cluster_dcenter_Middle{j});
      end
      
      if ~isempty(nodeStatus_IG_cluster_dcenter_Middle2{k})  & ~isempty(nodeStatus_IG_cluster_dcenter_Middle2{j})  
        ss_dCenter_middle_IG2(k,j) = kstest2(nodeStatus_IG_cluster_dcenter_Middle2{k}, nodeStatus_IG_cluster_dcenter_Middle2{j});
      end
      
      if ~isempty(nodeStatus_IG_cluster_dcenter_Back{k}) & ~isempty(nodeStatus_IG_cluster_dcenter_Back{j})
       ss_dCenter_back_IG(k,j) = kstest2(nodeStatus_IG_cluster_dcenter_Back{k}, nodeStatus_IG_cluster_dcenter_Back{j});
      end
      
      if ~isempty(nodeStatus_IG_cluster_dcenter_Back2{k}) & ~isempty(nodeStatus_IG_cluster_dcenter_Back2{j})
       ss_dCenter_back_IG2(k,j) = kstest2(nodeStatus_IG_cluster_dcenter_Back2{k}, nodeStatus_IG_cluster_dcenter_Back2{j});
      end
      
     end
 end

 %plot bins for d,D,|D|
  
 figure(9)
 bar(nodeStatus_SG_clusterAVG_d2)
 title('Average Node status SG for 10 bins of distance = absolute value of d to center-axis')
 xlabel('bins')
 ylabel('average node status SG')
 
 figure(10)
 bar(nodeStatus_SG_clusterAVG_D32)
 title('Average node status SG for 10 bins of distance = d to center-axis')
 xlabel('bins')
 ylabel('average node status SG')
 
 figure(11)
 bar(nodeStatus_SG_clusterAVG_dcenter2)
 title('Average node status SG for 10 bins of distance = d to center (point)')
 xlabel('bins')
 ylabel('average node status SG')
 
 color_spec = [0 1 0; 0 0 1; 1 0 0];
 figure(12)
 bar(1:10, [nodeStatus_SG_clusterAVG_dcenter_Front2',nodeStatus_SG_clusterAVG_dcenter_Middle2', nodeStatus_SG_clusterAVG_dcenter_Back2'], 0.5, 'stack') 
 legend('Front', 'Middle', 'Back')
 colormap(color_spec)
 
 figure(13)
 bar(nodeStatus_IG_clusterAVG_d2)
 title('Average Node status IG for 10 bins of distance = absolute value of d to center-axis')
 xlabel('bins')
 ylabel('average node status IG')
 
 figure(14)
 bar(nodeStatus_IG_clusterAVG_D32)
 title('Average node status IG for 10 bins of distance = d to center-axis')
 xlabel('bins')
 ylabel('avg node status IG')
 
 figure(15)
 bar(nodeStatus_IG_clusterAVG_dcenter2)
 title('Average node status IG for 10 bins of distance = d to center (point)')
 xlabel('bins')
 ylabel('average node status IG')
 
 color_spec = [0 1 0; 0 0 1; 1 0 0];
 figure(16)
 bar(1:10, [nodeStatus_IG_clusterAVG_dcenter_Front2',nodeStatus_IG_clusterAVG_dcenter_Middle2', nodeStatus_IG_clusterAVG_dcenter_Back2'], 0.5, 'stack') 
 legend('Front', 'Middle', 'Back')
 colormap(color_spec)
 
  
  %how many trials per bin
 figure(17)
 bar(nodeStatus_SG_clusterSize_d2)
 title('How many trials per 10 bins of distance = absolute value of d to center-axis')
 xlabel('bins')
 ylabel('number of trials')
 
 figure(18)
 bar(nodeStatus_SG_clusterSize_D32)
 title('How many trials per 10 bins of distance = d to center-axis')
 xlabel('bins')
 ylabel('number of trials')
 
 figure(19)
 bar(nodeStatus_SG_clusterSize_dcenter2)
 title('How many trials per 10 bins of distance = d to center (point)')
 xlabel('bins')
 ylabel('number of trials')
 
 color_spec = [0 1 0; 0 0 1; 1 0 0];
 figure(20)
 bar(1:10, [nodeStatus_SG_clusterSize_dcenter_Front2',nodeStatus_SG_clusterSize_dcenter_Middle2', nodeStatus_SG_clusterSize_dcenter_Back2'], 0.5, 'stack') 
 legend('Front', 'Middle', 'Back')
 colormap(color_spec)
 %}
  disp('statistical results on the bins')
  
  ss_D3_SG  
  ss_d_SG
  ss_dCenter_SG
  ss_dCenter_front_SG
  ss_dCenter_middle_SG
  ss_dCenter_back_SG
       
  ss_D3_IG 
  ss_d_IG
  ss_dCenter_IG
  ss_dCenter_front_IG
  ss_dCenter_middle_IG
  ss_dCenter_back_IG
 
  %PART II: bins contain data for only ONE node of choice (node)
  %the same thing is done here as for PART I
  
  max_D3_forNode = max(Distance3(:, node));
  min_D3_forNode = min(Distance3(:, node));
  max_d_forNode  = max(distance(:, node));
  min_d_forNode = min(distance(:, node));
  min_dist_to_center_forNode = min(dist_to_center2(:, node));
  max_dist_to_center_forNode = max(dist_to_center2(:, node));
  
  D3_interval_forNode = (max_D3_forNode-min_D3_forNode)/10;
  d_interval_forNode = (max_d_forNode-min_d_forNode)/10;
  dcenter_interval_forNode = (max_dist_to_center_forNode - min_dist_to_center_forNode)/10;
  
  ss_D3_SG_forNode = NaN*ones(10,10);
  ss_d_SG_forNode = NaN*ones(10,10);
  ss_dCenter_SG_forNode = NaN*ones(10,10);
  ss_dCenter_front_SG_forNode = NaN*ones(10,10);
  ss_dCenter_middle_SG_forNode = NaN*ones(10,10);
  ss_dCenter_back_SG_forNode = NaN*ones(10,10);
  
  ss_D3_IG_forNode = NaN*ones(10,10);
  ss_d_IG_forNode = NaN*ones(10,10);
  ss_dCenter_IG_forNode = NaN*ones(10,10);
  ss_dCenter_front_IG_forNode = NaN*ones(10,10);
  ss_dCenter_middle_IG_forNode = NaN*ones(10,10);
  ss_dCenter_back_IG_forNode = NaN*ones(10,10);
  
  dist_to_center_longV = dist_to_center2(:,node);
  distance_longV = distance(:,node);
  Distance3_longV = Distance3(:,node);
  
  [dist_to_center_longV ind_dist_to_center_longV]= sort(dist_to_center_longV);
  [distance_longV ind_distance_longV]= sort(distance_longV);
  [Distance3_longV ind_Distance3_longV] = sort(Distance3_longV);
 
  m = length(distance_longV);
  
 for k=1:10    %quantifying for one node, how the distance from center influences node status
     %k is the bin number, there are 10 bins
     %rows are frames, cols are nodes
     
  cols_d_sameBin{k} = find(cols_d{k}==node);
  rows_d_sameBin{k} = rows_d{k}(cols_d_sameBin{k});
  cols_D3_sameBin{k} = find(cols_D3{k}==node);
  rows_D3_sameBin{k} = rows_D3{k}(cols_D3_sameBin{k});
  cols_dcenter_sameBin{k} = find(cols_dcenter{k}==node);
  rows_dcenter_sameBin{k} = rows_dcenter{k}(cols_dcenter_sameBin{k});
  
  [rows_d{k}] = find(distance(:,node)>=min_d_forNode+(k-1)*d_interval_forNode & distance(:,node)<=min_d_forNode+k*d_interval_forNode);    
  [rows_D3{k}] = find(Distance3(:,node)>=min_D3_forNode+(k-1)*D3_interval_forNode & Distance3(:,node)<=min_D3_forNode+k*D3_interval_forNode);
  [rows_dcenter{k}] = find(dist_to_center2(:,node) >= min_dist_to_center_forNode + (k-1)*dcenter_interval_forNode & dist_to_center2(:,node) <= min_dist_to_center_forNode+k*dcenter_interval_forNode);
  
 if k~=10 
  index_d2{k} = ind_distance_longV((k-1)*floor(m/10)+1:k*floor(m/10));
  index_D32{k} = ind_Distance3_longV((k-1)*floor(m/10)+1:k*floor(m/10));
  index_dcenter2{k} = ind_dist_to_center_longV((k-1)*floor(m/10)+1:k*floor(m/10));
 else 
  index_d2{k} = ind_distance_longV((k-1)*floor(m/10)+1:m);
  index_D32{k} = ind_Distance3_longV((k-1)*floor(m/10)+1:m);
  index_dcenter2{k} = ind_dist_to_center_longV((k-1)*floor(m/10)+1:m);   
 end

  cols_d2_sameBin{k} = find(cols_d2{k}==node);
  rows_d2_sameBin{k}=rows_d2{k}(cols_d2_sameBin{k});
  cols_D32_sameBin{k} = find(cols_D32{k}==node);
  rows_D32_sameBin{k} = rows_D32{k}(cols_D32_sameBin{k});
  cols_dcenter2_sameBin{k} = find(cols_dcenter2{k}==node);
  rows_dcenter2_sameBin{k} = rows_d2{k}(cols_dcenter2_sameBin{k});  

  nodeStatus_SG_cluster_d_forNode{k}=[];
  nodeStatus_SG_cluster_D3_forNode{k}=[];
  nodeStatus_SG_cluster_dcenter_forNode{k} = [];
  nodeStatus_IG_cluster_d_forNode{k}=[];
  nodeStatus_IG_cluster_D3_forNode{k}=[];
  nodeStatus_IG_cluster_dcenter_forNode{k} = [];
  
  nodeStatus_SG_cluster_d_forNode_sameBin{k}=[];
  nodeStatus_SG_cluster_D3_forNode_sameBin{k}=[];
  nodeStatus_SG_cluster_dcenter_forNode_sameBin{k} = [];
  nodeStatus_IG_cluster_d_forNode_sameBin{k}=[];
  nodeStatus_IG_cluster_D3_forNode_sameBin{k}=[];
  nodeStatus_IG_cluster_dcenter_forNode_sameBin{k} = [];
  
  nodeStatus_SG_cluster_d_forNode2{k}=[];
  nodeStatus_SG_cluster_D3_forNode2{k}=[];
  nodeStatus_SG_cluster_dcenter_forNode2{k} = [];
  nodeStatus_IG_cluster_d_forNode2{k}=[];
  nodeStatus_IG_cluster_D3_forNode2{k}=[];
  nodeStatus_IG_cluster_dcenter_forNode2{k} = [];
  
  nodeStatus_SG_cluster_d_forNode2_sameBin{k}=[];
  nodeStatus_SG_cluster_D3_forNode2_sameBin{k}=[];
  nodeStatus_SG_cluster_dcenter_forNode2_sameBin{k} = [];
  nodeStatus_IG_cluster_d_forNode2_sameBin{k}=[];
  nodeStatus_IG_cluster_D3_forNode2_sameBin{k}=[];
  nodeStatus_IG_cluster_dcenter_forNode2_sameBin{k} = [];
  
  for j=1:length(rows_d{k})
      nodeStatus_SG_cluster_d_forNode{k} = [nodeStatus_SG_cluster_d_forNode{k}, nodeStatus_all(rows_d{k}(j), node)];
      nodeStatus_IG_cluster_d_forNode{k} = [nodeStatus_IG_cluster_d_forNode{k}, nodeStatus_IG_all(rows_d{k}(j)+start_frame-1, node)];
  end
  
  if isempty(rows_d{k})
      nodeStatus_SG_cluster_d_forNode{k} = [0];
      nodeStatus_IG_cluster_d_forNode{k} = [0];
  end
  
  for j=1:length(rows_d_sameBin{k})
      nodeStatus_SG_cluster_d_forNode_sameBin{k} = [nodeStatus_SG_cluster_d_forNode_sameBin{k}, nodeStatus_all(rows_d_sameBin{k}(j), node)];
      nodeStatus_IG_cluster_d_forNode_sameBin{k} = [nodeStatus_IG_cluster_d_forNode_sameBin{k}, nodeStatus_IG_all(rows_d_sameBin{k}(j)+start_frame-1, node)];
  end
  
  if isempty(rows_d_sameBin{k})
      nodeStatus_SG_cluster_d_forNode_sameBin{k} = [0];
      nodeStatus_IG_cluster_d_forNode_sameBin{k} = [0];
  end
  
  for j=1:length(index_d2{k})
      nodeStatus_SG_cluster_d_forNode2{k} = [nodeStatus_SG_cluster_d_forNode2{k}, nodeStatus_all(index_d2{k}(j), node)];
      nodeStatus_IG_cluster_d_forNode2{k} = [nodeStatus_IG_cluster_d_forNode2{k}, nodeStatus_IG_all(index_d2{k}(j)+start_frame-1, node)];
  end
  
  if isempty(index_d2{k})
      nodeStatus_SG_cluster_d_forNode2{k} = [0];
      nodeStatus_IG_cluster_d_forNode2{k} = [0];
  end
  
  for j=1:length(rows_d2_sameBin{k})
      nodeStatus_SG_cluster_d_forNode2_sameBin{k} = [nodeStatus_SG_cluster_d_forNode2_sameBin{k}, nodeStatus_all(rows_d2_sameBin{k}(j), node)];
      nodeStatus_IG_cluster_d_forNode2_sameBin{k} = [nodeStatus_IG_cluster_d_forNode2_sameBin{k}, nodeStatus_IG_all(rows_d2_sameBin{k}(j)+start_frame-1, node)];
  end
 
  if isempty(rows_d2_sameBin{k})
      nodeStatus_SG_cluster_d_forNode2_sameBin{k} = [0];
      nodeStatus_IG_cluster_d_forNode2_sameBin{k} = [0];
  end
  
  for j=1:length(rows_D3{k})
      nodeStatus_SG_cluster_D3_forNode{k} = [nodeStatus_SG_cluster_D3_forNode{k}, nodeStatus_all(rows_D3{k}(j), node)];
      nodeStatus_IG_cluster_D3_forNode{k} = [nodeStatus_IG_cluster_D3_forNode{k}, nodeStatus_IG_all(rows_D3{k}(j)+start_frame-1, node)];
  end
  
  if isempty(rows_D3{k})
      nodeStatus_SG_cluster_D3_forNode{k} = [0];
      nodeStatus_IG_cluster_D3_forNode{k} = [0];
  end
  
  for j=1:length(rows_D3_sameBin{k})
      nodeStatus_SG_cluster_D3_forNode_sameBin{k} = [nodeStatus_SG_cluster_D3_forNode_sameBin{k}, nodeStatus_all(rows_D3_sameBin{k}(j), node)];
      nodeStatus_IG_cluster_D3_forNode_sameBin{k} = [nodeStatus_IG_cluster_D3_forNode_sameBin{k}, nodeStatus_IG_all(rows_D3_sameBin{k}(j)+start_frame-1, node)];
  end
  
  if isempty(rows_D3_sameBin{k})
      nodeStatus_SG_cluster_D3_forNode_sameBin{k} = [0];
      nodeStatus_IG_cluster_D3_forNode_sameBin{k} = [0];
  end
  
  for j=1:length(index_D32{k})
      nodeStatus_SG_cluster_D3_forNode2{k} = [nodeStatus_SG_cluster_D3_forNode2{k}, nodeStatus_all(index_D32{k}(j), node)];
      nodeStatus_IG_cluster_D3_forNode2{k} = [nodeStatus_IG_cluster_D3_forNode2{k}, nodeStatus_IG_all(index_D32{k}(j)+start_frame-1, node)];
  end
  
  if isempty(index_D32{k})
      nodeStatus_SG_cluster_D3_forNode2{k} = [0];
      nodeStatus_IG_cluster_D3_forNode2{k} = [0];
  end
  
  for j=1:length(rows_D32_sameBin{k})
      nodeStatus_SG_cluster_D3_forNode2_sameBin{k} = [nodeStatus_SG_cluster_D3_forNode2_sameBin{k}, nodeStatus_all(rows_D32_sameBin{k}(j), node)];
      nodeStatus_IG_cluster_D3_forNode2_sameBin{k} = [nodeStatus_IG_cluster_D3_forNode2_sameBin{k}, nodeStatus_IG_all(rows_D32_sameBin{k}(j)+start_frame-1, node)];
  end
  
  if isempty(rows_D32_sameBin{k})
      nodeStatus_SG_cluster_D3_forNode2_sameBin{k} = [0];
      nodeStatus_IG_cluster_D3_forNode2_sameBin{k} = [0];
  end
  
  nodeStatus_SG_cluster_dcenter_Front_forNode{k} = [];
  nodeStatus_SG_cluster_dcenter_Middle_forNode{k} = [];
  nodeStatus_SG_cluster_dcenter_Back_forNode{k} = [];
  nodeStatus_IG_cluster_dcenter_Front_forNode{k} = [];
  nodeStatus_IG_cluster_dcenter_Middle_forNode{k} = [];
  nodeStatus_IG_cluster_dcenter_Back_forNode{k} = [];
  
  nodeStatus_SG_cluster_dcenter_Front_forNode_sameBin{k} = [];
  nodeStatus_SG_cluster_dcenter_Middle_forNode_sameBin{k} = [];
  nodeStatus_SG_cluster_dcenter_Back_forNode_sameBin{k} = [];
  nodeStatus_IG_cluster_dcenter_Front_forNode_sameBin{k} = [];
  nodeStatus_IG_cluster_dcenter_Middle_forNode_sameBin{k} = [];
  nodeStatus_IG_cluster_dcenter_Back_forNode_sameBin{k} = [];
 
  nodeStatus_SG_cluster_dcenter_Front_forNode2{k} = [];
  nodeStatus_SG_cluster_dcenter_Middle_forNode2{k} = [];
  nodeStatus_SG_cluster_dcenter_Back_forNode2{k} = [];
  nodeStatus_IG_cluster_dcenter_Front_forNode2{k} = [];
  nodeStatus_IG_cluster_dcenter_Middle_forNode2{k} = [];
  nodeStatus_IG_cluster_dcenter_Back_forNode2{k} = [];
 
  nodeStatus_SG_cluster_dcenter_Front_forNode2_sameBin{k} = [];
  nodeStatus_SG_cluster_dcenter_Middle_forNode2_sameBin{k} = [];
  nodeStatus_SG_cluster_dcenter_Back_forNode2_sameBin{k} = [];
  nodeStatus_IG_cluster_dcenter_Front_forNode2_sameBin{k} = [];
  nodeStatus_IG_cluster_dcenter_Middle_forNode2_sameBin{k} = [];
  nodeStatus_IG_cluster_dcenter_Back_forNode2_sameBin{k} = [];

  for j=1:length(rows_dcenter{k})
      
     nodeStatus_SG_cluster_dcenter_forNode{k} = [nodeStatus_SG_cluster_dcenter_forNode{k}, nodeStatus_all(rows_dcenter{k}(j), node)];
     nodeStatus_IG_cluster_dcenter_forNode{k} = [nodeStatus_IG_cluster_dcenter_forNode{k}, nodeStatus_IG_all(rows_dcenter{k}(j)+start_frame-1, node)];
     
     if dir_node2(rows_dcenter{k}(j), node)>=-pi/3 & dir_node2(rows_dcenter{k}(j), node)<=pi/3
         nodeStatus_SG_cluster_dcenter_Front_forNode{k} = [nodeStatus_SG_cluster_dcenter_Front_forNode{k}, nodeStatus_all(rows_dcenter{k}(j), node)];
         nodeStatus_IG_cluster_dcenter_Front_forNode{k} = [nodeStatus_IG_cluster_dcenter_Front_forNode{k}, nodeStatus_IG_all(rows_dcenter{k}(j)+start_frame-1, node)];
     elseif (dir_node2(rows_dcenter{k}(j), node)>=pi/3 & dir_node2(rows_dcenter{k}(j), node)<=2*pi/3) | (dir_node2(rows_dcenter{k}(j),node)<=-pi/3 & dir_node2(rows_dcenter{k}(j), node)>=-2*pi/3)
         nodeStatus_SG_cluster_dcenter_Middle_forNode{k} = [nodeStatus_SG_cluster_dcenter_Middle_forNode{k}, nodeStatus_all(rows_dcenter{k}(j), node)];
         nodeStatus_IG_cluster_dcenter_Middle_forNode{k} = [nodeStatus_IG_cluster_dcenter_Middle_forNode{k}, nodeStatus_IG_all(rows_dcenter{k}(j)+start_frame-1, node)];
     elseif dir_node2(rows_dcenter{k}(j), node)>=2*pi/3 | dir_node2(rows_dcenter{k}(j), node)<=-2*pi/3
         nodeStatus_SG_cluster_dcenter_Back_forNode{k} = [nodeStatus_SG_cluster_dcenter_Back_forNode{k}, nodeStatus_all(rows_dcenter{k}(j), node)];
         nodeStatus_IG_cluster_dcenter_Back_forNode{k} = [nodeStatus_IG_cluster_dcenter_Back_forNode{k}, nodeStatus_IG_all(rows_dcenter{k}(j)+start_frame-1, node)];
     end
  end
  
  if isempty(rows_dcenter{k})
      nodeStatus_SG_cluster_dcenter_forNode{k} = [0];
      nodeStatus_IG_cluster_dcenter_forNode{k} = [0];
  end
  
   nodeStatus_SG_clusterAVG_d_forNode(k) = mean(nodeStatus_SG_cluster_d_forNode{k});
   nodeStatus_SG_clusterAVG_D3_forNode(k) = mean(nodeStatus_SG_cluster_D3_forNode{k});
   nodeStatus_SG_clusterAVG_dcenter_forNode(k) = mean(nodeStatus_SG_cluster_dcenter_forNode{k});
   nodeStatus_SG_clusterAVG_dcenter_Front_forNode(k) = mean(nodeStatus_SG_cluster_dcenter_Front_forNode{k});
   nodeStatus_SG_clusterAVG_dcenter_Middle_forNode(k) = mean(nodeStatus_SG_cluster_dcenter_Middle_forNode{k});
   nodeStatus_SG_clusterAVG_dcenter_Back_forNode(k) = mean(nodeStatus_SG_cluster_dcenter_Back_forNode{k});
   
   nodeStatus_IG_clusterAVG_d_forNode(k) = mean(nodeStatus_IG_cluster_d_forNode{k});
   nodeStatus_IG_clusterAVG_D3_forNode(k) = mean(nodeStatus_IG_cluster_D3_forNode{k});
   nodeStatus_IG_clusterAVG_dcenter_forNode(k) = mean(nodeStatus_IG_cluster_dcenter_forNode{k});
   nodeStatus_IG_clusterAVG_dcenter_Front_forNode(k) = mean(nodeStatus_IG_cluster_dcenter_Front_forNode{k});
   nodeStatus_IG_clusterAVG_dcenter_Middle_forNode(k) = mean(nodeStatus_IG_cluster_dcenter_Middle_forNode{k});
   nodeStatus_IG_clusterAVG_dcenter_Back_forNode(k) = mean(nodeStatus_IG_cluster_dcenter_Back_forNode{k});
   
   nodeStatus_SG_clusterSize_d_forNode(k) = length(nodeStatus_SG_cluster_d_forNode{k});
   nodeStatus_SG_clusterSize_D3_forNode(k) = length(nodeStatus_SG_cluster_D3_forNode{k});
   nodeStatus_SG_clusterSize_dcenter_forNode(k) = length(nodeStatus_SG_cluster_dcenter_forNode{k});
   nodeStatus_SG_clusterSize_dcenter_Front_forNode(k) = length(nodeStatus_SG_cluster_dcenter_Front_forNode{k});
   nodeStatus_SG_clusterSize_dcenter_Middle_forNode(k) = length(nodeStatus_SG_cluster_dcenter_Middle_forNode{k});
   nodeStatus_SG_clusterSize_dcenter_Back_forNode(k) = length(nodeStatus_SG_cluster_dcenter_Back_forNode{k});
   
  for j=1:length(rows_dcenter_sameBin{k})
      
     nodeStatus_SG_cluster_dcenter_forNode_sameBin{k} = [nodeStatus_SG_cluster_dcenter_forNode_sameBin{k}, nodeStatus_all(rows_dcenter_sameBin{k}(j), node)];
     nodeStatus_IG_cluster_dcenter_forNode_sameBin{k} = [nodeStatus_IG_cluster_dcenter_forNode_sameBin{k}, nodeStatus_IG_all(rows_dcenter_sameBin{k}(j)+start_frame-1, node)];    
     
     if dir_node2(rows_dcenter_sameBin{k}(j), node)>=-pi/3 & dir_node2(rows_dcenter_sameBin{k}(j), node)<=pi/3
         nodeStatus_SG_cluster_dcenter_Front_forNode_sameBin{k} = [nodeStatus_SG_cluster_dcenter_Front_forNode_sameBin{k}, nodeStatus_all(rows_dcenter_sameBin{k}(j), node)];
         nodeStatus_IG_cluster_dcenter_Front_forNode_sameBin{k} = [nodeStatus_IG_cluster_dcenter_Front_forNode_sameBin{k}, nodeStatus_IG_all(rows_dcenter_sameBin{k}(j)+start_frame-1, node)];
            
     elseif (dir_node2(rows_dcenter_sameBin{k}(j), node)>=pi/3 & dir_node2(rows_dcenter_sameBin{k}(j), node)<=2*pi/3) | (dir_node2(rows_dcenter_sameBin{k}(j), node)<=-pi/3 & dir_node2(rows_dcenter_sameBin{k}(j), node)>=-2*pi/3)
         nodeStatus_SG_cluster_dcenter_Middle_forNode_sameBin{k} = [nodeStatus_SG_cluster_dcenter_Middle_forNode_sameBin{k}, nodeStatus_all(rows_dcenter_sameBin{k}(j),node)];
         nodeStatus_IG_cluster_dcenter_Middle_forNode_sameBin{k} = [nodeStatus_IG_cluster_dcenter_Middle_forNode_sameBin{k}, nodeStatus_IG_all(rows_dcenter_sameBin{k}(j)+start_frame-1, node)];
         
     elseif dir_node2(rows_dcenter_sameBin{k}(j), node)>=2*pi/3 | dir_node2(rows_dcenter_sameBin{k}(j), node)<=-2*pi/3
         nodeStatus_SG_cluster_dcenter_Back_forNode_sameBin{k} = [nodeStatus_SG_cluster_dcenter_Back_forNode_sameBin{k}, nodeStatus_all(rows_dcenter_sameBin{k}(j), node)];
         nodeStatus_IG_cluster_dcenter_Back_forNode_sameBin{k} = [nodeStatus_IG_cluster_dcenter_Back_forNode_sameBin{k}, nodeStatus_IG_all(rows_dcenter_sameBin{k}(j)+start_frame-1, node)];
         
     end 
  end
 
  if isempty(rows_dcenter_sameBin{k})
      nodeStatus_SG_cluster_dcenter_forNode_sameBin{k} = [0];
      nodeStatus_IG_cluster_dcenter_forNode_sameBin{k} = [0];
  end
  
   nodeStatus_SG_clusterAVG_d_sameBin(k) = mean(nodeStatus_SG_cluster_d_forNode_sameBin{k});
   nodeStatus_SG_clusterAVG_D3_sameBin(k) = mean(nodeStatus_SG_cluster_D3_forNode_sameBin{k});
   nodeStatus_SG_clusterAVG_dcenter_sameBin(k) = mean(nodeStatus_SG_cluster_dcenter_forNode_sameBin{k});
   nodeStatus_SG_clusterAVG_dcenter_Front_sameBin(k) = mean(nodeStatus_SG_cluster_dcenter_Front_forNode_sameBin{k});
   nodeStatus_SG_clusterAVG_dcenter_Middle_sameBin(k) = mean(nodeStatus_SG_cluster_dcenter_Middle_forNode_sameBin{k});
   nodeStatus_SG_clusterAVG_dcenter_Back_sameBin(k) = mean(nodeStatus_SG_cluster_dcenter_Back_forNode_sameBin{k});
  
   nodeStatus_IG_clusterAVG_d_sameBin(k) = mean(nodeStatus_IG_cluster_d_forNode_sameBin{k});
   nodeStatus_IG_clusterAVG_D3_sameBin(k) = mean(nodeStatus_IG_cluster_D3_forNode_sameBin{k});
   nodeStatus_IG_clusterAVG_dcenter_sameBin(k) = mean(nodeStatus_IG_cluster_dcenter_forNode_sameBin{k});
   nodeStatus_IG_clusterAVG_dcenter_Front_sameBin(k) = mean(nodeStatus_IG_cluster_dcenter_Front_forNode_sameBin{k});
   nodeStatus_IG_clusterAVG_dcenter_Middle_sameBin(k) = mean(nodeStatus_IG_cluster_dcenter_Middle_forNode_sameBin{k});
   nodeStatus_IG_clusterAVG_dcenter_Back_sameBin(k) = mean(nodeStatus_IG_cluster_dcenter_Back_forNode_sameBin{k});
  
   nodeStatus_SG_clusterSize_d_sameBin(k) = length(nodeStatus_SG_cluster_d_forNode_sameBin{k});
   nodeStatus_SG_clusterSize_D3_sameBin(k) = length(nodeStatus_SG_cluster_D3_forNode_sameBin{k});
   nodeStatus_SG_clusterSize_dcenter_sameBin(k) = length(nodeStatus_SG_cluster_dcenter_forNode_sameBin{k});
   nodeStatus_SG_clusterSize_dcenter_Front_sameBin(k) = length(nodeStatus_SG_cluster_dcenter_Front_forNode_sameBin{k});
   nodeStatus_SG_clusterSize_dcenter_Middle_sameBin(k) = length(nodeStatus_SG_cluster_dcenter_Middle_forNode_sameBin{k});
   nodeStatus_SG_clusterSize_dcenter_Back_sameBin(k) = length(nodeStatus_SG_cluster_dcenter_Back_forNode_sameBin{k});
   
 for j=1:length(index_d2{k})
     
     nodeStatus_SG_cluster_dcenter_forNode2{k} = [nodeStatus_SG_cluster_dcenter_forNode2{k}, nodeStatus_all(index_dcenter2{k}(j), node)];
     nodeStatus_IG_cluster_dcenter_forNode2{k} = [nodeStatus_IG_cluster_dcenter_forNode2{k}, nodeStatus_IG_all(index_dcenter2{k}(j)+start_frame-1, node)];
     
     if dir_node2(index_dcenter2{k}(j), node)>=-pi/3 & dir_node2(index_dcenter2{k}(j), node)<=pi/3
         nodeStatus_SG_cluster_dcenter_Front_forNode2{k} = [nodeStatus_SG_cluster_dcenter_Front_forNode2{k}, nodeStatus_all(index_dcenter2{k}(j), node)];
         nodeStatus_IG_cluster_dcenter_Front_forNode2{k} = [nodeStatus_IG_cluster_dcenter_Front_forNode2{k}, nodeStatus_IG_all(index_dcenter2{k}(j)+start_frame-1, node)];
     elseif (dir_node2(index_dcenter2{k}(j), node)>=pi/3 & dir_node2(index_dcenter2{k}(j), node)<=2*pi/3) | (dir_node2(index_dcenter2{k}(j),node)<=-pi/3 & dir_node2(index_dcenter2{k}(j), node)>=-2*pi/3)
         nodeStatus_SG_cluster_dcenter_Middle_forNode2{k} = [nodeStatus_SG_cluster_dcenter_Middle_forNode2{k}, nodeStatus_all(index_dcenter2{k}(j), node)];
         nodeStatus_IG_cluster_dcenter_Middle_forNode2{k} = [nodeStatus_IG_cluster_dcenter_Middle_forNode2{k}, nodeStatus_IG_all(index_dcenter2{k}(j)+start_frame-1, node)];
     elseif dir_node2(index_dcenter2{k}(j), node)>=2*pi/3 | dir_node2(index_dcenter2{k}(j), node)<=-2*pi/3
         nodeStatus_SG_cluster_dcenter_Back_forNode2{k} = [nodeStatus_SG_cluster_dcenter_Back_forNode2{k}, nodeStatus_all(index_dcenter2{k}(j), node)];
         nodeStatus_IG_cluster_dcenter_Back_forNode2{k} = [nodeStatus_IG_cluster_dcenter_Back_forNode2{k}, nodeStatus_IG_all(index_dcenter2{k}(j)+start_frame-1, node)];
     end
    
 end
   
 if isempty(index_d2{k})
     nodeStatus_SG_cluster_dcenter_forNode2{k} = [0];
     nodeStatus_IG_cluster_dcenter_forNode2{k} = [0];
 end
 
   nodeStatus_SG_clusterAVG_d_forNode2(k) = mean(nodeStatus_SG_cluster_d_forNode2{k});
   nodeStatus_SG_clusterAVG_D3_forNode2(k) = mean(nodeStatus_SG_cluster_D3_forNode2{k});
   nodeStatus_SG_clusterAVG_dcenter_forNode2(k) = mean(nodeStatus_SG_cluster_dcenter_forNode2{k});
   nodeStatus_SG_clusterAVG_dcenter_Front_forNode2(k) = mean(nodeStatus_SG_cluster_dcenter_Front_forNode2{k});
   nodeStatus_SG_clusterAVG_dcenter_Middle_forNode2(k) = mean(nodeStatus_SG_cluster_dcenter_Middle_forNode2{k});
   nodeStatus_SG_clusterAVG_dcenter_Back_forNode2(k) = mean(nodeStatus_SG_cluster_dcenter_Back_forNode2{k});
   
   nodeStatus_IG_clusterAVG_d_forNode2(k) = mean(nodeStatus_IG_cluster_d_forNode2{k});
   nodeStatus_IG_clusterAVG_D3_forNode2(k) = mean(nodeStatus_IG_cluster_D3_forNode2{k});
   nodeStatus_IG_clusterAVG_dcenter_forNode2(k) = mean(nodeStatus_IG_cluster_dcenter_forNode2{k});
   nodeStatus_IG_clusterAVG_dcenter_Front_forNode2(k) = mean(nodeStatus_IG_cluster_dcenter_Front_forNode2{k});
   nodeStatus_IG_clusterAVG_dcenter_Middle_forNode2(k) = mean(nodeStatus_IG_cluster_dcenter_Middle_forNode2{k});
   nodeStatus_IG_clusterAVG_dcenter_Back_forNode2(k) = mean(nodeStatus_IG_cluster_dcenter_Back_forNode2{k});
   
   nodeStatus_SG_clusterSize_d_forNode2(k) = length(nodeStatus_SG_cluster_d_forNode2{k});
   nodeStatus_SG_clusterSize_D3_forNode2(k) = length(nodeStatus_SG_cluster_D3_forNode2{k});
   nodeStatus_SG_clusterSize_dcenter_forNode2(k) = length(nodeStatus_SG_cluster_dcenter_forNode2{k});
   nodeStatus_SG_clusterSize_dcenter_Front_forNode2(k) = length(nodeStatus_SG_cluster_dcenter_Front_forNode2{k});
   nodeStatus_SG_clusterSize_dcenter_Middle_forNode2(k) = length(nodeStatus_SG_cluster_dcenter_Middle_forNode2{k});
   nodeStatus_SG_clusterSize_dcenter_Back_forNode2(k) = length(nodeStatus_SG_cluster_dcenter_Back_forNode2{k});
   
 for j=1:length(rows_dcenter2_sameBin{k})
      
     nodeStatus_SG_cluster_dcenter_forNode2_sameBin{k} = [nodeStatus_SG_cluster_dcenter_forNode2_sameBin{k}, nodeStatus_all(rows_dcenter2_sameBin{k}(j), node)];
     nodeStatus_IG_cluster_dcenter_forNode2_sameBin{k} = [nodeStatus_IG_cluster_dcenter_forNode2_sameBin{k}, nodeStatus_IG_all(rows_dcenter2_sameBin{k}(j)+start_frame-1, node)];    
     
     if dir_node2(rows_dcenter2_sameBin{k}(j), node)>=-pi/3 & dir_node2(rows_dcenter2_sameBin{k}(j), node)<=pi/3
         nodeStatus_SG_cluster_dcenter_Front_forNode2_sameBin{k} = [nodeStatus_SG_cluster_dcenter_Front_forNode2_sameBin{k}, nodeStatus_all(rows_dcenter2_sameBin{k}(j), node)];
         nodeStatus_IG_cluster_dcenter_Front_forNode2_sameBin{k} = [nodeStatus_IG_cluster_dcenter_Front_forNode2_sameBin{k}, nodeStatus_IG_all(rows_dcenter2_sameBin{k}(j)+start_frame-1, node)];
            
    elseif (dir_node2(rows_dcenter2_sameBin{k}(j), node)>=pi/3 & dir_node2(rows_dcenter2_sameBin{k}(j), node)<=2*pi/3) | (dir_node2(rows_dcenter2_sameBin{k}(j), node)<=-pi/3 & dir_node2(rows_dcenter2_sameBin{k}(j), node)>=-2*pi/3)
         nodeStatus_SG_cluster_dcenter_Middle_forNode2_sameBin{k} = [nodeStatus_SG_cluster_dcenter_Middle_forNode2_sameBin{k}, nodeStatus_all(rows_dcenter2_sameBin{k}(j), node)];
         nodeStatus_IG_cluster_dcenter_Middle_forNode2_sameBin{k} = [nodeStatus_IG_cluster_dcenter_Middle_forNode2_sameBin{k}, nodeStatus_IG_all(rows_dcenter2_sameBin{k}(j)+start_frame-1, node)];
         
     elseif dir_node2(rows_dcenter2_sameBin{k}(j), node)>=2*pi/3 | dir_node2(rows_dcenter2_sameBin{k}(j), node)<=-2*pi/3
         nodeStatus_SG_cluster_dcenter_Back_forNode2_sameBin{k} = [nodeStatus_SG_cluster_dcenter_Back_forNode2_sameBin{k}, nodeStatus_all(rows_dcenter2_sameBin{k}(j), node)];
         nodeStatus_IG_cluster_dcenter_Back_forNode2_sameBin{k} = [nodeStatus_IG_cluster_dcenter_Back_forNode2_sameBin{k}, nodeStatus_IG_all(rows_dcenter2_sameBin{k}(j)+start_frame-1, node)];       
     end
      
 end
 
 if isempty(rows_dcenter2_sameBin{k})
     nodeStatus_SG_cluster_dcenter_forNode2_sameBin{k} = [0];
     nodeStatus_IG_cluster_dcenter_forNode2_sameBin{k} = [0];
 end
 
   nodeStatus_SG_clusterAVG_d_forNode2_sameBin(k) = mean(nodeStatus_SG_cluster_d_forNode2_sameBin{k});
   nodeStatus_SG_clusterAVG_D3_forNode2_sameBin(k) = mean(nodeStatus_SG_cluster_D3_forNode2_sameBin{k});
   nodeStatus_SG_clusterAVG_dcenter_forNode2_sameBin(k) = mean(nodeStatus_SG_cluster_dcenter_forNode2_sameBin{k});
   nodeStatus_SG_clusterAVG_dcenter_Front_forNode2_sameBin(k) = mean(nodeStatus_SG_cluster_dcenter_Front_forNode2_sameBin{k});
   nodeStatus_SG_clusterAVG_dcenter_Middle_forNode2_sameBin(k) = mean(nodeStatus_SG_cluster_dcenter_Middle_forNode2_sameBin{k});
   nodeStatus_SG_clusterAVG_dcenter_Back_forNode2_sameBin(k) = mean(nodeStatus_SG_cluster_dcenter_Back_forNode2_sameBin{k});
  
   nodeStatus_IG_clusterAVG_d_forNode2_sameBin(k) = mean(nodeStatus_IG_cluster_d_forNode2_sameBin{k});
   nodeStatus_IG_clusterAVG_D3_forNode2_sameBin(k) = mean(nodeStatus_IG_cluster_D3_forNode2_sameBin{k});
   nodeStatus_IG_clusterAVG_dcenter_forNode2_sameBin(k) = mean(nodeStatus_IG_cluster_dcenter_forNode2_sameBin{k});
   nodeStatus_IG_clusterAVG_dcenter_Front_forNode2_sameBin(k) = mean(nodeStatus_IG_cluster_dcenter_Front_forNode2_sameBin{k});
   nodeStatus_IG_clusterAVG_dcenter_Middle_forNode2_sameBin(k) = mean(nodeStatus_IG_cluster_dcenter_Middle_forNode2_sameBin{k});
   nodeStatus_IG_clusterAVG_dcenter_Back_forNode2_sameBin(k) = mean(nodeStatus_IG_cluster_dcenter_Back_forNode2_sameBin{k});
  
   nodeStatus_SG_clusterSize_d_forNode2_sameBin(k) = length(nodeStatus_SG_cluster_d_forNode2_sameBin{k});
   nodeStatus_SG_clusterSize_D3_forNode2_sameBin(k) = length(nodeStatus_SG_cluster_D3_forNode2_sameBin{k});
   nodeStatus_SG_clusterSize_dcenter_forNode2_sameBin(k) = length(nodeStatus_SG_cluster_dcenter_forNode2_sameBin{k});
   nodeStatus_SG_clusterSize_dcenter_Front_forNode2_sameBin(k) = length(nodeStatus_SG_cluster_dcenter_Front_forNode2_sameBin{k});
   nodeStatus_SG_clusterSize_dcenter_Middle_forNode2_sameBin(k) = length(nodeStatus_SG_cluster_dcenter_Middle_forNode2_sameBin{k});
   nodeStatus_SG_clusterSize_dcenter_Back_forNode2_sameBin(k) = length(nodeStatus_SG_cluster_dcenter_Back_forNode2_sameBin{k});
   
 end
 
 for k=1:10
     for j=1:10
       ss_D3_SG_forNode(k,j) = kstest2(nodeStatus_SG_cluster_D3_forNode{k}, nodeStatus_SG_cluster_D3_forNode{j});  
       ss_d_SG_forNode(k,j) = kstest2(nodeStatus_SG_cluster_d_forNode{k}, nodeStatus_SG_cluster_d_forNode{j});
       ss_dCenter_SG_forNode(k,j) = kstest2(nodeStatus_SG_cluster_dcenter_forNode{k}, nodeStatus_SG_cluster_dcenter_forNode{j});
     if ~isempty(nodeStatus_SG_cluster_dcenter_Front_forNode{k})  & ~isempty(nodeStatus_SG_cluster_dcenter_Front_forNode{j})
       ss_dCenter_front_SG_forNode(k,j) = kstest2(nodeStatus_SG_cluster_dcenter_Front_forNode{k}, nodeStatus_SG_cluster_dcenter_Front_forNode{j});
     end
     if ~isempty(nodeStatus_SG_cluster_dcenter_Middle_forNode{k})  & ~isempty(nodeStatus_SG_cluster_dcenter_Middle_forNode{j})
       ss_dCenter_middle_SG_forNode(k,j) = kstest2(nodeStatus_SG_cluster_dcenter_Middle_forNode{k}, nodeStatus_SG_cluster_dcenter_Middle_forNode{j});
     end
     if ~isempty(nodeStatus_SG_cluster_dcenter_Back_forNode{k}) & ~isempty(nodeStatus_SG_cluster_dcenter_Back_forNode{j})
       ss_dCenter_back_SG_forNode(k,j) = kstest2(nodeStatus_SG_cluster_dcenter_Back_forNode{k}, nodeStatus_SG_cluster_dcenter_Back_forNode{j});
     end  
       
       ss_D3_IG_forNode(k,j) = kstest2(nodeStatus_IG_cluster_D3_forNode{k}, nodeStatus_IG_cluster_D3_forNode{j});  
       ss_d_IG_forNode(k,j) = kstest2(nodeStatus_IG_cluster_d_forNode{k}, nodeStatus_IG_cluster_d_forNode{j});
       ss_dCenter_IG_forNode(k,j) = kstest2(nodeStatus_IG_cluster_dcenter_forNode{k}, nodeStatus_IG_cluster_dcenter_forNode{j});
      if ~isempty(nodeStatus_SG_cluster_dcenter_Front_forNode{k})  & ~isempty(nodeStatus_SG_cluster_dcenter_Front_forNode{j})
       ss_dCenter_front_IG_forNode(k,j) = kstest2(nodeStatus_IG_cluster_dcenter_Front_forNode{k}, nodeStatus_IG_cluster_dcenter_Front_forNode{j});
      end
      if ~isempty(nodeStatus_IG_cluster_dcenter_Middle_forNode{k})  & ~isempty(nodeStatus_IG_cluster_dcenter_Middle_forNode{j})  
       ss_dCenter_middle_IG_forNode(k,j) = kstest2(nodeStatus_IG_cluster_dcenter_Middle_forNode{k}, nodeStatus_IG_cluster_dcenter_Middle_forNode{j});
      end
      if ~isempty(nodeStatus_IG_cluster_dcenter_Back_forNode{k}) & ~isempty(nodeStatus_IG_cluster_dcenter_Back_forNode{j})
       ss_dCenter_back_IG_forNode(k,j) = kstest2(nodeStatus_IG_cluster_dcenter_Back_forNode{k}, nodeStatus_IG_cluster_dcenter_Back_forNode{j});
      end
      
     end
 end
 
 for k=1:10
     for j=1:10
       ss_D3_SG_forNode2(k,j) = kstest2(nodeStatus_SG_cluster_D3_forNode2{k}, nodeStatus_SG_cluster_D3_forNode2{j});  
       ss_d_SG_forNode2(k,j) = kstest2(nodeStatus_SG_cluster_d_forNode2{k}, nodeStatus_SG_cluster_d_forNode2{j});
       ss_dCenter_SG_forNode2(k,j) = kstest2(nodeStatus_SG_cluster_dcenter_forNode2{k}, nodeStatus_SG_cluster_dcenter_forNode2{j});
     if ~isempty(nodeStatus_SG_cluster_dcenter_Front_forNode2{k})  & ~isempty(nodeStatus_SG_cluster_dcenter_Front_forNode2{j})
       ss_dCenter_front_SG_forNode2(k,j) = kstest2(nodeStatus_SG_cluster_dcenter_Front_forNode2{k}, nodeStatus_SG_cluster_dcenter_Front_forNode2{j});
     end
     if ~isempty(nodeStatus_SG_cluster_dcenter_Middle_forNode2{k})  & ~isempty(nodeStatus_SG_cluster_dcenter_Middle_forNode2{j})
       ss_dCenter_middle_SG_forNode2(k,j) = kstest2(nodeStatus_SG_cluster_dcenter_Middle_forNode2{k}, nodeStatus_SG_cluster_dcenter_Middle_forNode2{j});
     end
     if ~isempty(nodeStatus_SG_cluster_dcenter_Back_forNode2{k}) & ~isempty(nodeStatus_SG_cluster_dcenter_Back_forNode2{j})
       ss_dCenter_back_SG_forNode2(k,j) = kstest2(nodeStatus_SG_cluster_dcenter_Back_forNode2{k}, nodeStatus_SG_cluster_dcenter_Back_forNode2{j});
     end  
       
       ss_D3_IG_forNode2(k,j) = kstest2(nodeStatus_IG_cluster_D3_forNode2{k}, nodeStatus_IG_cluster_D3_forNode2{j});  
       ss_d_IG_forNode2(k,j) = kstest2(nodeStatus_IG_cluster_d_forNode2{k}, nodeStatus_IG_cluster_d_forNode2{j});
       ss_dCenter_IG_forNode2(k,j) = kstest2(nodeStatus_IG_cluster_dcenter_forNode2{k}, nodeStatus_IG_cluster_dcenter_forNode2{j});
      if ~isempty(nodeStatus_SG_cluster_dcenter_Front_forNode2{k})  & ~isempty(nodeStatus_SG_cluster_dcenter_Front_forNode2{j})
       ss_dCenter_front_IG_forNode2(k,j) = kstest2(nodeStatus_IG_cluster_dcenter_Front_forNode2{k}, nodeStatus_IG_cluster_dcenter_Front_forNode2{j});
      end
      if ~isempty(nodeStatus_IG_cluster_dcenter_Middle_forNode2{k})  & ~isempty(nodeStatus_IG_cluster_dcenter_Middle_forNode2{j})  
       ss_dCenter_middle_IG_forNode2(k,j) = kstest2(nodeStatus_IG_cluster_dcenter_Middle_forNode2{k}, nodeStatus_IG_cluster_dcenter_Middle_forNode2{j});
      end
      if ~isempty(nodeStatus_IG_cluster_dcenter_Back_forNode2{k}) & ~isempty(nodeStatus_IG_cluster_dcenter_Back_forNode2{j})
       ss_dCenter_back_IG_forNode2(k,j) = kstest2(nodeStatus_IG_cluster_dcenter_Back_forNode2{k}, nodeStatus_IG_cluster_dcenter_Back_forNode2{j});
      end
      
     end
 end
 
 % plot bins for d,D,|D|. There are two options here: using the bins with
 % equal number of trials (considered as in PART I), denoted by _sameBin
 % OR using the bins from PART I
 
 figure(21)
 %bar(nodeStatus_SG_clusterAVG_d_forNode2_sameBin)
 bar(nodeStatus_SG_clusterAVG_d_forNode2)
 title('Average Node status SG for 10 bins of distance = absolute value of d to center-axis')
 xlabel('bins')
 ylabel('average node status SG')
 
 figure(22)
 %bar(nodeStatus_SG_clusterAVG_D3_forNode2_sameBin)
 bar(nodeStatus_SG_clusterAVG_D3_forNode2)
 title('Average node status SG for 10 bins of distance = d to center-axis')
 xlabel('bins')
 ylabel('average node status SG')
 
 figure(23)
 %bar(nodeStatus_SG_clusterAVG_dcenter_forNode2_sameBin)
 bar(nodeStatus_SG_clusterAVG_dcenter_forNode2)
 title('Average node status SG for 10 bins of distance = d to center (point)')
 xlabel('bins')
 ylabel('average node status SG')
 
 color_spec = [0 1 0; 0 0 1; 1 0 0];
 figure(24)
 disp('figure 24')
 %bar(1:10, [nodeStatus_SG_clusterAVG_dcenter_Front_forNode2_sameBin', nodeStatus_SG_clusterAVG_dcenter_Middle_forNode2_sameBin', nodeStatus_SG_clusterAVG_dcenter_Back_forNode2_sameBin'], 0.5, 'stack') 
 bar(1:10, [nodeStatus_SG_clusterAVG_dcenter_Front_forNode2', nodeStatus_SG_clusterAVG_dcenter_Middle_forNode2', nodeStatus_SG_clusterAVG_dcenter_Back_forNode2'], 0.5, 'stack') 
 legend('Front', 'Middle', 'Back')
 colormap(color_spec)
 title('Average node status SG binning accoring to distance to center(point) - front, middle, back')
 
 figure(25)
 %bar(nodeStatus_IG_clusterAVG_d_forNode2_sameBin)
 bar(nodeStatus_IG_clusterAVG_d_forNode2)
 title('Average Node status IG for 10 bins of distance = absolute value of d to center-axis')
 xlabel('bins')
 ylabel('average node status IG')
 
 figure(26)
 %bar(nodeStatus_IG_clusterAVG_D3_forNode2_sameBin)
 bar(nodeStatus_IG_clusterAVG_D3_forNode2)
 title('Average node status IG for 10 bins of distance = d to center-axis')
 xlabel('bins')
 ylabel('avg node status IG')
 
 figure(27)
 %bar(nodeStatus_IG_clusterAVG_dcenter_forNode2_sameBin)
 bar(nodeStatus_IG_clusterAVG_dcenter_forNode2)
 title('Average node status IG for 10 bins of distance = d to center (point)')
 xlabel('bins')
 ylabel('average node status IG')
 
 color_spec = [0 1 0; 0 0 1; 1 0 0];
 figure(28)
 %bar(1:10, [nodeStatus_IG_clusterAVG_dcenter_Front_forNode2_sameBin',nodeStatus_IG_clusterAVG_dcenter_Middle_forNode2_sameBin', nodeStatus_IG_clusterAVG_dcenter_Back_forNode2_sameBin'], 0.5, 'stack') 
 bar(1:10, [nodeStatus_IG_clusterAVG_dcenter_Front_forNode2',nodeStatus_IG_clusterAVG_dcenter_Middle_forNode2', nodeStatus_IG_clusterAVG_dcenter_Back_forNode2'], 0.5, 'stack') 
 legend('Front', 'Middle', 'Back')
 colormap(color_spec)
 
 figure(29)
 %bar(nodeStatus_SG_clusterSize_d_forNode2_sameBin)
 bar(nodeStatus_SG_clusterSize_d_forNode2)
 title('How many trials per 10 bins of distance = absolute value of d to center-axis')
 xlabel('bins')
 ylabel('number of trials')
 
 figure(30)
 %bar(nodeStatus_SG_clusterSize_D3_forNode2_sameBin)
 bar(nodeStatus_SG_clusterSize_D3_forNode2)
 title('How many trials per 10 bins of distance = d to center-axis')
 xlabel('bins')
 ylabel('number of trials')
 
 figure(31)
 %bar(nodeStatus_SG_clusterSize_dcenter_forNode2_sameBin)
 bar(nodeStatus_SG_clusterSize_dcenter_forNode2)
 title('How many trials per 10 bins of distance = d to center (point)')
 xlabel('bins')
 ylabel('number of trials')
 
 color_spec = [0 1 0; 0 0 1; 1 0 0];
 figure(32)
 bar(1:10, [nodeStatus_SG_clusterSize_dcenter_Front_forNode2',nodeStatus_SG_clusterSize_dcenter_Middle_forNode2', nodeStatus_SG_clusterSize_dcenter_Back_forNode2'], 0.5, 'stack') 
 %bar(1:10, [nodeStatus_SG_clusterSize_dcenter_Front_forNode2_sameBin',nodeStatus_SG_clusterSize_dcenter_Middle_forNode2_sameBin', nodeStatus_SG_clusterSize_dcenter_Back_forNode2_sameBin'], 0.5, 'stack') 
 legend('Front', 'Middle', 'Back')
 colormap(color_spec)
 title('How many trials per bin for d to center(point) - front, middle, back')
 
  disp('statistical results on the bins')
  
  ss_D3_SG_forNode  
  ss_d_SG_forNode
  ss_dCenter_SG_forNode
  ss_dCenter_front_SG_forNode
  ss_dCenter_middle_SG_forNode
  ss_dCenter_back_SG_forNode
       
  ss_D3_IG_forNode 
  ss_d_IG_forNode
  ss_dCenter_IG_forNode
  ss_dCenter_front_IG_forNode
  ss_dCenter_middle_IG_forNode
  ss_dCenter_back_IG_forNode
  
  ss_D3_SG_forNode2  
  ss_d_SG_forNode2
  ss_dCenter_SG_forNode2
  ss_dCenter_front_SG_forNode2
  ss_dCenter_middle_SG_forNode2
  ss_dCenter_back_SG_forNode2
       
  ss_D3_IG_forNode2
  ss_d_IG_forNode2
  ss_dCenter_IG_forNode2
  ss_dCenter_front_IG_forNode2
  ss_dCenter_middle_IG_forNode2
  ss_dCenter_back_IG_forNode2
  %}
  %{
  
  disp('Does high number of trials mean high node status?')
  
  disp('ns SG and d')
  [c p] = corrcoef(nodeStatus_SG_clusterAVG_d_forNode2_sameBin, nodeStatus_SG_clusterSize_d_forNode2_sameBin)
  
  disp('ns SG and D3')
  [c p] = corrcoef(nodeStatus_SG_clusterAVG_D3_forNode2_sameBin, nodeStatus_SG_clusterSize_D3_forNode2_sameBin)
  
  disp('ns SG and distToCenter')
  [c p] = corrcoef(nodeStatus_SG_clusterAVG_dcenter_forNode2_sameBin, nodeStatus_SG_clusterSize_dcenter_forNode2_sameBin)
  
  disp('ns SG and distToCenter Front')
  [c p] = corrcoef(nodeStatus_SG_clusterAVG_dcenter_Front_forNode2_sameBin, nodeStatus_SG_clusterSize_dcenter_Front_forNode2_sameBin)
  
  disp('ns SG and distToCenter Middle')
  [c p] = corrcoef(nodeStatus_SG_clusterAVG_dcenter_Middle_forNode2_sameBin, nodeStatus_SG_clusterSize_dcenter_Middle_forNode2_sameBin)
  
  disp('ns SG and distToCenter Back')
  [c p] = corrcoef(nodeStatus_SG_clusterAVG_dcenter_Back_forNode2_sameBin, nodeStatus_SG_clusterSize_dcenter_Back_forNode2_sameBin)
  
  disp('ns IG and d')
  [c p] = corrcoef(nodeStatus_IG_clusterAVG_d_forNode2_sameBin, nodeStatus_SG_clusterSize_d_forNode2_sameBin)
  
  disp('ns IG and D3')
  [c p] = corrcoef(nodeStatus_IG_clusterAVG_D3_forNode2_sameBin, nodeStatus_SG_clusterSize_D3_forNode2_sameBin)
  
  disp('ns IG and dTocenter')
  [c p] = corrcoef(nodeStatus_IG_clusterAVG_dcenter_forNode2_sameBin, nodeStatus_SG_clusterSize_dcenter_forNode2_sameBin)

  disp('ns IG and dToCenter Front')
  [c p] = corrcoef(nodeStatus_IG_clusterAVG_dcenter_Front_forNode2_sameBin, nodeStatus_SG_clusterSize_dcenter_Front_forNode2_sameBin)

  disp('ns IG and dToCenter Middle')
  [c p] = corrcoef(nodeStatus_IG_clusterAVG_dcenter_Middle_forNode2_sameBin, nodeStatus_SG_clusterSize_dcenter_Middle_forNode2_sameBin)
  
  disp('ns IG and dToCenter Back')
  [c p] = corrcoef(nodeStatus_IG_clusterAVG_dcenter_Back_forNode2_sameBin, nodeStatus_SG_clusterSize_dcenter_Back_forNode2_sameBin)
  %}
  else
      
      % for every node average of node status SG and IG, and average of
      % d,D,|D| for all frames
      
      numberOfNodes = handles.numberOfNodes;
     
      %compute node status SG,IG and smooth out these functions
      [nodeStatus_IG, nodeStatus_avg2_IG] = nodeStatusIG(handles, hObject, start_frame, end_frame);
      nodeStatus_IG_normalized = normalize(nodeStatus_IG);  
      nodeStatus_avg2_IG_normalized = normalize(nodeStatus_avg2_IG);
      
        %[nodeStatus_sorted, Index_nodeStatus_sorted] = sort(nodeStatus_IG(handles.lb,:));
       %}%{
       
      [nodeStatus, nodeStatus_avg] = compute_nodeStatus(hObject, handles, adj, start_frame, end_frame);
      nodeStatus_normalized = normalize(nodeStatus);
      nodeStatus_avg_normalized = normalize(nodeStatus_avg);
      
      %avg_toCenter is the average d for each node
      [dist_to_center, avg_toCenter] = distToCenter(handles, hObject, start_frame, end_frame, numberOfNodes,1);
      avg_toCenter_normalized = normalize(avg_toCenter);
      [Distance2, Distance3, distance, leader, backer, N_leader, N_backer] = LeaderBacker(hObject, handles, start_frame, end_frame,1);
      
      %average of D and |D| for each node o
      for o=1:numberOfNodes
        Distance3_avg(o) = mean(Distance3(:,o));
        distance_avg(o) = mean(distance(:,o));
      end
      
      %normalize D, |D|
      Distance3_avg = Distance3_avg/max(Distance3_avg);
      distance_avg = distance_avg/max(distance_avg);
      nodeStatus_avg = nodeStatus_avg/max(nodeStatus_avg);
      
      figure(1)
      plot(1:numberOfNodes, nodeStatus_avg_normalized, 'x', 'Color', 'b', 'MarkerSize', 14);
      hold on
      plot(1:numberOfNodes, avg_toCenter_normalized, 'x', 'Color', 'r', 'MarkerSize', 14);
      title('For each node, average node status SG vs. average distance to center')
      xlabel('node')
      ylabel('Average node status SG and average distance to center')
      legend('node status SG', 'average dist from center')
      
      disp('average node status SG and average distance to center correlation, d')
      %[xcf, lags, bounds] = crosscorr(nodeStatus_avg, avg_toCenter);
      %xcf'
      %bounds'
      [c p] = corrcoef(nodeStatus_avg, avg_toCenter);
      c
      p
      
      figure(2)
      plot(1:numberOfNodes, nodeStatus_avg_normalized, 'x', 'Color', 'b', 'MarkerSize', 14);
      hold on
      plot(1:numberOfNodes, Distance3_avg, 'x', 'Color', 'r', 'MarkerSize', 14)
      title('For each node, average node status SG and average distance from the center-axis')
      xlabel('node')
      ylabel('Average node status SG and average distance from the center-axis')
      legend('average node status SG', 'average distance from center-axis')
      
      disp('average node status SG and average distance to center-axis correlation, D')
     % [xcf, lags, bounds] = crosscorr(nodeStatus_avg, Distance3_avg);
     % xcf'
     % bounds'
      [c p] = corrcoef(nodeStatus_avg, Distance3_avg);
      c
      p
      
      figure(3)
      plot(1:numberOfNodes, nodeStatus_avg_normalized, 'x', 'Color', 'b', 'MarkerSize', 14)
      hold on
      plot(1:numberOfNodes, distance_avg, 'x', 'Color', 'r', 'MarkerSize', 14)
      hold on
      plot(1:numberOfNodes, nodeStatus_avg2_IG_normalized, 'x', 'Color', 'g', 'MarkerSize', 14)
      title('For each node, average node status SG, node status IG, average distance from center')
      xlabel('node')
      legend('average node status SG', 'average distance from center', 'average node status IG')
      
      disp('average node status SG and |D|')
     % [xcf, lags, bounds] = crosscorr(nodeStatus_avg, nodeStatus_avg2_IG);
     % xcf'
     % bounds'
      [c p] = corrcoef(nodeStatus_avg, distance_avg);
      c
      p
      
      figure(4)
       plot(1:numberOfNodes, nodeStatus_avg_normalized, 'x', 'Color', 'b', 'MarkerSize', 14)
       hold on
       plot(1:numberOfNodes, nodeStatus_avg2_IG_normalized, 'x', 'Color', 'r', 'MarkerSize', 14)
       title('For each node, average node status SG and average node status IG') 
       xlabel('node') 
       legend('average node status SG', 'average node status IG')
      
       disp('correlation between average node status of sensing graph and average node status of influence graph')
       [c p] = corrcoef(nodeStatus_avg_normalized, nodeStatus_avg2_IG);
      c
      p
      
        figure(5)
        plot(1:numberOfNodes, avg_toCenter, 'x', 'Color', 'r', 'MarkerSize', 14)
       hold on
        plot(1:numberOfNodes, nodeStatus_avg2_IG_normalized, 'x', 'Color', 'b', 'MarkerSize', 14)
        title('For each node, average distance to center and node status IG')
        xlabel('node')
        legend('average distance to center', 'node status IG')
        
        disp('average distance to center, d,  and node status of interaction graph correlated')
     % [xcf, lags, bounds] = crosscorr(avg_toCenter, nodeStatus_avg2_IG);
     % xcf'
     % bounds'
      [c p] = corrcoef(avg_toCenter, nodeStatus_avg2_IG);
      c
      p
      
        figure(6)
        plot(1:numberOfNodes, Distance3_avg, 'x', 'Color', 'r', 'MarkerSize', 14)
        hold on
        plot(1:numberOfNodes, nodeStatus_avg2_IG_normalized, 'x', 'Color', 'b', 'MarkerSize', 14)
        title('For each node, average distance from center-axis and average node status IG')
        xlabel('node')
        legend('average distance to center', 'node status IG')
      
      disp('average distance from center-axis, D, and node status of interaction graph correlated')
      %[xcf, lags, bounds] = crosscorr(Distance3_avg, nodeStatus_avg2_IG);
      %xcf'
      %bounds'
      [c p] = corrcoef(Distance3_avg, nodeStatus_avg2_IG);
      c
      p
        
        figure(7)
        plot(1:numberOfNodes, distance_avg, 'x', 'Color', 'r', 'MarkerSize', 14)
       hold on
        plot(1:numberOfNodes, nodeStatus_avg2_IG_normalized, 'x', 'Color', 'b', 'MarkerSize', 14)
        title('For each node, average distance from center and average node status IG')
        legend('average distance to center', 'node status IG')
        
        disp('average absolute value of distance to center-axis, |D|, and average node status of the influence graph, ns IG')
        [c p] = corrcoef(distance_avg, nodeStatus_avg2_IG);
        c
        p
        
end

end

%normalizes any vector by dividing by the max element of it
function x = normalize(x)

[m n] = size(x);

for i=1:m
        m = max(x(i,:));
        x(i,:) = x(i,:)/m;
end

end

function group_properties(hObject, handles)
    
    lb = handles.lb;
    ub = handles.ub;

    prompt = {'Select start frame:','Select end frame:'};
    answer = inputdlg(prompt, 'Select range of frames');
    start_frame = str2num(answer{1});
    end_frame = str2num(answer{2});
  
    if start_frame<25
        error('Cannot have a starting frame <25, because nodeStatus_IG cannot be computed otherwise');
        return
    end
    
    prompt = {'What should the weight matrix include?                                                                                           Select 1). Distance2: the distance from the center-axis;                                                                                                                                                    Select 2). distance : The abs of the distance to the center-axis                                                                                                                                                    Select 3). distance_toCenter: the distance to the center of the group                                                                                                                                                       Select 4). All 3 position measures'};
    option = str2num(cell2mat(inputdlg(prompt)));
    
    numberOfNodes = handles.numberOfNodes;
    
    %[positionMatrix, distanceMatrix, avg_distance, avg_distance_inTime, avg_distance_inTime_closest1, avg_distance_inTime_closest2, avg_distanceAll, avg_dist_index, closest_neighbors, relativeAngles, HeadingDiff, angleToclosest, Avg_angleToclosest, Avg_angleOverall] = computeAll_pos_dist(handles, start_frame, end_frame, option);
    [dist_toCenter, avg_toCenter] = distToCenter(handles, hObject, start_frame, end_frame, 13,1);
    max_dist_toCenter = max(max(dist_toCenter));
    dist_toCenter = dist_toCenter/max_dist_toCenter;
    
    [Distance2, Distance3, distance, leader, backer, N_leader, N_backer] = LeaderBacker(hObject, handles, start_frame, end_frame,1);
    max_Distance2 = max(max(Distance2));
    max_distance = max(max(distance));
    Distance2 = Distance2/max_Distance2;
    distance = distance/max_distance;
    
    [nodeStatus_IG, nodeStatus_avg2_IG] = nodeStatusIG(handles, hObject, start_frame, end_frame);      
     %betweenness_c = betweenness(adj);
     adj = choose_graph(hObject, handles, 0);
     [nodeStatus, nodeStatus_avg] = compute_nodeStatus(hObject, handles, adj, start_frame, end_frame);
    %[distNeigh, distNeigh_avg] = dist_to_neighb(hObject, handles);
    
    cnt = 1; %entries in known data 
    cnt2 = 1;
   
    for j=1:numberOfNodes
     for i = start_frame:end_frame-1-20    
      %for i = start_frame:start_frame+100-1
       HugeMatrix(cnt,cnt2) = j;
       cnt2 = cnt2+1;
       HugeMatrix(cnt, cnt2) = Distance2(i-start_frame+1, j);
       cnt2 = cnt2+1;
       HugeMatrix(cnt, cnt2) = distance(i-start_frame+1, j);
       cnt2 = cnt2+1;
       HugeMatrix(cnt, cnt2) = dist_toCenter(i-start_frame+1, j);
       cnt2 = cnt2+1;
       HugeMatrix(cnt, cnt2) = nodeStatus_IG(i, j);
       cnt2 = cnt2+1;
       HugeMatrix(cnt, cnt2) = nodeStatus(i-start_frame+1, j);
     cnt = cnt+1;
     cnt2 = 1;
     end   
    end
    
    [M N] = size(HugeMatrix);
    %colors = jet(numberOfNodes);
    colors = jet(13);
    eps = 1000;
    
    for i=1:M
       for j=1:M 
           switch option
               case 1
              W(i,j) = exp(-norm([HugeMatrix(i,1)-HugeMatrix(j,1), HugeMatrix(i,4)-HugeMatrix(j,4), HugeMatrix(i,5)-HugeMatrix(j,5)])/eps);     
               case 2
              W(i,j) = exp(-norm([HugeMatrix(i,2)-HugeMatrix(j,2), HugeMatrix(i,4)-HugeMatrix(j,4), HugeMatrix(i,5)-HugeMatrix(j,5)])/eps);
               case 3
              W(i,j) = exp(-norm(HugeMatrix(i,3:5)-HugeMatrix(j,3:5))/eps);     
               case 4
              W(i,j) = exp(-norm(HugeMatrix(i,:)-HugeMatrix(j,:))/eps);
           end
       end
    end 
    
    for i=1:M
        D(i,i) = sum(W(i,:));
        A(i,:) = W(i,:)/sum(W(i,:));
    end
    
    close all
    
    figure(1)
    for e = 1:M
        scatter3([HugeMatrix(e,2)], [HugeMatrix(e,5)], [HugeMatrix(e,6)], 50, colors(HugeMatrix(e,1),:), 'fill')
        hold on;   
    end
    
    xlabel('Distance2')
    ylabel('node status IG')
    zlabel('node status')
    
     figure(2)
    for e = 1:M
        scatter3([HugeMatrix(e,3)], [HugeMatrix(e,5)], [HugeMatrix(e,6)], 50, colors(HugeMatrix(e,1),:), 'fill')
        hold on;
    end
    
    xlabel('distance')
    ylabel('node status IG')
    zlabel('node status')
    
     figure(3)
    for e = 1:M
        scatter3([HugeMatrix(e,4)], [HugeMatrix(e,5)], [HugeMatrix(e,6)], 50, colors(HugeMatrix(e,1),:), 'fill')
        hold on;
    
    end
    
    xlabel('distance_toCenter')
    ylabel('node status IG')
    zlabel('node status')
    
    heatm = HeatMap(W);
    
    figure(4)
    imagesc(W)
    
    diffusion_map(W, D, A, start_frame, end_frame)
    
end

%a way of figuring out the structure of the data - look up diffusion map
function diffusion_map(W, D, A, start_frame, end_frame)

colors = jet(13);

n = length(W);

L = D-A;
[V, E] = eigs(A,130);

t=1;
D_map = [];
cnt=1;
for i=1:13
  for j=1:10 
    D_map = [D_map; E(cnt,cnt)^t*V((i-1)*(end_frame-start_frame-1-20)+j,:)];
    cnt = cnt+1;
  end
end

disp('size of D_map')
size(D_map)
D_map(:,2:4)

cnt = 1;
figure(6)
for j=3:130
  plot([D_map(j,2)], [D_map(j,3)], 'x', 'linewidth', 4);
  hold on
end

cnt=1;
figure(7)
for j=3:130 
 scatter3([D_map(j,2)], [D_map(j,3)], [D_map(j,4)], 'x', 'linewidth', 4)   
 hold on   
end

end

function spectral_clustering

     
end

function max_prob_Callback(hObject, eventdata, handles)
end
function max_prob_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function Increment2_Callback(hObject, eventdata, handles)
end
function Increment2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function Increment_Callback(hObject, eventdata, handles)
end
function Increment_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function ProbCheck_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function savewhat_pushbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savewhat_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over plot_position.
function plot_position_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to plot_position (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over numNeigh.
function numNeigh_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to numNeigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end
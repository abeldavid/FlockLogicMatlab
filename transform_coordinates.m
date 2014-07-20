first=1;

answer = inputdlg('What is the last frame?');
last = str2num(cell2mat(answer));

location = uigetdir('','Please browse to Dance folder');
addpath(location);
cd(location);

for i=first:last 
    file =[num2str(i) '.dat']
    data = importdata(file);
    x = data(:,1);
    y = data(:,2);
    theta = data(:,3);
    cd('D:\Work stuff\FlockLogic\version 9')
    [x y theta] = transf_in_familiar_coord_video9(x, y, theta);
    data = [x y theta]
    cd(location)
    
    save(file, 'data', '-ascii')
end

cd('D:\Work stuff\FlockLogic\version 9')
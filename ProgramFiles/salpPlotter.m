%Launches salp plotter

%Find out the exact location this file is being launched from
stk = dbstack;
filePath = which(stk(1).file);
slashIndices = strfind(filePath,'\');
%Get program files directory
prgFiles = filePath(1:slashIndices(end));

%Make sure that the SalpPlotter directory is on the current path
addpath([prgFiles,'SalpPlotter']);

%Delete the information we generated here to cover our tracks
clear stk filePath slashIndices prgFiles;

%Let user know something happened
disp('Starting GUI');

%Launch GUI
salpPlotterGUI;
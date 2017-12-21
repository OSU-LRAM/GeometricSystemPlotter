% GUI data saver test

function GUI_property_scraper(handles)

% list of things we want to grab from the GUI. GUIDE property inspector
% lists all properties that can be set by the user. 
% Tag is read first, and everything else can be written using that tag
% props = {'Tag','BackgroundColor','Enable','FontName','FontSize','FontWeight',...
%          'InnerPosition','OuterPosition','Position','String','Value','Visible'}';
props = {'BackgroundColor','Enable','FontName','FontSize','FontWeight',...
         'Position','String','Value','Visible'}';

% get the names of all the objects in the GUI
tagNames = fieldnames(handles);
propertyList.Tag = tagNames;

% save the props list for later loading
propertyList.propSelections = props; 

% for each property selected and each tag name, gather data
for i = 1:numel(props)
    
    for j = 1:numel(tagNames) % skip the figure at the top
        if isprop(handles.(tagNames{j}),(props{i}))
            % gather the existing property for that object
            propertyList.(props{i}){j,1} = handles.(tagNames{j}).(props{i});
        else
            % if that property doesn't exist, fill cell with NaN to avoid errors
            propertyList.(props{i}){j,1} = NaN;
        end
    end
    
end


% Ask the user where to save the file and what to call it
[file,path] = uiputfile('propertyList.mat','Save Property File As');
path_file = fullfile(path,file);
save(path_file,'propertyList')

% save(fullfile([sysplotter_inputpath '\sysplotter_data'],[snake_sys,'_',lambda_func,'_',length_sys,'_runtimes.mat']),'opt_time');

% for accessing properties later: 
% get(allProps.handleProps{i},'Tag')

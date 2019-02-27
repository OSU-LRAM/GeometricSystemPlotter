function [] = sysf_config_gui_savevals(hObject)
%sysf_config_gui_savevals: This function will be the callback for the save
%button in the config gui. It will save all of the values from the gui back
%to the mat file and will then close the gui
% Inputs: 
% hObject: this is the object handle to the whole gui

handles = guidata(hObject);

%Get the path to the mat save file
matFilePath = handles.varargin{2};

%Get the system structure
s = load(matFilePath,'s');
s = s.s;

%%%%%% 
% Replace fields in struct s with array + numeric values from gui
systemfields = handles.text_inputs.systemfields;
guielements = handles.text_inputs.guielements;


%Loop through all the gui elements and either set it or empty it
for k = 1:size(systemfields,2)
    str = guielements{k}.String;
    if ~strcmpi(str,'None')
        if isempty(str)
            value = 0;
        else
            value = eval(str);
        end
        
        s = setfield(s,systemfields{k}{:},value);
    else
        s = rmfieldnested(s,systemfields{k});
    end
end

%Get the values from the drop down lists
dropdownfields = handles.text_inputs.dropdownfields;
dropdownguielements = handles.text_inputs.dropdownguielements;
dropdownoptions = handles.text_inputs.dropdownoptions;

%Loop through all the dropdown elements and either set it or empty it
for k = 1:size(dropdownfields,2)
    val = dropdownguielements{k}.Value;
    str = dropdownoptions{k}{val};
    if strcmpi(str,'None')
        s = rmfieldnested(s,dropdownfields{k});
    elseif strcmpi(str,'Use Number')
        numberbox = handles.(strjoin(dropdownfields{k},'_'));
        value = eval(numberbox.String);
        s = setfield(s,dropdownfields{k}{:},value);
    else
        value = str;
        s = setfield(s,dropdownfields{k}{:},value);
    end
end

if ~strcmpi(handles.visual_grid,'None') && isfieldnested(s,{'visual','grid'})
    ndgridnum = s.visual.grid;
    s.visual.grid = cell(s.visual.cellsize);
    [s.visual.grid{:}] = ndgrid(ndgridnum);
else
    s = rmfieldnested(s,{'visual','grid'});
end

if ~strcmpi(handles.geometry_link_number.String,'None')
    temparray = num2cell(eval(handles.geometry_link_number.String));
    [s.geometry.constraint_list.link_number] = temparray{:};
else
    s = rmfieldnested(s,{'geometry','constraint_list','link_number'});
end

if handles.geometry_constraint_direction.Data ~= -1
    tempcell = num2cell(handles.geometry_constraint_direction.Data,2);
   
    
    [s.geometry.constraint_list.constraint_direction] = tempcell{:};
    
else
    s = rmfieldnested(s,{'geometry','constraint_list','link_number'});
end


                       

save(matFilePath,'s');
close(hObject);

end


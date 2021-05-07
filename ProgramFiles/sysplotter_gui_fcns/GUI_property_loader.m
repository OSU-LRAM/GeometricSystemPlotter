% visual property editor function

function GUI_property_loader(handles,propertyFile)

load(propertyFile);
props = propertyList.propSelections;
tagNames = propertyList.Tag;

for i = 1:numel(props)
    
    for j = 1:numel(tagNames) % skip the figure at the top
        
        % if the property exists and it's not the tag name, rewrite it
        if isfield(handles,tagNames{j}) && isprop(handles.(tagNames{j}),(props{i})) && ~strcmp(tagNames{j},'Tag')
            
            currentProp = handles.(tagNames{j}).(props{i}); % current object property
            userProp = propertyList.(props{i}){j};          % user defined property
            
            %Stop refresh button from automatically resizing to cover
            %export button
            if strcmp(tagNames{j},'refresh_gui')
                continue
            end
            
            % I really need to make this cleaner. Too many if statements
            if ~isequal(currentProp,userProp)
                % don't change the position of the window
                if ~strcmp(tagNames{j},'figure1')
                    % don't edit the systemmenu names
                    if ~strcmp(tagNames{j},'systemmenu') && ~strcmp(props{i},'UserData')
                        % if ~strcmp(tagNames{j},'systemmenu') && ~strcmp(props{i},'String')
                        % don't edit the shape change menu names
                        if ~strcmp(tagNames{j},'shapechangemenu') && ~strcmp(props{i},'UserData')
                            % if ~strcmp(tagNames{j},'shapechangemenu') && ~strcmp(props{i},'String')
                            % don't edit the metric stretch menu
                            if ~strcmp(tagNames{j},'stretchmenu') && ~strcmp(props{i},'UserData')
                                % if ~strcmp(tagNames{j},'stretchmenu') && ~strcmp(props{i},'String')
                                set(handles.(tagNames{j}),(props{i}),userProp);
                                % end
                            end
                            % end
                        end
                        % end
                    end
                end
            end
        end
    end
end

end
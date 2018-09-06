function objects = draw_set(objects)
%draw all the objects in the input structure array

%ensure that the global coordinate vertices of all objects are updated
for i = 1:numel(objects)
    
    objects(i) = update_vertices(objects(i));
    
end

%draw all the objects
for i = 1:numel(objects)
    
    objects(i) = draw_object(objects(i));
    
end
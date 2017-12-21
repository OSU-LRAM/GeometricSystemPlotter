% test points within boundary
% dummy file. may turn into function later. 

clear; 

n = 200;
mn = 0;
mx = 10;
x = rand(n,1)*mx;
y = rand(n,1)*mx;

ng = 11; % number of grid spaces
agrid = linspace(mn,mx,ng); % generate grid along one dimension
srange = (agrid(2)-agrid(1))/2; % radius of boundary

[a1,a2] = ndgrid(agrid,agrid);

bound_type = 'square';

for i = 1:ng
    
    for j = 1:ng
        
        d1 = [];
        d2 = []; % initialize deltas for later use
        
        pt = [a1(i,j) a2(i,j)]; % 1x2 coordinate matrix for grid location
        figure(1); cla
        hold on
        grid on
        box on
        scatter(x,y,'b+');
        for k = 1:n 
       
            if strcmp(bound_type,'circle')
                viscircles([a1(i,j),a2(i,j)],srange);
            elseif strcmp(bound_type,'square')
                plot([pt(1)-srange,pt(1)+srange,pt(1)+srange,pt(1)-srange,pt(1)-srange],[pt(2)-srange,pt(2)-srange,pt(2)+srange,pt(2)+srange,pt(2)-srange],'r')
            end
            
            % if the data point is within the defined range, store the
            % index and the value at that point.
            if strcmp(bound_type,'circle')
                if ((pt(1) - x(k))^2 + (pt(2) - y(k))^2) <= srange
                    
                    % add [x_index,delta] to end matrix
                    d1 = [d1;k,(x(k)-pt(1))]; %#ok<AGROW>
                    d2 = [d2;k,(y(k)-pt(2))]; %#ok<AGROW>
                    
                    plot(x(k),y(k),'ko')
                    %                 pause(0.1)
                    
                end
            elseif strcmp(bound_type,'square')
                if abs(pt(1) - x(k)) <= srange && abs(pt(2) - y(k)) <= srange
                    
                    % add [x_index,delta] to end matrix
                    d1 = [d1;k,(x(k)-pt(1))]; %#ok<AGROW>
                    d2 = [d2;k,(y(k)-pt(2))]; %#ok<AGROW>
                    
                    plot(x(k),y(k),'ko')
                    %                 pause(0.1)
                    
                end
            end
        end
        
        axis equal
        axis([mn mx mn mx])
%         disp(['points within range: ',num2str(size(d1,1))])
        pause(0.05)
    end
end
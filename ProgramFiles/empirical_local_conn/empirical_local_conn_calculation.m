% Local connection matrix from alpha_dots and deltas

function output = empirical_local_conn_calculation(xd,ad,d1,d2)

ad_matrix = [];
xbd = [];

% if there aren't any points or there aren't enough points within range
if isempty(d1) || size(d1,1) < 6
    A_matrix = zeros(6,1);
else
    for k = 1:size(d1,1)
        
        idx = d1(k,1);                  % find index which xd array uses
        xbd(k,1) = xd(idx); %#ok<AGROW> % find associated body velocity value
        ad_matrix = [ad_matrix;...      % build ad matrix for pseudo-inverse
            ad(idx,1),ad(idx,2),ad(idx,1)*d1(k,2),ad(idx,1)*d2(k,2),ad(idx,2)*d1(k,2),ad(idx,2)*d2(k,2)]; %#ok<AGROW>
        %      ad1       ad2         ad1*d1            ad1*d2            ad2*d1            ad2*d2
        
    end
    
    if size(ad_matrix,1) ~= size(xbd,1)
        error('Alpha_dot and x_dot matrices not similar sizes.')
    end
    
    if any(isnan(xbd))  %~isempty(find(isnan(xbd),1)) 
     
        nanid = isnan(xbd);
        [nanrows,~] = find(isnan(xbd));
        
        xbd(nanid) = [];
        ad_matrix(unique(nanrows),:) = [];
                
    elseif any(isnan(ad_matrix))
        
        nanid = any(isnan(ad_matrix));
        [nanrows,~] = find(isnan(ad_matrix));
        
        xbd(nanid) = [];
        ad_matrix(unique(nanrows),:) = [];
    
    end
    
    A_matrix = pinv(ad_matrix)*xbd;

end

output = A_matrix;

end
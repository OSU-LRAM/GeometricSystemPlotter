function output = color_Red

% Spot color for output
output.spot = [234 14 30]/255;

% Secondary color (e.g., grey)
output.secondary = [100 100 118]/255;
    
% Colormap for surfaces and other "intense color" plots
load('BlackWhiteRedColormap','blackwhitered');
output.colormap = blackwhitered;
            
% Colormap for contour plots
load('BlackGreyRedColormap','blackgreyred');
output.colormap_contour = blackgreyred;
    
end


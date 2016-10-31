function output = color_OSU

% Spot color for output
output.spot = [216 90 26]/230;

% Secondary color (e.g., grey)
output.secondary = [171 175 166]/300;
    
% Colormap for surfaces and other "intense color" plots
load('OSUColormap','OSU');
output.colormap = OSU;
            
% Colormap for contour plots
load('OSU_contourColormap','OSU_contour');
output.colormap_contour = OSU_contour;
    
end
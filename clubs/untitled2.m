% Load the .mat file
data = load('clubs.mat');

% Extract the variables
sections = data.sections; % 7x4 matrix representing cutting lines
bbox_poly = data.bbox_poly; % 4x2 matrix representing bounding box vertices

% Write the 'sections' to 'cutting_lines.txt'
dlmwrite('cutting_lines.txt', sections, 'delimiter', ' ');

% Write the 'bbox_poly' to 'bbox_poly.txt'
dlmwrite('bbox_poly.txt', bbox_poly, 'delimiter', ' ');

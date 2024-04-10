% Download intersections.m from http://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections
% and put it in a folder on your MATLAB path. Put this function on your path too.

function [xi, yi] = polyxpoly(x1, y1, x2, y2, opt)
    [xi, yi] = intersections(x1, y1, x2, y2);
end
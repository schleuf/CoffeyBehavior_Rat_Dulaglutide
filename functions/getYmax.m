function [ymax] = getYmax(data)
    ymax = nanmax(data);
    ymax = ymax + .05 * ymax;
end
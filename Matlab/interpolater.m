function [plotdata] = interpolater(xvals, yvals, listdata)
    % turns a list of xyz values into an interpolated surface over xvals,
    % yvals
    [xq, yq] = meshgrid(xvals, yvals);
    plotdata = griddata(listdata(:, 1), listdata(:, 2), listdata(:, 3), xq, yq);
    
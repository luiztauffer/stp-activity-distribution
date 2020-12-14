
function mymap = my_cmap(minval,maxval,centerval)

bottomcolor = [0 0 1];      % color for minimum data value (blue = [0 0 1])
indexColor = [1 1 1];       % color for indexed data value (white = [1 1 1])
topcolor = [1 0 0];         % color for maximum data value (red = [1 0 0])

% Calculate where proportionally indexValue lies between minimum and
% maximum values
indexValue = centerval;             % value for which to set a particular color
largest = maxval;
smallest = minval;
index = 100*abs(indexValue-smallest)/(largest-smallest);

% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
            linspace(bottomcolor(2),indexColor(2),100*index)',...
            linspace(bottomcolor(3),indexColor(3),100*index)'];
        
% Create color map ranging from index color to top color
% Multipling number of points by 100 adds more resolution
customCMap2 = [linspace(indexColor(1),topcolor(1),100*(100-index))',...
            linspace(indexColor(2),topcolor(2),100*(100-index))',...
            linspace(indexColor(3),topcolor(3),100*(100-index))'];
        
mymap = [customCMap1;customCMap2];  % Combine colormaps
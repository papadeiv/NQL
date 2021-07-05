function color = cmap(ndata,flag)
    % Interpolation method on 128 points with multi-colors (8 colors)
    percentage = [100;85;70;55;40;25;10;0];
    rgb = [255,0,0; % color 1
           204,204,0; % color 2
           51,255,51; % color 3
           102,255,255; % color 4
           0,51,102; % color 5
           102,102,255; % color 6
           204,0,204; % color 7
           255,153,204]; % color 8
    normalised = rgb./255;
    map = interp1(percentage,normalised,linspace(100,0,128), 'pchip');
    % Continuous gradient method with 2 colors
    color1 = [255,51,255]/255;
    color2 = [0,255,255]/255;
    if flag==1
        color = map(ceil(linspace(1,128,ndata)),:);
    else
        color = [linspace(color1(1),color2(1),ndata)', linspace(color1(2),color2(2),ndata)', linspace(color1(3),color2(3),ndata)'];
    end
end
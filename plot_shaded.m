function plot_shaded(x, y1, y2, clr, sat)
    % Plot a shaded between line y1 and y2, assuming the x axis is the same
    % for both. 
    %
    % colour = clr
    % saturation = sat

    x2 = [x, fliplr(x)];
    inBetween = [y1, fliplr(y2)];
    
    % apply saturation
    pltclr = clr*sat; 
    pltclr(pltclr>1) = 1;
    
    %plot
    fill(x2, inBetween, pltclr,'LineStyle','none')
end
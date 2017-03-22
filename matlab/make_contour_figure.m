function make_contour_figure(vertices,u,v)
%MAKE_CONTOUR_FIGURE Summary of this function goes here
%   Detailed explanation goes here
    xmin = 0;
    xmax = 1;
    ymin = 0;
    ymax = 1;
    xlin = linspace(xmin,xmax,300);
    ylin = linspace(ymin,ymax,300);
    [X,Y] = meshgrid(xlin,ylin);
    U = griddata(vertices(:,1),vertices(:,2),u,X,Y,'linear');
    V = griddata(vertices(:,1),vertices(:,2),v,X,Y,'linear');
 

    figure
    subplot(1,2,1)
    contourf(X,Y,U,10)
    xlim([xmin xmax])
    ylim([ymin ymax])
    subplot(1,2,2)
    contourf(X,Y,V,10)
    xlim([xmin xmax])
    ylim([ymin ymax])
    figure
    subplot(1,2,1)
    mesh(X,Y,U)
    xlim([xmin xmax])
    ylim([ymin ymax])
    subplot(1,2,2)
    mesh(X,Y,V)
    xlim([xmin xmax])
    ylim([ymin ymax])

end


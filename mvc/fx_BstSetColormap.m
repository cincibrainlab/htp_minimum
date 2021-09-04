function cfig = fx_BstSetColormap(hFig, min, max)
% adjust color map for consistency
currentColormapInfo = getappdata(hFig, 'Colormap');
currentColorMap = bst_colormaps('GetColormap', hFig);
currentColorMap.MinValue = min;
currentColorMap.MaxValue = max;
bst_colormaps('SetColormap', ...
    currentColormapInfo.Type, currentColorMap);
bst_colormaps('FireColormapChanged');

% get separate colorbar image

%%
cfig = figure(965);
set(gcf,'color','w');
ax = axes;
c = colorbar(ax,'South');

c.Ticks = [-5 5];
c.FontSize = 70;
colormap(currentColorMap.CMap);
caxis([min max]);
ax.Visible = 'off';
x=get(c,'Position');
x(3)=x(3)*1; 
x(4)=x(4)*1.5; 
x1=get(gca,'position');

set(c,'Position',x)
set(gca,'position',x1)

F = getframe(gcf);
[currentColorBar2, ~] = frame2im(F);
%[~, b] = imcrop(currentColorBar2);
%figure; imshow(imcrop(currentColorBar2, [74.5100  318.5100  225.9800   48.9800]));
cfig = imcrop(currentColorBar2, [   42.5100  211.5100  481.9800  155.9800]);

close(figure(965));

end
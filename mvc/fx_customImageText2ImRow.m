function imgout = fx_customImageText2ImRow( textstr, view_size )

% put a small x in the center of a new figure
label_figure = figure('menubar','none') ;
set(gcf,'color','w');
ah = gca ;
th = text(1,1,textstr,'FontSize',75) ;
set(ah,'visible','off','xlim',[0 2],'ylim',[0 2],'Position',[0 0 1 1]) ;
set(th,'visible','on','HorizontalAlignment','center','VerticalAlignment','middle')
F = getframe(label_figure);
[X, ~] = frame2im(F);
imgout = fx_customImageCropImage(X, [150 view_size(2)]);
close(label_figure);

end
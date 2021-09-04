function freq_label2 = fx_customImage3dPlotRowLabels( row_labels, column_label )


whitesquare = 255 .* ones(size(column_label,1), ....
    size(row_labels,2),3, 'uint8');


currentColorBar = fx_BstSetColormap(hFig, colorMapMin, colorMapMax);
position_small_image = [25 0];
A = whitesquare;
B = imresize(currentColorBar,.9);
% B=padarray(B,[1 1]); %If you want borders
Y=position_small_image(1);
X=position_small_image(2);
A((1:size(B,1))+X,(1:size(B,2))+Y,:) = B;
freq_label2 = vertcat(freq_label, A);


end
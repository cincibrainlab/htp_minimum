function img = fx_customImageMakeWhiteBackground( img )
[L,W,D] = size(img);
for i = 1:L
    for j = 1:W
        for k = 1:D
            if A(i,j,k) == 0
                A(i,j,k) = 255;
            end
        end
    end
end
end
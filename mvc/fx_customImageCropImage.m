function croppedimg = fx_customImageCropImage(rawimg, imgsize)

    imgsize_original = size(rawimg);
    imgsize_original = imgsize_original(1:2);
    %imgsize_original = [imgsize_original(2) imgsize_original(1)]; % Width, Height

    win = centerCropWindow2d(imgsize_original, imgsize);
    croppedimg = imcrop(rawimg, win);

end
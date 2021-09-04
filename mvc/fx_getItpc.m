%% Optimized ITPC
function itpc = fx_getItpc( bstsub )

% sample data for cfg
sampleDat = reshape(bstsub.ImageGridAmp(1,:), 1626, []);

cfg.centerfreq = 40;
cfg.srate = 500;
cfg.pnts = size(sampleDat,1);
cfg.trials = size(sampleDat,2);

cfg.n_wavelet = cfg.pnts;
cfg.n_data = cfg.pnts * cfg.trials;
cfg.n_convolution = cfg.n_wavelet+cfg.n_data-1;
cfg.n_conv_pow2   = pow2(nextpow2(cfg.n_convolution));

cfg.wavelettime   = -cfg.pnts/cfg.srate/2:1/cfg.srate: ...
    cfg.pnts/cfg.srate/2-1/cfg.srate;
cfg.wavelet = exp(2*1i*pi*cfg.centerfreq.*cfg.wavelettime) .* ...
    exp(-cfg.wavelettime.^2./(2*((4/(2*pi*cfg.centerfreq))^2)))/...
    cfg.centerfreq;
cfg.timeVec = bstsub.Time(1:cfg.pnts) .* 1000;

itpc = zeros(size( bstsub.ImageGridAmp, 1 ),1);
for vert_i = 1 : size( bstsub.ImageGridAmp, 1 )
    
    dat = gpuArray( reshape(bstsub.ImageGridAmp(vert_i,:), 1626, []) );
    dat = fft(reshape(dat(:,:),1,[]),cfg.n_conv_pow2);
    eegconv = ifft(fft(cfg.wavelet,cfg.n_conv_pow2).*dat);
    %figure; plot(abs(eegconv(1:1000)))
    eegconv = eegconv(1:cfg.n_convolution);
    
    eegconv = reshape(eegconv(floor((cfg.pnts-1)/2):end-1-ceil((cfg.pnts-1)/2)), ...
        cfg.pnts,cfg.trials);
    
    
    idx1 = fx_itc_getTimeIdx(cfg.timeVec, 750);
    idx2 = fx_itc_getTimeIdx(cfg.timeVec, 1000);

    itpc(vert_i) = mean((round(1000*abs(mean(exp(1i*angle(eegconv(idx1:idx2,:))))))/1000)');
    
end
end


%%

function idx = fx_itc_getTimeIdx( timevector, searchTime )

[~, idx] = min(abs(timevector-searchTime));

end

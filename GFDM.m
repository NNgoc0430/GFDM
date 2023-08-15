close all;
clear;
if not(isfolder("figure"))
    mkdir("figure")
end
%% (.p1) parameter sets for GFDM 
gfdm = struct;
gfdm.K = 512; % Number of samples per sub-symbol
gfdm.M = 15; % Number of sub-symbols
gfdm.gfdmM = gfdm.M; 
gfdm.blockLength = gfdm.M*gfdm.K; 
gfdm.CP = 0.1; % Percentage of cyclic prefix
gfdm.pulse = 'rc'; % Pulse shaping filter ("rc" is raised cosine)
gfdm.a = 0.1; 
gfdm.mu = 4; % Modulation order of the QAM symbol
gfdm.subcarriers = 1:201; 
gfdm.subsymbols = 2:15; 
gfdm.blocks = 10; 
%% (.p2) parameter sets for OFDM 
ofdm = struct;
ofdm.K = 512;
ofdm.M = 1;
ofdm.gfdmM = gfdm.M;
ofdm.blockLength = ofdm.K;
ofdm.CP = 0.1;
ofdm.pulse = 'rc_td'; % ("rc_td" is a rectangular filter in frequency domain)
ofdm.a = 0.1;            
ofdm.mu = 4;
ofdm.subcarriers = 1:201;
ofdm.subsymbols = 1;
ofdm.blocks = gfdm.blocks;
%% (.s1) Simulate GFDM and OFDM transmission
[GFDMS, sym_gfdm_tx, cpx_gfdm_sym_tx, g1] = tx_simul(gfdm);
[OFDMS, sym_ofdm_tx, cpx_ofdm_sym_tx, g2] = tx_simul(ofdm);
%% (.s2) Simulate reception with different SNR values 
snr = 1:21; gfdm_ser = []; ofdm_ser = [];
% execution_rx_GFDM = []; execution_rx_OFDM = [];
[theory_ser] = berawgn(snr, "qam", gfdm.mu^2);
for sii=snr
%     tic
    [sym_rx1, cpx_sym_rx1, ser1] = rx_simul(gfdm, awgn(GFDMS, sii), sym_gfdm_tx, g1);
%     execution_rx_GFDM = [execution_rx_GFDM; toc];
%     tic
    [sym_rx2, cpx_sym_rx2, ser2] = rx_simul(ofdm, awgn(OFDMS, sii), sym_ofdm_tx, g2);
%     execution_rx_OFDM = [execution_rx_OFDM; toc];
    gfdm_ser = [gfdm_ser; ser1];
    ofdm_ser = [ofdm_ser; ser2];
end
%% (.f1) Plot the PSD of the GFDM/OFDM signals
    f = figure("Name", "PSD");
    length = min(length(OFDMS), length(GFDMS));
    t = linspace(-gfdm.K/2, gfdm.K/2, 2*length+1); 
    t = t(1:end-1)';
    plot(t, mag2db(fftshift(abs(fft(OFDMS(1:length), 2*length))))/2, 'b');
    hold on;
    plot(t, mag2db(fftshift(abs(fft(GFDMS(1:length), 2*length))))/2, 'r');
    hold off;
    ylim([-40, 30]);
    xlabel('f/F'); ylabel('PSD [dB]');
    grid();
    legend({'OFDM', 'GFDM'});
    print(['figure' filesep 'PSD_gfdm_ofdm'], "-dpng");
    set(f, 'visible', 'off'); 
    %% (.f2) Plot the Eye Diagrams/Scatter Plots for different SNR value 
    if mod(sii, 3) == 0 
        print_eye_scatter("gfdm", gfdm, cpx_sym_rx1, sii);
        print_eye_scatter("ofdm", ofdm, cpx_sym_rx2, sii);
    end
    %% (.f3) Plot the SER values of the GFDM/OFDM signals 
print_ser_q('SER', snr, theory_ser, gfdm_ser, ofdm_ser);
% (.es) Function that can print Scatter plot of different formats
function print_eye_scatter(name, opt, sym, sii)
    f2 = scatterplot(sym, 1, 0, 'r.');
    hold on;
    scatterplot(qammod(0:opt.mu^2-1, opt.mu^2, 'gray'), 1, 0, 'k*', f2);
    title(strcat('Scatter plot,', {' '}, name, {', '}, 'SNR=', string(sii), 'dB'));
    print(strcat('figure', filesep, name, '_scatter_snr_', string(sii)), "-dpng");
    set(f2, 'visible', 'off'); 
end

% (.bsq) Function that can print SER
function print_ser_q(name, snr, theory, gfdm, ofdm)
    f = figure('Name', name);
    if (size(theory) > 0)
        semilogy(snr, theory,'-x');
        hold on;
    end
    semilogy(snr, gfdm, 'o');
    hold on;
    semilogy(snr, ofdm, '*');
    if (size(theory) > 0)
        legend('theory', 'gfdm', 'ofdm');
    else 
        legend('gfdm', 'ofdm');
    end
    xlabel('SNR, dB');
    ylabel(name);
    xlim([1, 20]);
    print(['figure' filesep name], "-dpng");
    set(f, 'visible', 'off'); 
end
% (.t) Simulation of transmitter GFDM/OFDM
function [signal, sym_tx, cpx_sym_tx, g] = tx_simul(opt)
    signal = []; sym_tx = []; cpx_sym_tx = [];
    
    % (.t1) Generate the transmitter pulse g
    if strcmp(opt.pulse, 'rc') 
        t = linspace(-opt.M/2, opt.M/2, opt.M*opt.K + 1); 
        t = t(1:end-1); 
        t = t';
        g = (sinc(t) .* cos(pi*opt.a*t) ./ (1-4*opt.a*opt.a*t.*t));
        g = fftshift(g);
        g(opt.K+1:opt.K:end) = 0;
        g = g / sqrt(sum(g.*g));
    elseif strcmp(opt.pulse, 'rc_td')
        g = zeros(opt.M*opt.K, 1);
        g(1:opt.K) = 1;
        g = g / sqrt(sum(g.*g));
    end    
    for b = 1:opt.blocks  
        gfdmM = 1;
        
        % (.t2) OFDM case
        if opt.M == 1
            gfdmM = opt.gfdmM;
        end  
        for m = 1 : gfdmM 
            
            % (.t3) Generate random data and map to QAM constellation
            sym = randi(2 ^ opt.mu, length(opt.subsymbols) * length(opt.subcarriers), 1) - 1;
            cpx_sym = qammod(sym, 2^opt.mu, 'gray')/ sqrt(2/3 * (2^opt.mu - 1));
            sym_tx = [sym_tx; sym];
            cpx_sym_tx = [cpx_sym_tx; cpx_sym];
            
            % (.t4) Map data to correct dimensions or pad it with zeros
            if length(opt.subcarriers)== opt.K && length(opt.subsymbols) == opt.M
                D = reshape(cpx_sym, opt.K, opt.M);
            else
                Dm = reshape(cpx_sym, length(opt.subcarriers), length(opt.subsymbols));
                res1 = zeros(opt.K, length(opt.subsymbols));
                res1(opt.subcarriers, :) = Dm;
                res = zeros(opt.K, opt.M);
                res(:, opt.subsymbols) = res1;
                D = res;
            end
            
            % (.t5) Modulate in time
            DD = repmat(opt.K*ifft(D), opt.M, 1);
            x = zeros(opt.K*opt.M, 1);
            for m=1:opt.M
                symbol = DD(:, m) .* g;
                symbol = circshift(symbol, opt.K*(m-1));
                x = x + symbol;
            end
            
            % (.t6) Add CP
            cp = ceil(opt.CP*opt.blockLength);
            xcp = [x(end-cp + (1:cp), :); x];
            signal = [signal; xcp];
        end
    end
end

% (.r) Simulation of receiver GFDM/OFDM
function [sym_rx, cpx_sym_rx, ber, ser] = rx_simul(opt, signal, sym_tx, g)
    sym_rx = []; cpx_sym_rx = []; ber = 0; ser = 0;
    for b = 1:opt.blocks 
        cp = ceil(opt.CP*opt.blockLength);
        gfdmM = 1;
        dim1 = opt.blockLength + cp;
        dim2 = dim1;       
        
        % (.r1) OFDM case
        if opt.M == 1
            gfdmM = opt.gfdmM;
            dim1 = opt.gfdmM*opt.K + opt.gfdmM*cp;
            dim2 = opt.K + cp;
        end       
        for m = 1 : gfdmM   
            block = signal((b-1)*dim1 + (m-1)*dim2 + (1:dim2));
            
            % (.r2) Remove CP
            block = block(cp + 1:end);
            
            % (.r3) Calculate the transmitter pulse g (Matched Filter)
            g = g(round(1:opt.K/opt.K:end));
            G = conj(fft(g));
            
            % (.r4) Demodulate in frequency
            L = length(G) / opt.M; 
            Xhat = fft(block);
            Dhat = zeros(opt.K, opt.M);
            for k=1:opt.K
                carrier = circshift(Xhat, ceil(L*opt.M/2) - opt.M*(k-1));
                carrier = fftshift(carrier(1:L*opt.M));
                carrierMatched = carrier .* G;
                dhat = ifft(sum(reshape(carrierMatched, opt.M, L), 2)/L);
                Dhat(k,:) = dhat;
            end
            
            % (.r5) Unmap
            if length(opt.subcarriers) == opt.K && length(subsymbols) == opt.M
                dhat_mf = reshape(Dhat, numel(Dhat), 1);
            else
                Dm = Dhat(opt.subcarriers, opt.subsymbols);
                dhat_mf = reshape(Dm, numel(Dm), 1);
            end
            
            % (.r6) Map constellation points to data
            dhat_mf = dhat_mf * sqrt(2/3 * (2^opt.mu - 1));
            shm = qamdemod(dhat_mf, 2^opt.mu, 'gray');
            cpx_sym_rx = [cpx_sym_rx; dhat_mf];
            sym_rx = [sym_rx; shm];
        end
    end
    
    % (.r7) Calculate SER
    [dim1, dim2] = size(sym_tx);
    sym_b_tx = de2bi(sym_tx, opt.mu);
    sym_b_rx = de2bi(sym_rx, opt.mu);
    err = [];
    for ii = 1:size(sym_b_tx)
        err = [err; nnz(xor(sym_b_tx(ii, :), sym_b_rx(ii, :)))];
    end
    ser = double(nnz(sym_tx - sym_rx)/dim1);
    ber = sum(err)/(dim1*4);   
end    
    
    

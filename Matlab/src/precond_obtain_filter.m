%% precond_obtain_filter: obtain preconditioning filter
% Author: Yongjie Zhuang
% Fxx: estimated filter such that ref signal x = Fxx * v, v is white noise,
% dimension is N_xx*Nr*Nr
% h_Ge: estimated filter of secondary path, dimension is Nt_Ge*Ne*Ns
% N_Finv: length of filter to estimate the inverse of Fxx
% N_mininv: length of filter to estimate the inverse of Ge_min and Ge_all
% fig_flag: 1: plot, 0: do not plot
function [Fxx_inv,Ge_min_inv,Ge_all]= precond_obtain_filter(Fxx,h_Ge,N_Finv,N_mininv,plot_flag)
%% obtain dimension and check input data
[N_xx,Nr,Nr2] = size(Fxx);
[N_Ge,Ne,Ns] = size(h_Ge);
comp1 = [Nr];
comp2 = [Nr2];
if any(comp1~=comp2)
    error('The dimension of input data does not match');
end
%% Fit F_inv
den_fac = 8; % density factor, determine how many points being fitted
Fxx_freq = zeros(Nr,Nr,N_Finv*den_fac);
% obtain freq domain
for ii = 1:Nr
    for jj = 1:Nr
        [h,w] = freqz(Fxx(:,ii,jj),1,N_Finv*den_fac);
        Fxx_freq(ii,jj,:) = h;
    end
end
% get inverse
for ii = 1:N_Finv*den_fac
    Fxx_freq(:,:,ii) = inv(Fxx_freq(:,:,ii));
end
% fit Fxx_inv
Fxx_inv = zeros(N_Finv,Nr,Nr);
for ii = 1:Nr
    for jj = 1:Nr
        H_measure = squeeze(Fxx_freq(ii,jj,:));
        Fxx_inv(:,ii,jj) = invfreqz(H_measure,w,N_Finv-1,0);
    end
end

if plot_flag
    % frequency plot
    % get a large plot
    figure('Units','normalized','Position',[0.1 0.1 0.8 0.8])
    figureH = tiledlayout(Nr, Nr*2);
    title(figureH, 'Inverse of Fxx Frequency Response')
    for ii = 1:Nr
        for jj = 1:Nr
            idx = (ii-1)*Nr + jj;
            % Magnitude subplot
            nexttile((idx-1)*2+1)
            [H_fit, w] = freqz(Fxx_inv(:,ii,jj), 1, N_Finv*den_fac);
            H_measure = squeeze(Fxx_freq(ii,jj,:));
            plot(w, 20*log10(abs(H_measure)), 'r-', 'LineWidth', 2)
            hold on
            plot(w, 20*log10(abs(H_fit)), 'b--', 'LineWidth', 2)
            xlabel('Frequency (rad/sample)')
            ylabel('Magnitude (dB)')
            title(['Mag: (' num2str(ii) ',' num2str(jj) ')'])
            grid on

            % Phase subplot
            nexttile((idx-1)*2+2)
            plot(w, angle(H_measure), 'r-', 'LineWidth', 2)
            hold on
            plot(w, angle(H_fit), 'b--', 'LineWidth', 2)
            xlabel('Frequency (rad/sample)')
            ylabel('Phase (rad)')
            ylim([-pi,pi])
            title(['Phase: (' num2str(ii) ',' num2str(jj) ')'])
            grid on
        end
    end
    legend('Measured','Fitted')
    
    % impulse response plot
    figure
    for ii = 1:Nr
        for jj = 1:Nr
            idx = (ii-1)*Nr + jj;    % subplot index
            subplot(Nr, Nr, idx);    % create Nr x Nr grid of subplots
            plot(Fxx_inv(:,ii,jj),'LineWidth',2);
            title(['Inverse Fxx: ' num2str(ii) ' ' num2str(jj)]);
            xlabel('Samples');
            ylabel('Amplitude');
            grid on;
        end
    end
end


    
%% Fit G_min
GeGe = zeros(2*N_Ge-1,Ns,Ns);
for ii = 1:Ns
    for jj = ii:Ns
        for kk = 1:Ns
            GeGe(:,ii,jj) = GeGe(:,ii,jj)+ conv(flip(h_Ge(:,kk,ii)),h_Ge(:,kk,jj));
        end
    end
end
% only half of GeGe is needed, the other half is filled to lower triangle
GeGe_p = zeros(N_Ge,Ns,Ns);
for ii = 1:Ns
    for jj = ii:Ns
        GeGe_p(:,ii,jj) = GeGe(N_Ge:end,ii,jj);
        GeGe_p(:,jj,ii) = GeGe(N_Ge:-1:1,ii,jj);
    end
end
% get spectral factorization
[Ge_min_H,~,~]= mul_specFact(GeGe_p);
% note that we get Ge_min', we need to recover the Ge_min
Ge_min = zeros(N_Ge,Ns,Ns);
for ii = 1:Ns
    for jj = 1:Ns
        Ge_min(:,ii,jj) = Ge_min_H(:,jj,ii);
    end
end
%% Fit G_min_inv
Ge_min_freq = zeros(Ns,Ns,N_mininv*den_fac);
% obtain freq domain
for ii = 1:Ns
    for jj = 1:Ns
        [h,w] = freqz(Ge_min(:,ii,jj),1,N_mininv*den_fac);
        Ge_min_freq(ii,jj,:) = h;
    end
end
% get inverse
for ii = 1:N_mininv*den_fac
    Ge_min_freq(:,:,ii) = inv(Ge_min_freq(:,:,ii));
end
% fit G_min_inv
Ge_min_inv = zeros(N_mininv,Ns,Ns);
for ii = 1:Ns
    for jj = 1:Ns
        H_measure = squeeze(Ge_min_freq(ii,jj,:));
        Ge_min_inv(:,ii,jj) = invfreqz(H_measure,w,N_mininv-1,0);
    end
end
if plot_flag
    % frequency plot
    figure('Units','normalized','Position',[0.1 0.1 0.8 0.8])
    figureH = tiledlayout(Ns, Ns*2);
    title(figureH, 'Inverse of G_{min} Frequency Response')
    for ii = 1:Ns
        for jj = 1:Ns
            idx = (ii-1)*Ns + jj;
            % Magnitude subplot
            nexttile((idx-1)*2+1)
            [H_fit, w] = freqz(Ge_min_inv(:,ii,jj), 1, N_mininv*den_fac);
            H_measure = squeeze(Ge_min_freq(ii,jj,:));
            plot(w, 20*log10(abs(H_measure)), 'r-', 'LineWidth', 2)
            hold on
            plot(w, 20*log10(abs(H_fit)), 'b--', 'LineWidth', 2)
            xlabel('Frequency (rad/sample)')
            ylabel('Magnitude (dB)')
            title(['Mag: (' num2str(ii) ',' num2str(jj) ')'])
            grid on

            % Phase subplot
            nexttile((idx-1)*2+2)
            plot(w, angle(H_measure), 'r-', 'LineWidth', 2)
            hold on
            plot(w, angle(H_fit), 'b--', 'LineWidth', 2)
            xlabel('Frequency (rad/sample)')
            ylabel('Phase (rad)')
            ylim([-pi,pi])
            title(['Phase: (' num2str(ii) ',' num2str(jj) ')'])
            grid on
        end
    end
    legend('Measured','Fitted')

    % impulse response plot
    figure
    for ii = 1:Ns
        for jj = 1:Ns
            idx = (ii-1)*Ns + jj;    % subplot index
            subplot(Ns, Ns, idx);    % create Ns x Ns grid of subplots
            plot(Ge_min_inv(:,ii,jj),'LineWidth',2);
            title(['G_{min}^{-1}: ' num2str(ii) ' ' num2str(jj)]);
            xlabel('Samples');
            ylabel('Amplitude');
            grid on;
        end
    end
end


%% fit G_all
Ge_freq = zeros(Ne,Ns,N_mininv*den_fac);
% obtain freq domain
for ii = 1:Ne
    for jj = 1:Ns
        [h,w] = freqz(h_Ge(:,ii,jj),1,N_mininv*den_fac);
        Ge_freq(ii,jj,:) = h;
    end
end
% get G_all_freq
Ge_all_freq = zeros(Ne,Ns,N_mininv*den_fac);
for ii = 1:N_mininv*den_fac
    Ge_all_freq(:,:,ii) = Ge_freq(:,:,ii)*Ge_min_freq(:,:,ii);
end
Ge_all = zeros(N_mininv,Ne,Ns);
for ii = 1:Ne
    for jj = 1:Ns
        H_measure = squeeze(Ge_all_freq(ii,jj,:));
        Ge_all(:,ii,jj) = invfreqz(H_measure,w,N_mininv-1,0);
    end
end
if plot_flag
    % frequency plot
    figure('Units','normalized','Position',[0.1 0.1 0.8 0.8])
    figureH = tiledlayout(Ne, Ns*2);
    title(figureH, 'Ge_{all} Frequency Response')
    for ii = 1:Ne
        for jj = 1:Ns
            idx = (ii-1)*Ns + jj;
            % Magnitude subplot
            nexttile((idx-1)*2+1)
            [H_fit, w] = freqz(Ge_all(:,ii,jj), 1, N_mininv*den_fac);
            H_measure = squeeze(Ge_all_freq(ii,jj,:));
            plot(w, 20*log10(abs(H_measure)), 'r-', 'LineWidth', 2)
            hold on
            plot(w, 20*log10(abs(H_fit)), 'b--', 'LineWidth', 2)
            xlabel('Frequency (rad/sample)')
            ylabel('Magnitude (dB)')
            title(['Mag: (' num2str(ii) ',' num2str(jj) ')'])
            grid on

            % Phase subplot
            nexttile((idx-1)*2+2)
            plot(w, angle(H_measure), 'r-', 'LineWidth', 2)
            hold on
            plot(w, angle(H_fit), 'b--', 'LineWidth', 2)
            xlabel('Frequency (rad/sample)')
            ylabel('Phase (rad)')
            ylim([-pi,pi])
            title(['Phase: (' num2str(ii) ',' num2str(jj) ')'])
            grid on
        end
    end
    legend('Measured','Fitted')

    % impulse response plot
    figure
    for ii = 1:Ne
        for jj = 1:Ns
            idx = (ii-1)*Ns + jj;     % subplot index
            subplot(Ne, Ns, idx);     % create Ne x Ns grid of subplots
            plot(Ge_all(:,ii,jj), 'LineWidth', 2);
            title(['G_{all}: ' num2str(ii) ' ' num2str(jj)]);
            xlabel('Samples');
            ylabel('Amplitude');
            grid on;
        end
    end
end


end
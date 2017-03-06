clear all;
% x = [1:4096];
% yx = 30*cos(2*x*pi/76)';
% zx = 20*sin(2*x*pi/76)';
% yx = zeros(4096,1);
% zx =[1:4096]';
data_in = [];
data_snr = [];
data_snr_all = [];

data_all = [];
%% 10.5
% for mk=5:1:30
for bandbit = 9:1:64;
    for n= 1:20
        yx = rand(1,4096,'double')';
        zx = rand(1,4096,'double')';
        
        a = complex(yx,zx);
        ES = mean(abs((a).^2));
        a = a/sqrt(ES);
        
        data_in = [data_in,a];
        
        temp_m = 2^11;
%         bandbit = 12;
        databit = 18;
        data_out = matlab_ifft(round(a*temp_m),bandbit ,databit)/temp_m;
        data_ifft = ifft(a)*64;
        
        %  figure(1);
        %  plot([1:4096],real(data_out),'r',[1:4096],imag(data_out),'k');
        %  figure(2);
        %  plot([1:4096],real(data_fft),'r',[1:4096],imag(data_fft),'k');
        FFT_SNR = 10*log10(mean(abs(data_ifft).^2)/mean(abs(data_ifft-data_out).^2));
        data_snr = [data_snr,FFT_SNR];
    end
%             plot([1:100],data_snr);
      data_snr_all = [data_snr_all,mean(data_snr)];
end
       plot([9:64],data_snr_all);
%     data_snr_all = [data_snr_all,mean(data_snr)];
% end
%   plot([5:30],data_snr_all);
%%
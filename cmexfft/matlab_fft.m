function data_out = matlab_fft(data_in,wr,tempbit,databit)
% ����������complex����,�����������4096*1�����ľ���
% ���������complex����,�����������4096*1
data_real = real(data_in);
data_image = imag(data_in);
temp = [data_real;data_image];

temp = [wr;databit;tempbit;temp];

data_temp = matlabfft(temp);
c = reshape(data_temp,4096,2);
data_out = complex(c(:,1),c(:,2));
end
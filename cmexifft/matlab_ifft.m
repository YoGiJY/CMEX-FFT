function data_out = matlab_ifft(data_in,wr,databit)
% ����������complex����,�����������4096*1�����ľ���
% ���������complex����,�����������4096*1
data_real = real(data_in);
data_image = imag(data_in);
temp = [data_real;data_image];

temp = [wr;databit;temp];

data_temp = matlabifft(temp);
c = reshape(data_temp,4096,2);
data_out = complex(c(:,1),c(:,2));
end
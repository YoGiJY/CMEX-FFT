function data_out = matlab_ifft(data_in,wr,databit)
% 输入数据是complex数据,输入的数据是4096*1这样的矩阵
% 输出数据是complex数据,输出的数据是4096*1
data_real = real(data_in);
data_image = imag(data_in);
temp = [data_real;data_image];

temp = [wr;databit;temp];

data_temp = matlabifft(temp);
c = reshape(data_temp,4096,2);
data_out = complex(c(:,1),c(:,2));
end
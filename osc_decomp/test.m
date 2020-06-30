%load("/Users/matsukiakari/Documents/doctor/RA_work/result_osc_decomp/result_WT130604_L");
%load("/Users/matsukiakari/Documents/doctor/RA_work/result_osc_decomp/result_WT130604_L_decomp_cov.mat")
fileID = fopen("/Users/matsukiakari/Documents/doctor/RA_work/LR_v1603/WT130604matt005.txt", "r")
formatSpec = "%f %f %f";
sizeA = [3 Inf];
A = fscanf(fileID, formatSpec, sizeA);
y =A(3,:);

fs = 1e4;
%plot_decomp(y,fs,decomp_mu,decomp_cov,osc_param,K);
%plot_phase(y,fs,decomp_phase,decomp_mu,decomp_cov,K);
%y_lowpass = 
lowpass(y,2,fs)

%disp(size(y_lowpass))
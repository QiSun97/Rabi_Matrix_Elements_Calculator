function [fft_init,step,force_endtime,fft_T]=FFTforce(nground, n, eq)
ex_pop=zeros(length(eq.evTime),1);
for i=nground+1:n
    ex_pop=ex_pop+real(squeeze(eq.evolution(i,i,:)));
end
fft_init=find(eq.evTime(:)>0,1);
fft_endtime = length(eq.evTime(:));
fft_endtime = fft_endtime-1+mod(fft_endtime-fft_init,2);

%     actual_time_start = eq.evTime(fft_init);
%     actual_time_end = eq.evTime(fft_endtime);
%     interptime = linspace(actual_time_start, actual_time_end, fft_endtime-fft_init+1);

fft_data = ex_pop(fft_init:fft_endtime);
fft_datalength = fft_endtime-fft_init+1;
fft_Y = abs(fft(fft_data));
fft_cutoff = 2;
fft_Y = fft_Y(fft_cutoff+1:(fft_datalength/2)+1);
fft_freq = (fft_cutoff:(fft_datalength)/2)/fft_datalength;
[~,peak_freq] = max(fft_Y);
fft_T = 1/(fft_freq(peak_freq));

if fft_T>1000
    step = floor(fft_T/1000)+1;
    force_endtime = floor(fft_T)+fft_init+1;
elseif fft_T<200
    step = 1;
    force_endtime = floor(floor(400/fft_T)*fft_T)+fft_init+1;
else
    step = 1;
    force_endtime = floor(fft_T)+fft_init+1;
end
if force_endtime>=fft_endtime
    disp('something is wrong, your force end time is too long')
    force_endtime = fft_endtime;
end
for index=fft_init:step:force_endtime
end
force_endtime=index;
end
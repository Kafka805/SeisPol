% FILBUTT uses a Butterworth window to filter the data.  The corner
%	  frequencies are specified in the command line, as is the 
%	  sampling rate.  If the lower limit is set to zero, then the
%	  filter is lowpass.  Similarly, if the upper limit is zero,
%	  the filter is highpass.  If both limits are set, the filter
%	  will be bandpass.
%       
%	  y = FILBUTT(Data,Samp_Rate,Low_Limit,High_Limit)
%       (c) Em Schnorr, 2022
function [filtered] = filbutt(Data,Samp_Rate,Low_Limit,High_Limit)

Nyquist=Samp_Rate/2;
Low_Param=Low_Limit/Nyquist;
High_Param=High_Limit/Nyquist;
if Low_Limit == 0
   [filb,fila]=butter(1,High_Param);
elseif High_Limit == 0
   [filb,fila]=butter(1,Low_Param,'high');
else 
   [filb,fila]=butter(1,[Low_Param High_Param]);
end
filtered = filter(filb,fila,Data);


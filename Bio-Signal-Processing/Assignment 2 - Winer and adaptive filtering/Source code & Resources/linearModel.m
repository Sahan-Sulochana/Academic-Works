%% LINEAR MODEL

function [y] = linearModel(yi)
% length = n_stop - n_start;%length of the array
% n = 1:length;


y = zeros(1,84);
for i = 1: length(yi)

    if i > 79
        y(i) = 0;
    elseif i > 70 
        y(i) = ((0.211062+0.108938) - (0))*(i - 79)/(70 - 79) + (0);  
    elseif i > 68
        y(i) = (0.211062+0.108938);
    elseif i > 55
        y(i) = ((0.211062+0.108938) - (0))*(i - 55)/(68 - 55) + (0);  
    elseif i > 40
        y(i) = 0;
    elseif i > 37
        y(i) = (0 - (-0.548938))*(i - 37)/(40 - 37) + (-0.548938);
    elseif i > 34
        y(i) = (1.83106 - (-0.548938))*(i - 37)/(34 - 37) + (-0.548938);
    elseif i > 30 
        y(i) = (1.83106 - (-0.398938))*(i - 30)/(34 - 30) + (-0.398938);
    elseif i > 28
        y(i) = ((0) - (-0.398938))*(i - 30)/(28 - 30) + (-0.398938);
    elseif i > 16
        y(i) = 0;
    elseif i > 12
        y(i) = ((0.0110618+0.128938) - (0))*(i - 16)/(12 - 16) + (0);
    elseif i > 11
        y(i) = (0.0110618+0.128938);
    elseif i > 8 
        y(i) = ((0.0110618+0.128938) - (0))*(i - 8)/(11 - 8) + (0);
    else
        y(i) = 0;
    end
end
y = y + ones(1,84)*0.001;
% tMod = mod(t,0.22);
% y = zeros(size(n));
% 
% for i = tMod
%     j = floor(i*fs) + 1;
%     %Neutral
%     if (0<=i && i<=0.006)
%         y(j) = -0.115;
%     %P wave
%     elseif (0.006 < i && i<=0.012)
%         y(j) = 35*i-0.33;
%     elseif (0.012<i && i<=0.022)
%         y(j) = -18*i+0.311;
%     %Neutral
%     elseif (0.022<i && i<=0.044)
%         y(j) = -0.115;
%     %QRS comlpex
%     elseif (0.044<i && i<=0.048)%End at Q
%         y(j) = -52.50*i+2.19;    
%     elseif (0.048<i && i<=0.054)%Q--->R 
%         y(j) = 326.67*i-16.01;     
%     elseif (0.054<i && i<=0.06)%R--->S
%         y(j) = -350*i+20.54;    
%     elseif (0.06<i && i<=0.069)
%         y(j) = 35*i-2.56;   
%     elseif (0.069<i && i<=0.094)
%         y(j) = -0.115;
%     %T wave
%     elseif (0.094<i && i<=0.12)
%         y(j) = 13.85*i-1.42;    
%     elseif (0.12<i && i<=0.13)
%         y(j) = 0.245;    
%     elseif (0.13<i && i<=0.154)
%         y(j) = -15*i+2.2;
%     elseif (0.154<i && i<=0.168)
%         y(j) = -0.115;
%     end
% end
end
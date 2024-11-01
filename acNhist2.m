function [bins,p]=acNhist2(c,figFlag)
%accurate histogram using Scott's rule

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute bins width by Scott's rule 
w=3.49*std(c) * length(c)^(-1/3);
% establish bins vector
bins=(min(c)+w/2):w:max(c);
% use hist to compute absolute frequencies
freq=hist(c,bins);
% return relative frequencies
p=freq/(length(c)*w);
%p=freq/(length(c)); 
if figFlag
    figure,plot(bins,p,'k','LineWidth',2);
    grid minor;
    title('Distribution of the signal')
end

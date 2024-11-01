clc;
clear all;
close all;

t = 0;
endTime = 6000;

n_arrived = 0;
n_departure = 0;
n_cus = 0;
m_arrived = 0;
step_s = 1;
m_left = 0;

ar_t = [];
de_t = [];
soj_t = [];
n_in_ar = [];

mu_1 = 10;
lam_1 = 7;  % 3/10 sytem is stable state with no queue length and 

queue_limit = 1000;

event = zeros(1, 3);
event(1) = exprnd(1 / lam_1);
event(2) = inf;               
event(3) = exprnd(1 / step_s); 

measure = 0; 
m_value = [];
m_t = [];
m_in_ar = [];

no_of_chefs = 31;

while t <= endTime
    [value, index] = min(event);
    if index == 1 && (isempty(n_in_ar) || n_cus <= queue_limit)  %sum(n_in_ar(max(1, n_arrived - n_cus + 1):n_arrived)
        t = event(1);
        n_arrived = n_arrived + 1;
        n_cus = n_cus + 1;
        ar_t(n_arrived) = t; 
        n_in_ar(n_arrived) = generate_group_size_exp(); 
        if n_cus == 1
            event(2) = t + exprnd(n_in_ar(end) / mu_1);
        end
        event(1) = t + exprnd(1 / lam_1);
    elseif n_cus >= queue_limit && index == 1 %(sum(n_in_ar(max(1, n_arrived - n_cus + 1):n_arrived)) 
            t = event(1);
            m_left = m_left + 1;
            m_t(m_left) = t;
            m_in_ar(m_left) = generate_group_size_exp();
            event(1) = t + exprnd(1/lam_1);
    elseif index == 2
        t = event(2);     
        n_departure = n_departure + 1;
        n_cus = n_cus - 1;
        de_t(n_departure) = t;
        time_sent = t - ar_t(n_arrived - n_cus);
        soj_t(n_departure) = time_sent;
        if n_cus > 0
            event(2) = t + exprnd(n_in_ar(n_arrived - n_cus) / mu_1);
        else
            event(2) = inf; % No departure if no customers left
        end
        if no_of_chefs > 1
        for i=1:(no_of_chefs - 1)
            if n_cus > 0
                t = event(2);     
                n_departure = n_departure + 1;
                n_cus = n_cus - 1;
                de_t(n_departure) = t;
                time_sent = t - ar_t(n_arrived - n_cus);
                soj_t(n_departure) = time_sent;
                if n_cus > 0
                    event(2) = t + exprnd(n_in_ar(n_arrived - n_cus) / mu_1);
                else
                    event(2) = inf; % No departure if no customers left
                end
            end
        end
        end

    else
        measure = measure + 1;
        m_value(measure) = sum(n_in_ar(max(1, n_arrived - n_cus + 1):n_arrived));
        m_value_without(measure) = n_cus + 1;
        event(3) = event(3) + exprnd(1 / step_s); 
    end
end

event(1) = inf;
while n_cus > 0 
    [value, index] = min(event);
    if index == 2 
        t = event(2);   
        n_departure = n_departure + 1;
        n_cus = n_cus - 1;
        de_t(n_departure) = t;
        time_sent = t - ar_t(n_arrived - n_cus);
        soj_t(n_departure) = time_sent;
        if n_cus > 0
            event(2) = t + exprnd(n_in_ar(n_arrived - n_cus) / mu_1);  %change the mu according to the no of people arrived 
        else
            event(2) = inf; % No departure if no customers left
        end
        if no_of_chefs > 1
        for i=1:(no_of_chefs-1)
            if n_cus > 0
                n_departure = n_departure + 1;
                n_cus = n_cus - 1;
                de_t(n_departure) = t;
                time_sent = t - ar_t(n_arrived - n_cus);
                soj_t(n_departure) = time_sent;
                event(2) = t + exprnd(n_in_ar(n_arrived - n_cus) / mu_1);  %change the mu according to the no of people arrived 
            else
                event(2) = inf; % No departure if no customers left
            end
        end
        end

    else
        measure = measure + 1;
        m_value(measure) = sum(n_in_ar(max(1, n_arrived - n_cus + 1):n_arrived));
        event(3) = event(3) + exprnd(1 / step_s);
    end
end


tp = find(n_in_ar == 1);
[bins_s1,freq_s1] = acNhist2(soj_t(tp),0);
tp = find(n_in_ar == 2);
[bins_s2,freq_s2] = acNhist2(soj_t(tp),0);
tp = find(n_in_ar == 3);
[bins_s3,freq_s3] = acNhist2(soj_t(tp),0);
tp = find(n_in_ar == 4);
[bins_s4,freq_s4] = acNhist2(soj_t(tp),0);
tp = find(n_in_ar == 5);
[bins_s5,freq_s5] = acNhist2(soj_t(tp),0);
tp = find(n_in_ar == 6);
[bins_s6,freq_s6] = acNhist2(soj_t(tp),0);
[bins_s,freq_s] = acNhist2(soj_t,0);

figure;
midpoints = 1:7;  % Define midpoints
edges = [midpoints-0.5, midpoints(end)+0.5];  % Calculate edges from midpoints
subplot(3,2,1);
histogram(n_in_ar, 'BinEdges', edges,"Normalization","pdf");
xlabel("family size");
ylabel("probabilty");


subplot(3,2,2);
plot(bins_s1, freq_s1, 'r', 'DisplayName', 'Group Size 1','LineWidth',2); % Red for group size 1
hold on;
plot(bins_s2, freq_s2, 'g', 'DisplayName', 'Group Size 2','LineWidth',2); % Green for group size 2
plot(bins_s3, freq_s3, 'b', 'DisplayName', 'Group Size 3','LineWidth',2); % Blue for group size 3
plot(bins_s4, freq_s4, 'm', 'DisplayName', 'Group Size 4','LineWidth',2); % Magenta for group size 4
plot(bins_s5, freq_s5, 'c', 'DisplayName', 'Group Size 5','LineWidth',2); % Cyan for group size 5
plot(bins_s6, freq_s6, 'k', 'DisplayName', 'Group Size 6','LineWidth',2); % Black for group size 6
plot(bins_s, freq_s, 'y', 'DisplayName', 'All Group Sizes','LineWidth',2); % Yellow dashed line for all group sizes
xlabel("sojour time");
ylabel("probabilty");
legend("s1","s2","s3","s4","s5","6","all");

[bins_m_value, freq_m_value] = acNhist2(m_value,0);
[bins_m_value_without,freq_m_value_without] = acNhist2(m_value_without,0);

subplot(3,2,3);
%plot(bins_m_value, freq_m_value, 'r', 'DisplayName', 'Group Size 1','LineWidth',5); 
plot(bins_m_value_without,freq_m_value_without, 'g', 'DisplayName', 'Group Size 2','LineWidth',1); 

xlabel("queue length");
ylabel("probabitly");

[bins_mt,fre_mt] = acNhist2(m_t,0);
subplot(3,2,4);
plot(ar_t,'g','LineWidth',5);
hold on
plot(de_t,"r",'LineWidth',2);
ylabel("arrival or depature time");


legend("arrival time","depature time","coustomer left time without taking service");
subplot(3,2,5);
[bin_left,freq_left] = acNhist2(m_in_ar,0);
plot(bin_left,freq_left,'k');
xlabel("family size left");
ylabel("probabilty");


subplot(3,2,6);
plot(bins_mt,fre_mt,'k');
xlabel("Time at which people left");
ylabel("probabilty");


disp("mean is ");
disp(mean(n_in_ar/10)*lam_1);

% ------ Local function defined at the end of the script ------
function group_size = generate_group_size_exp()
    exp_sample = exprnd(1); 
    max_value = -log(1 - 0.95); 
    scaled_value = 1 + (exp_sample / max_value) * (7 - 1); 
    group_size = round(min(max(1, scaled_value), 7)); 
end

% Plot four-panel figure showing methane breakdown trajectories

emis = [-57, -284, 2.8, 5];   % estimated emitted composition Whitehill, CH2D2 is guess
atm = [-47.1, -86, 3.0, 5-1000*log(0.9266)];   % atmosheric (Rigby et al 2012,
                                    % and predicted, Whitehill et al 2017)
                                    % and predicted, Wang et al 2017 eqn 21
                                    %   assuming DD behaves the same and OH
                                    %   only

wo.ac = csvread('./whiticaroutlines/acetoclastic.csv');
wo.cr = csvread('./whiticaroutlines/CO2-reduction.csv');
wo.th = csvread('./whiticaroutlines/thermogenic.csv');
wo.ms = csvread('./whiticaroutlines/main_side.csv');
wo.mt = csvread('./whiticaroutlines/main_top2.csv');

wos = struct2cell(wo);
wof = fieldnames(wo);

%% Calculate 13D and DD clumping at equilibrium

Tc = [0:10:1000]';
Tk = Tc+273.15;

K13D = 1 + 0.0355502./Tk - 433.038./Tk.^2 + 1270210.0./Tk.^3 - 5.94804e8./Tk.^4 ...
    + 1.196630e11./Tk.^5 - 9.07230e12./Tk.^6;
KDD = (6/16) * (1 + 0.183798./Tk - 785.483./Tk.^2 +1056280.0./Tk.^3 + 9.37307e7./Tk.^4 ...
    - 8.919480e10./Tk.^5 + 9.901730e12./Tk.^6);     % Young et al., 2017

D64eq = 1000*(K13D - 1);
D65eq = 1000*(KDD*16/6 - 1);


% D/H fractionations
load('./alphas from 0 to 1000 C/alphas.mat')
Tc = Tk-273.15;

amliq = bsxfun(@times, a(:,5), bsxfun(@times, a(:,[2,3,6]).^-1, a(:,7).^-1));
emliq = 1000*(amliq - 1);  %  C, S , BW

% Equilibrium fractionations
eqs = csvread('./alphas from 0 to 1000 C/Equils.csv',3);
%1/T	T / K	T / C	CO2g/CH4	H2Ol/CH4	CH4/stoch	bicarb/CH4	CH4/CO2g	CH4/H2Ol	CH4/stoch	CH4/bicarb

Tic = [0:10:1000]';  % T to interpolate
C13i = interp1(eqs(:,3), eqs(:,[8,11]), Tic);
D13Di = interp1(eqs(:,3),eqs(:,10),Tic);


%% Calculate Clumping for Different Sinks

% Conditions
% 1-5:  diffusion, AeOM, AOM, OH, Cl

init = [C13i(Tic==20,1), emliq(Tc==20,1), D64eq(Tc==20), D65eq(Tc==20)];
init = emis;

figure(2); clf
hts = tight_subplot(2, 2, [0.02 0.02], [0.08], [0.08])


%% TOP LEFT 
axes(hts(1))

plot(emliq(Tc<=370,1), D65eq(Tc<=370), '-', 'Color', colors('Black'), 'LineWidth', 1); hold on;
labelat = [0, 100, 200, 300];
markat = [0:20:360];
plot(emliq(ismember(Tc,markat),1), repmat(D65eq(ismember(Tc,markat)),1,1), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('Black'), ...
    'Color', colors('Black'), 'LineWidth', 0.5); hold on;
text(emliq(ismember(Tc,labelat),1)+10, repmat(D65eq(ismember(Tc,labelat)),1,1), ...
    strcat(strread(num2str(labelat),'%s'), repmat({' ｰC'},length(labelat),1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7, 'Color', colors('Black'));


plot([-500 500], [0 0], ':', 'Color', colors('Dark Gray')); hold on;
plot([0 0], [-500 500], ':', 'Color', colors('Dark Gray')); hold on;

out = RunMoxModel5(1,init);   % diffusion
plot(out(:,3), out(:,5), 'k:'); hold on;
labelat = [];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),3), out(ismember(out(:,1),markat),5), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Black'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),3)-0.1, out(ismember(out(:,1),labelat),5), ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7,...
    'Color', colors('Black'));


out = RunMoxModel5(2,init);   % AeOM
plot(out(:,3), out(:,5), '-', 'Color', colors('Crimson')); hold on;
labelat = [70 60 50 40];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),3), out(ismember(out(:,1),markat),5), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Crimson'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),3)+10, out(ismember(out(:,1),labelat),5), ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7, ...
    'Color', colors('Crimson'));


out = RunMoxModel5(3,init);   % AOM
plot(out(:,3), out(:,5), ':', 'Color', colors('Copper')); hold on;
labelat = [70 60 50 40];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),3), out(ismember(out(:,1),markat),5), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Copper'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),3)+5, out(ismember(out(:,1),labelat),5), ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 7, ...
    'Color', colors('Copper'));


out = RunMoxModel5(4,init);   % OH
plot(out(:,3), out(:,5), '-', 'Color', colors('Tufts Blue')); hold on;
labelat = [100:-20:40];  % specifiy individual values because for some reason
                            % in R2012b there is a bug where 
                            % [1 0.8 0.6 0.4] == [1:-0.2:0.4] is not true!
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),3), out(ismember(out(:,1),markat),5), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Tufts Blue'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),3)+5, out(ismember(out(:,1),labelat),5), ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 7, ...
    'Color', colors('Tufts Blue'));


out = RunMoxModel5(5,init);   % Cl
plot(out(:,3), out(:,5), ':', 'Color', colors('Tufts Blue')); hold on;
labelat = [100 90 80 70];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),3), out(ismember(out(:,1),markat),5), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Tufts Blue'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),3)-2, out(ismember(out(:,1),labelat),5)+0.1, ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 7, ...
    'Color', colors('Tufts Blue'));

% atmospheric - predicted
plot(atm(3), atm(4), 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Tufts Blue'), 'LineWidth', 0.75); hold on;
plot(atm(3), atm(4), 'pentagram', ...
    'MarkerSize', 7, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Tufts Blue'), 'LineWidth', 0.5); hold on;


ylim([-40,+40])
xlim([-450,-50])

set(gca(),'TickLength',3*get(gca(),'TickLength'))
set(gca(),'XMinorTick','on','YMinorTick','on')
set(gca(),'XAxisLocation','top');
set(gca(),'YAxisLocation','left');

ylabel('{\Delta}^{12}CH_2D_2 [云')
xlabel('{\delta}D [云')

xl=xlim
yl=ylim  %panel label
text(xl(1)+0.05*diff(xl),diff(yl)*0.95+yl(1),'A', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

axis square




%% BOTTOM LEFT 
axes(hts(3))

for i = 1:length(wof)
    pts = wos{i}; ls = ':'
    if ~any(ismember(wof{i}, {'ms', 'mt'})), pts = [pts; pts(1,:)]; ls = 'none'; end     % close the loop
    plot(pts(:,1), pts(:,2), 'LineStyle', ls, 'Color', colors('Gray'), 'LineWidth', 1); hold on;
end

plot([-500 500], [0 0], ':', 'Color', colors('Dark Gray')); hold on;
plot([0 0], [-500 500], ':', 'Color', colors('Dark Gray')); hold on;

% % equilibrium fractioations
% plot(C13i(Tc<=370,2), repmat(D13Di(Tc<=370,1),1,1), 'Color', colors('Burnt Orange'), 'LineWidth', 1); hold on;
% markat = [0:20:400];
% plot(C13i(ismember(Tc,markat),2), repmat(D13Di(ismember(Tc,markat),1),1,1), 'o', ...
%     'MarkerSize', 2.5, 'MarkerFaceColor', colors('Burnt Orange'), ...
%     'Color', colors('Burnt Orange'), 'LineWidth', 0.5)

plot(repmat(emliq(:,1),1,1), C13i(:,1), 'Color', colors('Black'), 'LineWidth', 1)
labelat = [0, 100, 200, 300, 500];
markat = [0:20:1000];
plot(repmat(emliq(ismember(Tc,markat)),1,1), C13i(ismember(Tc,markat),1), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('Black'), ...
    'Color', colors('Black'), 'LineWidth', 0.5)
text(repmat(emliq(ismember(Tc,labelat)),1,1)+5, C13i(ismember(Tc,labelat),1), ...
    strcat(strread(num2str(labelat),'%s'), repmat({' ｰC'},length(labelat),1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 7, ...
    'Color', colors('Black'));

% kinetic fractionations
out = RunMoxModel5(1,init);   % diffusion
plot(out(:,3), out(:,2), 'k:'); hold on;
labelat = [60 40 20];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),3), out(ismember(out(:,1),markat),2), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Black'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),3)-5, out(ismember(out(:,1),labelat),2), ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 7,...
    'Color', colors('Black'));


out = RunMoxModel5(2,init);   % AeOM
plot(out(:,3), out(:,2), '-', 'Color', colors('Crimson')); hold on;
labelat = [60 50];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),3), out(ismember(out(:,1),markat),2), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Crimson'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),3), out(ismember(out(:,1),labelat),2)+5, ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7, ...
    'Color', colors('Crimson'));


out = RunMoxModel5(3,init);   % AOM
plot(out(:,3), out(:,2), ':', 'Color', colors('Copper')); hold on;
labelat = [60 50];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),3), out(ismember(out(:,1),markat),2), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Copper'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),3), out(ismember(out(:,1),labelat),2)-5, ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 7, ...
    'Color', colors('Copper'));


out = RunMoxModel5(4,init);   % OH
plot(out(:,3), out(:,2), '-', 'Color', colors('Tufts Blue')); hold on;
labelat = [100  80  60];  % specifiy individual values because for some reason
                            % in R2012b there is a bug where 
                            % [1 0.8 0.6 0.4] == [1:-0.2:0.4] is not true!
                            % ..rather this is just a precisioin issue when
                            % defining non-integral intervals
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),3), out(ismember(out(:,1),markat),2), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Tufts Blue'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),3), out(ismember(out(:,1),labelat),2)-1, ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 7, ...
    'Color', colors('Tufts Blue'));


out = RunMoxModel5(5,init);   % Cl
plot(out(:,3), out(:,2), ':', 'Color', colors('Tufts Blue')); hold on;
labelat = [80 70 60];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),3), out(ismember(out(:,1),markat),2), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Tufts Blue'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),3), out(ismember(out(:,1),labelat),2), ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 7, ...
    'Color', colors('Tufts Blue'));

% atmospheric 
plot(atm(2), atm(1), 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Tufts Blue'), 'LineWidth', 0.75); hold on;
plot(atm(2), atm(1), 'pentagram', ...
    'MarkerSize', 7, 'MarkerFaceColor', colors('Tufts Blue'), ...
    'Color', colors('Tufts Blue'), 'LineWidth', 0.5); hold on;

xlim([-450,-50])
ylim([-115,-5])

set(gca(),'TickLength',3*get(gca(),'TickLength'))
set(gca(),'XMinorTick','on','YMinorTick','on')
set(gca(),'XAxisLocation','bottom');
set(gca(),'YAxisLocation','left');

set(gca(),'YDir', 'reverse')

ylabel('{\delta}^{13}C [云')
xlabel('{\delta}D [云')

xl=xlim
yl=ylim  %panel label
text(xl(1)+0.05*diff(xl),diff(yl)*0.05+yl(1),'C', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

axis square

%% TOP RIGHT 
axes(hts(2))

plot(D64eq, D65eq, '-', 'Color', colors('Black'), 'LineWidth', 1); hold on;
labelat = [0, 100, 200, 300];
markat = [0:20:1000];
plot(D64eq(ismember(Tc,markat),1), repmat(D65eq(ismember(Tc,markat)),1,1), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('Black'), ...
    'Color', colors('Black'), 'LineWidth', 0.5); hold on;
text(D64eq(ismember(Tc,labelat),1)-0.1, repmat(D65eq(ismember(Tc,labelat)),1,1), ...
    strcat(strread(num2str(labelat),'%s'), repmat({' ｰC'},length(labelat),1)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 7, 'Color', colors('Black'));


plot([-500 500], [0 0], ':', 'Color', colors('Dark Gray')); hold on;
plot([0 0], [-500 500], ':', 'Color', colors('Dark Gray')); hold on;

out = RunMoxModel5(1,init);   % diffusion
plot(out(:,4), out(:,5), 'k:'); hold on;
labelat = [];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),4), out(ismember(out(:,1),markat),5), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Black'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),4), out(ismember(out(:,1),labelat),5)+1, ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 7,...
    'Color', colors('Black'));


out = RunMoxModel5(2,init);   % AeOM
plot(out(:,4), out(:,5), '-', 'Color', colors('Crimson')); hold on;
labelat = [60 50];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),4), out(ismember(out(:,1),markat),5), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Crimson'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),4)+0.3, out(ismember(out(:,1),labelat),5), ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7, ...
    'Color', colors('Crimson'));


out = RunMoxModel5(3,init);   % AOM
plot(out(:,4), out(:,5), ':', 'Color', colors('Copper')); hold on;
labelat = [60 50];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),4), out(ismember(out(:,1),markat),5), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Copper'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),4)-0.2, out(ismember(out(:,1),labelat),5), ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 7, ...
    'Color', colors('Copper'));


out = RunMoxModel5(4,init);   % OH
plot(out(:,4), out(:,5), '-', 'Color', colors('Tufts Blue')); hold on;
labelat = [100:-20:40];  % specifiy individual values because for some reason
                            % in R2012b there is a bug where 
                            % [1 0.8 0.6 0.4] == [1:-0.2:0.4] is not true!
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),4), out(ismember(out(:,1),markat),5), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Tufts Blue'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),4)+0.3, out(ismember(out(:,1),labelat),5), ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7, ...
    'Color', colors('Tufts Blue'));


out = RunMoxModel5(5,init);   % Cl
plot(out(:,4), out(:,5), ':', 'Color', colors('Tufts Blue')); hold on;
labelat = [80 70];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),4), out(ismember(out(:,1),markat),5), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Tufts Blue'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),4), out(ismember(out(:,1),labelat),5)+0.1, ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 7, ...
    'Color', colors('Tufts Blue'));

% atmospheric - predicted
% plot(atm(3), atm(4), 'o', ...
%     'MarkerSize', 8, 'MarkerFaceColor', colors('White'), ...
%     'Color', colors('Tufts Blue'), 'LineWidth', 0.75); hold on;
plot(atm(3), atm(4), 'pentagram', ...
    'MarkerSize', 7, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Tufts Blue'), 'LineWidth', 0.5); hold on;

ylim([-40,+40])
xlim([-11,+11])

set(gca(),'TickLength',3*get(gca(),'TickLength'))
set(gca(),'XMinorTick','on','YMinorTick','on')
set(gca(),'XAxisLocation','top');
set(gca(),'YAxisLocation','right');

ylabel('{\Delta}^{12}CH_2D_2 [云')
xlabel('{\Delta}^{13}CH_3D [云')

xl=xlim
yl=ylim  %panel label
text(xl(1)+0.05*diff(xl),diff(yl)*0.95+yl(1),'B', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

axis square


%% BOTTOM RIGHT 
axes(hts(4))


plot([-200 100], [0 0], ':', 'Color', colors('Dark Gray')); hold on;
plot([0 0], [-500 500], ':', 'Color', colors('Dark Gray')); hold on;

% % equilibrium fractioations
% plot(C13i(Tc<=370,2), repmat(D13Di(Tc<=370,1),1,1), 'Color', colors('Burnt Orange'), 'LineWidth', 1); hold on;
% markat = [0:20:400];
% plot(C13i(ismember(Tc,markat),2), repmat(D13Di(ismember(Tc,markat),1),1,1), 'o', ...
%     'MarkerSize', 2.5, 'MarkerFaceColor', colors('Burnt Orange'), ...
%     'Color', colors('Burnt Orange'), 'LineWidth', 0.5)

plot(repmat(D13Di(:,1),1,1), C13i(:,1), 'Color', colors('Black'), 'LineWidth', 1)
labelat = [0, 100, 200, 300, 500];
markat = [0:20:1000];
plot(repmat(D13Di(ismember(Tc,markat)),1,1), C13i(ismember(Tc,markat),1), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('Black'), ...
    'Color', colors('Black'), 'LineWidth', 0.5)
text(repmat(D13Di(ismember(Tc,labelat)),1,1)+0.1, C13i(ismember(Tc,labelat),1), ...
    strcat(strread(num2str(labelat),'%s'), repmat({' ｰC'},length(labelat),1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 7, ...
    'Color', colors('Black'));

% kinetic fractionations
out = RunMoxModel5(1,init);   % diffusion
plot(out(:,4), out(:,2), 'k:'); hold on;
labelat = [ 20 10];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),4), out(ismember(out(:,1),markat),2), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Black'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),4), out(ismember(out(:,1),labelat),2), ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 7,...
    'Color', colors('Black'));


out = RunMoxModel5(2,init);   % AeOM
plot(out(:,4), out(:,2), '-', 'Color', colors('Crimson')); hold on;
labelat = [20 10];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),4), out(ismember(out(:,1),markat),2), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Crimson'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),4), out(ismember(out(:,1),labelat),2), ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 7, ...
    'Color', colors('Crimson'));


out = RunMoxModel5(3,init);   % AOM
plot(out(:,4), out(:,2), ':', 'Color', colors('Copper')); hold on;
labelat = [20 10];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),4), out(ismember(out(:,1),markat),2), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Copper'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),4), out(ismember(out(:,1),labelat),2)+0.1, ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 7, ...
    'Color', colors('Copper'));


out = RunMoxModel5(4,init);   % OH
plot(out(:,4), out(:,2), '-', 'Color', colors('Tufts Blue')); hold on;
labelat = [20 10];  % specifiy individual values because for some reason
                            % in R2012b there is a bug where 
                            % [1 0.8 0.6 0.4] == [1:-0.2:0.4] is not true!
                            % ..rather this is just a precisioin issue when
                            % defining non-integral intervals
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),4), out(ismember(out(:,1),markat),2), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Tufts Blue'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),4), out(ismember(out(:,1),labelat),2), ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 7, ...
    'Color', colors('Tufts Blue'));


out = RunMoxModel5(5,init);   % Cl
plot(out(:,4), out(:,2), ':', 'Color', colors('Tufts Blue')); hold on;
labelat = [100:-10:60];
markat = [100:-10:10];
plot(out(ismember(out(:,1),markat),4), out(ismember(out(:,1),markat),2), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Tufts Blue'), 'LineWidth', 0.5); hold on;
text(out(ismember(out(:,1),labelat),4), out(ismember(out(:,1),labelat),2)-0.1, ...
    strcat(strread(num2str(labelat),'%s'), repmat({'%'},length(labelat),1)), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 7, ...
    'Color', colors('Tufts Blue'));

% atmospheric - predicted
plot(atm(3), atm(1), 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Tufts Blue'), 'LineWidth', 0.75); hold on;
plot(atm(3), atm(1), 'pentagram', ...
    'MarkerSize', 7, 'MarkerFaceColor', colors('White'), ...
    'Color', colors('Tufts Blue'), 'LineWidth', 0.5); hold on;

ylim([-115,-5])
xlim([-11,+11])

set(gca(),'TickLength',3*get(gca(),'TickLength'))
set(gca(),'XMinorTick','on','YMinorTick','on')
set(gca(),'XAxisLocation','bottom');
set(gca(),'YAxisLocation','right');

set(gca(),'YDir', 'reverse')

ylabel('{\delta}^{13}C [云')
xlabel('{\Delta}^{13}CH_3D [云')

xl=xlim
yl=ylim  %panel label
text(xl(1)+0.05*diff(xl),diff(yl)*0.05+yl(1),'D', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

axis square



%% save
set(gcf, 'Position', [  698    59   615   615]);
set(gcf,'PaperPositionMode','auto')
print(gcf(), '-depsc2', [mfilename '.eps']);



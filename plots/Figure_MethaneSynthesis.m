% Plot two panel figure showing methane synthesis trajectories

wo.ac = csvread('./whiticaroutlines/acetoclastic.csv');
wo.cr = csvread('./whiticaroutlines/CO2-reduction.csv');
wo.th = csvread('./whiticaroutlines/thermogenic.csv');
wo.ms = csvread('./whiticaroutlines/main_side.csv');
wo.mt = csvread('./whiticaroutlines/main_top2.csv');

wos = struct2cell(wo);
wof = fieldnames(wo);

figure(1); clf
hts = tight_subplot(2, 1, [0.02 0.02], [0.11 0.03], [0.11 0.03])


%% Main Panel (Whiticar)

axes(hts(2))

for i = 1:length(wof)
    pts = wos{i}; ls = ':'
    if ~any(ismember(wof{i}, {'ms', 'mt'})), pts = [pts; pts(1,:)]; ls = 'none'; end     % close the loop
    plot(pts(:,1), pts(:,2), 'LineStyle', ls, 'Color', colors('Gray'), 'LineWidth', 1); hold on;
end

for it = [0, 20, 40]
    load(['./Rees_out/D18_' sprintf('%0g',it) 'C.mat'])
    a2 = d13C5stor;
    a3 = dD5stor;
    a4 = D185stor;
    x = a3;
    y = a2;
    if it == 0
        plot(x(:,1), y(:,1), 'Color', colors('Dark Gray')); hold on;
        plot(x(:,1)', y(:,1)', 'Color', colors('Dark Gray')); hold on;
    elseif it == 20
        plot(x, y, '-', 'Color', colors('Light Gray')); hold on;
        plot(x', y', ':', 'Color', colors('Light Gray')); hold on;
        plot(x(1,:)', y(1,:)', '-', 'Color', colors('Dark Gray')); hold on;
    elseif it == 40
        plot(x(:,end), y(:,end), 'Color', colors('Dark Gray')); hold on;
        plot(x(:,end)', y(:,end)', 'Color', colors('Dark Gray')); hold on;  
    end
end

% D/H fractionations
load('./alphas from 0 to 1000 C/alphas.mat')
Tc = Tk-273.15;

amliq = bsxfun(@times, a(:,5), bsxfun(@times, a(:,[2,3,6]).^-1, a(:,7).^-1));
emliq = 1000*(amliq - 1);

% Equilibrium fractionations
eqs = csvread('./alphas from 0 to 1000 C/Equils.csv',3);
%1/T	T / K	T / C	CO2g/CH4	H2Ol/CH4	CH4/stoch	bicarb/CH4	CH4/CO2g	CH4/H2Ol	CH4/stoch	CH4/bicarb

Tic = [0:10:1000]';  % T to interpolate
C13i = interp1(eqs(:,3), eqs(:,[8,11]), Tic)
D13Di = interp1(eqs(:,3),eqs(:,10),Tic)

plot(emliq(Tc<=370,1), repmat(C13i(Tc<=370,2),1,1), 'Color', colors('Burnt Orange'), 'LineWidth', 1)
% plot(emliq(Tc<=370,2:3), repmat(C13i(Tc<=370,2),1,2), 'Color', colors('Burnt Orange'), 'LineWidth', 0.5)
markat = [0:20:400];
plot(emliq(ismember(Tc,markat),1), repmat(C13i(ismember(Tc,markat),2),1,1), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('Burnt Orange'), ...
    'Color', colors('Burnt Orange'), 'LineWidth', 0.5)

plot(emliq(Tc<=370,1), repmat(C13i(Tc<=370,1),1,1), 'Color', colors('Black'), 'LineWidth', 1)
% plot(emliq(Tc<=370,2:3), repmat(C13i(Tc<=370,1),1,2), 'Color', colors('Black'), 'LineWidth', 0.5)
labelat = [0, 40, 100, 200, 300, 500];
markat = [0:20:400];
plot(emliq(ismember(Tc,markat),1), repmat(C13i(ismember(Tc,markat)),1,1), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('Black'), ...
    'Color', colors('Black'), 'LineWidth', 0.5)
text(emliq(ismember(Tc,labelat),1)+10, repmat(C13i(ismember(Tc,labelat)),1,1), ...
    strcat(strread(num2str(labelat),'%s'), repmat({' ｰC'},length(labelat),1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7, 'Color', colors('Black'));

xlim([-450,-50])
ylim([-115,-5])

set(gca(),'YDir', 'reverse')

set(gca(),'TickLength',3*get(gca(),'TickLength'))
set(gca(),'XMinorTick','on','YMinorTick','on')

% xlabel('{\delta}D_{CH_4} - {\delta}D_{H_2O} [云')
% ylabel('{\delta}^{13}C_{CH_4} - {\delta}^{13}C_{CO_2 or DIC} [云')

xlabel('{\delta}D [云')
ylabel('{\delta}^{13}C [云')

xl=xlim
yl=ylim  %panel label
text(xl(1)+0.05*diff(xl),diff(yl)*0.05+yl(1),'B', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

axis square

%% Top Panel (Delta vs dD)

axes(hts(1))

plot([-450 -50], [0 0], 'k:'); hold on;

for it = [0, 20, 40]
    load(['./Rees_out/D18_' sprintf('%0g',it) 'C.mat'])
    a2 = d13C5stor;
    a3 = dD5stor;
    a4 = D185stor;
    x = a3;
    y = a4;
    if it == 0
        plot(x(:,1), y(:,1), 'Color', colors('Dark Gray')); hold on;
        plot(x(:,1)', y(:,1)', 'Color', colors('Dark Gray')); hold on;
    elseif it == 20
        plot(x, y, '-', 'Color', colors('Light Gray')); hold on;
        plot(x', y', ':', 'Color', colors('Light Gray')); hold on;
        plot(x(1,:)', y(1,:)', '-', 'Color', colors('Dark Gray')); hold on;
    elseif it == 40
        plot(x(:,end), y(:,end), 'Color', colors('Dark Gray')); hold on;
        plot(x(:,end)', y(:,end)', 'Color', colors('Dark Gray')); hold on;  
    end
end

plot(emliq(Tc<=370,1), repmat(D13Di(Tc<=370),1,1), 'Color', colors('Black'), 'LineWidth', 1)
plot(emliq(Tc<=370,2:3), repmat(D13Di(Tc<=370),1,2), 'Color', colors('Black'), 'LineWidth', 0.5)
labelat = [0, 40, 100, 200, 300, 500];
markat = [0:20:400];
plot(emliq(ismember(Tc,markat),1), repmat(D13Di(ismember(Tc,markat)),1,1), 'o', ...
    'MarkerSize', 2.5, 'MarkerFaceColor', colors('Black'), ...
    'Color', colors('Black'), 'LineWidth', 0.5)
text(emliq(ismember(Tc,labelat),1)+10, repmat(D13Di(ismember(Tc,labelat)),1,1), ...
    strcat(strread(num2str(labelat),'%s'), repmat({' ｰC'},length(labelat),1)), ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7, 'Color', colors('Black'));



xlim([-450,-50])
ylim([-7,+9])

set(gca(),'XAxisLocation','top');
set(gca(),'TickLength',3*get(gca(),'TickLength'))
set(gca(),'XMinorTick','on','YMinorTick','on')

xlabel('{\delta}D [云')
ylabel('{\Delta}^{13}CH_3D [云')

xl=xlim
yl=ylim  %panel label
text(xl(1)+0.05*diff(xl),diff(yl)*0.95+yl(1),'A', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

axis square

%% Save/Clean Up

set(gcf(), 'Position', [587    59   376   607]);
print(gcf(), '-depsc2', [mfilename '.eps']);

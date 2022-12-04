%% Load Catera's experimental data
load('data.mat')

% add path for functions and results
dir_opt = '../Optimization'
addpath(dir_opt)

%% Load klingmuller model to get original parameter values
model_ifn = sbmlimport('modelev1.xml');

param_org = cell2mat(get(sbioselect(model_ifn, 'Type', 'parameter'), 'Value'));
param_names = get(sbioselect(model_ifn, 'Type', 'parameter'), 'Name');
param_org([65:66, 69:72, 75, 78]) = []; % remove duplicated entries for the delays
param_names([65:66, 69:72, 75, 78]) = []; % remove duplicated entries for the delays

% add parameter for active recepter complex formation for SOCS3
param_names = [param_names(1:6); 'kinhBySOCS3'; param_names(7:end)];
param_org = [param_org(1:6); 1; param_org(7:end)];

% rename kinhBySOCS
param_names{strcmp('kinhBySOCS', param_names)} = 'kinhBySOCS1';

% uncouple SOCS1mRNA degradation by IRF2
param_names = [param_names(1:58); 'hlSOCS1mRNAByIRF2'; param_names(59:end)];
param_org = [param_org(1:58); param_org(58); param_org(59:end)];

% uncouple IRF2mRNA degradation and IRF2 protein synthesis
param_names = [param_names; 'hlIRF2mRNA'];
param_org = [param_org; param_org(find(strcmp('delayIRF2', param_names)))];

% remove degByUSP18 which is not used based on the paper supplementary
param_org(find(strcmp('degByUSP18', param_names))) = [];
param_names(find(strcmp('degByUSP18', param_names))) = [];

%% get max time for simulation
tmp = cellfun(@(x) data_beta.(x).t_exp, fieldnames(data_beta), 'UniformOutput',false)
t_max = max(vertcat(tmp{:}));

%% Set doses
IFNlambdadose=1.1e7;
IFNbetadose=5e4;

free = 1:74; 
load([dir_opt, '/Results/3_free/PSO15/PSO_60particles_n15_final.mat']);
sse_3 = res.fval;
chosen_3 = free(floor(res.x(1:3))); 
param_values_3 = res.x(4:end);
free = setdiff(1:74,chosen_3);

load([dir_opt, '/Results/6_free/PSO10/PSO_90particles_n10_final.mat']);
sse_6 = res.fval;
chosen_6 = [chosen_3, free(floor(res.x(1:3)))]; 
param_values_6 = res.x(4:end);
free = setdiff(1:74,chosen_6);

load([dir_opt, '/Results/9_free/PSO6/PSO_120particles_n6_final.mat']);
sse_9 = res.fval;
chosen_9 = [chosen_6, free(floor(res.x(1:3)))]; 
param_values_9 = res.x(4:end);
free = setdiff(1:74,chosen_9);

load([dir_opt, '/Results/12_free/PSO10/PSO_150particles_n10_final.mat']);
sse_12 = res.fval;
chosen_12 = [chosen_9, free(floor(res.x(1:3)))]; 
param_values_12 = res.x(4:end);
free = setdiff(1:74,chosen_12);

load([dir_opt, '/Results/15_free/PSO14/PSO_180particles_n14_final.mat']);
sse_15 = res.fval;
chosen_15 = [chosen_12, free(floor(res.x(1:3)))]; 
param_values_15 = res.x(4:end);
free = setdiff(1:74,chosen_15);

load([dir_opt, '/Results/18_free/PSO5/PSO_210particles_n5_final.mat']);
sse_18 = res.fval;
chosen_18 = [chosen_15, free(floor(res.x(1:3)))]; 
param_values_18 = res.x(4:end);
free = setdiff(1:74,chosen_18);

load([dir_opt, '/Results/21_free/PSO23/PSO_480particles_n23_final.mat']);
sse_21 = res.fval;
chosen_21 = [chosen_18, free(floor(res.x(1:3)))]; 
param_values_21 = res.x(4:end);
free = setdiff(1:74,chosen_21);

chosen_all = {chosen_3, chosen_6, chosen_9, chosen_12, chosen_15, chosen_18, chosen_21};
param_values_all = {param_values_3, param_values_6, param_values_9, param_values_12, param_values_15, param_values_18,param_values_21};

out = solveODE_withSS (param_org, param_names, [], {}, repmat(0,41,1), 0:t_max, 0);
y0 = out.y0;

outs = [];
outs_names=[];
SSE_all = [];
SSE_total = [];
for i = 1:7    
    new_param_values = @(x) [log10(param_org(chosen_all{i})) + x'; log10(param_org(setdiff(1:length(param_org), chosen_all{i})))] ;
    new_param_names = {param_names{chosen_all{i}}, param_names{setdiff(1:length(param_org), chosen_all{i})}}';
        
    out = solveODE_withSS(10.^new_param_values(param_values_all{i}), new_param_names,[],{}, y0, 0:t_max, IFNbetadose);
    SSE_all = [SSE_all, {SSE_eval(out,data_beta,false, true)}];
    SSE_total = [SSE_total, {SSE_eval(out,data_beta,true, true)}];
    outs = [outs, {out}];
    outs_names = [outs_names, {[num2str(3*i), ' param free']}];
end

% remove the non fitted points i.e. with no replicates for plotting
data_beta_sub = data_beta;

fn = fieldnames(data_beta_sub);
for i = 1:size(fn,1)
    data_beta_sub.(fn{i}).t_exp = data_beta_sub.(fn{i}).t_exp(~isnan(data_beta_sub.(fn{i}).y_exp_sem));
    data_beta_sub.(fn{i}).y_exp = data_beta_sub.(fn{i}).y_exp(~isnan(data_beta_sub.(fn{i}).y_exp_sem));
    data_beta_sub.(fn{i}).y_exp_sem = data_beta_sub.(fn{i}).y_exp_sem(~isnan(data_beta_sub.(fn{i}).y_exp_sem));
end

cm = flipud(parula(length(outs)+1));
cm = cm(2:end,:);
plot_beta_paper(outs, outs_names, data_beta_sub, plot_beta_metadata, ['plots_allfits_',date], 'colors',  cm);
plot_beta_paper(outs(7), outs_names(7), data_beta_sub, plot_beta_metadata, ['plots_bestfit',date], 'colors',  parula(7));

% export table of parameters
table_out = table(param_names,log10(param_org), ...
    repmat(nan, size(param_org)), ...
    repmat(nan, size(param_org)), ...
    repmat(nan, size(param_org)), ...
    repmat(nan, size(param_org)), ...
    repmat(nan, size(param_org)), ...
    repmat(nan, size(param_org)), ...
    repmat(nan, size(param_org)));
    
table_out.Properties.VariableNames = {'Name', 'Value_log10', 'Free_3_x_factor', 'Free_6_x_factor', 'Free_9_x_factor', 'Free_12_x_factor', 'Free_15_x_factor', 'Free_18_x_factor', 'Free_21_x_factor'}
 
table_out.Free_3_x_factor(chosen_3) = param_values_3;
table_out.Free_6_x_factor(chosen_6) = param_values_6;
table_out.Free_9_x_factor(chosen_9) = param_values_9;
table_out.Free_12_x_factor(chosen_12) = param_values_12;
table_out.Free_15_x_factor(chosen_15) = param_values_15;
table_out.Free_18_x_factor(chosen_18) = param_values_18;
table_out.Free_21_x_factor(chosen_21) = param_values_21;

writetable(table_out, ['fits.xlsx'])


% Make promoter mutants
new_param_values = @(x) [log10(param_org(chosen_all{7})) + x'; log10(param_org(setdiff(1:length(param_org), chosen_all{7})))] ;
new_param_names = {param_names{chosen_all{7}}, param_names{setdiff(1:length(param_org), chosen_all{7})}}';
    
param_STAT2mut = new_param_values(param_values_all{7});
find(strcmp('synthSTAT2mRNA', new_param_names))
param_STAT2mut(54) = -Inf;
outs_STAT2mut = solveODE_withSS(10.^param_STAT2mut, new_param_names,[],{}, outs{7}.y0, 0:t_max, IFNbetadose, 'run_ss', false);

param_IRF9mut = new_param_values(param_values_all{7});
find(strcmp('synthIRF9mRNA', new_param_names))
param_IRF9mut(55) = -Inf;
outs_IRF9mut = solveODE_withSS(10.^param_IRF9mut, new_param_names,[],{}, outs{7}.y0, 0:t_max, IFNbetadose, 'run_ss', false);

param_STAT1mut = new_param_values(param_values_all{7});
find(strcmp('synthSTAT1mRNA', new_param_names))
param_STAT1mut(52) = -Inf;
outs_STAT1mut = solveODE_withSS(10.^param_STAT1mut, new_param_names,[],{}, outs{7}.y0, 0:t_max, IFNbetadose, 'run_ss', false);


outs_mut = {outs{7}, outs_STAT2mut, outs_IRF9mut, outs_STAT1mut};
outs_mut_names = {'WT', 'STAT2 promoter mutant', 'IRF9 promoter mutant', 'STAT1 promoter mutant'};

cm=parula(8);
cm=[cm(1,:); 208/256, 70/256, 230/256; 60/256, 200/256, 44/256; 255/256, 164/256, 0/256;];
plot_beta_paper_WTnorm(outs_mut, outs_mut_names, 1, data_beta, plot_beta_metadata, ['plots_prom_mut',date], 'colors',  cm);


% Predict pulse
out_pulse_1 = solveODE_withSS(10.^new_param_values(param_values_all{7}), new_param_names,[],{}, outs{7}.y0, 0:15, IFNbetadose, 'run_ss', false);
out_pulse_2 = solveODE_withSS(10.^new_param_values(param_values_all{7}), new_param_names,[],{}, out_pulse_1.fit_beta.y_sim(end,:), 0:(t_max-15), 0, 'run_ss', false);

out_pulse = out_pulse_1;
out_pulse.fit_beta.t_sim = [out_pulse_1.fit_beta.t_sim;15+out_pulse_2.fit_beta.t_sim(2:end)];
out_pulse.fit_beta.y_sim = [out_pulse_1.fit_beta.y_sim;out_pulse_2.fit_beta.y_sim(2:end,:)];

out_pulse.fit_beta_chx.t_sim = [out_pulse_1.fit_beta_chx.t_sim; 15+out_pulse_2.fit_beta_chx.t_sim(2:end)];
out_pulse.fit_beta_chx.y_sim = [out_pulse_1.fit_beta_chx.y_sim; out_pulse_2.fit_beta_chx.y_sim(2:end,:)];

cm=parula(8);
plot_beta_paper_pulse({outs{7},out_pulse}, {'Sustained IFN', 'IFN pulse'}, 1, data_beta_pulse, plot_beta_pulse_metadata, 'pulse', 'colors',  cm([1,4],:), 'data_name', 'IFN pulse data');

% Predict low beta
out_low = solveODE_withSS(10.^new_param_values(param_values_all{7}), new_param_names,[],{}, outs{7}.y0, 0:max(data_beta_low.ISGF3n_low.t_exp), 0.1 * IFNbetadose, 'run_ss', false);

cm=parula(8);
plot_beta_paper_pulse({outs{7},out_low}, {'High IFN', 'Low IFN'}, 1, data_beta_low, plot_beta_low_metadata, 'low', 'colors',  cm([1,4],:), 'data_name', 'Low IFN data');



% Predict low beta
out_low_1 = solveODE_withSS(10.^new_param_values(param_values_all{7}), new_param_names,[],{}, outs{7}.y0, 0:t_max, 0.01 * IFNbetadose, 'run_ss', false);
out_low_2 = solveODE_withSS(10.^new_param_values(param_values_all{7}), new_param_names,[],{}, outs{7}.y0, 0:t_max, 0.025 * IFNbetadose, 'run_ss', false);
out_low_3 = solveODE_withSS(10.^new_param_values(param_values_all{7}), new_param_names,[],{}, outs{7}.y0, 0:t_max, 0.05 * IFNbetadose, 'run_ss', false);
out_low_4 = solveODE_withSS(10.^new_param_values(param_values_all{7}), new_param_names,[],{}, outs{7}.y0, 0:t_max, 0.1 * IFNbetadose, 'run_ss', false);
out_low_5 = solveODE_withSS(10.^new_param_values(param_values_all{7}), new_param_names,[],{}, outs{7}.y0, 0:t_max, 0.25 * IFNbetadose, 'run_ss', false);
out_low_6 = solveODE_withSS(10.^new_param_values(param_values_all{7}), new_param_names,[],{}, outs{7}.y0, 0:t_max, 0.5 * IFNbetadose, 'run_ss', false);

cm=winter(7);
plot_beta_paper_multidose({outs{7}, out_low_6, out_low_5, out_low_4, out_low_3, out_low_2, out_low_1}, {'High IFN', '0.5x', '0.25x', '0.1x', '0.05x', '0.025x', '0.01x'}, 1, data_beta, plot_beta_metadata, 'multidose', 'colors', cm);

cm=winter(5);
plot_beta_paper_multidose({outs{7}, out_low_6, out_low_5, out_low_4, out_low_3}, {'High IFN', '0.5x', '0.25x', '0.1x', '0.05x'}, 1, data_beta, plot_beta_metadata, 'multidose_200min', 'colors', cm, 'tmax', 200);

ISGF3_bar=[outs{7}.fit_beta.y_sim(outs{7}.fit_beta.t_sim==60,13), outs{7}.fit_beta.y_sim(outs{7}.fit_beta.t_sim==180,13); ...
           out_low_6.fit_beta.y_sim(out_low_6.fit_beta.t_sim == 60,13), out_low_6.fit_beta.y_sim(out_low_6.fit_beta.t_sim == 180,13); ...
           out_low_5.fit_beta.y_sim(out_low_5.fit_beta.t_sim == 60,13), out_low_5.fit_beta.y_sim(out_low_5.fit_beta.t_sim == 180,13); ...
           out_low_4.fit_beta.y_sim(out_low_4.fit_beta.t_sim == 60,13), out_low_4.fit_beta.y_sim(out_low_4.fit_beta.t_sim == 180,13); ...
           out_low_3.fit_beta.y_sim(out_low_3.fit_beta.t_sim == 60,13), out_low_3.fit_beta.y_sim(out_low_3.fit_beta.t_sim == 180,13) ...
] / outs{7}.fit_beta.y_sim(outs{7}.fit_beta.t_sim==20,13);

X = categorical({'High IFN', '0.5x', '0.25x', '0.1x', '0.05x'});
X = reordercats(X,{'High IFN', '0.5x', '0.25x', '0.1x', '0.05x'});

hf = figure(); 
hold on;
% grey bar for legend
b2 = bar(X, ISGF3_bar, 'FaceColor','flat');
b2(1).CData = repmat(0.5,size(cm)) + 0.4*(repmat(1, size(cm))-repmat(0.5,size(cm)));
b2(2).CData = repmat(0.5,size(cm)) + 0.4*(repmat(0, size(cm))-repmat(0.5,size(cm)));
set(b2, {'DisplayName'}, {'1h','3h'}')
% overlay grey bar with color bar
b = bar(X, ISGF3_bar, 'FaceColor','flat');
b(1).CData = cm + 0.4*(repmat(1, size(cm))-cm);
b(2).CData = cm + 0.4*(repmat(0, size(cm))-cm);
ylabel('Normalized ISGF3n')
% show grey legend
lgd=legend(b2);
lgd.Box = 'off';
lgd.FontSize = 14/sqrt(7);
ax = gca(hf); 
ax.Units = 'pixels';
axis_size = [140, 105];
sizeDiff = axis_size - ax.Position(3:4); 
hf.Position(3:4) = hf.Position(3:4) + sizeDiff; 
ax.Position(3:4) = axis_size;
print('multidose_bargraph', '-dpsc', '-append');
system(['ps2pdf ', 'multidose_bargraph.ps ', 'multidose_bargraph.pdf']);
system('rm multidose_bargraph.ps');


table_dose = table(outs{7}.fit_beta.t_sim(1:5:721), ...
                   outs{7}.fit_beta.y_sim(1:5:721,13), ...
                   out_low_6.fit_beta.y_sim(1:5:721,13), ...
                   out_low_5.fit_beta.y_sim(1:5:721,13), ...
                   out_low_4.fit_beta.y_sim(1:5:721,13), ...
                   out_low_3.fit_beta.y_sim(1:5:721,13), ...
                   out_low_2.fit_beta.y_sim(1:5:721,13), ...
                   out_low_1.fit_beta.y_sim(1:5:721,13) ...
);
table_dose.Properties.VariableNames = {'Time', 'High_IFN', 'x0_5xIFN', 'x0_25xIFN', 'x0_1xIFN', 'x0_05xIFN', 'x0_025xIFN', 'x0_01xIFN'}
writetable(table_dose, ['multidose.xlsx'])

table_dose_norm = table_dose;
table_dose_norm(:,2:end) = array2table(table2array(table_dose(:,2:end))/table2array(table_dose(table_dose.Time == 20,2)))
writetable(table_dose_norm, ['multidose_norm.xlsx'])



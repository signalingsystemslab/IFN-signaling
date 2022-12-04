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

%% Load Experimental data
t_exp_USP18_mRNA = xlsread('Data/20201112/USP18 mRNA.xlsx', 'Time Points')*60;
y_exp_USP18_mRNA = 2.^xlsread('Data/20201112/USP18 mRNA.xlsx', 'IFNBeta');
y_exp_USP18_mRNA_sem = 2.^xlsread('Data/20201112/USP18 mRNA.xlsx', 'IFNBetaSEM');

t_exp_STAT2_mRNA = xlsread('Data/20201112/STAT2 mRNA.xlsx', 'Time Points')*60;
y_exp_STAT2_mRNA = 2.^xlsread('Data/20201112/STAT2 mRNA.xlsx', 'IFNBeta');
y_exp_STAT2_mRNA_sem = 2.^xlsread('Data/20201112/STAT2 mRNA.xlsx', 'IFNBetaSEM');

t_exp_STAT1_mRNA = xlsread('Data/20201112/STAT1 mRNA.xlsx', 'Time Points')*60;
y_exp_STAT1_mRNA = 2.^xlsread('Data/20201112/STAT1 mRNA.xlsx', 'IFNBeta');
y_exp_STAT1_mRNA_sem = 2.^xlsread('Data/20201112/STAT1 mRNA.xlsx', 'IFNBetaSEM');

t_exp_SOCS1_mRNA = xlsread('Data/20201112/SOCS1 mRNA.xlsx', 'Time Points')*60;
y_exp_SOCS1_mRNA = 2.^xlsread('Data/20201112/SOCS1 mRNA.xlsx', 'IFNBeta');
y_exp_SOCS1_mRNA_sem = 2.^xlsread('Data/20201112/SOCS1 mRNA.xlsx', 'IFNBetaSEM');

t_exp_SOCS3_mRNA = xlsread('Data/20201112/SOCS3 mRNA.xlsx', 'Time Points')*60;
y_exp_SOCS3_mRNA = 2.^xlsread('Data/20201112/SOCS3 mRNA.xlsx', 'IFNBeta');
y_exp_SOCS3_mRNA_sem = 2.^xlsread('Data/20201112/SOCS3 mRNA.xlsx', 'IFNBetaSEM');

t_exp_IRF9_mRNA = xlsread('Data/20201112/IRF9 mRNA.xlsx', 'Time Points')*60;
y_exp_IRF9_mRNA = 2.^xlsread('Data/20201112/IRF9 mRNA.xlsx', 'IFNBeta');
y_exp_IRF9_mRNA_sem = 2.^xlsread('Data/20201112/IRF9 mRNA.xlsx', 'IFNBetaSEM');

t_exp_STAT1n = xlsread('Data/20201112/Nuclear total STAT1 protein.xlsx', 'Time Points');
y_exp_STAT1n = xlsread('Data/20201112/Nuclear total STAT1 protein.xlsx', 'IFNbeta');
y_exp_STAT1n_sem = xlsread('Data/20201112/Nuclear total STAT1 protein.xlsx', 'IFNbeta SEM');

t_exp_STAT2n = xlsread('Data/20201112/Nuclear total STAT2 protein.xlsx', 'Time Points');
y_exp_STAT2n = xlsread('Data/20201112/Nuclear total STAT2 protein.xlsx', 'IFNbeta');
y_exp_STAT2n_sem = xlsread('Data/20201112/Nuclear total STAT2 protein.xlsx', 'IFNbeta SEM');

t_exp_IRF9n = xlsread('Data/20201112/Nuclear total IRF9 protein.xlsx', 'Time Points');
y_exp_IRF9n = xlsread('Data/20201112/Nuclear total IRF9 protein.xlsx', 'IFNbeta');
y_exp_IRF9n_sem = xlsread('Data/20201112/Nuclear total IRF9 protein.xlsx', 'IFNbeta SEM');

t_exp_pSTAT1n = xlsread('Data/20201112/Nuclear pSTAT1 protein.xlsx', 'Time Points');
y_exp_pSTAT1n = xlsread('Data/20201112/Nuclear pSTAT1 protein.xlsx', 'IFNbeta');
y_exp_pSTAT1n_sem = xlsread('Data/20201112/Nuclear pSTAT1 protein.xlsx', 'IFNbeta SEM');

t_exp_pSTAT2n = xlsread('Data/20201112/Nuclear pSTAT2 protein.xlsx', 'Time Points');
y_exp_pSTAT2n = xlsread('Data/20201112/Nuclear pSTAT2 protein.xlsx', 'IFNbeta');
y_exp_pSTAT2n_sem = xlsread('Data/20201112/Nuclear pSTAT2 protein.xlsx', 'IFNbeta SEM');

t_exp_ISGF3n = xlsread('Data/20201112/Nuclear active ISGF3 complex AllNormBeta20min.xlsx', 'Time Points');
y_exp_ISGF3n = xlsread('Data/20201112/Nuclear active ISGF3 complex AllNormBeta20min.xlsx', 'IFNbeta');
y_exp_ISGF3n_sem = xlsread('Data/20201112/Nuclear active ISGF3 complex AllNormBeta20min.xlsx', 'IFNbeta SEM');

t_exp_ISGF3n_CHX = xlsread('Data/20201112/Cycloheximide Nuclear active ISGF3 complex.xlsx', 'Time Points');
y_exp_ISGF3n_CHX = xlsread('Data/20201112/Cycloheximide Nuclear active ISGF3 complex.xlsx', 'IFNbetaCyclo');
y_exp_ISGF3n_sem_CHX = xlsread('Data/20201112/Cycloheximide Nuclear active ISGF3 complex.xlsx', 'IFNbetaCyclo SEM');

t_exp_STAT1c = xlsread('Data/20201112/Cytoplasmic total STAT1 protein.xlsx', 'Time Points');
y_exp_STAT1c = xlsread('Data/20201112/Cytoplasmic total STAT1 protein.xlsx', 'IFNbeta');
y_exp_STAT1c_sem = xlsread('Data/20201112/Cytoplasmic total STAT1 protein.xlsx', 'IFNbeta SEM');

t_exp_STAT2c = xlsread('Data/20201112/Cytoplasmic total STAT2 protein.xlsx', 'Time Points');
y_exp_STAT2c = xlsread('Data/20201112/Cytoplasmic total STAT2 protein.xlsx', 'IFNbeta');
y_exp_STAT2c_sem = xlsread('Data/20201112/Cytoplasmic total STAT2 protein.xlsx', 'IFNbeta SEM');

t_exp_IRF9c = xlsread('Data/20201112/Cytoplasmic total IRF9 protein.xlsx', 'Time Points');
y_exp_IRF9c = xlsread('Data/20201112/Cytoplasmic total IRF9 protein.xlsx', 'IFNbeta');
y_exp_IRF9c_sem = xlsread('Data/20201112/Cytoplasmic total IRF9 protein.xlsx', 'IFNbeta SEM');

t_exp_pSTAT1c = xlsread('Data/20201112/Cytoplasmic pSTAT1 protein.xlsx', 'Time Points');
y_exp_pSTAT1c = xlsread('Data/20201112/Cytoplasmic pSTAT1 protein.xlsx', 'IFNbeta');
y_exp_pSTAT1c_sem = xlsread('Data/20201112/Cytoplasmic pSTAT1 protein.xlsx', 'IFNbeta SEM');

t_exp_pSTAT2c = xlsread('Data/20201112/Cytoplasmic pSTAT2 protein.xlsx', 'Time Points');
y_exp_pSTAT2c = xlsread('Data/20201112/Cytoplasmic pSTAT2 protein.xlsx', 'IFNbeta');
y_exp_pSTAT2c_sem = xlsread('Data/20201112/Cytoplasmic pSTAT2 protein.xlsx', 'IFNbeta SEM');

% rescale all mRNA such that FC at t=0 is 1
y_exp_USP18_mRNA_sem = y_exp_USP18_mRNA_sem/y_exp_USP18_mRNA(1);
y_exp_USP18_mRNA = y_exp_USP18_mRNA/y_exp_USP18_mRNA(1); 

y_exp_STAT2_mRNA_sem = y_exp_STAT2_mRNA_sem/y_exp_STAT2_mRNA(1);
y_exp_STAT2_mRNA = y_exp_STAT2_mRNA/y_exp_STAT2_mRNA(1); 

y_exp_STAT1_mRNA_sem = y_exp_STAT1_mRNA_sem/y_exp_STAT1_mRNA(1);
y_exp_STAT1_mRNA = y_exp_STAT1_mRNA/y_exp_STAT1_mRNA(1); 

y_exp_SOCS1_mRNA_sem = y_exp_SOCS1_mRNA_sem/y_exp_SOCS1_mRNA(1);
y_exp_SOCS1_mRNA = y_exp_SOCS1_mRNA/y_exp_SOCS1_mRNA(1); 

y_exp_SOCS3_mRNA_sem = y_exp_SOCS3_mRNA_sem/y_exp_SOCS3_mRNA(1);
y_exp_SOCS3_mRNA = y_exp_SOCS3_mRNA/y_exp_SOCS3_mRNA(1); 

y_exp_IRF9_mRNA_sem = y_exp_IRF9_mRNA_sem/y_exp_IRF9_mRNA(1);
y_exp_IRF9_mRNA = y_exp_IRF9_mRNA/y_exp_IRF9_mRNA(1); 

% add timepoint at 180min as the average of 120 and 240min
t_exp_ISGF3n = [t_exp_ISGF3n(1:9); mean(t_exp_ISGF3n(9:10)); t_exp_ISGF3n(10:14)];
y_exp_ISGF3n = [y_exp_ISGF3n(1:9); mean(y_exp_ISGF3n(9:10)); y_exp_ISGF3n(10:14)];
y_exp_ISGF3n_sem = [y_exp_ISGF3n_sem(1:9); mean(y_exp_ISGF3n_sem(9:10)); y_exp_ISGF3n_sem(10:14)];

%% Prepare data structure for simulations 
data_beta = struct( ...
    'ISGF3n',       struct('t_exp', t_exp_ISGF3n,       'y_exp', y_exp_ISGF3n,      'y_exp_sem', y_exp_ISGF3n_sem,      'species', [13],            'coeff', [1],       'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'ISGF3n_CHX',   struct('t_exp', t_exp_ISGF3n_CHX,   'y_exp', y_exp_ISGF3n_CHX,  'y_exp_sem', y_exp_ISGF3n_sem_CHX,  'species', [13],            'coeff', [1],       'CHX', true,     'fit_coeff', 1, 'mRNA', false), ...
    'STAT1c',       struct('t_exp', t_exp_STAT1c,       'y_exp', y_exp_STAT1c,      'y_exp_sem', y_exp_STAT1c_sem,      'species', [6,8,12,15],     'coeff', [1,1,1,2], 'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'pSTAT1c',      struct('t_exp', t_exp_pSTAT1c,      'y_exp', y_exp_pSTAT1c,     'y_exp_sem', y_exp_pSTAT1c_sem,      'species', [8,12,15],       'coeff', [1,1,2],   'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'STAT1n',       struct('t_exp', t_exp_STAT1n,       'y_exp', y_exp_STAT1n,      'y_exp_sem', y_exp_STAT1n_sem,     'species', [10,9,13,16],    'coeff', [1,1,1,2], 'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'pSTAT1n',      struct('t_exp', t_exp_pSTAT1n,      'y_exp', y_exp_pSTAT1n,     'y_exp_sem', y_exp_pSTAT1n_sem,     'species', [9,13,16],       'coeff', [1,1,2],   'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'STAT2c',       struct('t_exp', t_exp_STAT2c,       'y_exp', y_exp_STAT2c,      'y_exp_sem', y_exp_STAT2c_sem,      'species', [7,8,12],        'coeff', [1,1,1],   'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'pSTAT2c',      struct('t_exp', t_exp_pSTAT2c,      'y_exp', y_exp_pSTAT2c,     'y_exp_sem', y_exp_pSTAT2c_sem,     'species', [8,12],          'coeff', [1,1],     'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'STAT2n',       struct('t_exp', t_exp_STAT2n,       'y_exp', y_exp_STAT2n,      'y_exp_sem', y_exp_STAT2n_sem,      'species', [11,9,13],       'coeff', [1,1,1],   'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'pSTAT2n',      struct('t_exp', t_exp_pSTAT2n,      'y_exp', y_exp_pSTAT2n,     'y_exp_sem', y_exp_pSTAT2n_sem,     'species', [9,13],          'coeff', [1,1],     'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'IRF9c',        struct('t_exp', t_exp_IRF9c,        'y_exp', y_exp_IRF9c,       'y_exp_sem', y_exp_IRF9c_sem,       'species', [19,12],         'coeff', [1,1],     'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'IRF9n',        struct('t_exp', t_exp_IRF9n,        'y_exp', y_exp_IRF9n,       'y_exp_sem', y_exp_IRF9n_sem,       'species', [13,14],         'coeff', [1,1],     'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'STAT1mRNA',    struct('t_exp', t_exp_STAT1_mRNA,   'y_exp', y_exp_STAT1_mRNA,  'y_exp_sem', y_exp_STAT1_mRNA_sem,  'species', [24],            'coeff', [1],       'CHX', false,    'fit_coeff', 1, 'mRNA', true), ...
    'STAT2mRNA',    struct('t_exp', t_exp_STAT2_mRNA,   'y_exp', y_exp_STAT2_mRNA,  'y_exp_sem', y_exp_STAT2_mRNA_sem,  'species', [30],            'coeff', [1],       'CHX', false,    'fit_coeff', 1, 'mRNA', true), ...
    'IRF9mRNA',     struct('t_exp', t_exp_IRF9_mRNA,    'y_exp', y_exp_IRF9_mRNA,   'y_exp_sem', y_exp_IRF9_mRNA_sem,   'species', [33],            'coeff', [1],       'CHX', false,    'fit_coeff', 1, 'mRNA', true), ...
    'USP18mRNA',    struct('t_exp', t_exp_USP18_mRNA,   'y_exp', y_exp_USP18_mRNA,  'y_exp_sem', y_exp_USP18_mRNA_sem,  'species', [38],            'coeff', [1],       'CHX', false,    'fit_coeff', 1, 'mRNA', true), ...
    'SOCS1mRNA',    struct('t_exp', t_exp_SOCS1_mRNA,   'y_exp', y_exp_SOCS1_mRNA,  'y_exp_sem', y_exp_SOCS1_mRNA_sem,  'species', [40],            'coeff', [1],       'CHX', false,    'fit_coeff', 1, 'mRNA', true), ...
    'SOCS3mRNA',    struct('t_exp', t_exp_SOCS3_mRNA,   'y_exp', y_exp_SOCS3_mRNA,  'y_exp_sem', y_exp_SOCS3_mRNA_sem,  'species', [41],            'coeff', [1],       'CHX', false,    'fit_coeff', 1, 'mRNA', true));
     
plot_beta_metadata = struct( ...
    'ISGF3n',       struct('title', 'ISGF3n',          'ylabel', 'Normalized ISGF3n'), ...
    'ISGF3n_CHX',   struct('title', 'ISGF3n - CHX',     'ylabel', 'Normalized ISGF3n'), ...
    'STAT1c',       struct('title', 'Total STAT1c',     'ylabel', 'Normalized Total STAT1c'), ...
    'pSTAT1c',      struct('title', 'Total pSTAT1c',    'ylabel', 'Normalized Total pSTAT1c'), ...
    'STAT1n',       struct('title', 'Total STAT1n',     'ylabel', 'Normalized Total STAT1n'), ...
    'pSTAT1n',      struct('title', 'Total pSTAT1n',    'ylabel', 'Normalized Total pSTAT1n'), ...
    'STAT2c',       struct('title', 'Total STAT2c',     'ylabel', 'Normalized Total STAT2c'), ...
    'pSTAT2c',      struct('title', 'Total pSTAT2c',    'ylabel', 'Normalized Total pSTAT2c'), ...
    'STAT2n',       struct('title', 'Total STAT2n',     'ylabel', 'Normalized Total STAT2n'), ...
    'pSTAT2n',      struct('title', 'Total pSTAT2n',    'ylabel', 'Normalized Total pSTAT2n'), ...
    'IRF9c',        struct('title', 'Total IRF9c',      'ylabel', 'Normalized Total IRF9c'), ...
    'IRF9n',        struct('title', 'Total IRF9n',      'ylabel', 'Normalized Total IRF9n'), ...
    'STAT1mRNA',    struct('title', 'STAT1_{mRNA}',     'ylabel', 'STAT1_{mRNA} log2FC'), ...
    'STAT2mRNA',    struct('title', 'STAT2_{mRNA}',     'ylabel', 'STAT2_{mRNA} log2FC'), ...
    'IRF9mRNA',     struct('title', 'IRF9_{mRNA}',      'ylabel', 'IRF9_{mRNA} log2FC'), ...
    'USP18mRNA',    struct('title', 'USP18_{mRNA}',     'ylabel', 'USP18_{mRNA} log2FC'), ...
    'SOCS1mRNA',    struct('title', 'SOCS1_{mRNA}',     'ylabel', 'SOCS1_{mRNA} log2FC'), ...
    'SOCS3mRNA',    struct('title', 'SOCS3_{mRNA}',     'ylabel', 'SOCS3_{mRNA} log2FC'));

% Search 3 by three
p.obj = parpool(16);
timeout = seconds(5);

% calculate time max of experiments
tmp = cellfun(@(x) data_beta.(x).t_exp, fieldnames(data_beta), 'UniformOutput',false)
t_max = max(vertcat(tmp{:}));

IFNbeta_dose = 5e4; 

% get first steady state
out = solveODE_withSS (param_org, param_names, [], {}, repmat(0,41,1), 0:t_max, 0);
y0 = out.y0;
    
%% run optimisation 3 by 3
chosen = [];
best_values_overall = []; % list best test for each number of free param
best_pso_overall = []; % list best pso index for each number of free param

for i = 1:7
    n_free = 3*i; 

    lb = [repmat(1,1,3) repmat(-5,1,n_free)];
    ub = [repmat(size(param_org,1)+1-(n_free -3), 1,3),repmat(5,1,n_free)];

    free_idx = setdiff(1:size(param_org,1),chosen);
    
    swarmsize = min(500,10*(n_free+3));

    SSE_beta=@(x) SSE_eval(solveODE_withSS(param_org([chosen, free_idx(floor(x(1:3)))]).*10.^x(4:end)', param_names([chosen, free_idx(floor(x(1:3)))]), param_org(setdiff(1:size(param_org,1),[chosen,free_idx(floor(x(1:3)))])), param_names(setdiff(1:size(param_org,1),[chosen, free_idx(floor(x(1:3)))])), y0, 0:t_max, IFNbeta_dose), data_beta, true, true);
    SSE_beta_vec=@(x) vecfcn(x, SSE_beta, p.obj, timeout) ;

    mkdir(['Results/', num2str(n_free), '_free']);
    save(['Results/', num2str(n_free), '_free/previously_chosen.mat'],'chosen', '-v7.3');
    
    best_value=nan;
    best_pso=nan;
    best_chosen=nan;
    for n_pso=1:20
        rng('shuffle'); % initialize random generator based on time

        dir = ['Results/', num2str(n_free), '_free/PSO', num2str(n_pso)];
        mkdir(dir);

        outfun = @(optimValues,state) save_pso(optimValues, state, [dir, '/PSO_',num2str(swarmsize),'_n', num2str(n_pso), '_part_iter']);
        options_pso = optimoptions('particleswarm', 'Display', 'iter', 'DisplayInterval', 25, 'SwarmSize', swarmsize, 'UseVectorized', true, 'OutputFcn', outfun);
        [x,fval,exitflag,output] = particleswarm(SSE_beta_vec, 3+n_free, lb, ub,  options_pso);
        res = struct('x', x, 'fval', fval, 'exitflag', exitflag, 'output', output); 
        save([dir,'/PSO_', num2str(swarmsize), 'particles_n', num2str(n_pso), '_final.mat'], 'res', '-v7.3');
        
        if isnan(best_value) || res.fval < best_value
            best_value = res.fval;
            best_pso = n_pso;
            best_chosen = free_idx(floor(res.x(1:3)));
        end
        
        disp(['PSO : ', num2str(n_pso), ' (best ', num2str(best_pso), ')']);
        disp(['Value : ', num2str(res.fval), ' (best ', num2str(best_value), ')']);
        disp(['Chosen parameters : ', strjoin({param_names{free_idx(floor(res.x(1:3)))}}, ', '), ' (best ',  strjoin({param_names{best_chosen}}, ', '), ')']) ;
    end
    
    if ~isempty(best_values_overall) 
        if best_value >= min(best_values_overall) % if not better than before try again
            swarmsize = min(500,20*(n_free+3));
            for n_pso=21:40
                rng('shuffle'); % initialize random generator based on time

                dir = ['Results/', num2str(n_free), '_free/PSO', num2str(n_pso)];
                mkdir(dir);

                outfun = @(optimValues,state) save_pso(optimValues, state, [dir, '/PSO_',num2str(swarmsize),'_n', num2str(n_pso), '_part_iter']);
                options_pso = optimoptions('particleswarm', 'Display', 'iter', 'DisplayInterval', 25, 'SwarmSize', swarmsize, 'UseVectorized', true, 'OutputFcn', outfun);
                [x,fval,exitflag,output] = particleswarm(SSE_beta_vec, 3+n_free, lb, ub,  options_pso);
                res = struct('x', x, 'fval', fval, 'exitflag', exitflag, 'output', output); 
                save([dir,'/PSO_', num2str(swarmsize), 'particles_n', num2str(n_pso), '_final.mat'], 'res', '-v7.3');
                
                if isnan(best_value) || res.fval < best_value
                    best_value = res.fval;
                    best_pso = n_pso;
                    best_chosen = free_idx(floor(res.x(1:3)));
                end
                
                disp(['PSO : ', num2str(n_pso), ' (best ', num2str(best_pso), ')']);
                disp(['Value : ', num2str(res.fval), ' (best ', num2str(best_value), ')']);
                disp(['Chosen parameters : ', strjoin({param_names{free_idx(floor(res.x(1:3)))}}, ', '), ' (best ',  strjoin({param_names{best_chosen}}, ', '), ')']) ;
            end
        end
    end
    
    best_values_overall = [best_values_overall, best_value];
    best_pso_overall = [best_pso_overall, best_pso];
    
    % load best 
    chosen = [chosen, best_chosen];
    free_idx = setdiff(1:size(param_org,1),chosen);
    n_free = n_free + 3;
    
    disp(['Best Values : ', strjoin(arrayfun(@(x) num2str(x),best_values_overall, 'UniformOutput', false), ', ')]);
    disp(['Best PSOs : ', strjoin(arrayfun(@(x) num2str(x),best_pso_overall, 'UniformOutput', false), ', ')]);
    disp(['Chosen parameters : ', strjoin({param_names{chosen}}, ', ')]);
end

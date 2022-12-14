% load all beta data
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

t_exp_ISGF3n = xlsread('Data/20210817/Nuclear active ISGF3 complex AllNormBeta20min.xlsx', 'Time');
y_exp_ISGF3n = xlsread('Data/20210817/Nuclear active ISGF3 complex AllNormBeta20min.xlsx', 'IFNbeta');
y_exp_ISGF3n_sem = xlsread('Data/20210817/Nuclear active ISGF3 complex AllNormBeta20min.xlsx', 'IFNbeta SEM');

t_exp_ISGF3n_CHX = xlsread('Data/20201112/Cycloheximide Nuclear active ISGF3 complex.xlsx', 'Time Points');
y_exp_ISGF3n_CHX = xlsread('Data/20201112/Cycloheximide Nuclear active ISGF3 complex.xlsx', 'IFNbetaCyclo');
y_exp_ISGF3n_sem_CHX = xlsread('Data/20201112/Cycloheximide Nuclear active ISGF3 complex.xlsx', 'IFNbetaCyclo SEM');

t_exp_STAT1c = xlsread('Data/20210513/Cytoplasmic total STAT1 protein.xlsx', 'Time Points');
y_exp_STAT1c = xlsread('Data/20210513/Cytoplasmic total STAT1 protein.xlsx', 'IFNbeta');
y_exp_STAT1c_sem = xlsread('Data/20210513/Cytoplasmic total STAT1 protein.xlsx', 'IFNbeta SEM');

t_exp_STAT2c = xlsread('Data/20210513/Cytoplasmic total STAT2 protein.xlsx', 'Time Points');
y_exp_STAT2c = xlsread('Data/20210513/Cytoplasmic total STAT2 protein.xlsx', 'IFNbeta');
y_exp_STAT2c_sem = xlsread('Data/20210513/Cytoplasmic total STAT2 protein.xlsx', 'IFNbeta SEM');

t_exp_IRF9c = xlsread('Data/20210513/Cytoplasmic total IRF9 protein.xlsx', 'Time Points');
y_exp_IRF9c = xlsread('Data/20210513/Cytoplasmic total IRF9 protein.xlsx', 'IFNbeta');
y_exp_IRF9c_sem = xlsread('Data/20210513/Cytoplasmic total IRF9 protein.xlsx', 'IFNbeta SEM');

t_exp_pSTAT1c = xlsread('Data/20210513/Cytoplasmic pSTAT1 protein.xlsx', 'Time Points');
y_exp_pSTAT1c = xlsread('Data/20210513/Cytoplasmic pSTAT1 protein.xlsx', 'IFNbeta');
y_exp_pSTAT1c_sem = xlsread('Data/20210513/Cytoplasmic pSTAT1 protein.xlsx', 'IFNbeta SEM');

t_exp_pSTAT2c = xlsread('Data/20210513/Cytoplasmic pSTAT2 protein.xlsx', 'Time Points');
y_exp_pSTAT2c = xlsread('Data/20210513/Cytoplasmic pSTAT2 protein.xlsx', 'IFNbeta');
y_exp_pSTAT2c_sem = xlsread('Data/20210513/Cytoplasmic pSTAT2 protein.xlsx', 'IFNbeta SEM');

% rescale all mRNA to FC
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

%% Prepare data structure for simulations 
%% Species with delay before mRNA
data_beta = struct( ...
    'ISGF3n',       struct('t_exp', t_exp_ISGF3n,       'y_exp', y_exp_ISGF3n,      'y_exp_sem', y_exp_ISGF3n_sem,      'species', [13],            'coeff', [1],       'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'ISGF3n_CHX',   struct('t_exp', t_exp_ISGF3n_CHX,   'y_exp', y_exp_ISGF3n_CHX,  'y_exp_sem', y_exp_ISGF3n_sem_CHX,  'species', [13],            'coeff', [1],       'CHX', true,     'fit_coeff', 1, 'mRNA', false), ...
    'STAT1c',       struct('t_exp', t_exp_STAT1c,       'y_exp', y_exp_STAT1c,      'y_exp_sem', y_exp_STAT1c_sem,      'species', [6,8,12,15],     'coeff', [1,1,1,2], 'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'pSTAT1c',      struct('t_exp', t_exp_pSTAT1c,      'y_exp', y_exp_pSTAT1c,     'y_exp_sem', y_exp_pSTAT1c_sem,      'species', [8,12,15],       'coeff', [1,1,2],   'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'STAT1n',       struct('t_exp', t_exp_STAT1n,       'y_exp', y_exp_STAT1n,      'y_exp_sem', y_exp_STAT1n_sem,     'species', [10,9,13,16],    'coeff', [1,1,1,2], 'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'pSTAT1n',      struct('t_exp', t_exp_pSTAT1n,      'y_exp', y_exp_pSTAT1n,     'y_exp_sem', y_exp_pSTAT1n_sem,     'species', [9,13,16],       'coeff', [1,1,2],   'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'STAT2c',       struct('t_exp', t_exp_STAT2c,       'y_exp', y_exp_STAT2c,      'y_exp_sem', y_exp_STAT2c_sem,      'species', [7,8,12],        'coeff', [1,1,1],   'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'pSTAT2c',      struct('t_exp', t_exp_pSTAT2c,      'y_exp', y_exp_pSTAT2c,     'y_exp_sem', y_exp_pSTAT2c_sem,     'species', [8,12],          'coeff', [1,1],     'CHX', false,    'fit_coeff', 0, 'mRNA', false), ...
    'STAT2n',       struct('t_exp', t_exp_STAT2n,       'y_exp', y_exp_STAT2n,      'y_exp_sem', y_exp_STAT2n_sem,      'species', [11,9,13],       'coeff', [1,1,1],   'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'pSTAT2n',      struct('t_exp', t_exp_pSTAT2n,      'y_exp', y_exp_pSTAT2n,     'y_exp_sem', y_exp_pSTAT2n_sem,     'species', [9,13],          'coeff', [1,1],     'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'IRF9c',        struct('t_exp', t_exp_IRF9c,        'y_exp', y_exp_IRF9c,       'y_exp_sem', y_exp_IRF9c_sem,       'species', [19,12],         'coeff', [1,1],     'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'IRF9n',        struct('t_exp', t_exp_IRF9n,        'y_exp', y_exp_IRF9n,       'y_exp_sem', y_exp_IRF9n_sem,       'species', [13,14],         'coeff', [1,1],     'CHX', false,    'fit_coeff', 1, 'mRNA', false), ...
    'STAT1mRNA',    struct('t_exp', t_exp_STAT1_mRNA,   'y_exp', y_exp_STAT1_mRNA,  'y_exp_sem', y_exp_STAT1_mRNA_sem,  'species', [24],            'coeff', [1],       'CHX', false,    'fit_coeff', 1, 'mRNA', true), ...
    'STAT2mRNA',    struct('t_exp', t_exp_STAT2_mRNA,   'y_exp', y_exp_STAT2_mRNA,  'y_exp_sem', y_exp_STAT2_mRNA_sem,  'species', [30],            'coeff', [1],       'CHX', false,    'fit_coeff', 1, 'mRNA', true), ...
    'IRF9mRNA',     struct('t_exp', t_exp_IRF9_mRNA,    'y_exp', y_exp_IRF9_mRNA,   'y_exp_sem', y_exp_IRF9_mRNA_sem,   'species', [33],            'coeff', [1],       'CHX', false,    'fit_coeff', 1, 'mRNA', true), ...
    'USP18mRNA',    struct('t_exp', t_exp_USP18_mRNA,   'y_exp', y_exp_USP18_mRNA,  'y_exp_sem', y_exp_USP18_mRNA_sem,  'species', [38],            'coeff', [1],       'CHX', false,    'fit_coeff', 1, 'mRNA', true), ...
    'SOCS1mRNA',    struct('t_exp', t_exp_SOCS1_mRNA,   'y_exp', y_exp_SOCS1_mRNA,  'y_exp_sem', y_exp_SOCS1_mRNA_sem,  'species', [40],            'coeff', [1],       'CHX', false,    'fit_coeff', 1, 'mRNA', true), ...
    'SOCS3mRNA',    struct('t_exp', t_exp_SOCS3_mRNA,   'y_exp', y_exp_SOCS3_mRNA,  'y_exp_sem', y_exp_SOCS3_mRNA_sem,  'species', [41],            'coeff', [1],       'CHX', false,    'fit_coeff', 1, 'mRNA', true), ...
    'Rec',          struct('t_exp', [],                 'y_exp', [],                'y_exp_sem', [],                    'species', [1,4],           'coeff', [1,1],     'CHX', false,    'fit_coeff', 0, 'mRNA', false), ...
    'IFN',          struct('t_exp', [],                 'y_exp', [],                'y_exp_sem', [],                    'species', [3,4],           'coeff', [1,1],     'CHX', false,    'fit_coeff', 0, 'mRNA', false), ...
    'SOCS1',        struct('t_exp', [],                 'y_exp', [],                'y_exp_sem', [],                    'species', [2],             'coeff', [1],       'CHX', false,    'fit_coeff', 0, 'mRNA', false), ...
    'SOCS3',        struct('t_exp', [],                 'y_exp', [],                'y_exp_sem', [],                    'species', [20],            'coeff', [1],       'CHX', false,    'fit_coeff', 0, 'mRNA', false), ...
    'USP18',        struct('t_exp', [],                 'y_exp', [],                'y_exp_sem', [],                    'species', [5],             'coeff', [1],       'CHX', false,    'fit_coeff', 0, 'mRNA', false), ...
    'IRF2',         struct('t_exp', [],                 'y_exp', [],                'y_exp_sem', [],                    'species', [35],            'coeff', [1],       'CHX', false,    'fit_coeff', 0, 'mRNA', false) ...
);


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
    'SOCS3mRNA',    struct('title', 'SOCS3_{mRNA}',     'ylabel', 'SOCS3_{mRNA} log2FC'), ...
    'Rec',          struct('title', 'Total Rec',        'ylabel', 'Total Rec'), ...
    'IFN',          struct('title', 'Total IFN',        'ylabel', 'Total IFN'), ...
    'SOCS1',        struct('title', 'SOCS1',            'ylabel', 'SOCS1'), ...
    'SOCS3',        struct('title', 'SOCS3',            'ylabel', 'SOCS3'), ...
    'USP18',        struct('title', 'USP18',            'ylabel', 'USP18'), ...
    'IRF2',         struct('title', 'IRF2',             'ylabel', 'IRF2') ...
);
   

% Pulse data 
t_exp_ISGF3n_pulse = xlsread('Data/20210817/Nuclear active ISGF3 complex AllNormBeta20min.xlsx', 'Time');
y_exp_ISGF3n_pulse = xlsread('Data/20210817/Nuclear active ISGF3 complex AllNormBeta20min.xlsx', 'Pulsebeta');

t_exp_ISGF3n_pulse = t_exp_ISGF3n_pulse(~isnan(y_exp_ISGF3n_pulse));
y_exp_ISGF3n_pulse = y_exp_ISGF3n_pulse(~isnan(y_exp_ISGF3n_pulse));

data_beta_pulse = struct(...
    'ISGF3n_pulse', struct('t_exp', t_exp_ISGF3n_pulse, 'y_exp', y_exp_ISGF3n_pulse, 'y_exp_sem', [], ...
                      'species', 13, 'coeff', [1], 'CHX', false, 'fit_coeff', 0, 'mRNA', false, 't_norm', data_beta.ISGF3n.t_exp(data_beta.ISGF3n.y_exp_sem == 0)) ...
);
plot_beta_pulse_metadata = struct( ...
    'ISGF3n_pulse',       struct('title', 'ISGF3n',          'ylabel', 'Normalized ISGF3n') ...
);
  
% Low IFN data 
t_exp_ISGF3n_low = xlsread('Data/20210817/Nuclear active ISGF3 complex AllNormBeta20min.xlsx', 'Time');
y_exp_ISGF3n_low = xlsread('Data/20210817/Nuclear active ISGF3 complex AllNormBeta20min.xlsx', 'Lowbeta');
y_exp_ISGF3n_low_sem = xlsread('Data/20210817/Nuclear active ISGF3 complex AllNormBeta20min.xlsx', 'Lowbeta SEM');

t_exp_ISGF3n_low = t_exp_ISGF3n_low(~isnan(y_exp_ISGF3n_low));
y_exp_ISGF3n_low_sem = y_exp_ISGF3n_low_sem(~isnan(y_exp_ISGF3n_low));
y_exp_ISGF3n_low = y_exp_ISGF3n_low(~isnan(y_exp_ISGF3n_low));

data_beta_low = struct(...
    'ISGF3n_low', struct('t_exp', t_exp_ISGF3n_low, 'y_exp', y_exp_ISGF3n_low, 'y_exp_sem', y_exp_ISGF3n_low_sem, ...
                      'species', 13, 'coeff', [1], 'CHX', false, 'fit_coeff', 0, 'mRNA', false, 't_norm', data_beta.ISGF3n.t_exp(data_beta.ISGF3n.y_exp_sem == 0)) ...
);
plot_beta_low_metadata = struct( ...
    'ISGF3n_low',       struct('title', 'ISGF3n - low dose',          'ylabel', 'Normalized ISGF3n') ...
);

save('data.mat', 'data_beta', 'plot_beta_metadata', 'data_beta_pulse', 'plot_beta_pulse_metadata', 'data_beta_low', 'plot_beta_low_metadata')

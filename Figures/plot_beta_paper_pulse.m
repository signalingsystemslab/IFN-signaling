function [] = plot_beta_paper_pulse(outs, outs_names, WT_idx, data, metadata, filename, varargin)

SSE = cellfun(@(x) SSE_eval(x, data, true, true), outs);
SSE_ind = cellfun(@(x) SSE_eval(x, data, false, true), outs, 'UniformOutput', false);

ordered=false;
for k=1:length(varargin)
    if strcmpi(varargin{k},'sort')
        ordered = varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
        break
    end
end

cm = flipud(parula(length(outs))); 
for k=1:length(varargin)
    if strcmpi(varargin{k},'colors')
        cm = varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
        break
    end
end


data_name = 'data';
for k=1:length(varargin)
    if strcmpi(varargin{k},'data_name')
        data_name = varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
        break
    end
end

if ordered
    [~, order_idx] = sort(SSE);
else
    order_idx=1:length(outs);
end

outs = {outs{order_idx}};
SSE = SSE(order_idx);
SSE_ind = {SSE_ind{order_idx}};

t_tmp = cellfun(@(x) data.(x).t_exp, fieldnames(data), 'UniformOutput',false);
t_all = 0:max(vertcat(t_tmp{:}));


% fs=14/sqrt(length(outs));
fs = 14/sqrt(7); % best only

fn = fieldnames(data);
for i = 1:size(fn,1)
    d = fn{i};

    hf = figure;
    hold on;
    if length(data.(d).t_exp) ~= 0
        t_exp = data.(d).t_exp;
        y_exp = data.(d).y_exp;
        y_exp_sd = data.(d).y_exp_sem;
    
        errorbar(t_exp, y_exp, y_exp_sd , '--ok');
        
        % get norm factor:
         if(data.(d).CHX)
            y_sim_all = outs{WT_idx}.fit_beta_chx.y_sim;
        else
            y_sim_all = outs{WT_idx}.fit_beta.y_sim;
        end
        
        y_sim_out_all = 0;
        for s = 1:size(data.(d).species,2)
            y_sim_out_all = y_sim_out_all + data.(d).coeff(1,s) * y_sim_all(:,data.(d).species(1,s));
        end
            
        alpha = 1/y_sim_out_all(1+data.(d).t_norm);
            
        if alpha < 1/max(y_sim_out_all)
            alpha = 1/max(y_sim_out_all);
        elseif alpha > 1/max([min(y_sim_out_all),0])
            alpha = 1/min(y_sim_out_all);
        end
        
        for j = 1:length(outs)
            if(data.(d).CHX)
                y_sim_all = outs{j}.fit_beta_chx.y_sim;
                y_sim = y_sim_all(1+t_exp,:);
            else
                y_sim_all = outs{j}.fit_beta.y_sim;
                y_sim = y_sim_all(1+t_exp,:);
            end

            y_sim_out=0;
            y_sim_out_all = 0;
            for s = 1:size(data.(d).species,2)
                y_sim_out = y_sim_out + data.(d).coeff(1,s) * y_sim(:,data.(d).species(1,s));
                y_sim_out_all = y_sim_out_all + data.(d).coeff(1,s) * y_sim_all(:,data.(d).species(1,s));
            end
            
            y_sim_out_norm = y_sim_out_all * alpha;
            if j == WT_idx
                plot(t_all(t_all <= max(data.(d).t_exp)), y_sim_out_norm(t_all <= max(data.(d).t_exp)), 'color', cm(j,:));
            else
                plot(t_all(t_all <= max(data.(d).t_exp)), y_sim_out_norm(t_all <= max(data.(d).t_exp)), 'color', cm(j,:), 'LineWidth', 1.5);
            end
        end
    else
        for j = 1:length(outs)
            if(data.(d).CHX)
                y_sim_all = outs{j}.fit_beta_chx.y_sim;
                y_sim = y_sim_all(1+t_exp,:);
            else
                y_sim_all = outs{j}.fit_beta.y_sim;
                y_sim = y_sim_all(1+t_exp,:);
            end
            
            y_sim_out_all = 0;
            for s = 1:size(data.(d).species,2)
                y_sim_out_all = y_sim_out_all + data.(d).coeff(1,s) * y_sim_all(:,data.(d).species(1,s));
            end
            plot(t_all, y_sim_out_all, 'color', cm(j,:));
        end
    end
    
    fm = fieldnames(metadata.(d));
    for k = 1:size(fm,1)
        f = fm{k};
        eval([f, '(metadata.(d).(f))']);
    end
    xlabel('Time (min)');
    if data.(d).mRNA
        set(gca, 'YScale', 'log')
    end

    lgd_names = outs_names;
    if length(data.(d).t_exp) ~= 0
        lgd_names = [data_name, outs_names];
    end
    legend(lgd_names, 'Location', 'southoutside', 'FontSize', fs, 'Box', 'off');
    pbaspect([1 1 1]);
    set(gca,'FontSize',6);

    temp = hf.Position;
    
    ax = gca(hf); 
    ax.Units = 'pixels';
    axis_size = [140, 105];
    sizeDiff = axis_size - ax.Position(3:4); 
    hf.Position(3:4) = hf.Position(3:4) + sizeDiff; 
    ax.Position(3:4) = axis_size;
    
    print(filename, '-dpsc', '-append');
end
system(['ps2pdf ', regexprep(filename, "([^\\]) ", "$1\\ "), '.ps ', regexprep(filename, "([^\\]) ", "$1\\ "), '.pdf']);
system(['rm ', regexprep(filename, "([^\\]) ", "$1\\ "), '.ps']);

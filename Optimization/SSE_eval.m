function SSE = SSE_eval(out, data, calc_sum, rmsd)
    
    if out.failed || ~out.ss
        SSE = NaN;
        %disp('failed');
        return
    end    
    
    SSE=zeros(size(fieldnames(data),1),1);
    fn = fieldnames(data);
    for i = 1:size(fieldnames(data),1)
        d = fn{i};
        t_exp = data.(d).t_exp;
        y_exp = data.(d).y_exp;
        y_exp_sem = data.(d).y_exp_sem;
        
        % don't fit if no replicates
        if(all(isnan(y_exp_sem)) || (calc_sum && data.(d).fit_coeff == 0))
            continue
        end
        
        t_exp = t_exp(~isnan(y_exp_sem));
        y_exp = y_exp(~isnan(y_exp_sem));
        y_exp_sem = y_exp_sem(~isnan(y_exp_sem));
        
        if ( data.(d).CHX )
            y_sim_all = out.fit_beta_chx.y_sim;
            y_sim = out.fit_beta_chx.y_sim(1+t_exp,:);
        else
            y_sim_all = out.fit_beta.y_sim;
            y_sim = out.fit_beta.y_sim(1+t_exp,:);
        end
              
        y_sim_out=0;
        y_sim_out_all = 0;
        for s = 1:size(data.(d).species,2)
            y_sim_out = y_sim_out + data.(d).coeff(1,s) * y_sim(:,data.(d).species(1,s));
            y_sim_out_all = y_sim_out_all + data.(d).coeff(1,s) * y_sim_all(:,data.(d).species(1,s));
        end
        
        if any(y_exp_sem == 0)
            alpha = 1/y_sim_out(y_exp_sem == 0);
            y_sim_out = y_sim_out(y_exp_sem ~= 0);
            y_exp = y_exp(y_exp_sem ~= 0);
            y_exp_sem = y_exp_sem(y_exp_sem ~= 0);
        elseif any(y_exp == 1)
            alpha = 1/y_sim_out(y_exp == 1);
            y_sim_out = y_sim_out(y_exp ~= 1);
            y_exp_sem = y_exp_sem(y_exp ~= 1);            
            y_exp = y_exp(y_exp ~= 1);
        else
            alpha = sum((y_sim_out .* y_exp)./y_exp_sem.^2)/sum(y_sim_out.^2 ./ y_exp_sem.^2); 
        end
        
        if alpha < 1/max(y_sim_out_all)
            alpha = 1/max(y_sim_out_all);
        elseif alpha > 1/max([min(y_sim_out_all),0])
            alpha = 1/min(y_sim_out_all);
        end
        
        y_sim_out_norm = y_sim_out * alpha;
        
        try
            SSE(i) = sum( (y_sim_out_norm - y_exp).^2 ./ (y_exp_sem.^2) );
            if rmsd
                SSE(i) = sqrt(SSE(i)/size(y_exp,1));
            end
        catch
            SSE(i) = NaN;
        end
    end
    
    if calc_sum
        SSE = sum(cellfun(@(fn) data.(fn).fit_coeff, fieldnames(data)) .* SSE);
    end
    % disp(SSE);
end

function out = solveODE_withSS (param, param_names, param_fixed, param_names_fixed, y0, t_out, dose, varargin )    
    % delay before mRNA
    
    args=[param_names, num2cell(param); param_names_fixed, num2cell(param_fixed)]';
    par = struct(args{:}); 
     
    Stoechio = sparse(size(y0,1),75) ;
    Stoechio(1,1:4) = [1, -1, -1, -1];                      % Rec
    Stoechio(2,[32,75]) = [-1, 1];                          % SOCS1
    Stoechio(3,4) = -1;                                     % IFN
    Stoechio(4,4:7) = [1, -1, -1, -1];                      % aRecIFN
    Stoechio(5,[8,74]) = [-1, 1];                           % USP18
    Stoechio(6,[9,15,23:25,71]) = [-1, -2, -1, 1, -1, 1];   % STAT1c
    Stoechio(7,[9,26:28,72]) = [-1, -1, 1, -1, 1];          % STAT2c
    Stoechio(8,[9,10,12]) = [1, -1, -1];                    % pSTAT1pSTAT2c 
    Stoechio(9,10:11) = [1, -1];                            % pSTAT1pSTAT2n
    Stoechio(10,[11,14,17,24:25]) = [1, 1, 2, -1, 1];       % STAT1n
    Stoechio(11,[11,14,27:28]) = [1, 1, -1, 1];             % STAT2n
    Stoechio(12,12:13) = [1, -1];                           % ISGF3c
    Stoechio(13,13:14) = [1, -1];                           % ISGF3n
    Stoechio(14,[14,29:30]) = [1, -1, 1] ;                  % IRF9n
    Stoechio(15,15:16) = [1, -1];                           % pSTAT1dimc
    Stoechio(16,16:17) = [1, -1];                           % pSTAT1dimn
    Stoechio(17,18:20) = [1, 1, -1];                        % OccGAS_ISREbs
    Stoechio(18,21:22) = [1, -1];                           % OccGASbs
    Stoechio(19,[12,29:31,73]) = [-1, 1, -1, -1, 1];        % IRF9c
    Stoechio(20,33:34) = [1, -1];                           % SOCS3
    Stoechio(21,35:37) = [1, 1, -1];                        % STAT1mRNA_LC1
    Stoechio(22,37:38) = [1, -1];                           % STAT1mRNA_LC2
    Stoechio(23,38:39) = [1, -1];                           % STAT1mRNA_LC3
    Stoechio(24,39:40) = [1, -1];                           % STAT1mRNA
    Stoechio(25,41:43) = [1, 1, -1];                        % STAT2mRNA_LC1
    Stoechio(26,43:44) = [1, -1];                           % STAT2mRNA_LC2
    Stoechio(27,44:45) = [1, -1];                           % STAT2mRNA_LC3
    Stoechio(28,45:46) = [1, -1];                           % STAT2mRNA_LC4
    Stoechio(29,46:47) = [1, -1];                           % STAT2mRNA_LC5
    Stoechio(30,47:48) = [1, -1];                           % STAT2mRNA
    Stoechio(31,49:51) = [1, 1, -1];                        % IRF9mRNA_LC1
    Stoechio(32,51:52) = [1, -1];                           % IRF9mRNA_LC2
    Stoechio(33,52:53) = [1, -1];                           % IRF9mRNA
    Stoechio(34,54:56) = [1, 1, -1];                        % IRF2mRNA
    Stoechio(35,[57,70]) = [1, -1];                         % IRF2
    Stoechio(36,58:60) = [1, 1, -1];                        % USP18mRNA_LC1
    Stoechio(37,60:61) = [1, -1];                           % USP18mRNA_LC2
    Stoechio(38,61:62) = [1, -1];                           % USP18mRNA
    Stoechio(39,63:65) = [1, 1, -1];                        % SOCS1mRNA_LC1
    Stoechio(40,65:66) = [1, -1];                           % SOCS1mRNA
    Stoechio(41,67:69) = [1, 1, -1];                        % SOCS3mRNA

    run_ss = true;
    for k = 1:length(varargin)
        if strcmpi(varargin{k},'run_ss')
            run_ss = varargin{k+1}; 
            varargin{k+1}=[]; 
            varargin{k}=[];
        end
    end
    
    failed = false;
    
    % run to steady state and update initial value variant
    if run_ss
        % run to steady state and update initial value variant
        ss_found = false;

        % options = odeset('NonNegative',ones(1,41), 'JPattern', JacPat);
        out.param = param;

        try         
            par.kchx = 1;
            J = @(t,y) 24 * jac_ifn_beta(t,y,par, Stoechio);
            options = odeset('NonNegative',ones(1,41), 'Jacobian', J); 
            
            tstart = 0; 
            tend = 10;
            sol = ode15s(@(t,y) 24 * ode_ifn_beta(t,y, par, Stoechio), [tstart (tstart + tend)/2 tend], y0, options); % run in day
            dy_h = ode_ifn_beta(sol.x(end),sol.y(:,end), par, Stoechio);

            while ~all(abs(dy_h) < 10^(-6))
                clear tstart tend dy_h;
                tstart = sol.x(end);
                tend = tstart + 10;
                sol = odextend(sol, @(t,y) 24 * ode_ifn_beta(t,y,par, Stoechio), tend); 
                
                % Remove earlier time points to free memory
                sol.y = sol.y(:,sol.x >= tstart); 
                sol.x = sol.x(sol.x >= tstart);
                dy_h = ode_ifn_beta(sol.x(end),sol.y(:,end),par, Stoechio);   

                % disp(max(dy))
                % if sol.x(end) >= 10*1000 || sol.x(end) ~= tend % 1000*10 day or failed 
                %     break
                % end
                
                if sol.x(end) ~= tend % failed 
                    break
                end
            end
            
            out.y0 = sol.y(:,end);
            
            if all(abs(dy_h) < 10^(-6))
                ss_found = true;
            end
        catch
            % disp(param);
        end
    else
        ss_found=true; 
        out.y0 = y0;
    end
    
    
    % run main simulation
    try      
        y1_v2 = out.y0;
        y1_v2(3) = dose; % set IFN dose
    
        par.kchx = 1;
        J = @(t,y) 1/60 * jac_ifn_beta(t,y,par,Stoechio);
        options = odeset('NonNegative',ones(1,41), 'Jacobian', J);

        [t_sim_out, y_sim_out] = ode15s(@(t,y) 1/60 * ode_ifn_beta(t,y,par,Stoechio), t_out, y1_v2, options); % in minutes

        if size(t_sim_out,1) ~= size(t_out,2)
            failed = true; 
        end
     
        % set for CHX run
        par.kchx = 0;
        J = @(t,y) 1/60 * jac_ifn_beta(t,y,par,Stoechio);
        options = odeset('NonNegative',ones(1,41), 'Jacobian', J);

        [t_sim_out_chx, y_sim_out_chx] = ode15s(@(t,y) 1/60 * ode_ifn_beta(t,y,par,Stoechio), t_out, y1_v2, options); % in minutes

        if size(t_sim_out_chx,1) ~= size(t_out,2)
            failed = true; 
        end
        
        out.fit_beta = struct( 't_sim', t_sim_out, 'y_sim', y_sim_out);
        out.fit_beta_chx = struct( 't_sim', t_sim_out_chx, 'y_sim', y_sim_out_chx);
    catch
        failed = true; 
    end
    
    % stop(timer_ode15s)
    
    out.ss = ss_found;
    out.failed = failed;
end

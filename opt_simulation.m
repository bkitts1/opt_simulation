
% nonlinear winrate
% vr distribution bimodal - how this changes the trajectory lines

% stats will remain persistent and will keep track of results from
% different parameters
if exist('stats','var')==0,
    stats = sparse([]);
end;

% ----------------------------------------------------
% Algorithm
usedistrib        = 0;
constraint        = 0;
dynamic           = 0;
nash              = 0;
pid               = 0;
exact             = 0; % exact target mode
acceptall         = 0;
adserver_shutdown = 0; % if 1 then adserver will shutdown if more imps than desired are bought
speedup           = 0; % when below pacing and kpi and new impression is high quality then can eject pacing
basecut           = 0;
buyaboveorig      = 0;

error_squared_power = 1;

lookahead_iterations = 0; % attempt to achieve convergence within this amount of time

zero_last_step = 1; % the last step when we exceed imps and budget, if zero_last_step then we force budgets to zero
% This is useful for final endpoint analysis
% However if we're more interested in trajectories then zero_last_step=0 -
% this will avoid a "kink" at the very end

buy_imps_overage = 1;

% graphs
onetime = 0;            % run with full animation and only one set of params
exterior_plot_only = 0; % only plot the exterior points in the parameter search

% pid controller
kprop = -.2;
kint  = -.1;
kdiff = 0;

% sub-periods
periods = 0; % 0 if don't use periods

% traffic volume and auction
impressions_per_cycle = 300; % how much traffic is generated each iteration
winrate_for_highbid   = 0.2; % 1 means all bids submitted win (non-competitive auction). 0.2 means 1 in 5 win
prediction_noise      = 0.0; % 0.1 means that actual kpi is additively offset to predicted with a standard deviation of 0.1 (ie. 10%)
prediction_bias       = -0.0;

% perturbation: percentage of the way through
starting_it = 6;
above_imps_list = [.5:.1:1.5];
above_kpi_list  = [.5:.1:1.5];

% Original targets
iterations                = 50; % number of cycles

impressions_budget_orig = 5000;
spend_budget_orig       = 50;
kpi_budget_orig         = 3500; % .7

% --------------------------------------------------------

step = 0.1; %0.001; % basecut variable
cutoff = 0; % basecut variable
figure(9);clf;hold on;

ecpm_calculated_orig      = 1000.*spend_budget_orig./impressions_budget_orig;
kpitargorig               = kpi_budget_orig ./ impressions_budget_orig;
spend_budget              = spend_budget_orig;
impressions_budget        = impressions_budget_orig;
kpi_budget                = kpi_budget_orig;


%vglobal = [0.1322    1.8678    0.5000    1.5000];
vglobal = [.5    1.5    0.5000    1.5000];


% error function
scaledown = 0.10;
kpis = 1;

desired_imps_this_cycle   = 0;
winrate = 0;

% this cycle
kpitarg      = [0.9]; % current targ
kpiweight    = ones(1,kpis);
pacingweight = 0;
impstargorig = impressions_budget ./ iterations;
impstarg     = impressions_budget ./ iterations;

% distrib bins
bins = 100;
distrib = zeros(bins,1);


% kpierr = vector of kpis across, iterations down

fprintf('cycle\tecpm_calculated\twon\tspend\timps\tspend_budget\timps_budget\tkpi\tkpi_budget\n');

if periods>0,
  periodlen = iterations ./ periods; % how long a period lasts
  period = (floor([1:iterations]' ./ periodlen)+1);
else
  periodlen = iterations;
  period = ones(iterations,1);
end;

impressions_desired_per_cycle = impressions_budget ./ (impressions_per_cycle .* iterations);

kpi_save = [];
impressions_save = [];
spend_save = [];

showgraphs = 0;
if onetime,
    showgraphs = 1; % ***** 0 for monte carlo
end;

impid=1;
kpiid=1;

% Run kpi and imp offsets
for impid = 1:length(above_imps_list),
    for kpiid = 1:length(above_kpi_list),

% perturbation for ecpm
above_ecpm = 1;

% perturbation for imps and kpi
if onetime, % **** hard override 0 for monte carlo
    above_imps = 1.5; %1.5; % 1.5
    above_kpi  = 0.05; %1.0; %0.05; %0.8;
else
    % read out the current above_imps and above_kpi
    above_imps = above_imps_list(impid);
    above_kpi  = above_kpi_list(kpiid);
    
    % only plot exterior
    if exterior_plot_only,
        if above_imps >= max(above_imps_list) || above_kpi <=min(above_kpi_list) || ...
           above_imps <= min(above_imps_list) || above_kpi >=max(above_kpi_list),
           % exterior
        else
          continue;
        end;
        % do all of them
    end;
end;

% reinitialize the its
ecpm_calculated = ecpm_calculated_orig;
spend           = 0;
impressions     = 0; % impressions this iteration
kpi             = zeros(1,kpis); % kpi this iteration

% reinitialize totals
spend_budget       = 0;
impressions_budget = 0;
kpi_budget         = 0;


% roll forward to starting_it
% ---------------------------

if starting_it > 1,
    % perturbation is implemented by assuming perturbed performance for some
    % number of itartions up to starting_it
    impressions(     1:starting_it,     1) = above_imps .* ones(starting_it,1)    * (impressions_budget_orig./iterations);
    ecpm_calculated( 1:starting_it,     1) = above_ecpm .* ones(starting_it,1)    * 1000.*spend_budget_orig./impressions_budget_orig;
    kpi(             1:starting_it,1:kpis) = above_kpi  .* impressions .* kpitargorig;
    spend(           1:starting_it,     1) = ecpm_calculated .* impressions ./ 1000;

    % adjust budgets for the roll forward
    kpi_budget         = kpi_budget_orig         - sum(kpi(         1:starting_it, 1:kpis));
    spend_budget       = spend_budget_orig       - sum(spend(       1:starting_it, 1));
    impressions_budget = impressions_budget_orig - sum(impressions( 1:starting_it, 1));
    kpitarg            = nan.*ones(starting_it,1);
end;

w=0.02; % model mapping bid to impressions
bid_for_desired_imps = [];


% Main loop
% -------------------------------------------------------------------------

for it = starting_it:iterations,

    % lookahead_iterations
    its_during_equal_delivery = iterations - (lookahead_iterations + it); % number of iterations that will need to run with original delivery
    imps_during_equal_delivery = its_during_equal_delivery .* (impressions_budget_orig ./ iterations); % imps to be served over that timeperiod

    if lookahead_iterations == 0 || it+lookahead_iterations >= iterations,
        desired_imps_this_cycle(it,1) = impressions_budget ./ (iterations-it+1); % iterations remaining
    else
        desired_imps_this_cycle(it,1) = (impressions_budget - imps_during_equal_delivery) ./ lookahead_iterations;
    end;
    
    bid_for_desired_imps(it) = desired_imps_this_cycle(it)./(w.*impressions_per_cycle);

    % sub-periods
    % ------------
    
    % if new period then restart budget
    if periods>0,        
        % Periods
        impressions_sofar = 0;
        spend_sofar       = 0;
        kpi_sofar         = 0;
        if it>1, % (otherwise these values are 0)
            impressions_sofar = sum(impressions( 1:it-1,1));
            spend_sofar       = sum(spend(       1:it-1,1));
            kpi_sofar         = sum(kpi(         1:it-1,1:kpis));
        end;
        
        % calc impressions still needed this period
        impressions_period = (impressions_budget_orig ./ iterations).*period(it).*periodlen - impressions_sofar;
        spend_period       = (spend_budget_orig       ./ iterations).*period(it).*periodlen - spend_sofar;
        kpi_period         = (kpi_budget_orig         ./ iterations).*period(it).*periodlen - kpi_sofar;

        % update ecpm_calculated
        ecpm_calculated(it,1) = 1000.*spend_period ./ impressions_period;    

        % update kpitarg
        kpitarg(it,1:kpis) = kpi_period ./ impressions_period; % kpis remaining / imps remaining
    else
        
        % Standard (no periods)
        % --------------------
        
       % keyboard;
        % No periods
        % update kpitarg
        kpitarg(it,1:kpis) = kpi_budget ./ impressions_budget; % kpis remaining / imps remaining
        % update ecpm_calculated
        ecpm_calculated(it,1) = 1000.*spend_budget ./ impressions_budget;            
        desired_kpi_this_cycle(it,1) = kpi_budget ./ iterations;
        desired_spend_this_cycle(it,1) = spend_budget ./ iterations;

        if lookahead_iterations>0 && it+lookahead_iterations < iterations,
            kpi_during_equal_delivery      = its_during_equal_delivery .* (kpi_budget_orig   ./ iterations); % imps to be served over that timeperiod
            spend_during_equal_delivery    = its_during_equal_delivery .* (spend_budget_orig ./ iterations); % imps to be served over that timeperiod
            desired_kpi_this_cycle(it,1)   = (kpi_budget   - kpi_during_equal_delivery)   ./ lookahead_iterations;
            desired_spend_this_cycle(it,1) = (spend_budget - spend_during_equal_delivery) ./ lookahead_iterations;
            kpitarg(it,1:kpis)    = desired_kpi_this_cycle(it,1)./desired_imps_this_cycle(it,1);
            ecpm_calculated(it,1) = 1000.*desired_spend_this_cycle(it,1)./desired_imps_this_cycle(it,1);
        end;
        % can we come up with a bid price that would lead to a particular
        % impressions pacing outcome?
        % model: bidprice -> winrate; then we get desired_imps_this_cycle;
        % invert to get bidprice; or equivalently, test all bid prices and
        % pick the one that produces the imps delivery desired
        % desired_imps this cycle impressions_ideal(it,1) = impressions_budget ./ 
    end;

    
    % pid controller
    % --------------
    
    if pid==1,
        
        % err
        if it>1,
            spend_err          = spend(it-1)        - (spend_budget_orig       ./ iterations);
            impressions_err    = impressions(it-1)  - (impressions_budget_orig ./ iterations);
            kpi_err            = kpi(it-1)          - (kpi_budget_orig         ./ iterations);
        end;

        % interr
        spend_interr       = (spend_budget_orig - spend_budget)             - (spend_budget_orig       .* (it./iterations)); 
        impressions_interr = (impressions_budget_orig - impressions_budget) - (impressions_budget_orig .* (it./iterations)); 
        kpi_interr         = (kpi_budget_orig - kpi_budget)                 - (kpi_budget_orig         .* (it./iterations));

        % deriverr
        if it>2,
            spend_differr        = diff( [ [spend(it-1)       spend(it-2)      ] - [(spend_budget_orig       .* ([it-1 it-2]./iterations))]] );
            impressions_differr  = diff( [ [impressions(it-1) impressions(it-2)] - [(impressions_budget_orig .* ([it-1 it-2]./iterations))]] );
            kpi_differr          = diff( [ [kpi(it-1)         kpi(it-2)        ] - [(kpi_budget_orig         .* ([it-1 it-2]./iterations))]] );
        end;
        
        % controller output plus spend the reference line amount
        spend_output       = kprop .*       spend_err + kint .*       spend_interr + kdiff .*       spend_differr + (spend_budget_orig       ./ iterations);
        impressions_output = kprop .* impressions_err + kint .* impressions_interr + kdiff .* impressions_differr + (impressions_budget_orig ./ iterations);
        kpi_output         = kprop .*         kpi_err + kint .*         kpi_interr + kdiff .*         kpi_differr + (kpi_budget_orig         ./ iterations);
        
        % These go into the actuator
        ecpm_calculated(it,1) = 1000.*spend_output ./ impressions_output;
        kpitarg(it,1:kpis)    = kpi_output ./ impressions_output;
        
        % The one last thing for the actuator, is to take the above
        % ecpm_calculated and kpitarg, and then calculate the bid price
        % The bid price then encapsulates our attempt to control
    end;

    if 0, %it==40,
        impid=length(above_imps_list);
        kpiid=length(above_kpi_list);
        break;
    end;

    ii = find(kpitarg>1);if ~isempty(ii), kpitarg(ii) = ones(size(ii));end;


    % Calc KPI error
    % --------------
    
    err = (kpitarg(it,:) ./ kpitargorig(1,:)).^error_squared_power;
    %if err(1) >= 0.8 && err(1) < 1.2, err(1) = err(1) .* scaledown; end;

    
    % KPI Exact error or unidirectional error
    % ----------------------------------------
    
    i = find(err<1.0); % any that are below error %i = setdiff(i,1); % remove the first kpi
    if ~isempty(i),
        if exact==1,
            % Over-performance gets turned into error
            err(i) = 1./err(i); % turn this into error
        else
             % Over-performance (already low error) gets scaled down
             % further
             err(i) = err(i).*scaledown; % de-scale the error
        end;
    end;

    
    % Calc pacing error
    % -----------------
    
    err_pacing = desired_imps_this_cycle(it,1)./(impressions_budget_orig./iterations); % required imps as a ratio of original
    i = find(err_pacing<1);
    if ~isempty(i), err_pacing(i) = 1./err_pacing(i);end;

    
    % Normalize the errors to sum to 1
    % --------------------------------
    
    if nash==1,
        % if nash then all kpis are weighted independently including pacing
        kpiweight(it,:)     = err        ./ (sum(err) + err_pacing);
        pacingweight(it,:)  = err_pacing ./ (sum(err) + err_pacing);
    else
        % pbase; calc weights only within the kpis. Pacing will be applied
        % as a separate term later
        kpiweight(it,:)    = err ./ sum(err);
        pacingweight(it,:) = 0;
    end;
    
    
    % Generate impressions & kpi predictions for this cycle
    % ------------------------------------------------------
    
    kpipred = rand(impressions_per_cycle,kpis);

    % Determine which bin each impression belongs
    for bin=1:bins,
        binid = find(floor(kpipred.*100-eps)+1 == bin);
        distrib(bin) = distrib(bin) + length(binid);
    end;

    
    % Calc Bid price
    % --------------
    
    % perfratio: ratio of impression perf to target perf
    kpiratio = (kpipred ./ (ones(impressions_per_cycle,1)*kpitarg(it,:)) );
    if exact==1,
        % in exact mode, traffic that has a perf greater than target
        % actually results in bids that are proportionately lower
        i = find(kpiratio>1);
        if ~isempty(i), kpiratio(i) = 1./kpiratio(i);end;
    end;
    
    if nash==1,
        kpibid = kpiweight(it,:) .* kpiratio .* ecpm_calculated(it) + pacingweight(it,1) .* bid_for_desired_imps(it); % was ecpm_calculated_orig
%        kpibid = kpiweight(it,:) .* (kpiratio .* ecpm_calculated_orig) + pacingweight(it,1) .* ecpm_calculated(it); % was ecpm_calculated_orig

            if speedup==1,
                if desired_imps_this_cycle(it,1) ./ (impressions_budget_orig./(iterations)) > 1,
                    % catching up
                    bidkpi    = kpiratio .* ecpm_calculated(it); %kpiratio .* bid_for_desired_imps(it);
                    i = find(bidkpi > bid_for_desired_imps(it)); % ie. kpiratio > 1
                else
                    % slowing down
                    bidkpi = kpiratio .* ecpm_calculated(it);
                    i = find(bidkpi < bid_for_desired_imps(it));
                end;
                if ~isempty(i),
                    kpibid(i) = bidkpi(i);
                end;
            end;

    else
%        kpibid = kpiratio .* ecpm_calculated(it); % pbase = kpiratio .* ecpm_calculated(it) as the scalar
        kpibid = bid_for_desired_imps(it) .* kpiratio; % pbase = kpiratio .* ecpm_calculated(it) as the scalar
    end;

    
    %keyboard;
    bid = kpibid; % sum( (ones(impressions_per_cycle,1)*kpiweight(it,:)) .* kpibid ,2 );

    
    % Apply Constraints
    % ------------------
    
    % usedistrib
    if usedistrib,
        prob = distrib./sum(distrib);
        vr = flipud( cumsum( flipud( [prob .* [1:100]'./100] )  ) ./ cumsum( flipud( prob ) ) ); % should be 100 bins
        bin = min(find(vr>=kpitarg(it)));
        if buyaboveorig==1, bin = min(find(vr>=min([kpitarg(it);kpitargorig])));end;
        if isempty(bin),
            % nothing is above the target; the distribution is wholly below
            bid = zeros(size(bid));
        else
            distrib_targ = bin ./ 100;
            ii = find(kpipred < distrib_targ);
            if ~isempty(ii), bid(ii) = zeros(size(ii));end;
        end;
    else

        if basecut==1,
            % basecut
            is_pacing = (desired_imps_this_cycle(it,1) ./ (impressions_budget_orig./(iterations)) < 1);

            if (is_pacing == 1 && step >= 0) || (is_pacing == 0 && step <= 0), 
                step = 1.2 * step;    % go faster in same direction
            else
                step = -0.3 * step;    % change direction and slow down
            end;
            if abs(step) < 0.001,
                step = 0.001 * step / abs(step); % clip step magnitude to 0.001
                % if step is small, step becomes -0.001 or +0.001
            end;
            if cutoff + step > 1,    % clip cutoff to 1
                step = 1 - cutoff;
                cutoff = 1;
            elseif cutoff + step < 0,    % clip cutoff to 0
                step = 0 - cutoff;
                cutoff = 0;
            else
                cutoff = cutoff + step;
            end;
            
            % Apply the constraint
            ii = find(kpipred < cutoff);
            if ~isempty(ii), bid(ii) = zeros(size(ii));end;                        
        end;

        % constraint bid price
        if constraint==1,
            if dynamic == 1,
                ii = find(kpipred < kpitarg(it));
                if buyaboveorig==1, ii = min(find(kpipred < min([kpitarg(it);kpitargorig])));end;
            else
                ii = find(kpipred < kpitargorig);
            end;
            if ~isempty(ii), bid(ii) = zeros(size(ii));end;
        end;
    end;
    

    % Hold the auction and determine which we win
    % --------------------------------------------
     
    % win = find(bid > 5);
    % low winrate (1 in 10 winrate)

    if acceptall==1,
        win = [1:size(bid,1)]';
    else
        win = find(rand(size(bid)).*bid > (1./winrate_for_highbid));    
    end;
    
    if impressions_budget<=0 || spend_budget <= 0,
        win = [];
    end;

    if isempty(win), fprintf('nowin\n');
        spend(it,1)=0;
        impressions(it,1)=0;
        kpi(it,1)=0;
        winrate(it,1) = inf;
        continue;
    end;

    % Induce bid->win model
    % ---------------------

    %keyboard;
    win10 = zeros(size(bid));win10(win) = ones(size(win));
    w = bid\win10;
    %figure(1);clf;plot(bid,win10,'o');hold on;plot(bid,bid.*w,'rx');
    % bid required to produce impressions at the desired delivery rate
    bid_for_desired_imps(it) = desired_imps_this_cycle(it)./(w.*impressions_per_cycle);
    % bid.*w = probwin; bid.*w.*impsavailable = desiredimps; bid = desiredimps/(w.*impsavailable)
    
    
    
    winlen = length(win); % number of winning impressions
    winrate(it) = desired_imps_this_cycle(it)./winlen;

    % Hard adserver shut down
    % -----------------------
    % Simulate a "hard adserver shut-down" - if more impressions are won
    % than the system wants based on pacing, it discards the remainder in this cycle.
    % Hard shut-down is /easy/ for ad-servers since they simply stop buying
    % after the count is reached. So rather than estimating impressions all
    % the time, at least in terms of buying too many, this mechanism is
    % usually enabled in practice.
    if adserver_shutdown,
        if length(win) > desired_imps_this_cycle(it) .* buy_imps_overage,
    %        fprintf('exceeded\n');
            if desired_imps_this_cycle(it)<1,
                win = [];
            else
    %            [a,i] = sort(rand(win));win = win(floor(desired_imps_this_cycle(it)));
                [a,i] = sort(rand(length(win),1));
                win = win(i(1:min([floor(desired_imps_this_cycle(it).*buy_imps_overage) impressions_per_cycle])));
            end;
        end;
    end;
        
    % Update spend and KPIs
    % ---------------------
    
    spend(it,1) = sum(bid(win)./1000).*0.1; % ****** WHY THE 0.1????
    impressions(it,1) = length(win);
    kpiactual = kpipred(win,1:kpis) + rand(length(win),1).*prediction_bias + randn(length(win),1).*prediction_noise;
    ii = find(kpiactual<0);if ~isempty(ii),kpiactual(ii)=zeros(length(ii),1);end;
    ii = find(kpiactual>1);if ~isempty(ii),kpiactual(ii) =ones(length(ii),1);end;
    kpi(it,1:kpis) = sum(kpiactual,1); % put some noise on here
    
    % **** REGRESSION DOWNWARDS NEEDED HERE

    kpi_budget         = kpi_budget_orig - sum(kpi(1:it,1:kpis));
    spend_budget       = spend_budget_orig - sum(spend(1:it,1));
    impressions_budget = impressions_budget_orig - sum(impressions(1:it,1));


    % Hard budget shut-down
    % ----------------------
    % if we reach zero this iteration, stop it at zero
    % This is a "visual aid" that prevents negative kpis and impressions - after running out of budget, on the last step
    % we can go negative on impressions etc. However if we zero that out, then it creates a "kink"
    % in the trajectory line, and so visually can make it harder to read.
    if zero_last_step,
%        if impressions_budget<0 | spend_budget<0 | kpi_budget<0, keyboard;end;
        if impressions_budget < 0, impressions(it,1) = impressions(it,1) + impressions_budget; impressions_budget = impressions_budget + -impressions_budget; end;
        if spend_budget       < 0, spend(it,1)       = spend(it,1)+        spend_budget;       spend_budget       = 0; end;
        if kpi_budget         < 0, kpi(it,1)         = kpi(it,1)+          kpi_budget;         kpi_budget         = 0; end;
    end;
    
    if showgraphs,
        fprintf('%d\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t\n',it,ecpm_calculated(it),length(win),spend(it,1),impressions(it,1),spend_budget,impressions_budget, kpi(it), kpi_budget,winrate(it),winlen);

        figure(1);clf;plot(ecpm_calculated);title('ecpm calculated');
        figure(2);clf;plot(kpitarg);title('kpitarg');
        figure(3);clf;plot(cumsum(spend));title('cum spend');hold on;plot(it,spend_budget_orig,'o');
        figure(4);clf;plot(cumsum(impressions));title('cum impressions');hold on;plot(it,impressions_budget_orig,'o');
        figure(5);clf;plot(cumsum(kpi));title('cum kpi');hold on;plot(it,kpi_budget_orig,'o');
        %figure(6);clf;plot(cumsum(kpi)./cumsum(impressions),cumsum(impressions)./( ([1:it]./iterations) .* impressions_budget),'-');
        figure(6);clf;plot(cumsum(kpi)./cumsum(impressions),cumsum(impressions)./( ([1:it]'.*impressions_desired_per_cycle)),'-');
        %figure(7);clf;plot(cumsum(impressions),'-');hold on;plot(cumsum(desired_imps_this_cycle),'r-');hold on;plot(([1:it]./iterations).*impressions_buget_orig,'k:');title('imps');legend('imps','imps desired');

        figure(7);clf;plot(cumsum(impressions),'-');hold on;plot(([1:it]./iterations).*impressions_budget_orig,'x:');title('imps');%legend('imps realized','imps orig desired');grid on;
%        grid on; %axis('equal');
%        v = axis;
%        axis([  min([v(1);v(3)]) max([v(2);v(4)]) min([v(1);v(3)]) max([v(2);v(4)]) ] );
        
        figure(8);clf;plot((cumsum(kpi)./cumsum(impressions)) ./ kpitargorig, cumsum(impressions)./cumsum(ones(it,1).*(impressions_budget_orig./iterations)),'r-');xlabel('kpi');ylabel('pacing');
        hold on;plot( (sum(kpi)./sum(impressions)) ./ kpitargorig, sum(impressions(1:it))./((impressions_budget_orig./iterations) .* it),'o'); grid on;
        axis('equal');
        v = axis;
        %axis([  min([v(1);v(3)]) max([v(2);v(4)]) min([v(1);v(3)]) max([v(2);v(4)]) ] );
        axis([min([vglobal(1);v(1)]) max([vglobal(2);v(2)]) min([vglobal(3);v(3)]) max([vglobal(4);v(4)])]);
        
        figure(9);hold on;plot(it,cutoff,'o');
        figure(10);hold on;plot(it,step,'x');
    end;

    
    % If we ran out of budget then stop
    % ----------------------------------
    
    %if impressions_budget<=0 || spend_budget <= 0, break; end;

    %    if it == 67, break; end;
%    if mod(it,100) == 0, pause; end;
end;

if ~isempty(kpi_save) && length(kpi)~= size(kpi_save,1), fprintf('uneven kpi len\t%d\tkpi_save len\t%d\n',length(kpi),size(kpi_save,1)); keyboard; continue; end;

kpi_save         = [kpi_save kpi];
impressions_save = [impressions_save impressions];
spend_save       = [spend_save spend];

if onetime, return; end;

end;end;


% Phase portraits
% -------------------------------------------------------------------------

% Thoughts: Adserver shut-down: This creates a convergence to goal - de-scales each
% cycle; so this is equivalent to the standard controller
% Without the above, shouldn't pbase throttle down? reduce the ecpm?

%i  = find(kpi_save==0);kpi_save(i)=nan.*ones(size(i));
%i  = find(impressions_save == 0); impressions_save(i)=nan.*ones(size(i));
%i  = find(spend_save       == 0); spend_save(i) = nan.*ones(size(i));
%nz = sum(isnan(kpi_save));

% Extract trajectories
% --------------------

% nan out the initial perturbation; remaining trajectories are under control system
kpitrajectory  = (cumsum(kpi_save,1)  ./cumsum(impressions_save,1)) ./ kpitargorig;
impstrajectory = cumsum(impressions_save,1) ./ cumsum(  ones(size(impressions_save)) .*(impressions_budget_orig./iterations),1);
if starting_it>1,
    kpitrajectory( 1:starting_it,:) = nan.*ones(starting_it,size(kpitrajectory,2));
    impstrajectory(1:starting_it,:) = nan.*ones(starting_it,size(kpitrajectory,2));
end;    

% Extract endpoints
% ------------------

kpiendpoint = (nansum(kpi_save,1)  ./nansum(impressions_save,1)) ./ kpitargorig;
impendpoint = nansum(impressions_save,1) ./ nansum(  ones(size(impressions_save)) .*(impressions_budget_orig./iterations),1);


% Summary stats on performance: Means and stdevs
% -----------------------------------------------

meanx    = mean([kpiendpoint' impendpoint']);
stdx     = std([kpiendpoint' impendpoint']);
medianx  = median([kpiendpoint' impendpoint']);

% Store the stats
paramindex = sum([usedistrib.*2.^0   constraint.*2.^1   dynamic.*2.^2   nash.*2.^3   pid.*2.^4   exact.*2.^5   acceptall.*2.^6   adserver_shutdown.*2.^7]) + 1;
paramstats = [meanx medianx stdx];
stats(paramindex,1:length(paramstats)) = paramstats;

% Report on the stats (all previous runs)
% Note that this will change usedistrib etc

stats_list = [];
paramstr_list = [];
i = find(full(stats(:,1))~=0); % find all stats entries
fprintf('usedistrib\tconstraint\tdynamic\tnash\tpid\texact\tacceptall\tadserver_shutdown\tkpimean\timpmean\tkpimed\timpmed\tkpistd\timpstd\n');
for ii=1:length(i),

        % Extract the stats
        paramstats = full(stats(i(ii),1:length(paramstats)));

        % Extract the key
        paramindex = i(ii)-1; % remove the offset
        
        b = floor(paramindex / 2.^7);
        adserver_shutdown = b;
        paramindex = paramindex-b.*2.^7;
        
        b = floor(paramindex / 2.^6);
        acceptall = b;
        paramindex = paramindex-b.*2.^6;

        b = floor(paramindex / 2.^5);
        exact = b;
        paramindex = paramindex-b.*2.^5;

        b = floor(paramindex / 2.^4);
        pid = b;
        paramindex = paramindex-b.*2.^4;

        b = floor(paramindex / 2.^3);
        nash = b;
        paramindex = paramindex-b.*2.^3;

        b = floor(paramindex / 2.^2);
        dynamic = b;
        paramindex = paramindex-b.*2.^2;

        b = floor(paramindex / 2.^1);
        constraint = b;
        paramindex = paramindex-b.*2.^1;

        b = floor(paramindex / 2.^0);
        usedistrib = b;
        paramindex = paramindex-b.*2.^0;

        stats_list = [stats_list;paramstats];
        paramstr = '';
        if usedistrib, paramstr = strcat(paramstr,'distr ');    end;
        %if constraint, paramstr = strcat(paramstr,'constr ');     end;
        %if dynamic,    paramstr = strcat(paramstr,'dynamic ');    end;

        if constraint==1 & dynamic==1, paramstr = 'dyn'; end;
        if constraint==1 & dynamic==0, paramstr = 'hard'; end;

        if nash,       paramstr = strcat(paramstr,'nash ');       end;
        if pid,        paramstr = strcat(paramstr,'pid ');        end;
        if exact,      paramstr = strcat(paramstr,'exact ');      end;
        if acceptall,  paramstr = strcat(paramstr,'all ');  end;
        if adserver_shutdown, paramstr = strcat(paramstr,'shtdwn ');end;
        if sum([usedistrib constraint dynamic,nash,pid,exact,acceptall,adserver_shutdown]')==0,
             paramstr = strcat(paramstr,'px ');
        end;
        paramstr_list{ii} = paramstr;        
        fprintf('%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t',usedistrib,constraint,dynamic,nash,pid,exact,acceptall,adserver_shutdown);
        fprintf('%f\t',paramstats);
        fprintf('\n');
end;


% Show targetplot (bullseye plot)
% ----------------------------------

figure(9);clf;hold on;
% higher version of matlab for these
%ax = gca;
%ax.YRuler.FirstCrossoverValue = 1;
%ax.XRuler.FirstCrossoverValue =1;
%ax.XAxisLocation = 'origin';
%ax.YAxisLocation = 'origin';

% plot circles on axis
radii = 50;
for i=1:6,
    plot(1,1,'ko','markersize',radii.*i);
end;
grid on;
axis(vglobal);
plot(stats_list(:,1),stats_list(:,2),'r.','markersize',20);
text(stats_list(:,1)+.01,stats_list(:,2)-.05,paramstr_list);
plot([vglobal(1) vglobal(2)]',[1 1]','k-');
plot([1 1]',[vglobal(3) vglobal(4)]','k-');
axis('square');
xlabel('kpi');
ylabel('pacing');

%plot(1,1,'o','markersize',100);plot(1,1,'o','markersize',150);grid on


% Phase portrait
% --------------

        figure(10);clf;plot(kpitrajectory, impstrajectory ,'r-');xlabel('kpi');ylabel('pacing');
        hold on;
        plot(kpiendpoint,   impendpoint    ,'o'); grid on; % impressions_budget_orig
        axis('equal');
        axis(vglobal);
xlabel('kpi');
ylabel('pacing');
        %v = axis;
        %axis([  min([v(1);v(3)]) max([v(2);v(4)]) min([v(1);v(3)]) max([v(2);v(4)]) ] )  

        % Points
        % --------
                
        figure(11);clf;
        hold on;
        plot(kpiendpoint,   impendpoint    ,'o'); xlabel('kpi');ylabel('pacing'); grid on; % impressions_budget_orig
        axis('equal');
        axis(vglobal);
xlabel('kpi');
ylabel('pacing');

        % Convhull
        % --------

        k = convhull(kpiendpoint,impendpoint,{'QJ','Pp'}); %%{'Qt','Pp'});
        figure(12);clf;
        hold on;%plot(kpiendpoint,   impendpoint    ,'o');
        xlabel('kpi');ylabel('pacing'); grid on; % impressions_budget_orig
        axis('equal');
        axis(vglobal);
        plot(kpiendpoint(k(:,1)),impendpoint(k(:,1)),'b-');
        plot(meanx(:,1),meanx(:,2),'ro','markersize',10);
xlabel('kpi');
ylabel('pacing');
        
        
        % Heatmap
        % --------
%        figure(10);clf;%plot((cumsum(kpi_save)./cumsum(impressions_save)) ./ (ones(size(kpi_save)).*kpitargorig), cumsum(impressions_save)./cumsum(ones(size(kpi_save,1),size(kpi_save,2)).*(impressions_budget_orig./iterations)),'r-');xlabel('kpi');ylabel('pacing');
 %       hold on;plot( (sum(kpi_save)./sum(impressions_save)) ./ kpitargorig, sum(impressions_save)./((impressions_budget_orig./iterations) .* it),'o'); grid on;
  %      axis('equal');
   %     axis(vglobal);

   x = [];
   y = [];
   z = [];
   
above_imps_list = [.5:.1:1.5];
above_kpi_list  = [.5:.1:1.5];

 for impid = 1:length(above_imps_list),
    for kpiid = 1:length(above_kpi_list),
        above_imps = above_imps_list(impid);
        above_kpi  = above_kpi_list(kpiid);
    
        % only plot exterior
        %if exterior_plot_only,
        x = [x;ones(10,1).*above_kpi ones(10,1).*above_imps]; z = [z;ones(10,1).*0];
            if above_imps >= max(above_imps_list) || above_kpi <=min(above_kpi_list) || ...
               above_imps <= min(above_imps_list) || above_kpi >=max(above_kpi_list),
               % exterior; keep going
        %    else
        %      continue;
        %      % ignore this
        %    end;
             end;
        
        x = [x;above_kpi above_imps]; z = [z;0];

    end; % kpiid
 end; % impid

   % add the endpoints
   x = [x; kpiendpoint' impendpoint'];
   z = [z; ones(length(kpiendpoint),1) ./ length(kpiendpoint)];
   
   figure(13);clf;
   % [minx,stepsize,zi,maxx] = showresults([kpiendpoint' impendpoint'],ones(length(kpiendpoint),1),10);
   [minx,stepsize,zi,maxx] = showresults_smooth(x,z,13);
   xlabel('kpi');ylabel('pacing');zlabel('prob');
   % we want final ending points
%   figure(10);clf;
   

        % x,z with the above
        
        % velocity vectors
        x = [(cumsum(kpi_save)./cumsum(impressions_save)) ./ (ones(size(kpi_save)).*kpitargorig)];
        y = [cumsum(impressions_save)./cumsum(ones(size(kpi_save,1),size(kpi_save,2)).*(impressions_budget_orig./iterations))];
        xstart = x(1,:);
        xend   = x(size(x,1),:);
        ystart = y(1,:);
        yend   = y(size(y,1),:);
        
        figure(14);clf;quiver(xstart,ystart,(xend-xstart),(yend-ystart));
        hold on;%plot(kpiendpoint,   impendpoint    ,'o');
        xlabel('kpi');ylabel('pacing'); %grid on; % impressions_budget_orig
        axis('equal');
        axis(vglobal);
        plot(kpiendpoint(k(:,1)),impendpoint(k(:,1)),'r-')
        plot(kpiendpoint,impendpoint,'r.');
        plot(meanx(:,1),meanx(:,2),'ks','markersize',20);
        plot(1,1,'k^','markersize',20);
        %plot([1 meanx(:,1)]',[1 meanx(:,2)]','k-');
        %plot([vglobal(1) vglobal(2)]',[1 1]','k:');
        %plot([1 1]',[vglobal(3) vglobal(4)]','k:');
        
        if 0,
            figure(15);clf;quiver(xstart,ystart,(mean(x(2:10,:),1)-xstart),(mean(y(2:10,:),1)-ystart));
            hold on;%plot(kpiendpoint,   impendpoint    ,'o');
            xlabel('kpi');ylabel('pacing'); %grid on; % impressions_budget_orig
            axis('equal');
            axis(vglobal);
            plot(kpiendpoint(k(:,1)),impendpoint(k(:,1)),'b-')
            plot(kpiendpoint,   impendpoint    ,'r.');
        end;
        
[meanx meanx-1;medianx medianx-1;stdx stdx]
        

        if 0,
          budget=100;targ=.1;pred=1;imps=1000;figure(2);clf;plot([1:imps],budget - budget.*(pred./targ).*log([1:imps]));

          
          bud = [];budget=100;imps = 1000;targ = .1; pred=1.0; while 1,bid = (pred./targ) .* (budget./imps); imps = imps - 1; budget=budget - bid; bud = [bud;budget]; figure(1);clf;plot(bud);if imps < 0,break;end;end; %end;
          
          figure(2);plot([0:1000],100.*exp(-0.1.*[0.001].*[0:1000]))
          figure(3);clf;plot([0:999],100.*exp(-1.*[0:999]./[1000-[0:999]]),'rx-');grid on;v=axis;axis([v(1) v(2) 0 100]);
          
        end;
        


% PID Controller experiment
% -------------------------

if 0,

kp = .5;ki = .1;kd = .1;
previous_error=0;integral=0;
output=0;
dt=1;
T=50;B=100;
c = 1.5;
integral = 0; derivative=0;error=0;
figure(1);clf;hold on;
for t=1:T,
    setpoint(t) = (B/T); %+rand(1,1).*0.1; % output from the plant that we'd like to see this t
    measured_value(t) = c .* (B/T) + output(t); % this is the actual progress
    %budget(t) = c.*(B/T) + 
    error(t+1)      = setpoint(t) - measured_value(t);
    integral(t+1)   = integral(t) + error(t) .* dt;
    derivative(t+1) = (error(t+1)-previous_error)./dt;
    output(t+1) = kp .* error(t) + ki .* integral(t+1) + kd .* derivative(t+1); % controller output
    previous_error = error(t);
    %fprintf('plant output\t%f\n',measured_value);
    %plot(t,setpoint,'rx');plot(t,measured_value,'bo');
    %pause;
end;
figure(2);clf;hold on;plot([1:T],setpoint([1:T]),'rx-');hold on;plot([1:T],measured_value([1:T]),'bo-');
legend('setpoint/reference signal','plant output');

end; % if 0



if 0,
    % symmetric error on pacing
    err_squared_power=1;r=0.1;elow=0.1;ehigh = 1./(1-elow);elow = 1-elow;kpiabove = [0.5:0.01:2]';kpiaboveorig = kpiabove;i=find(kpiabove<1);kpiabove(i)=1./kpiabove(i);i=find(kpiabove >= elow & kpiabove<=ehigh);kpiabove(i)=kpiabove(i).*r;kpiabove = kpiabove.^err_squared_power;figure(1);clf;plot(kpiaboveorig,kpiabove,'-');set(gca,'xscale','log');grid on

    % asymmetric error on kpi
end;


% Moat advertisers 6 months
%figure(1);clf;[c,h]=contourf(dat');colorbar;colormap(cool);set(gca,'xticklabel',num2str([0:0.2:2]'));set(gca,'yticklabel',num2str([0:0.2:2]'));grid on;xlabel('kpi/targ');ylabel('pacing/targ');hold on;plot(.25+means(2,1).*6,.25+means(1,1).*6,'s','markersize',20,'markerfacecolor','auto');plot(6,6,'^','markersize',20,'markerfacecolor','auto');%for i=100:150:300,plot(6,6,'ko','markersize',i);end;%clabel(c,h);
%figure(2);clf;[c,h]=contourf(dat2');colorbar;colormap(cool);set(gca,'xticklabel',num2str([0:0.2:2]'));set(gca,'yticklabel',num2str([0:0.2:2]'));grid on;xlabel('kpi/targ');ylabel('pacing/targ');hold on;plot(means(2,2).*6,.25+means(1,2).*6,'s','markersize',20,'markerfacecolor','auto');plot(6,6,'^','markersize',20,'markerfacecolor','auto');% for i=100:150:300,plot(6,6,'ko','markersize',i);end; %clabel(c,h);

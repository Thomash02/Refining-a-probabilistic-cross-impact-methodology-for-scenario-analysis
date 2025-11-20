% By running this module one can obtain results for the different
% simulation runs

%% Plot CDFs for Odds and Original
figure;
x = (1:length(CDF_Odds));
plot(x, CDF_Odds, 'b', x, CDF_Original, 'r');
hold on;
    
% Fill the area between the two curves
fill([x, fliplr(x)], [CDF_Odds, fliplr(CDF_Original)], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    
title('ECDFs of both methods for run No. 100');
legend('Odds-ratio interpretation', 'Original interpretation','Area between',Location='Best');
xlabel('Scenario Index');
ylabel('Cumulative Probability');
xlim([1,177147])
grid minor

%% Statistical tests (Table 5)

% Example 1: 10%, 20%, ..., 100%  (Table 5)
percScenarios = 0.10:0.10:1.00;
[statsTable, meanET] = compute_test_result_table( ...
    CDF_Original_all, CDF_Odds_all, ET_run, percScenarios, ...
    Joint_pd_Original, Joint_pd_Odds);

disp(statsTable)
fprintf('Mean ET (speed ratio): %.4f\n', meanET);
%% Statistical tests (Table 6)
% Example 2: 1%, 2%, ..., 10%  (Table 6)
percScenarios_small = 0.01:0.01:0.10;
[statsTable_small, meanET_small] = compute_test_result_table( ...
    CDF_Original_all, CDF_Odds_all, ET_run, percScenarios_small, ...
    Joint_pd_Original, Joint_pd_Odds);

disp(statsTable_small)
fprintf('Mean ET (speed ratio): %.4f\n', meanET);
%%
function [statsTable, meanET] = compute_test_result_table( ...
                            CDF_Original_all, CDF_Odds_all, ...
                            ET_run, percScenarios, Joint_pd_Original, ...
                            Joint_pd_Odds)
% COMPUTE_TEST_RESULT_TABLE
%   Build the test-result table for different percentages of scenarios.
%
%   Inputs
%     CDF_Original_all : [N x m] matrix, ECDF of original method for each run
%     CDF_Odds_all     : [N x m] matrix, ECDF of odds method for each run
%     ET_run           : [1 x m] mean runtime ratio per run (ETOD./ETO)
%     percScenarios    : vector of fractions, e.g. 0.10:0.10:1.00
%
%   Outputs
%     statsTable : MATLAB table with mean JS, TV, KS p, KS rejection %
%     meanET     : overall mean performance ratio across runs

    [nScen, m] = size(CDF_Original_all);
    numThresh  = numel(percScenarios);

    % Preallocate: (threshold x runs)
    JS   = zeros(numThresh, m);
    TV   = zeros(numThresh, m);
    KSp  = zeros(numThresh, m);
    KSre = zeros(numThresh, m);   % 0/1 reject indicator from kstest2

    for H = 1:m
        Original_cdf = CDF_Original_all(:, H);
        Odds_cdf     = CDF_Odds_all(:, H);

        Original = Joint_pd_Original(:, H);
        Odds = Joint_pd_Odds(:, H);

        for t = 1:numThresh
            frac        = percScenarios(t);               % e.g. 0.10
            cutoffIndex = max(1, floor(frac * nScen));    % #scenarios

            % Truncate ECDFs to the first frac% most probable scenarios
            Org_cut  = Original(1:cutoffIndex);
            Odds_cut = Odds(1:cutoffIndex);

            % Jensen–Shannon divergence on truncated ECDFs
            M = 0.5 * (Org_cut + Odds_cut);
            %epsVal = 1e-15;               % avoid log(0)
            %P = max(Org_cut,  epsVal);
            %Q = max(Odds_cut, epsVal);
            %M = max(M,       epsVal);
            P = Org_cut;
            Q = Odds_cut;

            JS(t, H) = 0.5 * (sum(P .* log(P ./ M)) + ...
                               sum(Q .* log(Q ./ M)));

            % Total variation distance
            TV(t, H) = 0.5 * sum(abs(Org_cut - Odds_cut));

            % KS test on truncated ECDFs
            [h, p]   = kstest2(Original_cdf, Odds_cdf);
            KSre(t,H) = h;
            KSp(t,H)  = p;
        end
    end

    % Mean over runs
    mean_JS   = mean(JS,  2);
    mean_TV   = mean(TV,  2);
    mean_KSp  = mean(KSp, 2);
    rejRate   = mean(KSre,2) * 100;   % % of runs rejecting H0

    meanET    = mean(ET_run);         % overall mean speed ratio

    % Build table (percentages instead of fractions)
    pct = percScenarios(:) * 100;

    statsTable = table(pct, mean_JS, mean_TV, mean_KSp, rejRate, ...
        'VariableNames', { ...
            'Percentage_of_scenarios', ...
            'JS_Divergence', ...
            'TV_Distance', ...
            'KS_test_p_value', ...
            'KS_test_rejections_percent' ...
        });
end

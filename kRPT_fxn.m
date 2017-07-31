% This function runs an ensemble of 1-D particle tracking simulations for
% the chemical reaction A + B --> C, where particles of C are ignored
% This function is driven by the script run_kRPT.m
function [totmass, runtimes] = kRPT_fxn(omega, k, D, NA, NB, C0, ...
    Avar, Bvar, maxtime, ens)

    % initialize varaibles and storage arrays
    initmass  =  C0 * omega;
    A  =  zeros(NA, 3);
    B  =  zeros(NB, 3);
    mindx  =  min(omega/NA, omega/NB);

    % array for plotting--only store values for certain time steps
    totmass  =  zeros(50,3);
    totmass(:, 1) = logspace(log10(0.1/k/C0), log10(maxtime/k/C0), 50);
    dt = min([.1/k/C0  4 * mindx^2/D]);

    % determine max colocation probability and print max distance a particle
    % can move in a time step for checking purposes
    c1max = 1/sqrt(2 * pi * (Avar + Bvar + 4 * D * dt));
    maxdm = (C0 * omega/NA) * k * dt * c1max

    % loop counter--prints to screen on each iteration
    loop = 0;

    % run time tracker
    runtimes = zeros(ens, 1);

    %  ensemble loop
    for ll = 1 : ens
        tic

        % print current loop to screen
        loop = loop + 1

        % particle positions    
        A(:, 1) = omega * rand(NA, 1);
        B(:, 1) = omega * rand(NB, 1);

        % initialize particle masses
        A(:, 2) = (C0 * omega/NA);
        B(:, 2) = (C0 * omega/NB);

        time = 0;
        kk = 1;

        while time < maxtime/k/C0

            % calculate cutoff distance for rangesearch
            var = Avar + Bvar + 4 * D * dt; 
            dist = 4 * sqrt(var);

            BS = max(50, NA/100);

            % perform rangesearch for nearest neighbors
            [idx, r] = rangesearch(B(:, 1), A(:, 1), dist, 'BucketSize', BS);

            for jj = 1 : NA     % A particle loop 

                ma = A(jj, 2);

                closeguys = idx{jj};     % This is the indices list of nearby particles
                closedists = r{jj};      % Associated radii

                if(length(closeguys) < 1) 
                    continue 
                end 

                for ii = 1 : length(closeguys)  % B particle loop
                    pidx = closeguys(ii);
                    % encounter probability
                    v_s = (1/sqrt(2 * pi * var)) *exp(closedists(ii)^2/(-2 * var));
                    % mass differential
                    dm = k * dt * ma * B(pidx, 2) * v_s;
                    % adjust masses
                    B(pidx, 2) = B(pidx, 2) - dm;
                    ma = ma - dm;
                end

                A(jj,2) = ma; 

            end   %  A particle loop end

            % Move the particles
            A(:, 1) = A(:, 1) + sqrt(2 * D * dt) * randn(NA, 1);
            B(:, 1) = B(:, 1) + sqrt(2 * D * dt) * randn(NB, 1);
            % bring particles that exit through boundary back into domain
            apple = find(A(:, 1) < 0);
            A(apple, 1) = -1 * A(apple, 1);
            apple = find(A(:, 1) > omega);
            A(apple, 1) = 2 * omega - A(apple, 1);
            apple = find(B(:, 1) < 0);
            B(apple, 1) = -1 * B(apple, 1);
            apple = find(B(:, 1) > omega);
            B(apple, 1) = 2 * omega - B(apple, 1);

            % update time
            time = time + dt;

            % store total mass of each species (normalized by initial concentration) and
            % current time near selected time steps
            % if ensemble is larger than one, calculate standard deviation for errorbar plot
            if(time >= totmass(kk, 1))
                totmass(kk, 1) = time; 
                fracmass = sum(A(:, 2))/initmass;
                incstore = fracmass - totmass(kk, 2);
                totmass(kk, 2) = ((ll - 1) * totmass(kk, 2) + fracmass)/ll;  % for the ensemble average
                if ll > 1
                    totmass(kk, 3) =  totmass(kk, 3) + incstore * (fracmass - totmass(kk, 2));
                end          
                kk = kk + 1;
            end

        end  % time loop

        % analytic, well-mixed solution, for reference
        analytic = 1./(1 + C0 * k * totmass(:, 1));
        figure(1) 
        loglog(totmass(:, 1), analytic);
        axis([0.1/k/C0 5000/k/C0 10^-3 1]);
        hold on
        errorbar(totmass(:, 1), totmass(:, 2), sqrt(totmass(:, 3)/ll), 'Marker', 'd', 'Color',...
            'r', 'LineStyle', 'none')
        hold off
        drawnow

        runtimes(loop) = toc

    end    % end the ensemble loop
end

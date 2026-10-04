M = 10;
c = 11;
d = 0.5 + 0.05*(0:M-1)';
phi0 = [-1.2;-.9;-.6;-.3;0;.3;.6;.9;1.2;1.5];
omega0 = [.6;.4;.2;0;-.2;-.4;-.6;.5;.1;-.3];
Nstep = round(Tsim/hEM);
Ncase = length(case_names);
all_metrics = zeros(Nmc,9,Ncase);
schedule_stats = zeros(Ncase,7);
case_elapsed = zeros(Ncase,1);
brownian_block_steps = 2000;
for icase = 1:Ncase
    pars = case_par(icase,:);
    strategy = pars(1);
    schedule_kind = pars(2);
    hs = pars(3);
    sigma1 = pars(4);
    sigma0 = pars(5);
    duty = pars(6);
    noise = pars(7);
    np = pars(8);
    gain = pars(9);
    topology = pars(10);
    bias_amp = pars(11);
    aperAmp = pars(12);
    sample_steps = max(1,round(hs/hEM));
    if abs(sample_steps*hEM-hs) > 1e-12
        error('Case %s: hs must be an integer multiple of hEM.',case_names{icase});
    end
    W = ones(M)-M*eye(M);
    pins = 1:np;
    if topology >= 2
        W = zeros(M);
        W(1,2:M) = 1;
        W(2:M,1) = 1;
        for inode = 2:M-1
            W(inode,inode+1) = 1;
        end
        W(M,2) = 1;
        W = W - diag(sum(W,2));
        if topology == 2
            pins = 1:np;
        elseif topology == 3
            if np+1 > M
                error('Leaf-first pinning requires np <= M-1.');
            end
            pins = 2:(np+1);
        end
    end
    kappa = zeros(M,1);
    kappa(pins) = gain;
    kappa_matrix = repmat(kappa,1,Nmc);
    damping_matrix = repmat(d,1,Nmc);
    bias = bias_amp*(.2+.8*(0:M-1)'/(M-1));
    bias_matrix = repmat(bias,1,Nmc);
    active_mask = false(1,Nstep);
    cycle_start = zeros(1,Nstep);
    boundaries = zeros(0,3);
    left = 0;
    mcycle = 0;
    while left < Nstep
        if schedule_kind == 1
            width = round(.05/hEM);
        else
            width = 4*round(.05*(1+aperAmp*sin(sqrt(2)*mcycle))/(4*hEM));
            width = max(4,width);
        end
        width = min(width,Nstep-left);
        right = left + width;
        on_width = round(duty*width);
        on_width = max(0,min(on_width,width));
        switch_off = left + on_width;
        if on_width > 0
            active_mask(left+1:switch_off) = true;
        end
        if width > 0
            cycle_start(left+1:right) = left;
        end
        boundaries(end+1,:) = [left*hEM switch_off*hEM right*hEM];
        left = right;
        mcycle = mcycle + 1;
    end
    if strategy == 1
        active_mask(:) = true;
        cycle_start(:) = 0;
    end
    achieved_active_fraction = mean(active_mask);
    tgrid = (0:Nstep)*hEM;
    active_time = [0 cumsum(active_mask)*hEM];
    D = duty*tgrid - active_time;
    running_min = inf;
    max_deficit = 0;
    for kk = 1:length(D)
        running_min = min(running_min,D(kk));
        max_deficit = max(max_deficit,D(kk)-running_min);
    end
    active_widths = boundaries(:,2)-boundaries(:,1);
    rest_widths = boundaries(:,3)-boundaries(:,2);
    schedule_stats(icase,:) = [ ...
        achieved_active_fraction, ...
        max_deficit, ...
        size(boundaries,1), ...
        min(active_widths), max(active_widths), ...
        min(rest_widths), max(rest_widths)];
    phi = repmat(phi0,1,Nmc);
    omega = repmat(omega0,1,Nmc);
    phiEvent = phi;
    omegaEvent = omega;
    zero_input = zeros(M,Nmc);
    Nu = zeros(1,Nmc);
    integral_error = zeros(1,Nmc);
    integral_energy = zeros(1,Nmc);
    integral_tail = zeros(1,Nmc);
    maximum_error = sum(phi.^2+omega.^2,1);
    Nsample = 0;
    tail_first = round(.8*Nstep);
    try
        rng(seed,'twister');
    catch
        randn('state',seed);
    end
    tic;
    block_start = 1;
    while block_start <= Nstep
        nb = min(brownian_block_steps,Nstep-block_start+1);
        dB_block = sqrt(hEM)*randn(Nmc,nb);
        for jb = 1:nb
            k = block_start + jb - 1;
            is_active = active_mask(k);
            is_start = (k==1);
            if k > 1
                is_start = is_active && ~active_mask(k-1);
            end
            if is_active
                is_sample = mod((k-1)-cycle_start(k),sample_steps)==0;
                if is_start || is_sample
                    Nsample = Nsample + 1;
                    rho2 = sum(phi.^2+omega.^2,1);
                    eta2 = sum((phi-phiEvent).^2+(omega-omegaEvent).^2,1);
                    if is_start || strategy==2
                        update = true(1,Nmc);
                    else
                        update = eta2 > sigma1*rho2 + sigma0;
                    end
                    phiEvent(:,update) = phi(:,update);
                    omegaEvent(:,update) = omega(:,update);
                    Nu = Nu + update;
                end
                uPhi = -kappa_matrix.*phiEvent;
                uOmega = -kappa_matrix.*omegaEvent;
            else
                uPhi = zero_input;
                uOmega = zero_input;
            end
            old_error = sum(phi.^2+omega.^2,1);
            energy = sum(uPhi.^2+uOmega.^2,1);
            driftPhi = omega + c*(W*phi) + uPhi;
            driftOmega = -sin(phi) - damping_matrix.*omega + ...
                c*(W*omega) + bias_matrix + uOmega;
            dB_step = dB_block(:,jb)';
            phi = phi + hEM*driftPhi;
            omega = omega + hEM*driftOmega + ...
                noise*omega.*repmat(dB_step,M,1);
            new_error = sum(phi.^2+omega.^2,1);
            integral_error = integral_error + .5*hEM*(old_error+new_error);
            integral_energy = integral_energy + hEM*energy;
            if k > tail_first
                integral_tail = integral_tail + .5*hEM*(old_error+new_error);
            end
            maximum_error = max(maximum_error,new_error);
        end
        clear dB_block;
        block_start = block_start + nb;
    end
    active_fraction = mean(active_mask);
    all_metrics(:,:,icase) = [ ...
        Nu', ...
        integral_error'/Tsim, ...
        integral_energy'/Tsim, ...
        integral_tail'/((Nstep-tail_first)*hEM), ...
        new_error', ...
        maximum_error', ...
        np*Nu', ...
        active_fraction*ones(Nmc,1), ...
        Nsample*ones(Nmc,1)];
    case_elapsed(icase) = toc;
    fprintf('%s: Nu=%.3f, J=%.9g, E=%.9g, tail=%.9g, maxPathMax=%.9g, active=%.6f, elapsed=%.1fs\n', ...
        case_names{icase}, ...
        mean(Nu), ...
        mean(integral_error)/Tsim, ...
        mean(integral_energy)/Tsim, ...
        mean(integral_tail)/((Nstep-tail_first)*hEM), ...
        max(maximum_error), ...
        active_fraction, ...
        case_elapsed(icase));
end

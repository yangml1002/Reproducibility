M = 10;
c = 11;
d = 0.5 + 0.05*(0:M-1)';
phi0 = [-1.2;-.9;-.6;-.3;0;.3;.6;.9;1.2;1.5];
omega0 = [.6;.4;.2;0;-.2;-.4;-.6;.5;.1;-.3];
Nstep = round(Tsim/hEM);
if abs(Nstep*hEM-Tsim)>1e-12
    error('Tsim must be an integer multiple of hEM.');
end
stride = max(1,round(.001/hEM));
record_steps = unique([0:stride:Nstep Nstep]);
curve_time = record_steps*hEM;
curves = zeros(length(case_names),length(record_steps));
all_metrics = zeros(Nmc,9,length(case_names));
schedule_records = cell(length(case_names),1);
state_example = cell(length(case_names),1);
event_example = cell(length(case_names),1);
try
    rng(seed,'twister');
catch
    randn('state',seed);
end
if exist('refine_bridge','var') && refine_bridge
    if mod(Nstep,2)~=0
        error('Nstep must be even for the bridge refinement check.');
    end
    dB_coarse = sqrt(2*hEM)*randn(Nstep/2,Nmc)';
    dB_bridge = sqrt(hEM/2)*randn(Nstep/2,Nmc)';
    dB_all = zeros(Nmc,Nstep);
    dB_all(:,1:2:end) = .5*dB_coarse + dB_bridge;
    dB_all(:,2:2:end) = .5*dB_coarse - dB_bridge;
    clear dB_coarse dB_bridge;
else
    dB_all = sqrt(hEM)*randn(Nstep,Nmc)';
end
for icase = 1:length(case_names)
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
    if size(case_par,2)>=12
        aper_amp = pars(12);
    else
        aper_amp = .2;
    end
    if duty<=0 || duty>1
        error('Active fraction must be in (0,1].');
    end
    if aper_amp<0 || aper_amp>=1
        error('Aperiodic amplitude must satisfy 0 <= A < 1.');
    end
    sample_steps = round(hs/hEM);
    if sample_steps<1 || abs(sample_steps*hEM-hs)>1e-12
        error('h_s must be a positive integer multiple of hEM.');
    end
    W = ones(M)-M*eye(M);
    pins = 1:np;
    if topology >= 2
        W = zeros(M);
        W(1,2:M) = 1;
        W(2:M,1) = 1;
        for inode=2:M-1
            W(inode,inode+1)=1;
        end
        W(M,2)=1;
        W = W-diag(sum(W,2));
        if topology==3
            pins=2:(np+1);
        end
    end
    kappa = zeros(M,1);
    kappa(pins)=gain;
    kappa_matrix = repmat(kappa,1,Nmc);
    damping_matrix = repmat(d,1,Nmc);
    bias = bias_amp*(.2+.8*(0:M-1)'/(M-1));
    bias_matrix = repmat(bias,1,Nmc);
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
    tail_first = round(.8*Nstep);
    active_mask = false(1,Nstep);
    cycle_start = zeros(1,Nstep);
    boundaries = [];
    left=0;
    mcycle=0;
    while left<Nstep
        if schedule_kind==1
            width = round(.05/hEM);
        else
            width = 4*round(.05*(1+aper_amp*sin(sqrt(2)*mcycle))/(4*hEM));
            width = max(4,width);
        end
        width = min(width,Nstep-left);
        right = left+width;
        on_width = round(duty*width);
        switch_off = min(left+on_width,Nstep);
        if switch_off>left
            active_mask(left+1:switch_off)=true;
        end
        cycle_start(left+1:right)=left;
        boundaries=[boundaries; left*hEM switch_off*hEM right*hEM];
        left=right;
        mcycle=mcycle+1;
    end
    if strategy==1
        active_mask(:)=true;
        cycle_start(:)=0;
    end
    schedule_records{icase}=boundaries;
    state_curve=zeros(length(record_steps),2*M);
    state_curve(1,:)=[phi(:,1)' omega(:,1)'];
    curves(icase,1)=mean(maximum_error);
    event_times=[];
    record_index=2;
    Nsample=0;
    tic;
    for k=1:Nstep
        is_active=active_mask(k);
        is_start=(k==1);
        if k>1
            is_start=is_active && ~active_mask(k-1);
        end
        if is_active
            is_sample = mod((k-1)-cycle_start(k),sample_steps)==0;
            if is_start || is_sample
                Nsample=Nsample+1;
                rho2=sum(phi.^2+omega.^2,1);
                eta2=sum((phi-phiEvent).^2+(omega-omegaEvent).^2,1);
                if is_start || strategy==2
                    update=true(1,Nmc);
                else
                    update=eta2>sigma1*rho2+sigma0;
                end
                phiEvent(:,update)=phi(:,update);
                omegaEvent(:,update)=omega(:,update);
                Nu=Nu+update;
                if update(1)
                    event_times=[event_times (k-1)*hEM];
                end
            end
            uPhi=-kappa_matrix.*phiEvent;
            uOmega=-kappa_matrix.*omegaEvent;
        else
            uPhi=zero_input;
            uOmega=zero_input;
        end
        old_error=sum(phi.^2+omega.^2,1);
        energy=sum(uPhi.^2+uOmega.^2,1);
        driftPhi=omega+c*(W*phi)+uPhi;
        driftOmega=-sin(phi)-damping_matrix.*omega+c*(W*omega)+bias_matrix+uOmega;
        phi=phi+hEM*driftPhi;
        omega=omega+hEM*driftOmega+noise*omega.*repmat(dB_all(:,k)',M,1);
        new_error=sum(phi.^2+omega.^2,1);
        integral_error=integral_error+.5*hEM*(old_error+new_error);
        integral_energy=integral_energy+hEM*energy;
        if k>tail_first
            integral_tail=integral_tail+.5*hEM*(old_error+new_error);
        end
        maximum_error=max(maximum_error,new_error);
        if record_index<=length(record_steps) && k==record_steps(record_index)
            curves(icase,record_index)=mean(new_error);
            state_curve(record_index,:)=[phi(:,1)' omega(:,1)'];
            record_index=record_index+1;
        end
    end
    active_fraction=sum(active_mask)/Nstep;
    tail_length=(Nstep-tail_first)*hEM;
    all_metrics(:,:,icase)=[Nu' ...
        integral_error'/Tsim ...
        integral_energy'/Tsim ...
        integral_tail'/tail_length ...
        new_error' ...
        maximum_error' ...
        np*Nu' ...
        active_fraction*ones(Nmc,1) ...
        Nsample*ones(Nmc,1)];
    state_example{icase}=state_curve;
    event_example{icase}=event_times;
    fprintf('%s: Nu=%.3f, J=%.9g, E=%.9g, tail=%.9g, active=%.6f, elapsed=%.1fs\n', ...
        case_names{icase},mean(Nu),mean(integral_error)/Tsim, ...
        mean(integral_energy)/Tsim,mean(integral_tail)/tail_length, ...
        active_fraction,toc);
end
clear dB_all;

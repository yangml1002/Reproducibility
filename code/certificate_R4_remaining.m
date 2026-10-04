alpha1 = 3;
eps1 = .5;
eps2 = .05;
muTheta = 10;
Nc = .05;
Lambda1 = 3;
Lambda2 = .2;
M = 10;
c = 11;
Ncase = length(case_names);
cert_values = NaN(Ncase,12);
schedule_check_pass = false(Ncase,1);
for icase = 1:Ncase
    hs = case_par(icase,3);
    sig1 = case_par(icase,4);
    sig0 = case_par(icase,5);
    duty = case_par(icase,6);
    alpha2 = case_par(icase,7);
    np = case_par(icase,8);
    gain = case_par(icase,9);
    topology = case_par(icase,10);
    qh = case_par(icase,11);
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
        W = W-diag(sum(W,2));
        if topology == 3
            pins = 2:(np+1);
        end
    end
    K = zeros(M);
    for ip = pins
        K(ip,ip) = gain;
    end
    lambda1 = max(real(eig(W'*W)));
    lambda2 = gain^2;
    Theta_max = muTheta*Nc;
    mu_e = gain*exp(Theta_max)/eps2;
    P1 = eps1 + alpha1 + alpha2^2 + gain*eps2 + ...
        .5*Lambda1 + .5*muTheta*(1-duty);
    P2 = eps1 + alpha1 + alpha2^2 - .5*muTheta*duty;
    eig1 = max(real(eig(P1*eye(M)+c*(W+W')/2-K)));
    eig2 = max(real(eig((P2+.5*Lambda2)*eye(M)+c*(W+W')/2)));
    A = 6*hs*(2*alpha1^2+c^2*lambda1)+4*alpha2^2;
    Esharp = (2*hs*A+12*hs^2*lambda2*(1+sig1))*exp(2*A*hs);
    Hsharp = (4*M*(3*hs+1)*hs*qh^2+12*hs^2*lambda2*sig0)*exp(2*A*hs);
    Ehat = 2*(Esharp+sig1);
    Hhat = 2*(Hsharp+sig0);
    b = 2*mu_e*Ehat;
    Lambda4 = mu_e*Hhat + (1/eps1+1)*M*qh^2*exp(Theta_max);
    r = NaN;
    rmin = NaN;
    R = NaN;
    if b < Lambda1
        lo = 0;
        hi = Lambda1;
        for it = 1:80
            mid = (lo+hi)/2;
            if mid-Lambda1+b*exp(mid*hs) > 0
                hi = mid;
            else
                lo = mid;
            end
        end
        r = (lo+hi)/2;
        rmin = min(r,Lambda2);
        R = 2*Lambda4/rmin;
    end
    if 4*mu_e*sig1 >= Lambda1
        hlimit = 0;
    else
        lo = 0;
        hi = .001;
        bracketed = false;
        for ib = 1:40
            Ab = 6*hi*(2*alpha1^2+c^2*lambda1)+4*alpha2^2;
            Eb = (2*hi*Ab+12*hi^2*lambda2*(1+sig1))*exp(2*Ab*hi);
            testval = 4*mu_e*(Eb+sig1);
            if ~isfinite(testval) || testval >= Lambda1
                bracketed = true;
                break;
            end
            hi = 2*hi;
        end
        if ~bracketed
            hlimit = NaN;
        else
            for it = 1:80
                mid = (lo+hi)/2;
                Am = 6*mid*(2*alpha1^2+c^2*lambda1)+4*alpha2^2;
                Em = (2*mid*Am+12*mid^2*lambda2*(1+sig1))*exp(2*Am*mid);
                testval = 4*mu_e*(Em+sig1);
                if ~isfinite(testval) || testval >= Lambda1
                    hi = mid;
                else
                    lo = mid;
                end
            end
            hlimit = (lo+hi)/2;
        end
    end
    active_ok = eig1 <= 0;
    rest_ok = eig2 <= 0;
    sampling_ok = b < Lambda1;
    cert_values(icase,:) = [ ...
        eig1, eig2, mu_e*Ehat, b, r, rmin, Lambda4, R, hlimit, ...
        active_ok, rest_ok, sampling_ok];
    schedule_check_pass(icase) = schedule_stats(icase,2) <= Nc + 1e-12;
end

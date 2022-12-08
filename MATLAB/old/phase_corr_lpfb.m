function symbsRec_ = phase_corr_lpfb(symbsRec, c)
    
    quad = imag(symbsRec);
    in_ph = real(symbsRec);
    N = length(quad);
    %{
    freq_log = zeros(1, N+1);
    freq=0;
    %err(k) does not need to be a list - this is just
    %here for tracking purposes.
    err =zeros(1,N);
    A = 0.575*100;
    B = 0.255*100;
    phase = 0;
    for k = 1:N
      %Adjust phase
      symbsRec(k) = symbsRec(k) * exp(-1i*phase);
      in_ph(k) = real(symbsRec(k));
      quad(k) = imag(symbsRec(k));
      err(k) = in_ph(k)*quad(k);
      freq = freq + B*err(k);
      freq_log(k) = freq*Fs/(2*pi);
      phase = phase + freq + A*err(k);
    end
    %}
    % Get a good estimate of initial phase
    phase_acc = zeros(1, N+1);
    phase_acc(1) = atan2(quad(1), in_ph(1));
    % err(k) does not need to be a list - this is just
    % here for tracking purposes.
    err =zeros(1,N);
    % our loop filter has transfer function:
    % (0.000439x^2 + 0.000878z + 0.000439)
    % ------------------------------------------------------
    % (z^2 - 1.951z + 0.9512)
    px = 0;
    px2 = 0;
    py =phase_acc(1) ;
    py2 = phase_acc(1) ;
    for k = 1:N
        % Adjust phase
        in_ph(k) = real(symbsRec(k) * exp(-1i*phase_acc(k)));
        quad(k) = imag(symbsRec(k) * exp(-1i*phase_acc(k)));
        err(k) = in_ph(k)*quad(k);
        % Implement loop filter
        % phase_acc(k+1) = 10*(0.00086024*err(k) + 0.0017*px +0.00086024*px2) +1.9785*py -0.9785*py2 ;
        phase_acc(k+1) = c(1)*(c(2)*err(k) + c(3)*px + c(4)*px2) + c(5)*py + c(6)*py2 ;
        px2 = px;
        px = err(k);
        py2 = py;
        py = phase_acc(k+1);
    end
    
    symbsRec_ = in_ph + 1i*quad;
end
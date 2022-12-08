function symbsRec_ = phase_corr_ls(symbsRec)
    in_ph = real(symbsRec);
    quad = imag(symbsRec);
    slope = in_ph'\quad';
    theta = atan(slope);
    symbsRec_ = symbsRec * exp(1i*(-theta));
end
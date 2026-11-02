function R = opk2R(omega,phi, kappa)

    a1=cos(phi)*cos(kappa);
    a2=-cos(phi)*sin(kappa);
    a3=sin(phi);
    b1=sin(omega)*sin(phi)*cos(kappa)+cos(omega)*sin(kappa);
    b2=-sin(omega)*sin(phi)*sin(kappa)+cos(omega)*cos(kappa);
    b3=-sin(omega)*cos(phi);
    c1=-cos(omega)*sin(phi)*cos(kappa)+sin(omega)*sin(kappa);
    c2=cos(omega)*sin(phi)*sin(kappa)+sin(omega)*cos(kappa);
    c3=cos(omega)*cos(phi);
    R=[a1 b1 c1;
       a2 b2 c2;
       a3 b3 c3];
end


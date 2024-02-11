%% Definition of the function - Load on the back due to the barbell 
function Rhand = f_rhand(pho,W,Alpha,Bheta,Gamma)
if pho==0
    Rhand=0;
else
    A=cos(Gamma).^2*(1-1/(pho^2));
    B=-2*W*cos(Alpha)*sin(Bheta).*cos(Gamma)-2*W*sin(Alpha)*cos(Bheta).*cos(Gamma);
    C=W^2;
    Delta=B.^2-4*A*C;
    Rhand=(-B-sqrt(Delta))./(2*A);
end
end
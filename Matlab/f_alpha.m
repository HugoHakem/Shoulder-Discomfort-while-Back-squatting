%% Definition of the function - Upper Body / shoulder inclination bis
function Alpha = f_alpha(DeltaH,Gamma,Theta,c,H1,H2,Bheta)
Alpha = atan((c-DeltaH -H1*cos(Gamma)*cos(Bheta)+H2*cos(Theta))./(c-H1*cos(Gamma)*sin(Bheta)));
Alpha(Alpha<0)=nan;
end

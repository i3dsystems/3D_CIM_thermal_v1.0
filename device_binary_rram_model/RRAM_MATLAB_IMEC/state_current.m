function current= state_current( gap, voltage ) 
    I0=6.14E-5;  % Initial current (A)
    a= 2.75e-10;
    b=0.43; 

    %gap0= 0.4E-9;
    %SL= 2e-11;
    %gap_min= 0.25E-9;

    %G0= I0*exp(-gap/gap_min)*10;
    current_tunneling = I0*exp(-gap/a)*sinh(voltage/b);
%     current_ohmic= G0*voltage;
%     F= 1/(1+exp((gap-gap0)/SL));
%     current= current_tunneling*(1-F) + current_ohmic*F;
    
    current = current_tunneling;
end
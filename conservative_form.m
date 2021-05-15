function [v,pressure,T,Mach_no,mass_flow,rho,rho_throat,P_throat,v_throat,T_throat,Mach_no_throat, mass_flow_throat] = conservative_form(n,x,dx,gamma,CFL,nt);
    
    %initial profile
    %all flow variables are in their non-dimensional form

    for i = 1:n
        if x(i) <= 0.5
            rho(i) = 1;                                  %density
            T(i)   = 1;                                  %temperature
        elseif (x(i) > 0.5 && x(i) <= 1.5)
            rho(i) = 1 - 0.366*(x(i) - 0.5);
            T(i)   = 1 - 0.167*(x(i) - 0.5);
        elseif (x(i) > 1.5 && x(i) <= 3)
            rho(i) = 0.634 - 0.3879*(x(i) - 1.5);
            T(i)   = 0.833 - 0.3507*(x(i) - 1.5);
            P(i)   = rho(i)*T(i);
        end

        A(i)         = 1 + 2.2*(x(i) - 1.5).^2;               %converging diverging area
        v(i)         = 0.59/(rho(i)*A(i));                    %velocity
        mass_flow(i) = rho(i).*A(i).*v(i);                    %mass flow rate
    end

    mass_flow_initial = rho.*A.*v;                    %initial mass flow rate

    throat = find(A == 1);                            %find throat location


    %defining inital solution vector

    U1 = rho.*A;
    U2 = rho.*A.*v;
    U3 = rho.*A.*((T./(gamma - 1)) + ((gamma/2).*v.^2));


    total_sim_time = 0;                               %initial flow time

    %Solving the equations (Conservative Form)
    %Maccormack Method


    for z = 1:nt;                   %time marching

        % Storing the initial data

        U1_old   = U1;
        U2_old   = U2;
        U3_old   = U3;

        %calculating the time step

        dt = min((CFL.*dx)./(sqrt(T) + v));

        %calculating the intial flux vector

        F1 = U2;
        F2 = ((U2.^2)./U1) + ((gamma - 1)/gamma)*(U3 - (gamma/2)*(U2.^2)./U1);
        F3 = ((gamma.*U2.*U3)./U1) - (gamma*(gamma-1)/2)*((U2.^3)./(U1.^2));

        %PRIDECTOR STEP

        for j = 2:n-1
            
            %source term
            J(j)  = (1/gamma)*rho(j)*T(j)*((A(j+1) - A(j))/dx);    

            %continuity equation
            dU1_dt_p(j) = -((F1(j+1) - F1(j))/dx);

            %momemtum equation
            dU2_dt_p(j) = -((F2(j+1) - F2(j))/dx) + J(j);

            %energy equation
            dU3_dt_p(j) = -((F3(j+1) - F3(j))/dx);

            %updating the solution vector
            U1(j) = U1(j) + dU1_dt_p(j)*dt;
            U2(j) = U2(j) + dU2_dt_p(j)*dt;
            U3(j) = U3(j) + dU3_dt_p(j)*dt;

        end

        %correcting flux vector

        F1 = U2;
        F2 = ((U2.^2)./U1) + ((gamma - 1)/gamma)*(U3 - (gamma/2)*(U2.^2)./U1);
        F3 = ((gamma.*U2.*U3)./U1) - (gamma*(gamma-1)/2)*((U2.^3)./(U1.^2));

        %calculating premitive variables

        rho = U1./A;
        v   = U2./U1;
        T   = (gamma - 1)*((U3./U1) - ((gamma/2)*(v.^2)));
        P   = rho.*T;

        %CORRECTOR STEP

        for k = 2:n-1

            J(k)  = (1/gamma)*rho(k)*T(k)*((A(k) - A(k-1))/dx);    %source term

            %continuity equation
            dU1_dt_c(k) = -((F1(k) - F1(k-1))/dx);

            %momemtum equation
            dU2_dt_c(k) = -((F2(k) - F2(k-1))/dx) + J(k);

            %energy equation
            dU3_dt_c(k) = -((F3(k) - F3(k-1))/dx);

        end

        %calculating the averaged time derivatives

        dU1_dt = 0.5*(dU1_dt_p + dU1_dt_c);
        dU2_dt = 0.5*(dU2_dt_p + dU2_dt_c);
        dU3_dt = 0.5*(dU3_dt_p + dU3_dt_c);


        %calculating corrected solution vector
        for m = 2:n-1
            U1(m) = U1_old(m) + dU1_dt(m)*dt;
            U2(m) = U2_old(m) + dU2_dt(m)*dt;
            U3(m) = U3_old(m) + dU3_dt(m)*dt;
        end

        %apply BC

        %inlet
        U1(1) = rho(1)*A(1);             %fixed value
        U2(1) = 2*U2(2) - U2(3);         %extrapolated              
        U3(1) = U1(1)*((T(1)/(gamma - 1)) + ((gamma/2)*(v(1)^2)));

        %outlet
        %extrapolated from the inside domain
        U1(n) = 2*U1(n-1) - U1(n-2);
        U2(n) = 2*U2(n-1) - U2(n-2);
        U3(n) = 2*U3(n-1) - U3(n-2);

        %calculating corrected primitive non-dimensional variable for next time step
        rho = U1./A;
        v   = U2./U1;
        T   = (gamma - 1)*((U3./U1) - ((gamma/2)*(v.^2)));
        pressure   = rho.*T;
        mass_flow = rho.*A.*v;
        Mach_no   = (v./sqrt(T));


        %calculating non-dimensioanl variables at the throat
        rho_throat(z) = rho(throat);
        v_throat(z)   = v(throat);
        T_throat(z)   = T(throat);
        P_throat(z)   = P(throat);
        mass_flow_throat(z) = mass_flow(throat);
        Mach_no_throat(z)   = Mach_no(throat);


        total_sim_time = total_sim_time + dt;  %updating the total flow time

        %plotting the numerical results
        
        %non-dimensional mass flow rate at the different time steps

        figure(6)
        if z == 1
            plot(x,mass_flow_initial,'color','K','LineWidth', 1.3, 'LineStyle','--');
            hold on
        elseif z == 50
            plot(x,mass_flow, 'm', 'LineWidth', 1.3);
            hold on
        elseif z == 100
            plot(x,mass_flow, 'c', 'LineWidth', 1.3);
            hold on
        elseif z == 200
            plot(x,mass_flow, 'g', 'LineWidth', 1.3);
            hold on
        elseif z == 400
            plot(x,mass_flow, 'y', 'LineWidth', 1.3);
            hold on
        elseif z == 700
            plot(x,mass_flow, 'k', 'LineWidth', 1.3);
            hold on
        end

        title("Mass Flow Rate Distribution across the Nozzle [Conservative Form]")
        xlabel("Nozzle X-Distance")
        ylabel("Non-Dimensional Mass Flow Rate")
        legend('0^t^h Timestep','50^t^h Timestep','100^t^h Timestep','200^t^h Timestep','400^t^h Timestep','700^t^h Timestep')
        axis([0 3 0.4 1])
        grid minor;

    end
    grid minor;

    %printing total flow time

    sprintf("Total Flow time for Conservative form = %.2f seconds", total_sim_time)
end
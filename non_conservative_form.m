function [v,pressure,T,Mach_no,mass_flow,rho,rho_throat,P_throat,v_throat,T_throat,Mach_no_throat, mass_flow_throat] = non_conservative_form(n,x,dx,gamma,CFL,nt)

    %initial profile
    %all flow variables are in their non-dimensional form
    
    rho = 1 - 0.3146*x;                   %density
    T   = 1 - 0.2314*x;                   %temperature
    v   = (0.1 + 1.09*x).*T.^0.5;         %velocity
    A   = 1 + 2.2*(x - 1.5).^2;           %converging diverging area

    mass_flow_initial = rho.*A.*v;        %mass flow rate

    throat = find(A == 1);                %find throat location

    total_sim_time = 0;                   %initial flow time

    %solving the discretized equations (Non-Conservative Form)
    %Maccormack Method
    
    for k = 1:nt;                            %time marching

        % Storing initial data

        T_old   = T;
        v_old   = v;
        rho_old = rho;

        % Calculating the time step

        dt = min((CFL.*dx)./(sqrt(T) + v));

        % Predictor method

        for j = 2:n-1;

            dvdx   = (v(j+1) - v(j))/dx;
            dlnAdx = (log(A(j+1)) - log(A(j)))/dx;
            drhodx = (rho(j+1) - rho(j))/dx;
            dTdx   = (T(j+1) - T(j))/dx;

            %continutiy equation
            drhodt_p(j) = -rho(j)*dvdx - rho(j)*v(j)*dlnAdx - v(j)*drhodx;

            %momemtum equation
            dvdt_p(j)   = -v(j)*dvdx - (1/gamma)*(dTdx + (T(j)/rho(j))*drhodx);

            %energy equation
            dTdt_p(j)   = -v(j)*dTdx - (gamma -1)*T(j)*(dvdx + v(j)*dlnAdx);

            %solution update
            v(j)   = v(j) + dvdt_p(j)*dt;
            rho(j) = rho(j) + drhodt_p(j)*dt;
            T(j)   = T(j)+ dTdt_p(j)*dt;

        end

        %corrector method

        for j = 2:n-1;

            dvdx   = (v(j) - v(j-1))/dx;
            dlnAdx = (log(A(j)) - log(A(j-1)))/dx;
            drhodx = (rho(j) - rho(j-1))/dx;
            dTdx   = (T(j) - T(j-1))/dx;

            %continutiy equation
            drhodt_c(j) = -rho(j)*dvdx - rho(j)*v(j)*dlnAdx - v(j)*drhodx;

            %momemtum equation
            dvdt_c(j)   = -v(j)*dvdx - (1/gamma)*(dTdx + (T(j)/rho(j))*drhodx);

            %energy equation
            dTdt_c(j)   = -v(j)*dTdx - (gamma -1)*T(j)*(dvdx + v(j)*dlnAdx);

        end

        % compute the average time derivative

        dvdt = 0.5*(dvdt_p + dvdt_c);
        drhodt = 0.5*(drhodt_p + drhodt_c);
        dTdt = 0.5*(dTdt_p + dTdt_c);

        % Solution update

        for i = 2:n-1
            v(i) = v_old(i) + dvdt(i)*dt;
            rho(i) = rho_old(i) + drhodt(i)*dt;
            T(i) = T_old(i) + dTdt(i)*dt;
        end

        %boundary condition

        %inlet
        rho(1) = 1;                 %fixed value
        T(1)   = 1;                 %fixed value
        v(1)   = 2*v(2) - v(3);     %extrapolated

        %outlet
        %extrapolated from the inside domain
        v(n)   = 2*v(n-1) - v(n-2);
        rho(n) = 2*rho(n-1) - rho(n-2);
        T(n)   = 2*T(n-1) - T(n-2);

        %calculating the non-dimensional Mass Flow Rate, Pressure & Mach Number

        mass_flow = rho.*A.*v;
        pressure  = rho.*T;
        Mach_no   = v./sqrt(T);

        %calculating non-dimensional variables at the nozzle throat
        
        rho_throat(k) = rho(throat);
        v_throat(k)   = v(throat);
        T_throat(k)   = T(throat);
        P_throat(k)   = pressure(throat);
        mass_flow_throat(k) = mass_flow(throat);
        Mach_no_throat(k)   = Mach_no(throat);

        total_sim_time = total_sim_time + dt;        %updating total flow time

        dvdt_throat(k) = dvdt(throat);
        drhodt_throat(k) = drhodt(throat);

        %plotting the numerical results
        
        %non-dimensional mass flow rate at the different time steps

        figure(1)
        
        if k == 1
            plot(x,mass_flow_initial,'color','K','LineWidth', 1.3, 'LineStyle','--');
            hold on
        elseif k == 50
            plot(x,mass_flow, 'm', 'LineWidth', 1.3);
            hold on
        elseif k == 100
            plot(x,mass_flow, 'c', 'LineWidth', 1.3);
            hold on
        elseif k == 150
            plot(x,mass_flow, 'g', 'LineWidth', 1.3);
            hold on
        elseif k == 200
            plot(x,mass_flow, 'y', 'LineWidth', 1.3);
            hold on
        elseif k == 700
            plot(x,mass_flow, 'k', 'LineWidth', 1.3);
            hold on
        end

        title("Mass Flow Rate Distribution across the Nozzle [Non-Conservative Form]")
        xlabel("Nozzle X-Distance")
        ylabel("Non-Dimensional Mass Flow Rate")
        legend('0^t^h Timestep', '50^t^h Timestep','100^t^h Timestep','150^t^h Timestep','200^t^h Timestep','700^t^h Timestep')
        axis([0 3 0.1 2])
        grid minor;

    end
    grid minor;
    
    %plotting the density & velocity time derivative variation with time at the nozzle throat
    
    figure(2)
    
    plot(linspace(1,nt,nt),dvdt_throat, 'k', 'LineWidth', 1.3);
    hold on
    plot(linspace(1,nt,nt),drhodt_throat,'g', 'LineWidth', 1.3);

    title("Time Derivatives of Velocity & Density [Non-Conservative Form]")
    xlabel("Number of Time Steps")
    ylabel("Residuals")
    legend({'$dv/dt$' , '$d\rho/dt$'},'Interpreter', 'latex');
    grid minor;
    legend({'$dv/dt$' , '$\partial{\rho}/\partial{t}$'},'Interpreter', 'latex');
    
    %plotting steady-state Mach No. & Density along with the Nozzle

    figure(3)
    plot(x,rho,'k', 'LineWidth', 1.3);
    xlabel('Nozzle X - Direction')
    ylabel('Non-Dimensional Density')

    yyaxis right
    plot(x,Mach_no,'g', 'LineWidth', 1.3);
    ylabel('Non-Dimensional Mach No.')
    grid minor;

    title("Steady State Density & Mach No. across Nozzle [Non-Conservative Form]")
    legend('Density','Mach Number');
    
    %printing total flow time

    sprintf("Total Flow time for non conservative form = %.2f seconds", total_sim_time)
end
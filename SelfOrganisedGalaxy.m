
tic
Rings=49;                           % Define the number of rings in the galaxy
Iterations=80000;                  % Define the nubmer of times the simulation repeats
x=[];                               % Initialise an empty matrix for which to generate x values of the galaxy matrix elements
y=[];                               % Initialise an empty matrix for which to generate y values of the galaxy matrix elements
P=[];                               % Initialise an empty matrix for which to record the number of critical sites each timestep
PROB=[];                            % Initialise an empty matrix for which to record the probability of star formation each iteration
AvalancheSize=[];                   % Initialise an empty matrix for which to record the size of avalanches each iteration
AvalancheTime=[];                   % Initialise an empty matrix for which to record the time of avalanches each iteration
EventSize=0;                        % Initialise EventSize as 0
EventTime=0;                        % Initialise EventTime as 0
%%
for i=1:Rings                                                               % Loop to convert the linear matrix into a polar matrix
    Cells=i*6;                                                              % Defining the number of cells in the current ring
    r=i*ones(1,Cells);                                                      % Defining the radius of each cell in the current ring to the centre
    theta=linspace(0+((2*pi)/Cells)/2, (2*pi)+((2*pi)/Cells)/2, Cells+1);   % Define angle of centre of cells ((2*pi)/Cells) to base in radians
    theta(end)=[];                                                          % Eliminating the final angle as it is the angle between the final and first cell, which is unneccesary
    [xx yy]=pol2cart(theta, r);                                             % Convert polar coordinates of cells into cartesian matrix
    x=[x,xx];                                                               % Assigning the polar coordinate of x values in the system to matrix x
    y=[y,yy];                                                               % Assigning the polar coordinate of y values in the system to matrix y
    % xx and yy are stored as the outer-most/boundary ring %

end

density=randi([0,5],1,length(x));                  % Population of randomly generated mass density values throughout an initialised matrix
Mass=sum(density);                                 % Calculate the mass of the system in the total density as a sum of the density values
MassDensity=Mass./length(x);                       % Calculate the average mass density of the system

%%
for n=1:Iterations                                              % Loop to perform the operation of the simulation equal to a pre-defined number of timesteps
    SFC=find(density>=6);                                       % Storing the active star-forming regions in matrix SFC
    %% When measuring Event Sizes remove %'s on this code block %%
    % If generating galaxy scatter plot keep %'s on code block %

    EventSize=EventSize+length(SFC);                          % Count to increase measure of event size of current perturbation iteration
    EventTime=EventTime+1;                                    % Count to increase measure of event time of current perturbation iteration
    if isempty(SFC)==1                                    % Condition which no active star-forming regions remain and the system is perturbed again
        if EventSize>0                                    % Filters out interations where event sizes are 0
            AvalancheSize=[AvalancheSize, EventSize];     % Recording the size of independent events once they conclude
            AvalancheTime=[AvalancheTime, EventTime];     % Recording the time of independent events once they conclude
        end
        EventSize=0;                                  % Once measure of event is finished, reset EventSize parameter
        EventTime=0;                                  % Once measure of event is finished, reset EventTime parameter

        for per=1:1                                        % Mechanism of perturbing the system
            t=randi(numel(density));                        % Generate a random location in density matrix for which to perturb
            density(t)=density(t)+1;                        % Perturb the location in question
        end
    end
    %% storing different stages of the simulation for graphical plotting
    if n==Iterations-1
        PrevStep=SFC;
    end

    if n==Iterations-2
        PrevStep2=SFC;
    end
    if n==Iterations-3
        PrevStep3=SFC;
    end
    if n==Iterations-4
        PrevStep4=SFC;
    end
    if n==Iterations-5
        PrevStep5=SFC;
    end
    if n==Iterations-6
        PrevStep6=SFC;
    end
    if n==Iterations-7
        PrevStep7=SFC;
    end
    if n==Iterations-8
        PrevStep8=SFC;
    end
    if n==Iterations-9
        PrevStep9=SFC;
    end
    if n==Iterations-10
        PrevStep10=SFC;
    end
    if n==Iterations-11
        PrevStep11=SFC;
    end

    %%
    for m=1:length(SFC)                                 % Loop to find relevant information for memeber of SFC
        XC=x(SFC(m));                                   % Member of the matrix x which is star-forming
        YC=y(SFC(m));                                   % Memeber of matrix y which is star-forming
        D=((x-XC).^2+(y-YC).^2).^0.5;                   % Calculating distance between site in question and other sites

        if ismember(XC,xx)==1 && ismember(YC,yy)==1     % If the site in question is in the outermost ring, set 4 neighbours
            [~,Neighbours]=mink(D,5);
        else                                            % If the site is not in the outer most ring, set 6 neighbours
            [~,Neighbours]=mink(D,7);
        end

        for k=2:length(Neighbours)                              % Loop to determine in Neighbours of the site in question become star-forming

            density(Neighbours(k))=density(Neighbours(k))+1;    % Mechanism of relaxation of the region that is no longer star-forming


        end
        density(Neighbours(1))=0;                       % Regions that were intially star forming at the beginning of the timestep undergo star death

    end
    C=0;                % Initialise a count at 0
    XstepTot=[];        % Initialise an empty array to introduce new rotated x coordinates
    YstepTot=[];        % Initialise an empty array to introduce new rotated y coordinates
    steps=10;           % Define the number of locations which the system will rotate, 200km/s is 10 steps, 250km/s is 13, 300km/s is 15 steps
    %%

    for s=1:Rings
        AngRot=steps*(2*pi/(6*s));              %Angle of rotation

        Xsr=(x(C+1:C+ s*6))*cos(AngRot) + (y(C+1:C+ s*6))*sin(AngRot);    %Rotating the cells of the rings in x
        Ysr=(y(C+1:C+ s*6))*cos(AngRot) - (x(C+1:C+ s*6))*sin(AngRot);    %Rotating the cells of the rings in y

        C=C+s*6;                                % Updating count so as the move to the next ring
        XstepTot=[XstepTot, Xsr];               % Update the array of rotated x values
        YstepTot=[YstepTot, Ysr];               % Update the array of rotated y values
    end


    x=XstepTot;             % Define the new coordinates of the x values
    y=YstepTot;             % Define the new coordinates of the y values

    p=sum(density>=6);                                  % Define p as the sum of all star forming regions
    P=[P, p];                                           % Record # of star forming regions each iteration in matrix P
    Mass=[Mass, sum(density)];                          % Store Mass parameter of each iteration as a matrix
    MassDensity=[MassDensity, Mass(n+1)./length(x)];    % Store MassDensity parameter of each iteration as a matrix
end

for t=1:(length(P)-1)                       % Calculating the probability of star formation throughout the simulation
    PROB=[PROB, (P(t+1)./(5.92*P(t)))];     % Formula for calculating probability as per Finn's final report
end
PROB=PROB(~isnan(PROB));                    % Removing values of PROB that are Inf or NaN
PROB=PROB(~isinf(PROB));                    % Removing values of PROB that are Inf or NaN
PROBavg=mean(PROB);                         % Average probability of the simulation



%% If measuring events must uncomment %%
muSize=min(AvalancheSize);                % Find smallest event size
nuSize=max(AvalancheSize);                % Find greatest event size
BSize=linspace(muSize,nuSize,nuSize);     % Create a linspace ranging between smallest and greatest event sizes
COUNTSize=hist(AvalancheSize,BSize);      % Count the number of events of various sizes in matrix AvalancheSize
DsSize=COUNTSize./(sum(COUNTSize));       % Calculate the probability density of each of the events of different sizes occurring

muTime=min(AvalancheTime);                % Find smallest event size
nuTime=max(AvalancheTime);                % Find greatest event size
BTime=linspace(muTime,nuTime,nuTime);     % Create a linspace ranging between smallest and greatest event sizes
COUNTTime=hist(AvalancheTime,BTime);      % Count the number of events of various sizes in matrix AvalancheSize
DsTime=COUNTTime./(sum(COUNTTime));       % Calculate the probability density of each of the events of different sizes occurring



%% non-linear regression
F=@(beta,BSize) beta(1).*(BSize.^(beta(2)).*exp(-BSize./(beta(3))));     % Function of line of best fit
initials=[1,0.9,10];                                                     % Initial guesses of beta coefficients
new_coeffs=nlinfit(BSize(1,2:200),DsSize(1,2:200),F,initials);           % Apply nlinfit function to fun F by applying assigned values of BSize and DsSize with generated initials
new_y=F(new_coeffs, BSize(1,2:200));                                     % Create non-linear regression y values using new_coeffs and BSize x values

%%
% Writing desired data to a seperate file %
filename_str='StarFormationProbability vs Perturbation Rate SOC';
fid=fopen(filename_str,'a');
fprintf(fid,'%15.10f %15.10f \n', per, PROBavg);
fclose(fid);

%%
% Graphical Scatter Plot of the galaxy %
figure
scatter(x,y,'white','.');
axis equal;
SFC=find(density>=6);
for q=1:length(SFC)
    XXc=x(SFC(q));
    YYc=y(SFC(q));
    hold on
    scatter(XXc,YYc,24,"black","filled");
end

for a=1:length(PrevStep)
    XXc=XstepTot(PrevStep(a));
    YYc=YstepTot(PrevStep(a));
    hold on
    scatter(XXc,YYc,22,'black','filled');
end

for b=1:length(PrevStep2)
    XXc=XstepTot(PrevStep2(b));
    YYc=YstepTot(PrevStep2(b));
    hold on
    scatter(XXc,YYc,20,'black','filled');
end

for c=1:length(PrevStep3)
    XXc=XstepTot(PrevStep3(c));
    YYc=YstepTot(PrevStep3(c));
    hold on
    scatter(XXc,YYc,18,'black','filled');
end

for d=1:length(PrevStep4)
    XXc=XstepTot(PrevStep4(d));
    YYc=YstepTot(PrevStep4(d));
    hold on
    scatter(XXc,YYc,14,'black','filled');
end

for e=1:length(PrevStep5)
    XXc=XstepTot(PrevStep5(e));
    YYc=YstepTot(PrevStep5(e));
    hold on
    scatter(XXc,YYc,12,'black','filled');
end
for ei=1:length(PrevStep6)
    XXc=XstepTot(PrevStep6(ei));
    YYc=YstepTot(PrevStep6(ei));
    hold on
    scatter(XXc,YYc,10,'black','filled');
end
for eii=1:length(PrevStep7)
    XXc=XstepTot(PrevStep7(eii));
    YYc=YstepTot(PrevStep7(eii));
    hold on
    scatter(XXc,YYc,8,'black','filled');
end
for eiii=1:length(PrevStep8)
    XXc=XstepTot(PrevStep8(eiii));
    YYc=YstepTot(PrevStep8(eiii));
    hold on
    scatter(XXc,YYc,6,'black','filled');
end
for eiv=1:length(PrevStep9)
    XXc=XstepTot(PrevStep9(eiv));
    YYc=YstepTot(PrevStep9(eiv));
    hold on
    scatter(XXc,YYc,4,'black','filled');
end
for ev=1:length(PrevStep10)
    XXc=XstepTot(PrevStep10(ev));
    YYc=YstepTot(PrevStep10(ev));
    hold on
    scatter(XXc,YYc,2,'black','filled');
end
for evi=1:length(PrevStep11)
    XXc=XstepTot(PrevStep11(evi));
    YYc=YstepTot(PrevStep11(evi));
    hold on
    scatter(XXc,YYc,1,'black','filled');
end

% Plot of Trend data of system %
% If events are not measured relevant plots must be commented out %
figure
plot([1:Iterations],P,'MarkerEdgeColor','red')
legend('20','50','80')
xlabel('Iterations')
ylabel('Number of Stars Present')
hold on

figure
plot(PROB,'MarkerEdgeColor','red')


figure
plot(BSize,DsSize)
title(['Linear Plot Prob Dens Event Size; Perturbation=',num2str(per)])

figure
loglog(BSize,DsSize)
title(['Log Plot Prob Dens Event Size; Perturbation=',num2str(per)])
hold on
plot(BSize(1,2:200),new_y, 'r')
text(20,0.01, new_coeffs(1)+"*(x\^("+new_coeffs(2)+"))*e\^(-x/("+new_coeffs(3)+")))")
text(20,0.02, 'Function of Line of Best Fit;')
title('Distribution of Event Sizes')
xlabel('Size of Event S')
ylabel('Probability Density of Event Size S')
legend('Real Data', 'Fitting Line')

figure
semilogy(BSize,DsSize)
title(['Lin Log Plot Prob Dens Event Size; Perturbation=',num2str(per)])

figure
semilogx(BSize,DsSize)
title(['Log Lin Plot Prob Dens Event Size; Perturbation=',num2str(per)])

figure
plot(BTime,DsTime)
title(['Linear Plot Prob Density Event Time; Perturbation=',num2str(per)])

figure
loglog(BTime,DsTime)
title(['Log Plot Prob Dens Event Time; Perturbation=',num2str(per)])

figure
plot(1:(Iterations+1), MassDensity)
title(['Mass Density over time; Perturbation=',num2str(per)])

toc
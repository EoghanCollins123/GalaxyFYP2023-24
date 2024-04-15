tic
Rings=50;                           % Define the number of rings in the galaxy
Iterations=800;                     % Define the nubmer of times the simulation repeats
x=[];                               % Initialise an empty matrix for which to generate x values of the galaxy matrix elements
y=[];                               % Initialise an empty matrix for which to generate y values of the galaxy matrix elements
N_s_T=[];                           % Initialise an empty matrix for which to record # of stars formed each timestep
NewStar=0;                          % Initialise NewStar parameter as 0
StarFormationFluctuations=[];
TotStarForm=[];                     % Initialise an empty matrix for which to count all stars formed throught simulation at each timestep
%%
for i=1:Rings                                                               % Loop to convert the linear matrix into a polar matrix
    Cells=i*6;                                                              % Defining the number of cells in the current ring
    r=i*ones(1,Cells);                                                      % Defining the radius of each cell in the current ring to the centre
    theta=linspace(0+((2*pi)/Cells)/2, (2*pi)+((2*pi)/Cells)/2, Cells+1);   % Define angle of centre of cells ((2*pi)/Cells) to base in radians
    theta(end)=[];                                                          % Eliminating the final angle as it is the angle between the final and first cell, which is unneccesary
    [xx yy]=pol2cart(theta, r);                                             % Convert polar coordinates of cells into cartesian matrix
    x=[x,xx];                                                               % Assigning the polar coordinate of x values in the system to matrix x
    y=[y,yy];                                                               % Assigning the polar coordinate of y values in the system to matrix y
    % xx and yy are stored as the outer-most/boundary ring

end

POP=0.01;                   % Initial population of randomly generated regions to be starforming
Z=rand(1,length(x));        % Generate Z, a matrix of random values between 0 and 1, which is the same size as matrix x
PROB=0.17;                  % Define probability of star-formation propagation
SFC=find(Z<POP);            % Storing the initially randomly generated starry regions for operation
%%
for n=1:Iterations          % Loop to perform the operation of the simulation equal to a pre-defined number of timesteps
    if n>1                  % Define the active star-forming regions each iteration after the first
        SFC=find(Z==1);     % Storing the active star-forming regions in matrix SFC
    end
    % storing different stages of the simulation for graphical plotting %
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
    %%
    Z=zeros(1,length(x));               % Regions that were intially star forming at the beginning of the timestep undergo star death
    for m=1:length(SFC)                 % Loop to find relevant information for memeber of SFC
        XC=x(SFC(m));                   % Member of the matrix x which is star-forming
        YC=y(SFC(m));                   % Memeber of matrix y which is star-forming
        D=((x-XC).^2+(y-YC).^2).^0.5;   % Calculating distance between site in question and other sites

        if ismember(XC,xx)==1 && ismember(YC,yy)==1     %If the site in question is in the outermost ring, set 4 neighbours
            [~,Neighbours]=mink(D,5);
        else                                            %If the site is not in the outer most ring, set 6 neighbours
            [~,Neighbours]=mink(D,7);
        end

        for k=2:length(Neighbours)                      % Loop to determine in Neighbours of the site in question become star-forming
            if rand<PROB                                % Neighbour becomes star-forming is a random number between 0 and 1 is less than the defined probability
                if Z(Neighbours(k))~=1                  % Set the neighbour as star forming is it is not already so
                    Z(Neighbours(k))=1;
                    NewStar=NewStar+1;                  % Increase the value of NewStar parameter to count the star forming regions created each iteration
                end
            end
        end
    end
    TotStarForm=[TotStarForm, NewStar];         %Total Stars formed through simulation
    C=0;                                        % Initialise a count at 0
    XstepTot=[];                                % Initialise an empty array to introduce new rotated x coordinates
    YstepTot=[];                                % Initialise an empty array to introduce new rotated y coordinates
    steps=10;                                   % Define the number of locations which the system will rotate, 10 steps correlates to 200km/s
    %%
    for s=1:Rings                               % Rotate each ring individually
        AngRot=steps*(2*pi/(6*s));              % Angle of rotation

        Xsr=(x(C+1:C+ s*6))*cos(AngRot) + (y(C+1:C+ s*6))*sin(AngRot);    %Rotating the cells of the rings in x
        Ysr=(y(C+1:C+ s*6))*cos(AngRot) - (x(C+1:C+ s*6))*sin(AngRot);    %Rotating the cells of the rings in y

        C=C+s*6;                                %Updating count so as the move to the next ring
        XstepTot=[XstepTot, Xsr];               %Update the array of rotated x values
        YstepTot=[YstepTot, Ysr];               %Update the array of rotated y values
    end


    x=XstepTot;                 % Define the new coordinates of the x values
    y=YstepTot;                 % Define the new coordinates of the y values
    p=sum(Z==1);                % Define p as the sum of all star forming regions
    N_s_T=[N_s_T, p];           % Define N_s_T as the # of stars formed during each interval



end

SFF=N_s_T./length(x);          % # Star regions per interval vs Total sites (StarFormation rate per Seiden and Schulman)
StarFormationRate=mean(SFF);   % Define mean star formation rate of the simulation for plotting at different probabilities


for nova=2:length(SFF)                                                          % Loop to calculate the fluctuations in star formation
    StarFFluctuations=(SFF(nova)-SFF(nova-1));                                  % Finding the fluctuation/difference in star formation rate between timesteps
    StarFormationFluctuations=[StarFormationFluctuations, StarFFluctuations];   % Storing the starformationfluctuations in a matrix
end

MeanStarFormationFluctuations=mean(StarFormationFluctuations);                  % Calculating the mean of fluctuations in star formation rate
RMSStarFormationFluctuations=rms(StarFormationFluctuations);                    % Calculating the root mean squared of the starformation fluctuations

% Writing desired data to a seperate file %
filename_str='Star Formation Rate vs Probability';
fid=fopen(filename_str,'a');
fprintf(fid,'%15.10f %15.10f \n', PROB, StarFormationRate);
fclose(fid);

%%
% Graphical Scatter Plot of the galaxy %
figure
scatter(x,y,'white','.');
axis equal;
title(['Probability=',num2str(PROB)])
SFC=find(Z==1);
for q=1:length(SFC)
    XXc=x(SFC(q));
    YYc=y(SFC(q));
    hold on
    scatter(XXc,YYc,18,"black","filled");
end

for a=1:length(PrevStep)
    XXc=XstepTot(PrevStep(a));
    YYc=YstepTot(PrevStep(a));
    hold on
    scatter(XXc,YYc,15,'black','filled');
end

for b=1:length(PrevStep2)
    XXc=XstepTot(PrevStep2(b));
    YYc=YstepTot(PrevStep2(b));
    hold on
    scatter(XXc,YYc,12,'black','filled');
end

for c=1:length(PrevStep3)
    XXc=XstepTot(PrevStep3(c));
    YYc=YstepTot(PrevStep3(c));
    hold on
    scatter(XXc,YYc,10,'black','filled');
end

for d=1:length(PrevStep4)
    XXc=XstepTot(PrevStep4(d));
    YYc=YstepTot(PrevStep4(d));
    hold on
    scatter(XXc,YYc,7,'black','filled');
end

for e=1:length(PrevStep5)
    XXc=XstepTot(PrevStep5(e));
    YYc=YstepTot(PrevStep5(e));
    hold on
    scatter(XXc,YYc,5,'black','filled');
end

% Plot of Trend data of system %
figure
loglog([1:n],SFF)
xlabel('Timesteps')
ylabel('Fractional Star Formation')
title('Star Formation per timestep as a fraction of Total Star Regions Formed')

figure
plot([1:n],SFF)
title('Star Formation per timestep normalised as a fraction of total system cells')
xlabel('Timesteps')
ylabel('Star Formation')


figure
plot(StarFormationFluctuations)
title('Fluctuations in Star Formation per timestep')
xlabel('Timestep')
ylabel('Fluctuations in Star Formation')



toc
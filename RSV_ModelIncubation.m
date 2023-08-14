% Originally created by Luke Riexinger March 2020
% Edited and modified by Sara Arena April 1 2022 
% Edited and modified by Anika Thatavarthy on April 6 2023

clear 
close all
clc

%% Simulation Time and Size
areaSize =  15;
total_time =  150*8; %frames; minutes
dt =  1; %time step

speed = 0.075;

%% Healthy People Properties
healthyNum =  28; %number of healthy cells at start of simulation

for h = 1:healthyNum
    healthyLocation(h, :) = PlaceCell(areaSize); %places each healthy cell randomly within the space
    healthyAngle(h,:) = randi(360);
    
    score = 100;
    vaccinated = randi(1000) < 578;
    handWash = randi(100) < 34;
    
    if vaccinated
        score = score * 0.25;
    end
    if handWash
        score = score * 0.8;
    end
    
    healthyScore(h) = score;
end


healthyNums = [healthyNum];
%% Infected People Properties 
infectedNum =  2; % number of immune cells at start of simulation
incubationOptions = [2; 3; 4; 4; 5; 5; 5; 6; 6; 7; 8];

for i = 1:infectedNum
    infectedLocation(i, :) = PlaceCell(areaSize); %places each immune cell randomly within the space
    infectedAngle(i,:) = randi(360);
    incubationPeriod(i,:) = incubationOptions(randi(11))*8;
end

infectedNums = [infectedNum];
%% Recovered People Properties

recoveredNum = 0;
recoveredLocation = [];
recoveredAngle = [];

recoveredNums = [recoveredNum];
%% Plot the initial locations of all agents
figure(1)
plot(healthyLocation(:, 1), healthyLocation(:, 2), 'ok')
hold on
plot(infectedLocation(:, 1), infectedLocation(:,2), 'or')
xlim([-areaSize, areaSize])
ylim([-areaSize, areaSize])
hold off


months = ["November" "December" "January" "February" "March"];
%% Simulation 
times = [1];
for t = 1:dt:total_time %loop through each instant in time based on tie step and simulation time
	% Update each virus in simulation
    for h = 1:healthyNum
        if(h > healthyNum)
            break;
        end
        
        angleIncrement = randi(30);
        angleDirection = randi(2);
        
        if(angleDirection == 1)
            angleIncrement = angleIncrement * -1;
        end
        
        healthyAngle(h,:) = healthyAngle(h,:) + angleIncrement;
        
        %make sure child does not go out of bounds
        if(OutOfBounds(healthyLocation(h,:) + [cosd(healthyAngle(h,:)), sind(healthyAngle(h,:))] .* speed, areaSize))
            healthyLocation(h,:) = PlaceCell(areaSize);
        end
            
        healthyLocation(h,:) = healthyLocation(h,:) + [cosd(healthyAngle(h,:)), sind(healthyAngle(h,:))] .* speed;
        
        for i = 1:infectedNum
            if(h > healthyNum)
                break;
            end
            distance = ComputeDistance(healthyLocation(h,:), infectedLocation(i,:));    % computing distance between healthy people and infected people
            if distance <= 2 && randi(100) < healthyScore(h)  % if infected cell is within a distance of 5
                infectedLocation = [infectedLocation; healthyLocation(h,:)];     % mark person as infected, remove from healthy
                infectedAngle = [infectedAngle; healthyAngle(h,:)];     % mark person as infected, remove from healthy
                incubationPeriod = [incubationPeriod; incubationOptions(randi(11))*8]; % generate an incubation period
                infectedNum = infectedNum + 1;
                healthyLocation(h,:) = [];
                healthyScore(h) = [];
                healthyAngle(h,:) = [];
                healthyNum = healthyNum - 1;
            end
        end
        
        % random chance an individual gets infected
        randomChance = randi(5000);
        if randomChance == 1
            fprintf("randomly infected!\n");
            infectedLocation = [infectedLocation; healthyLocation(h,:)];     % mark person as infected, remove from healthy
            infectedAngle = [infectedAngle; healthyAngle(h,:)];     % mark person as infected, remove from healthy
            incubationPeriod = [incubationPeriod; incubationOptions(randi(11))*16]; % generate an incubation period
            infectedNum = infectedNum + 1;
            healthyLocation(h,:) = [];
            healthyScore(h) = [];
            healthyAngle(h,:) = [];
            healthyNum = healthyNum - 1;
        end       
        
    end
    
    for i = 1:infectedNum
        
        if(i > infectedNum)
            break;
        end
        
        incubationPeriod(i,:) = incubationPeriod(i,:) - 1; % counting down the incubation period
        
        if(incubationPeriod(i,:) == 0) %if incubation period is over
            recoveredLocation = [recoveredLocation; infectedLocation(i,:)]; % add person to recovered list
            recoveredAngle = [recoveredAngle; infectedAngle(i,:)];
            recoveredNum = recoveredNum + 1;
            infectedLocation(i,:) = []; % remove from infected list
            infectedAngle(i,:) = []; % remove from infected list
            incubationPeriod(i,:) = [];   % remove from incubation period list
            infectedNum = infectedNum - 1;
        end
        
        if(i > infectedNum)
            break;
        end
        
        angleIncrement = randi(30);
        angleDirection = randi(2);
        
        if(angleDirection == 1)
            angleIncrement = angleIncrement * -1;
        end
        
        infectedAngle(i,:) = infectedAngle(i,:) + angleIncrement;
        
        %make sure child does not go out of bounds
        if(OutOfBounds(infectedLocation(i,:) + [cosd(infectedAngle(i,:)), sind(infectedAngle(i,:))] .* speed, areaSize))
            infectedLocation(i,:) = PlaceCell(areaSize);
        end
        
        infectedLocation(i,:) = infectedLocation(i,:) + [cosd(infectedAngle(i,:)), sind(infectedAngle(i,:))] .* speed;
    end
    
    for r = 1:recoveredNum
        angleIncrement = randi(30);
        angleDirection = randi(2);
        
        if(angleDirection == 1)
            angleIncrement = angleIncrement * -1;
        end
        
        recoveredAngle(r,:) = recoveredAngle(r,:) + angleIncrement;
        
        %make sure child does not go out of bounds
        if(OutOfBounds(recoveredLocation(r,:) + [cosd(recoveredAngle(r,:)), sind(recoveredAngle(r,:))] .* speed, areaSize))
            recoveredLocation(r,:) = PlaceCell(areaSize);
        end
        
        recoveredLocation(r,:) = recoveredLocation(r,:) + [cosd(recoveredAngle(r,:)), sind(recoveredAngle(r,:))] .* speed;
    end
    
    % pause(0.075) %pause gives a delay in updating our plot so we can see each time step

    % Plot the new locations and state of our system
    
    healthyNums(t,:) = healthyNum;
    infectedNums(t,:) = infectedNum;
    recoveredNums(t,:) = recoveredNum;  
    times(t,:) = t;
    
    figure(1)
    if healthyNum >=1
        plot(healthyLocation(:, 1), healthyLocation(:, 2), 'ok')
    end
    hold on
    
    if infectedNum >= 1
        plot(infectedLocation(:, 1), infectedLocation(:,2), 'or')
    end
    
    if recoveredNum >=1
        plot(recoveredLocation(:, 1), recoveredLocation(:,2), 'og')
    end
        
    if(t > total_time)
        fprintf("Simulation has ended!");
        break; 
    end
    
    text(0.5*areaSize,0.9*areaSize,['infected children = ' num2str(infectedNum)])
    text(0.5*areaSize,0.8*areaSize,['healthy children = ' num2str(healthyNum)])
    text(0.5*areaSize,0.7*areaSize,['recovered children = ' num2str(recoveredNum)])
    text(0.5*areaSize,0.6*areaSize,['t = ' num2str(t)])
    text(0.5*areaSize,0.5*areaSize,['month = ' num2str(months(floor(t/(30*8))+1))])
    xlim([-areaSize, areaSize])
    ylim([-areaSize, areaSize])
    
    hold off
    
    figure(2)
    
    hold on
    
    plot(times, healthyNums, 'ok');
    plot(times, infectedNums, 'or');
    plot(times, recoveredNums, 'og');
    
    xlabel('Time, hours');
    ylabel('Number, people');
    ylim([0, 30])
    
    legend("Healthy", "Infected", "Recovered");
    
    hold off
end

%% Function to Place Agents Randomly within the Simulation Area
function location = PlaceCell(areaSize)
    location = [rand(),rand()] .* 2 .* areaSize - areaSize;
end

%% Function to compute the distance between two agents
function distance = ComputeDistance(point1, point2)
   distance = sqrt( (point1(1) - point2(1))^2 + (point1(2) - point2(2)).^2 );
end

function oob = OutOfBounds(point, areaSize)
    oob = ~(point(1) > -1*areaSize & point(1) < areaSize & point(2) > -1*areaSize & point(2) < areaSize);
end

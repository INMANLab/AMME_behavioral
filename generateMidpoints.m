
function [midpointtableable] = generateMidpoints(table) 
% for a list of xyz coordinates of anodes and cathods, calculate the xyz
% coordinates of midpoint given a vector of labels [-1,0,1]
%the midpoint will contain an x,y,z coordinate, and a name

%electrode label will concatonate the two labels of the anode and cathode.

sessions = unique(table.sessionID);

midpointtableable = cell([0 5]);

figure;

for i = 1:length(sessions)

currentSession = sessions(i);

pair = table(ismember(table.sessionID,currentSession),:);

%v is the closest voxel of HPC here
v(1) = pair.closest_PRC_coord_RAS_X(1);
v(2) = pair.closest_PRC_coord_RAS_Y(1);
v(3) = pair.closest_PRC_coord_RAS_Z(1);

%this is for monopolar stim
if size(pair,1) == 1
        %add x,y,z, coordinates
midpoint = table2array(pair(:,3:5));
label = pair.stim_electrode_label;

%for bipolar stim
elseif size(pair,1) == 2

anodeIDX = find(pair.label == -1);
cathodeIDX = find(pair.label == 1);

anode = pair(anodeIDX,:);
anodeXYZ = table2array(anode(:,3:5));
cathode = pair(cathodeIDX,:);
cathodeXYZ = table2array(cathode(:,3:5));
midpoint = (anodeXYZ + cathodeXYZ)/2;
label = strcat(pair.stim_electrode_label(1), ',', pair.stim_electrode_label(2));

end



%plot the midpoint
plot3(midpoint(1),midpoint(2),midpoint(3),'o','MarkerSize',10,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor','none')%Midpoint will be black
hold on
%now plot voxe 1
plot3(v(1),v(2),v(3),'o','MarkerSize',10,'MarkerFaceColor',[0,0,1],'MarkerEdgeColor','none')%we can just make it blue for now

%midpoint to voxel 1
plot3([midpoint(1) v(1)],[midpoint(2) v(2)],[midpoint(3) v(3)],'-','LineWidth',2,'color','k')

%calculate distance

distance = sqrt((midpoint(1)-v(1))^2 + (midpoint(2)-v(2))^2 + (midpoint(3)-v(3))^2);

midpointtableable{i,1} = currentSession;
midpointtableable{i,2} = label;
midpointtableable{i,3} = midpoint(1);
midpointtableable{i,4} = midpoint(2);
midpointtableable{i,5} = midpoint(3);
midpointtableable{i,6} = distance;

end

midpointtableable=array2table(midpointtableable,'VariableNames',{'sessionID','electrodeLabel','midPointX','midPointY','midPointZ', 'distance'});
writetable(midpointtableable, 'PRC_midpointtableable.csv');
end



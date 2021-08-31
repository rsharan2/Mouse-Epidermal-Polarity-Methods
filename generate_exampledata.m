% Simulate a data set with x,y points generated randomly from squares on a grid. 
% Choose between several options for angle data.

% Set values for grid size (for x and y coordinates)
% grid = 490;


function Data=generate_exampledata(grid,angledist) 

Rad3Over2 = sqrt(3) / 2;
[X Y] = meshgrid(0:1:41);
n = size(X,1);

X = Rad3Over2 * X;
Y = Y + repmat([0 0.5],[size(X,1) size(X,2)/2]);

X_big = X*grid;
Y_big = Y*grid;

x_og = X_big(X_big <= 6000 & Y_big <= 2000);
y_og = Y_big(X_big <= 6000 & Y_big <= 2000);

x = x_og + 100.*rand(length(x_og),1) + 1;
y = y_og + 100.*rand(length(y_og),1) + 1;

%{
% if desired, initiate random number geneator to keep consistent results
% rng(1234,'twister');

% create semirandomized values of x and y with one point per square of a
% 'grid'x'grid' grid
% first generate modifiers to adjust random values 0-1 to the correct x
% and y ranges
modx = [];
mody = [];
modx(:,2) = 0:grid:6000;
modx(:,1) = grid;
mody(:,2)= 0:grid:2000;
mody(:,1) = grid;
% create a matrix of x and y values such that each corresponds to a random
% point within one square of the grid
x=rand(length(modx),length(mody)).*modx(:,1) + modx(:,2);
y=(rand(length(mody),length(modx)).*mody(:,1) + mody(:,2))';
% plot
% scatter(x(:),y(:), 'black','.')
% generate random angles within the desired range
%angle = mod(randn([length(modx),length(mody)]).*stdev+mean, 360);

% If regional differences are desired, select the third option below and then 
% choose a split (x or y, and the proportion of the anterior region)
%}
% If a normal distribution of angles is desired, select the mean and stdev to 
% describe this distribution and choose the frequency of downgrown follicles.

if angledist == 'Normal'
    % normal distribution of angles shaped by mean and stdev
    meanangle = 0;
    stdev = 5;
    f = 0.01;
    angle = mod(randn([length(modx),length(mody)]).*stdev+meanangle, 360);
    down = round(rand(length(modx),length(mody)).*(1+f)+.5);
    % concatenate all into one data table 'Data'
    Data = [reshape(x,[],1),reshape(y,[],1),reshape(angle,[],1), reshape(down,[],1)];


% If a random distribution of angles is desired, select the frequency of downgrown 
% follicles.

elseif angledist == 'Randomized'
    % random distribution of angles
    f = 0.2;
    angle = mod(rand([length(modx),length(mody)]).*360, 360);
    down = round(rand(length(modx),length(mody)).*(1+f)+.5);
    % concatenate all into one data table 'Data'
    Data = [reshape(x,[],1),reshape(y,[],1),reshape(angle,[],1), reshape(down,[],1)];

 
% If a regional distribution of angles is desired, select the 'split' option 
% that describes this distribution. This split can be either left-right, top-bottom, 
% or both (4-way), and can fit several proportions. Choose the distribution of 
% angles (select '1' if you do not want a split in a given direction) and the 
% frequency of downgrown follicles in each region. For the 4-way option, these 
% parameters will be averaged (i.e. the top left corner will be the average of 
% the top (y1) and left (x1) parameters.
% 
% 
% 
% Left-right split:

elseif angledist == "Left-right split"
    % split distribution of angles
    splitx = 2/3;
    meanx1 = 0;
    stdevx1 =5;
    fx1 = 0;
    meanx2 = 180;
    stdevx2 =40;
    fx2 = 0.25;
    
    angle = [(mod(randn([ceil(length(modx).*splitx), length(mody)]) .* stdevx1 + meanx1, 360));
        (mod(randn([floor(length(modx).*(1-splitx)), length(mody)]) .* stdevx2 + meanx2, 360))];
    
    down = round([ (rand([ceil(length(modx).*splitx), length(mody)]) .* (1+fx1) + .5);
        (rand([floor(length(modx).*(1-splitx)), length(mody)]) .* (1+ fx2) + .5) ]);
    % concatenate all into one data table 'Data'
    Data = [reshape(x,[],1),reshape(y,[],1),reshape(angle,[],1), reshape(down,[],1)];
   
% Top-bottom split:

elseif angledist == "Top-bottom split"
    % split distribution of angles
    splity = 2/3;
    meany1 = 0;
    stdevy1 =5;
    fy1 = 0.01;
    meany2 = 0;
    stdevy2 =5;
    fy2 = 0.01;
    
    angle = [ (mod(randn([length(modx), ceil(length(mody).*splity)]) .* stdevy1 + meany1, 360)) (mod(randn([length(modx), floor(length(mody).*(1-splity))]) .* stdevy2 + meany2, 360))];
    
    down = round([ (rand([length(modx), ceil(length(mody).*splity)]) .* (1+fy1) + .5) (rand([length(modx), floor(length(mody).*(1-splity))]) .* (1+fy2) + .5) ]);
    % concatenate all into one data table 'Data'
    Data = [reshape(x,[],1),reshape(y,[],1),reshape(angle,[],1), reshape(down,[],1)];
   

% 4-way split:

elseif angledist == '4-way split'
    % split distribution of angles
    splitx = 2/3;
    splity = 1/2;
    meanx1 = 0;
    stdevx1 =10;
    fx1 = 0.01;
    meanx2 = 180;
    stdevx2 =30;
    fx2 = 0.25;
    meany1 = 30;
    stdevy1 =10;
    fy1 = 0.01;
    meany2 = 330;
    stdevy2 =10;
    fy2 = 0.01;
    
    angle = [ (mod(randn([ceil(length(modx).*splitx), ceil(length(mody).*splity)]) .* mean([stdevx1 stdevy1]) + mean([meanx1 meany1]), 360)) (mod(randn([ceil(length(modx).*splitx), floor(length(mody).*(1-splity))]) .* mean([stdevx1 stdevy2]) + mean([meanx1 meany2]), 360));
        (mod(randn([floor(length(modx).*(1-splitx)), ceil(length(mody).*splity)]) .* mean([stdevx2 stdevy1]) + mean([meanx2 meany1]), 360)) (mod(randn([floor(length(modx).*(1-splitx)), floor(length(mody).*(1-splity))]) .* mean([stdevx2 stdevy2]) + mean([meanx2 meany2]), 360))];
    
    down = round([ (rand([ceil(length(modx).*splitx), ceil(length(mody).*splity)]) .* (1+mean(fx1+fy1)) + .5) (rand([ceil(length(modx).*splitx), floor(length(mody).*(1-splity))]) .* (1+mean(fx1+fy2)) + .5);
        (rand([floor(length(modx).*(1-splitx)), ceil(length(mody).*splity)]) .* (1+ mean(fx2+fy1)) + .5) (rand([floor(length(modx).*(1-splitx)), floor(length(mody).*(1-splity))]) .* (1+mean(fx2+fy2)) + .5) ]);
    % concatenate all into one data table 'Data'
    Data = [reshape(x,[],1),reshape(y,[],1),reshape(angle,[],1), reshape(down,[],1)];

% Wildtype angle distribution:

elseif angledist == "Wildtype"
    % concatenate all into one data table 'Data'
    angle = load('sim_data.mat','wt_angle');
    angle = angle.wt_angle;
    down = load('sim_data.mat','wt_down');
    down = down.wt_down;
    Data = [reshape(x,[],1),reshape(y,[],1),reshape(angle,[],1), reshape(down,[],1)];

% Angle distribution where there are local zones of reversal:

elseif angledist == "local_order"
    angle = load('sim_data.mat','local_order_angle');
    angle = angle.local_order_angle;
    down = load('sim_data.mat','local_order_down');
    down = down.local_order_down;
    Data = [reshape(x,[],1),reshape(y,[],1),reshape(angle,[],1), reshape(down,[],1)];

% Angle distribution in a mutant with disrupted PCP activity:

elseif angledist == "pcp_mut"
    angle = load('sim_data.mat','pcp_mut_angle');
    angle = angle.pcp_mut_angle;
    down = load('sim_data.mat','pcp_mut_down');
    down = down.pcp_mut_down;
    Data = [reshape(x,[],1),reshape(y,[],1),reshape(angle,[],1), reshape(down,[],1)];
end

% histogram(reshape(angle,[],1), 'NumBins', 80)
% assign some subset of points with frequence 'f' as downgrown
% concatenate all into one data table 'Data'
% reassign angle of downgrown points as -1
Data(Data(:,4)==2,3)= -1;

end
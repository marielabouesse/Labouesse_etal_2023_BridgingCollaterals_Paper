function  [angle_deg_final] = rotations_function_forpaper(neck_x,neck_y,tail_x,tail_y,time_vect,fiber_side,dt_ds)

% using head snout and tail basis

% inputs: 
% mouse body position list of body center (near neck) and tail base at all times t
% column 1: neck_x position neck 
% column 2: neck_y position neck
% column 3: tail_x position tail
% column 4: tail_y position tail

% EXAMPLES
% neck_x = [3;0;-3;0;-2]
% neck_y = [0;3;0;-3;-2]
% tail_x = [1;0;-1;0;-0.5]
% tail_y = [0;1;0;-1;-0.5]

% they were already corrected against the reference (0,0) lower left corner of open field
% they are in cm

% outputs: 
% column vector with angule at all times t, in degrees

% first we calculate the vectors 

% bodyvector is the vector of the mouse body position at time t
% we calculate it by substracting neck x and tail x; and neck y and tail y
% bodyvector has two columns, one is x and other is y coordinates of the vector
% each row corresponds to different t
% eg row 1: x1, y1

% convert to columns
neck_x = neck_x';
neck_y = neck_y';
tail_x = tail_x';
tail_y = tail_y';

bodyvector = [neck_x - tail_x neck_y - tail_y]; % on normalize par rapport a la tail
% we need to normalize the vector so that it can calculate the angle; cos you need unit vectors
bodyvector_normalized = bodyvector./sqrt(sum((bodyvector).*(bodyvector),2));

        
%% formula to calculate angle
	% angle = arctan2(x1y2 - y1x2, x1x2 + y1y2) where x1y1 is a multiplication and you have 2 arguments, separated by a comma
    % these angles are in rad 
    % The four-quadrant inverse tangent, atan2(Y,X), returns values in the closed interval [-pi,pi] based on the values of Y and X, as shown in the graphic.
    % in trig convention, anticlockwise is positive 
    
angle_rad=NaN(length(bodyvector_normalized)-1,1);
for i=1:length(bodyvector_normalized)-1
    angle_rad(i) = atan2(bodyvector_normalized(i,1)*bodyvector_normalized(i+1,2) - bodyvector_normalized(i,2)*bodyvector_normalized(i+1,1),bodyvector_normalized(i,1)*bodyvector_normalized(i+1,1) + bodyvector_normalized(i,2)*bodyvector_normalized(i+1,2));
%     angle_rad(i) = bodyvector_normalized(i,1)*bodyvector_normalized(i+1,2) - bodyvector_normalized(i,2)*bodyvector_normalized(i+1,1)
end
   
%% to convert from rad to deg:
angle_deg = rad2deg(angle_rad);

% plot
figure; clf; hold on;
plot(time_vect(1:end-1),angle_deg');
% end
xlim([-1 time_vect(end)])
ylim([-180 180])
title('Angle at time t')
xlabel('Time (sec)')
ylabel('Angle (deg)')
    



%%%long camber results guide


%cycle through pitch angles
%get camber from idealized/actual curve for a given pitch angle (based off
%wheel displacement)
%calculate max accel at those given cambers. compare to max accel normally
%with graphs
%(assume neutral camber)

pitch_vector = -5:.5:5;
for i = 1:numel(pitch_vector)
    pitch_angle = pitch_vector(i);
    if pitch<0
        STC_brake
    else
        STC_accel
    end
    calcWheelForcesAndDisplacements_long(pitch_angle,C)
    %if pitch =>braking, all wheel forces included, else only back 2
    singleTireCamberEvaluation(zero_camber)
    %store as ideal vector for comparison
    %-> ideal camber curve generated
    singleTireCamberEvaluation(ideal_camber)

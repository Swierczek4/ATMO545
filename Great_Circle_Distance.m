function dist = Great_Circle_Distance(A,B)
%% Gives great circle distance
% A and B must be 2 number arrays where the
% first number is latitude in degrees and the
% second number is longitude in degrees.
%%

A = A.*pi./180;
B = B.*pi./180;

X = cos(A(1))*cos(B(1))*cos(A(2)-B(2)) + sin(A(1))*sin(B(1));

dist = 6378*acos(X);
end


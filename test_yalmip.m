import casadi.*

% Declare variables
x = SX.sym('x');
y = SX.sym('y');
z = SX.sym('z');

p = SX.sym('p');

% Formulate the NLP
f = x^2 + p*100*z^2;
g = [z + (1-x)^2 - y;x-y];
nlp = struct('x',[x;y;z],'p',p,'f',f, 'g',g);


data = struct('x0',[2.5;3.0;0.75],'lbg',0,'ubg',0,'p',3);

casadi2yalmip(nlp,data,[false, false, false]);
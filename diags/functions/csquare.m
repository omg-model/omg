function rgb = csquare(n)

x=[0 0 0.5 1 1]'; % x cordinates of corners and centre
y=[0 1 0.5 0 1]'; % y cordinates of corners and centre
r=[0 1 1 0 1]';   % corresponding r values
g=[1 1 1 0 0]';   % corresponding g values
b=[0 0 1 1 1]';   % corresponding b values

xl=linspace(0,1,n); % unique values of trait x
yl=linspace(0,1,n); % unique values of trait y

[xx yy]=meshgrid(xl,yl); % generate trait grid

Fr=scatteredInterpolant(x,y,r); % r values across grid
Fg=scatteredInterpolant(x,y,g); % g values across grid
Fb=scatteredInterpolant(x,y,b); % b values across grid

rgb=[Fr(xx(:),yy(:)) Fg(xx(:),yy(:)) Fb(xx(:),yy(:))]; % rgb array

end
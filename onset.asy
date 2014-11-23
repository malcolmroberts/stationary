import graph;
size(40cm,10cm,false);


file fin = input("output.csv").line();
real[][] a = fin.dimension(0, 0);
a = transpose(a);
//write(a[0].length);

real[] ay = a[1];
real[] ax = sequence(ay.length);

draw(graph(ax, ay), invisible);

void drawbar(real x, real y, pen p=defaultpen, pen f=invisible) {
  real w = 1.0/3.0;
  pair p0 = (x - w, 0);
  pair p1 = (x + w, 0);
  pair p2 = (x + w, y);
  pair p3 = (x - w, y);

  filldraw(p0--p1--p2--p3--cycle,f,p);
}
  
for(int i=0; i < ay.length; ++i) {
  drawbar(ax[i],ay[i],red,red+opacity(0.25));
  if(ay[i] >= 0) label("s",(ax[i], ay[i]), NW, red);
}


file fref = input("ref.csv").line();
real[][] b = fref.dimension(0, 0);
b = transpose(b);
real[] by = b[0];

real[] bx = sequence(b[0].length);
draw(graph(bx, by), invisible);

for(int i=0; i < by.length; ++i) {
  drawbar(bx[i], by[i], blue, blue + opacity(0.25));
  if(by[i] >= 0)
    label("s",(bx[i], by[i]), NE, blue);
}

xaxis(BottomTop);

for(int i=0; i < bx.length; ++i) {
  label(string(i+1), bx[i], S);
}


yaxis("Onset of stationarity", LeftRight, RightTicks);

real xmean = 0.5*(min(ax[0],bx[0]) + max(ax[ax.length -1],bx[bx.length -1]));
real ymax = max(max(ay), max(by));
label("Comparison between detected (red) and reference (blue) onset of stationarity.",(xmean, -0.1 * ymax));
label("The label ``s\" indicated stationarity was detected.",(xmean, -0.1 * ymax), 4*S);


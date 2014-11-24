import graph;
size(40cm,10cm,false);

real[] ay;
{
  file fin = input("output.csv").line();
  real[][] a = fin.dimension(0, 0);
  a = transpose(a);
  //write(a[0].length);
  ay = a[1];
}

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
  

real offset = 0.005;

real[] pobs;
{
  file fpref = input("np.csv").line();
  real[][] prr = fpref.dimension(0, 0);
  prr = transpose(prr);
  pobs = prr[0];
}

for(int i = 0; i < ay.length; ++i) {
  drawbar(ax[i], ay[i], red, red + opacity(0.25));
  if(ay[i] >= 0) {
    if(pobs[i] > 0)
      label("\small{sp}",(ax[i], ay[i]), offset * NW, red);
    else
      label("\small{s}",(ax[i], ay[i]), offset * NW, red);
  }
}

real[] by;
{
  file fref = input("ref.csv").line();
  real[][] b = fref.dimension(0, 0);
  b = transpose(b);
  by = b[0];
}
  
real[] bx = sequence(by.length);
draw(graph(bx, by), invisible);


real[] pref;
{
  file fpref = input("pref.csv").line();
  real[][] prr = fpref.dimension(0, 0);
  prr = transpose(prr);
  pref = prr[0];
}

for(int i = 0; i < by.length; ++i) {
  drawbar(bx[i], by[i], blue, blue + opacity(0.25));
  if(by[i] >= 0) {
    if(pref[i] == 1)
      label("\small{sp}",(bx[i], by[i]), offset * NE, blue);
    else 
      label("\small{s}", (bx[i], by[i]), offset * NE, blue);
  }
}

xaxis(BottomTop);

for(int i=0; i < bx.length; ++i) {
  label(string(i+1), bx[i], S);
}


yaxis("Onset of stationarity", LeftRight, RightTicks);

real xmean = 0.5*(min(ax[0],bx[0]) + max(ax[ax.length -1],bx[bx.length -1]));
real ymax = max(max(ay), max(by));
label("Comparison between detected (red) and reference (blue) onset of stationarity.",(xmean, -0.1 * ymax));
label("The label ``s\" indicated stationarity was detected, ``sp\" indicates a periodic stationary state.",(xmean, -0.1 * ymax), 4*S);


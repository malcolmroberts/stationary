import graph;
size(40cm, 10cm, false);

void drawbar(real x, real y, pen p = defaultpen, pen f = invisible) {
  real w = 1.0 / 3.0;
  pair p0 = (x - w, 0);
  pair p1 = (x + w, 0);
  pair p2 = (x + w, y);
  pair p3 = (x - w, y);
  filldraw(p0--p1--p2--p3--cycle, f, p);
}

// offset for placement of sp text
real offset = 0.005;

// end of warmup from the observations
real[] sobs;
{
  file f = input("output.csv").line();
  real[][] a = f.dimension(0, 0);
  a = transpose(a);
  //write(a[0].length);
  sobs = a[1];
}

// presence of periodicity as per observations
real[] pobs;
{
  file f = input("np.csv").line();
  real[][] prr = f.dimension(0, 0);
  prr = transpose(prr);
  pobs = prr[0];
}

real[] x = sequence(sobs.length);

// draw an invisible graph to set up the axes
draw(graph(x, sobs), invisible);

for(int i = 0; i < sobs.length; ++i) {
  drawbar(x[i], sobs[i], red, red + opacity(0.25));
  if(sobs[i] >= 0) {
    if(pobs[i] > 0)
      label("\small{sp}",(x[i], sobs[i]), offset * NW, red);
    else
      label("\small{s}",(x[i], sobs[i]), offset * NW, red);
  }
}


// reference end of warmup
real[] wref;
{
  file f = input("ref/ref.csv").line();
  real[][] b = f.dimension(0, 0);
  b = transpose(b);
  wref = b[0];
}
// invisible graph to set up axes
real[] rx = sequence(wref.length);
draw(graph(rx, wref), invisible);


// presence of periodicity (reference)
real[] pref;
{
  file f = input("ref/pref.csv").line();
  real[][] prr = f.dimension(0, 0);
  prr = transpose(prr);
  pref = prr[0];
}

// presence of stationarity (reference)
real[] sref;
{
  file f = input("ref/sref.csv").line();
  real[][] a = f.dimension(0, 0);
  a = transpose(a);
  sref = a[0];
}

for(int i = 0; i < wref.length; ++i) {
  drawbar(rx[i], wref[i], blue, blue + opacity(0.25));
  if(sref[i] == 1) {
    if(pref[i] == 1)
      label("\small{sp}",(rx[i], wref[i]), offset * NE, blue);
    else 
      label("\small{s}", (rx[i], wref[i]), offset * NE, blue);
  }
}

xaxis(BottomTop);

for(int i=0; i < rx.length; ++i) {
  label(string(i+1), rx[i], S);
}

yaxis("End of warmpu period", LeftRight, RightTicks);
real xmean = 0.5 * (min(x[0],rx[0]) + max(x[x.length -1], rx[rx.length -1]));
real ymax = max(max(sobs), max(wref));
label("Comparison between detected (red) and reference (blue) onset of stationarity.", (xmean, -0.1 * ymax));
label("The label ``s\" indicated stationarity was detected, ``sp\" indicates a periodic stationary state.", (xmean, -0.1 * ymax), 4 * S);

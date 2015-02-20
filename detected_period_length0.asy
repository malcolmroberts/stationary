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

real[] efficiency;
{
  string filename;
  filename = getstring("FSfilename");
  file f = input(filename).line();
  real[][] a = f.dimension(0, 0);
  a = transpose(a);
  for(int i = 0; i < 62; ++i) {
    efficiency.push(a[17][i]); // R, #CycleLength
  }
}

real[] x = sequence(efficiency.length);

// draw an invisible graph to set up the axes
draw(graph(x, efficiency), invisible);

for(int i = 0; i < efficiency.length; ++i) {
  drawbar(x[i], efficiency[i], red, red + opacity(0.25));
}

xaxis(BottomTop);


real[] rx = x;

for(int i=0; i < rx.length; ++i) {
  label(string(i+1), rx[i], S);
}

label("Run number", 0.5 * (min(rx) + max(rx)), 4*S);

yaxis("\#CycleLength", LeftRight, RightTicks);

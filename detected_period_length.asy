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

real[] period_lengths;
{
  string filename = getstring("period lengths file"); 
  file f = input(filename).line();
  real[][] a = f.dimension(0, 0);
  a = transpose(a);
  period_lengths = a[0];
}

real[] x = sequence(period_lengths.length);

// draw an invisible graph to set up the axes
draw(graph(x, period_lengths), invisible);

for(int i = 0; i < period_lengths.length; ++i) {
  drawbar(x[i], period_lengths[i], red, red + opacity(0.25));
}

xaxis(BottomTop);


real[] rx = x;

for(int i=0; i < rx.length; ++i) {
  label(string(i+1), rx[i], S);
}

label("Run number", 0.5 * (min(rx) + max(rx)), 4*S);

yaxis("Detected Period Length", LeftRight, RightTicks);

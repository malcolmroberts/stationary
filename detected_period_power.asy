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

string period_powers_filename = "";
string nonperiod_powers_filename = "";

usersetting();

real[] period_powers;
{
  if(period_powers_filename == "")
    period_powers_filename = getstring("period powers file"); 
  file f = input(period_powers_filename).line();
  real[][] a = f.dimension(0, 0);
  a = transpose(a);
  period_powers = a[0];
}

real[] nonperiod_powers;
{
  if(nonperiod_powers_filename == "") 
    nonperiod_powers_filename = getstring("nonperiod powers file"); 
  file f = input(nonperiod_powers_filename).line();
  real[][] a = f.dimension(0, 0);
  a = transpose(a);
  nonperiod_powers = a[0];
}

real[] x = sequence(period_powers.length);

real[] relative_power;
for(int i = 0; i < period_powers.length; ++i) {
  real pp = period_powers[i];
  real npp = nonperiod_powers[i];
  relative_power.push(pp / (pp + npp + 1e-16));
}

// draw an invisible graph to set up the axes
draw(graph(x, relative_power), invisible);

for(int i = 0; i < relative_power.length; ++i) {
  drawbar(x[i], relative_power[i], red, red + opacity(0.25));
}

xaxis(BottomTop);


real[] rx = x;

for(int i=0; i < rx.length; ++i) {
  label(string(i+1), rx[i], S);
}

label("Run number", 0.5 * (min(rx) + max(rx)), 4*S);

yaxis("Relative power of the periodic component", LeftRight, RightTicks);

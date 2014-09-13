import graph;

size(300,200,IgnoreAspect);

string filenames="";//=getstring("filenames");
string sscale="";
string xlabel="";
string ylabel="";
bool logscale=false;

bool x95=false;
usersetting();

if(filenames == "")
  filenames=getstring("filenames");
if(sscale == "")
  sscale=getstring("scale linlin or loglin or linlog or loglog");

if(sscale == "loglog") {
  scale(Log,Log);
  logscale=true;
}
if(sscale == "loglin") {
  scale(Log,Linear);
  logscale=true;
}

if(sscale == "linlog") {
  scale(Linear,Log);
  logscale=true;
}


string filename;
int n=-1;
bool flag=true;
int lastpos;
while(flag) {
  ++n;
  int pos=find(filenames,",",lastpos);
  if(lastpos == -1) {filename=""; flag=false;}
  filename=substr(filenames,lastpos,pos-lastpos);
  if(flag) {
    write(filename);
    lastpos=pos > 0 ? pos+1 : -1;
    file fin=input(filename).line();
    
    real[][] a=fin.dimension(0,0);

    write(a[0].length);

    a=transpose(a);
    real[] x=a[0];
    real[] y=a[1];

    real xmin=min(x);
    //write(xmin);
    
    if(logscale) {
      draw(graph(x,y,x>0),Pen(n),texify("data file"));
    } else {
      real nf=1.96/sqrt(x.length);
      write();

      //-1.96/sqrt(length(x));
      draw(graph(x,y,x<9017),Pen(n),texify("data file"));

      if(x95) {
	// draw 95% confidence intervals
	yequals(nf,grey);
	yequals(0,grey);
	yequals(-nf,grey);
      }

      /*
      real period;
      //period=570.393251702;//595.0344806;
      period=566.502167343;
      period=533.531741204;
      period=575.541778777;
      //period=559.023729726;
      // period=555.052491453;
      //period=568.750385788;
      int xeq=1;
      while(xeq * period < max(x)) {
	xequals(period*xeq,grey);
	xeq += 1;
      }
      */
      
    }
  }
}

if(xlabel == "")
  xlabel=getstring("xlabel");
if(ylabel == "")
  ylabel=getstring("ylabel");

xaxis(xlabel,BottomTop,LeftTicks);
yaxis(ylabel,LeftRight,RightTicks);

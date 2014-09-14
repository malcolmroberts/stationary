import graph;

size(300,200,IgnoreAspect);

string filenames="";//=getstring("filenames");
string sscale="";
string xlabel="";
string ylabel="";
bool logscale=false;
int end=0;
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


bool end_set=false;

if(end != 0)
  end_set=true;

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

    a=transpose(a);
    real[] x=a[0];
    real[] y=a[1];

    if(!end_set)
      end=y.length;
    else
      end=min(end,y.length);
    real xmin=min(x);
    //write(xmin);
    
    if(logscale) {
      draw(graph(x,y,x>0),Pen(n),texify("data file"));
    } else {
      real nf=1.96/sqrt(x.length);
      
      //-1.96/sqrt(length(x));
      draw(graph(x,y,x<end),Pen(n),texify("data file"));

      if(x95) {
	// draw 95% confidence intervals
	yequals(nf,grey);
	yequals(0,grey);
	yequals(-nf,grey);
      }
      
    }
  }
}

if(xlabel == "")
  xlabel=getstring("xlabel");
if(ylabel == "")
  ylabel=getstring("ylabel");

xaxis(xlabel,BottomTop,LeftTicks);
yaxis(ylabel,LeftRight,RightTicks);

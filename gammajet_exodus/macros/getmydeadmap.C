#include <fstream>
#include <TH1.h>

void getmydeadmap()
{
  invertfile = new TFile("/direct/phenix+workarea/hgong/cvs/macros/allinvert.root","READ");
  TH2D *invert_west = (TH2D*)invertfile->Get("totalinvert_west");
  TH2D *invert_east = (TH2D*)invertfile->Get("totalinvert_east");

  ofstream fout("mydeadmap.lst");

  int gbiny = 36;
  int gbinx = 72;
  int sec = 0;
  
  double threshold = 0.0;
  int arm = 0;
  int ly;

  for( int i = 1; i <= 4*gbiny; i++ )
    {
      for( int j = 1; j <= gbinx; j++ )
	{
	  if( invert_west->GetBinContent(j,i) > threshold )
	    {
	      if( i <= gbiny ) 
		{
		  sec = 3;
		  ly = i;
		}
	      
	      if( (i > gbiny) && (i <= 2*gbiny) )
		{
		  sec = 2;
		  ly = i - gbiny;
		}
	      
	      if( (i > 2*gbiny) && (i <= 3*gbiny) )
		{
		  sec = 1;
		  ly = i - 2*gbiny;
		}
	      
	      if( i > 3*gbiny )      
		{
		  sec = 0;
		  ly = i - 3*gbiny;
		}

	      fout << arm << " ";
	      fout << sec << " ";
	      fout << j-1 << " ";
	      fout << ly-1 << " " << endl;
	    }
	}
    }
      
  arm = 1;

  for( i = 1; i <= 4*gbiny; i++ )
    {
      for( j = 1; j <= gbinx; j++ )
	{
	  if(invert_east->GetBinContent(j,i) > threshold)
	    {
	      if(i>2*gbiny && i<=3*gbiny)
		{
		  sec = 2;
		  ly = i - 2*gbiny;
		}
	      if(i>3*gbiny)
		{
		  sec = 3;
		  ly = i - 3*gbiny;
		}
	      if(sec == 2 || sec == 3)
		{
		  fout << arm << " ";
		  fout << sec << " ";
		  fout << j-1 << " ";
		  fout << ly-1 << " " << endl;
		}
	    }
	}
    }

  fout.close();

}
